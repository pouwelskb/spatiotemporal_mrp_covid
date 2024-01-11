# compare AIC

library(ggplot2)
library(bayesplot)
theme_set(bayesplot::theme_default())
library(dplyr)
library(tidyr)
library(sp)
library(rgeos)
library(rgdal)
library(sf)
library(class)
library(KernSmooth)
library(INLA)
library(arm)
library(spdep)
library(httr)

# code below works with data grouped by covariates of interest and positive counts
# stored in the variable called result and number of tests stored in the variable named ntrials.

# furthermore, it requires a post-stratification table (poststratI) to be already loaded containing the number of people
# living in the target population, taking into account the conditional distribution of all covariates of interest.
#####################################################################################
# data grouped by 

##################################################################################

download.file("https://github.com/pouwelskb/shapefiles/raw/main/Covid_Infection_Survey__December_2020__UK_BUC-shp%20(3).zip", 
              destfile = "shapefile.zip" , mode='wb')
outdir <- getwd()
zipF <- "shapefile.zip"
unzip(zipF, exdir = outdir)
file.remove(zipF)

shape.file.cis <- readOGR("Covid_Infection_Survey__December_2020__UK_BUC.shp")

unique(shape.file.cis@data$CIS20CD)

## ONLY ENGLAND, for which we have vaccination data:
shape.file.cis <- shape.file.cis[!shape.file.cis@data$CIS20CD %in% 
                                   c("J06000217",
                                     "J06000218",
                                     "J06000219",
                                     "J06000220",
                                     "J06000221",
                                     "J06000222",
                                     "J06000223",
                                     "J06000224",
                                     "J06000225",
                                     "J06000226",
                                     "J06000227",
                                     "J06000228",
                                     "J06000229",
                                     "J06000230",
                                     "J06000231",
                                     "J06000232",
                                     "J06000233"),]

# check for any data:
names(shape.file.cis@data)[c(2)] <- c("CIS.name")

# shape.file.cis <- sp::merge(shape.file.cis,check,by="CIS.name")
# get objectid into numeric values:
shape.file.cis@data$objectid <- shape.file.cis@data$OBJECTID
shape.file.cis@data$cis_num <- shape.file.cis@data$objectid


# check coordinates of centroids of polygons:
coordinates(shape.file.cis)

maxXY <- pmax(bbox(shape.file.cis)[,2])
minXY <- pmin(bbox(shape.file.cis)[,1])

# plot(shape.file.cis, xlim=c(minXY[1],maxXY[1]), ylim=c(minXY[2],maxXY[2]))


# construct the neighbour list:
lanb <- poly2nb(shape.file.cis, row.names=1:length(shape.file.cis))
summary(shape.file.cis)
# make shape file with boundaries compatible with INLA:
nb2INLA("CISdistrictEngland_BUC.graph",lanb)
CISdistrictUK.adj <- "CISdistrictEngland_BUC.graph"

inla.setOption(scale.model.default=F)
H <- inla.read.graph(filename="CISdistrictEngland_BUC.graph")

############################################################################
# binary adjacency matrix: 
cis_num <- shape.file.cis@data$cis_num
carto <- sf::st_as_sf(shape.file.cis)

#################################################################
# code needed to generate constraints
# spatial neighbourhoud matrix Rs
Rs <- matrix(0, H$n, H$n)
for (i in 1:H$n) {
  Rs[i,i]=H$nnbs[[i]]
  Rs[i,H$nbs[[i]]]=-1
}


# number of spatial units
S <- length(lanb)

# number of age groups 
A <- length(unique(data_grouped$ageg_small))
V <- length(unique(data_grouped$vac_stat_any))
E <- length(unique(data_grouped$ethnicityg))

# temporal structure matrix:
T <- length(unique(data_grouped$weeknr))
dif1 <- 1 # for rw1
dif2 <- 2 # for rw2
D1 <- diff(diag(T), differences=dif1)
D2 <- diff(diag(T), differences=dif2)
Rt1 <- t(D1)%*%D1 #rw1
Rt2 <- t(D2)%*%D2 #rw2
dim(Rt1)
dim(Rt2)


# define constrain matrices:
# type I interaction
id_nt <- diag(1, nrow=(S*T))
dim(id_nt)

# note: I work here with PC priors
# type I interaction, no additional constrains

# type II interaction + RW1:
R_1_2 <- kronecker(Rt1, Diagonal(S))
r_def_1_2 <- S
A_constr_1_2 <- kronecker(matrix(1,1,T),Diagonal(S))
A_constr_1_2 <- as(A_constr_1_2[-1,],"matrix")
R_1_2_scaled <- R_1_2*exp(mean(log(diag(MASS::ginv(as.matrix(Rt1)))))) 

# II interaction + RW2:
R_2_2 <- kronecker(Rt2, Diagonal(S))
r_def_2_2 <- 2*S
A_constr_2_2 <- kronecker(matrix(1,1,T),Diagonal(S))
A_constr_2_2 <- A_constr_2_2[-1,]
R_2_2_scaled <- R_2_2*exp(mean(log(diag(MASS::ginv(as.matrix(Rt2))))))

# Type III + RW1/RW2 (same regardless of RW or RW2 or idd on time)
R_1_3 <- kronecker(diag(T),Rs)
r_def_1_3 <- T
A_constr_1_3 <- kronecker(diag(T),matrix(1,1,S))   # LCAR, DCAR, ICAR
A_constr_1_3 <- A_constr_1_3[-1,]
R_1_3_scaled<- R_1_3*exp(mean(log(diag(INLA:::inla.ginv(Rs))))) # BYM2                                   

# Type IV + RW1
R_1_4 <- kronecker(Rt1,Rs)
r_def_1_4 <- S+T-1
A.1.1 <- kronecker(matrix(1,1,T),diag(S))
A.1.2 <- kronecker(diag(T),matrix(1,1,S))
A_constr_1_4 <- rbind(A.1.1[-1,],A.1.2[-1,])
R_1_4_scaled<- R_1_4*exp(mean(log(diag(INLA:::inla.ginv(Rs)))))*
  exp(mean(log(diag(INLA:::inla.ginv(Rt1)))))

## Type IV + RW1
R_2_4 <- kronecker(Rt2,Rs)
r_def_2_4 <- 2*S+T-2
A.2.1 <- kronecker(matrix(1,1,T),diag(S))
A.2.2 <- kronecker(diag(T),matrix(1,1,S))
A_constr_2_4 <- rbind(A.2.1[-1,],A.2.2[-1,])
R_2_4_scaled<- R_2_4*exp(mean(log(diag(INLA:::inla.ginv(Rs)))))*
  exp(mean(log(diag(INLA:::inla.ginv(Rt2)))))

summary(data_grouped$weeknr.cis_num)


#rw1 + intI
form1_1 <- paste("f(cis_num, model='bym2', graph=Rs, constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01)),phi=list(prior='pc', param=c(0.5,0.5))))", sep="")
form1_1 <-  paste(form1_1,"+ f(weeknr, model='rw1', constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))", sep="")
form1_1 <- paste(form1_1,"+ f(weeknr_cisnum, model='iid', constr=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))", sep="")

#rw1 + intII:
form1_2 <- paste("f(cis_num, model='bym2', graph=Rs, constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01)),phi=list(prior='pc', param=c(0.5,0.5))))", sep="")
form1_2 <- paste(form1_2,"+ f(weeknr, model='rw1', constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))", sep="")
form1_2 <- paste(form1_2,"+ f(weeknr_cisnum, model='generic0', Cmatrix=R, rankdef=r.def, constr=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))), extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))", sep="")

# 
icar_formula_intII_rw1 = result ~ 1 + sex + vac_stat_any + 
  f(ethnicityg,
    model="iid",
    hyper = list(prec = list(prior = "pc.prec", param=c(1,0.1)))
  ) +
  f(ageg_small,
    model="iid",
    hyper= list(prec = list(prior="pc.prec", param=c(1,0.1))) # penalised complexicity (PC) prior for hyperprios
  ) +
  f(cis_num,
    model="bym2",
    graph=Rs,
    scale.model = TRUE,
    constr = TRUE,
    hyper = list(prec = list(prior="pc.prec", param=c(1,0.1)))
  ) +
  f(weeknr, model='rw1', constr=TRUE, scale.model=TRUE, 
    hyper=list(prec=list(prior='pc.prec', param=c(1,0.1)))) +
  f(weeknr.cis_num, model='generic0', Cmatrix=as.matrix(R_1_2_scaled,nrow=dim(R_1_2_scaled)[1]), 
    rankdef=r_def_1_2, constr=TRUE, 
    hyper=list(prec=list(prior='pc.prec', param=c(1,0.1))), 
    extraconstr=list(A=as.matrix(A_constr_1_2,nrow=dim(A_constr_1_2)[1]), e=rep(0,nrow(A_constr_1_2)))) + 
  
  f(weeknr.vac_stat_any, model="iid", 
    hyper = list(prec = list(prior = "pc.prec", param=c(1,0.01)))) +
  f(weeknr.age, model="iid",
    hyper = list(prec = list(prior = "pc.prec", param=c(1,0.01)))) + 
  f(weeknr.ethnicityg, model="iid", 
    hyper = list(prec = list(prior = "pc.prec", param=c(1,0.01))))  


mod_icar_intII_rw1 = inla(icar_formula_intII_rw1,
                      family="binomial",
                      data=data, Ntrials = data$ntrials,
                      control.predictor = list(compute=TRUE),
                      control.compute = list(waic=TRUE,
                                             config=TRUE),
                      control.inla = list(cmin = 0),
                      control.fixed = list(prec.intercept=0.01, prec=list(default=0.01)),
                      num.threads = 15,
                      verbose = TRUE
)

summary(mod_icar_intII_rw1)
mod_icar_intII_rw1$waic$waic 

# intII_rw2
icar_formula_intII_rw2 = result ~ 1 + sex + vac_stat_any + 
  f(ethnicityg,
    model="iid",
    hyper = list(prec = list(prior = "pc.prec", param=c(1,0.1)))
  ) +
  f(ageg_small,
    model="iid",
    hyper= list(prec = list(prior="pc.prec", param=c(1,0.1))) # penalised complexicity (PC) prior for hyperprios
  ) +
  f(cis_num,
    model="bym2",
    graph=Rs,
    scale.model = TRUE,
    constr = TRUE,
    hyper = list(prec = list(prior="pc.prec", param=c(1,0.1)))
  ) +
  f(weeknr, model='rw2', constr=TRUE, scale.model=TRUE, 
    hyper=list(prec=list(prior='pc.prec', param=c(1,0.1)))) +
  f(weeknr.cis_num, model='generic0', Cmatrix=as.matrix(R_2_2_scaled,nrow=dim(R_2_2_scaled)[1]), 
    rankdef=r_def_2_2, constr=TRUE, 
    hyper=list(prec=list(prior='pc.prec', param=c(1,0.1))), 
    extraconstr=list(A=as.matrix(A_constr_2_2,nrow=dim(A_constr_2_2)[1]), e=rep(0,nrow(A_constr_2_2)))) + 
  
  f(weeknr.vac_stat_any, model="iid", 
    hyper = list(prec = list(prior = "pc.prec", param=c(1,0.01)))) +
  f(weeknr.age, model="iid",
    hyper = list(prec = list(prior = "pc.prec", param=c(1,0.01)))) + 
  f(weeknr.ethnicityg, model="iid", 
    hyper = list(prec = list(prior = "pc.prec", param=c(1,0.01))))  


mod_icar_intII_rw2 = inla(icar_formula_intII_rw2,
                          family="binomial",
                          data=data, Ntrials = data$ntrials,
                          control.predictor = list(compute=TRUE),
                          control.compute = list(waic=TRUE,
                                                 config=TRUE),
                          control.inla = list(cmin = 0),
                          control.fixed = list(prec.intercept=0.01, prec=list(default=0.01)),
                          num.threads = 15,
                          verbose = TRUE
)

summary(mod_icar_intII_rw2)
mod_icar_intII_rw2$waic$waic

# Int-IV and rw1:
icar_formula_intIV_rw1 = result ~ 1 + sex + vac_stat_any + 
  f(ethnicityg,
    model="iid",
    hyper = list(prec = list(prior = "pc.prec", param=c(1,0.1)))
  ) +
  f(ageg_small,
    model="iid",
    hyper= list(prec = list(prior="pc.prec", param=c(1,0.1))) # penalised complexicity (PC) prior for hyperprios
  ) +
  f(cis_num,
    model="bym2",
    graph=Rs,
    scale.model = TRUE,
    constr = TRUE,
    hyper = list(prec = list(prior="pc.prec", param=c(1,0.1)))
  ) +
  f(weeknr, model='rw1', constr=TRUE, scale.model=TRUE, 
    hyper=list(prec=list(prior='pc.prec', param=c(1,0.1)))) +
  f(weeknr.cis_num, model='generic0', Cmatrix=as.matrix(R_1_4_scaled,nrow=dim(R_1_4_scaled)[1]), 
    rankdef=r_def_1_4, constr=TRUE, 
    hyper=list(prec=list(prior='pc.prec', param=c(1,0.1))), 
    extraconstr=list(A=as.matrix(A_constr_1_4,nrow=dim(A_constr_1_4)[1]), e=rep(0,nrow(A_constr_1_4)))) + 
  
  f(weeknr.vac_stat_any, model="iid", 
    hyper = list(prec = list(prior = "pc.prec", param=c(1,0.01)))) +
  f(weeknr.age, model="iid",
    hyper = list(prec = list(prior = "pc.prec", param=c(1,0.01)))) + 
  f(weeknr.ethnicityg, model="iid", 
    hyper = list(prec = list(prior = "pc.prec", param=c(1,0.01))))  


mod_icar_intIV_rw1 = inla(icar_formula_intIV_rw1,
                          family="binomial",
                          data=data, Ntrials = data$ntrials,
                          control.predictor = list(compute=TRUE),
                          control.compute = list(waic=TRUE,
                                                 config=TRUE),
                          control.inla = list(cmin = 0),
                          control.fixed = list(prec.intercept=0.01, prec=list(default=0.01)),
                          num.threads = 15,
                          verbose = TRUE
)

summary(mod_icar_intIV_rw1)
mod_icar_intIV_rw1$waic$waic 


# Int-I and rw1:
icar_formula_intIV_rw2 = result ~ 1 + sex + vac_stat_any + 
  f(ethnicityg,
    model="iid",
    hyper = list(prec = list(prior = "pc.prec", param=c(1,0.1)))
  ) +
  f(ageg_small,
    model="iid",
    hyper= list(prec = list(prior="pc.prec", param=c(1,0.1))) # penalised complexicity (PC) prior for hyperprios
  ) +
  f(cis_num,
    model="bym2",
    graph=Rs,
    scale.model = TRUE,
    constr = TRUE,
    hyper = list(prec = list(prior="pc.prec", param=c(1,0.1)))
  ) +
  f(weeknr, model='rw2', constr=TRUE, scale.model=TRUE, 
    hyper=list(prec=list(prior='pc.prec', param=c(1,0.1)))) +
  f(weeknr.cis_num, model='generic0', Cmatrix=as.matrix(R_2_4_scaled,nrow=dim(R_2_4_scaled)[1]), 
    rankdef=r_def_2_4, constr=TRUE, 
    hyper=list(prec=list(prior='pc.prec', param=c(1,0.1))), 
    extraconstr=list(A=as.matrix(A_constr_2_4,nrow=dim(A_constr_2_4)[1]), e=rep(0,nrow(A_constr_2_4)))) + 
  
  f(weeknr.vac_stat_any, model="iid", 
    hyper = list(prec = list(prior = "pc.prec", param=c(1,0.01)))) +
  f(weeknr.age, model="iid",
    hyper = list(prec = list(prior = "pc.prec", param=c(1,0.01)))) + 
  f(weeknr.ethnicityg, model="iid", 
    hyper = list(prec = list(prior = "pc.prec", param=c(1,0.01))))  


mod_icar_intIV_rw2 = inla(icar_formula_intIV_rw2,
                          family="binomial",
                          data=data, Ntrials = data$ntrials,
                          control.predictor = list(compute=TRUE),
                          control.compute = list(waic=TRUE,
                                                 config=TRUE),
                          control.inla = list(cmin = 0),
                          control.fixed = list(prec.intercept=0.01, prec=list(default=0.01)),
                          num.threads = 15,
                          verbose = TRUE
)

summary(mod_icar_intIV_rw2)
mod_icar_intIV_rw2$waic$waic 

########################################################################################################
#
num_sex <- length(unique(poststrat$sex))
num_sex.weeknr <- length(unique(poststrat$sex))*length(unique(poststrat$weeknr))
num_vac_stat_num <- length(unique(poststrat$vac_stat_num))
num_ethnicityg <- length(unique(poststrat$ethnicityg))
num_ageg_small <- length(unique(poststrat$ageg_small))
num_cis_num <- length(unique(poststrat$cis_num))
num_region_num <- length(unique(poststrat$region_num))
num_region_num.ageg_small <- length(unique(poststrat$region_num))*length(unique(poststrat$ageg_small))
num_region_num.ethnicityg <- length(unique(poststrat$region_num))*length(unique(poststrat$ethnicityg))

num_ageg_small.weeknr <- length(unique(poststrat$ageg_small))*length(unique(poststrat$weeknr))
num_vac_stat_num.weeknr <- length(unique(poststrat$vac_stat_num))*length(unique(poststrat$weeknr))
num_weeknr <- length(unique(poststrat$weeknr))
num_cis_num.weeknr <- length(unique(poststrat$cis_num))*length(unique(poststrat$weeknr))
num_region_num.weeknr <- length(unique(poststrat$region_num))*length(unique(poststrat$weeknr))
num_ethnicityg.weeknr <- length(unique(poststrat$ethnicityg))*length(unique(poststrat$weeknr))

summary(data_grouped$weeknr)

unique(data_grouped$region_num)
#############################################################################################
#####################################################################################

# with int II, best fitting model: 
summary(mod_icar_intIV_rw1)

# take posterior samples:
icar_intIV_samples = inla.posterior.sample(n=6000, result=mod_icar_intIV_rw1,
                                           num.threads = 20)


check <- rownames(icar_intIV_samples[[1]]$latent)
#View(data.frame(check))

#View(check)

# intercept --------------------------------------
icar_intercept_index_icar = grep("Intercept", rownames(icar_intIV_samples[[1]]$latent))

icar_intercept_samples_f = function(x) { #x=samples[[i]]
  return(x$latent[icar_intercept_index_icar])
}

intercept_samples_icar = unlist(lapply(icar_intIV_samples, icar_intercept_samples_f))


# sex --------------------------------------------------
sex_indices_icar = grep("^sex", rownames(icar_intIV_samples[[1]]$latent))
icar_sex_samples_f = function(x) {
  return(x$latent[sex_indices_icar])
}

U_sex_samples_icar = matrix(unlist(lapply(icar_intIV_samples, icar_sex_samples_f)))


# vac_stat_any--------------------------------------------------------
vac_stat_any_indices_icar = grep("^vac_stat_any1:", rownames(icar_intIV_samples[[1]]$latent))

vac_stat_any_samples_f = function(x) {
  return(x$latent[vac_stat_any_indices_icar])
}

U_vac_stat_any_samples_icar = matrix(unlist(lapply(icar_intIV_samples, vac_stat_any_samples_f)))

# ethnicityg ---------------------------------------------------------------
ethnicity_indices_icar = grep("^ethnicityg:", rownames(icar_intIV_samples[[1]]$latent))
icar_ethnicity_samples_f = function(x) {
  return(x$latent[ethnicity_indices_icar])
}

U_ethnicity_samples_icar = matrix(unlist(lapply(icar_intIV_samples, icar_ethnicity_samples_f)),
                                  ncol=num_ethnicityg, 
                                  byrow=TRUE)


# ageg_small --------------------------------------------------------
age_indices_icar = grep("^ageg_small:", rownames(icar_intIV_samples[[1]]$latent))

icar_age_samples_f = function(x) {
  return(x$latent[age_indices_icar])
}

U_age_samples_icar = matrix(unlist(lapply(icar_intIV_samples, icar_age_samples_f)),
                            ncol=num_ageg_small,
                            byrow=TRUE)

# cis area ---------------------------------------------------------
cisnum_indices_icar = grep("^cis_num:", rownames(icar_intIV_samples[[1]]$latent))[(1):(num_cis_num)] # see example from Gelman's group

icar_cisnum_samples_f = function(x) {
  return(x$latent[cisnum_indices_icar])
}

U_cisnum_samples_icar = matrix(unlist(lapply(icar_intIV_samples, icar_cisnum_samples_f)),
                               ncol=num_cis_num,
                               byrow=TRUE)
dim(U_cisnum_samples_icar)



# weeknr.age idd ------------------------------------------------------------
ageg_small.int_indices_icar = grep("weeknr.age:",rownames(icar_intIV_samples[[1]]$latent))

icar_ageg_small.int_samples_f = function(x) {
  return(x$latent[ageg_small.int_indices_icar])
}

length(unique(data_grouped$ageg_small))*length(unique(data_grouped$weeknr))

U_ageg_small.int_samples_icar = matrix(unlist(lapply(icar_intIV_samples, icar_ageg_small.int_samples_f)),
                                       ncol=num_ageg_small.weeknr,
                                       byrow=TRUE)
dim(U_ageg_small.int_samples_icar)



# weeknr.vac_stat_any ---------------------------------------------------------------------------------
vac_stat_any.int_indices_icar = grep("^weeknr.vac_stat_any:",rownames(icar_intIV_samples[[1]]$latent))

icar_vac_stat_any.int_samples_f = function(x) {
  return(x$latent[vac_stat_any.int_indices_icar])
}

length(unique(data_grouped$vac_stat_any))*length(unique(data_grouped$weeknr))

U_vac_stat_any.int_samples_icar = matrix(unlist(lapply(icar_intIV_samples, icar_vac_stat_any.int_samples_f)),
                                         ncol=num_vac_stat_num.weeknr,
                                         byrow=TRUE)
dim(U_vac_stat_any.int_samples_icar)


# weeknr.ethnicityg iid ---------------------------------------------------------------------------------
ethnicityg.int_indices_icar = grep("^weeknr.ethnicityg:",rownames(icar_intIV_samples[[1]]$latent))

icar_ethnicityg.int_samples_f = function(x) {
  return(x$latent[ethnicityg.int_indices_icar])
}

U_ethnicityg.int_samples_icar = matrix(unlist(lapply(icar_intIV_samples, icar_ethnicityg.int_samples_f)),
                                       ncol=num_ethnicityg.weeknr,
                                       byrow=TRUE)
dim(U_ethnicityg.int_samples_icar)



# first order random walk weeknr --------------------------------------------------------
weeknr_indices_icar = grep("^weeknr:",rownames(icar_intIV_samples[[1]]$latent))

icar_weeknr_samples_f = function(x) {
  return(x$latent[weeknr_indices_icar])
}

U_weeknr_samples_icar = matrix(unlist(lapply(icar_intIV_samples, icar_weeknr_samples_f)),
                               ncol=num_weeknr,
                               byrow=TRUE)
dim(U_weeknr_samples_icar)



# cis_num.int idd ------------------------------------------------------------

cisnum.int_indices_icar = grep("^weeknr.cis_num:",rownames(icar_intIV_samples[[1]]$latent))

icar_cisnum.int_samples_f = function(x) {
  return(x$latent[cisnum.int_indices_icar])
}

U_cisnum.int_samples_icar = matrix(unlist(lapply(icar_intIV_samples, icar_cisnum.int_samples_f)),
                                   ncol=num_cis_num.weeknr,
                                   byrow=TRUE)
dim(U_cisnum.int_samples_icar)



#############################################################################################
# matrix below stores the posterior linear predictors for each poststrat cell in poststrat
postpred_sim_icar = matrix(0,
                           length(intercept_samples_icar),
                           dim(poststratI)[1])


dim(postpred_sim_icar)
class(postpred_sim_icar)

library(foreach)
dim(postpred_sim_icar)

parallel::detectCores()
n.cores <- 40
# create the cluster: 
my.cluster <- parallel::makeCluster(
  n.cores,
  type="PSOCK"
)
print(my.cluster)
# register it to be used by %dopar%
doParallel::registerDoParallel(cl=my.cluster)
# check if it is registered
foreach::getDoParRegistered()
foreach::getDoParWorkers()



clusterExport(cl=my.cluster, varlist=c('intercept_samples_icar',
                                       'poststratI',
                                       'U_sex_samples_icar',
                                       'U_vac_stat_any_samples_icar',
                                       'U_ethnicity_samples_icar',
                                       'U_age_samples_icar',
                                       'U_cisnum_samples_icar',
                                       'U_weeknr_samples_icar',
                                       'U_ageg_small.int_samples_icar',
                                       'U_vac_stat_any.int_samples_icar',
                                       'U_ethnicityg.int_samples_icar',
                                       'U_cisnum.int_samples_icar') , envir=environment())


postpred_sim_icar <- parSapply(my.cluster, 1:length(poststratI$weeknr), function(i)
  arm::invlogit(intercept_samples_icar +
                  
                  `if`(unlist(poststratI[i,c("sex")], use.names=FALSE)==2,U_sex_samples_icar[,1],
                       rep(0,6000)) +
                  `if`(unlist(poststratI[i,c("vac_stat_num")], use.names=FALSE)==1,U_vac_stat_any_samples_icar[,1],
                       rep(0,6000)) +
                  
                  U_ethnicity_samples_icar[,unlist(poststratI[i,c("ethnicityg")], use.names=FALSE)] +
                  U_age_samples_icar[,unlist(poststratI[i,c("ageg_small")], use.names=FALSE)] +
                  U_cisnum_samples_icar[,as.numeric(poststratI[i,c("cis_num")])] +
                
                  U_weeknr_samples_icar[,as.numeric(poststratI[i,c("weeknr")])] +
                 
                  
                  U_ageg_small.int_samples_icar[,as.numeric(poststratI[i,c("weeknr.age")])] +
                 
                  U_vac_stat_any.int_samples_icar[,as.numeric(poststratI[i,c("weeknr.vac_stat_num")])] +
                  U_ethnicityg.int_samples_icar[,as.numeric(poststratI[i,c("weeknr.ethnicityg")])] +
                  
                  U_cisnum.int_samples_icar[,as.numeric(poststratI[i,c("weeknr.cis_num")])] 
               
  )
)


dim(postpred_sim_icar)



################################################################################################
# post-stratify for each CIS area:
poststratI <- poststratI%>%arrange(weeknr.cis_num, sex,
                                   weeknr.age, weeknr.ethnicityg,
                                   weeknr.vac_stat_any
)
cisnum_time_mrp_samples_icar = matrix(0, 6000, num_cis_num.weeknr) 

# sort levels so I know what relates to what:
levels_weeknr.cis_num = sort(unique(poststratI$weeknr.cis_num))


N_sub = rep(0, num_cis_num) # the N for each time-region based off poststrat matrix
g_counter = 1
for (g in levels_weeknr.cis_num) {
  N_sub[g_counter] = sum(poststratI[which(poststratI$weeknr.cis_num==g),]$N)
  g_counter = g_counter + 1
}

g_counter = 1
for (g in levels_weeknr.cis_num) {
  # print(g)
  for (i in which(poststratI$weeknr.cis_num==g)) {
    cisnum_time_mrp_samples_icar[,g_counter] = cisnum_time_mrp_samples_icar[, g_counter] + 
      (postpred_sim_icar[,i] * poststratI$N[i] / N_sub[g_counter])
  }
  g_counter = g_counter + 1
}

names(poststratI)

cistime_results <- poststratI%>%group_by(cis_num, CIS.code, Region_Name, region_num, weeknr, weeknr.cis_num)%>%
  summarise(n=n())%>%dplyr::select(-n)%>%arrange(weeknr.cis_num)

dim(cisnum_time_mrp_samples_icar)
# so there are 8584 combinations, each 116 is a week: 
check <- cistime_results%>%filter(weeknr==1)

check2 <- data.frame(apply(cisnum_time_mrp_samples_icar[,1:116],
                           MARGIN=1,
                           FUN=rank,
))
check$rank_median <- apply(check2,
                           MARGIN=1,
                           FUN=quantile,
                           probs=0.5)
check$rank_liqr <- apply(check2,
                         MARGIN=1,
                         FUN=quantile,
                         probs=0.25)
check$rank_uiqr <- apply(check2,
                         MARGIN=1,
                         FUN=quantile,
                         probs=0.75)
check$rank_ll95 <- apply(check2,
                         MARGIN=1,
                         FUN=quantile,
                         probs=0.025)
check$rank_ul95 <- apply(check2,
                         MARGIN=1,
                         FUN=quantile,
                         probs=0.975)

# do a for loop, for each week. 

check <- cistime_results
check$rank_median <- NA
check$rank_liqr <- NA
check$rank_uiqr <- NA
check$rank_ll95 <- NA
check$rank_ul95 <- NA
check$prob_top10 <- NA
check$prob_bottom10 <- NA

for (i in min(cistime_results$weeknr):max(cistime_results$weeknr)) {
  check2 <- data.frame(apply(cisnum_time_mrp_samples_icar[,(1+(116*(i-1))):(116*i)],
                             MARGIN=1,
                             FUN=rank,
  ))
  check$rank_median[(1+(116*(i-1))):(116*i)] <- apply(check2,
                                                      MARGIN=1,
                                                      FUN=quantile,
                                                      probs=0.5)
  check$rank_liqr[(1+(116*(i-1))):(116*i)] <- apply(check2,
                                                    MARGIN=1,
                                                    FUN=quantile,
                                                    probs=0.25)
  check$rank_uiqr[(1+(116*(i-1))):(116*i)] <- apply(check2,
                                                    MARGIN=1,
                                                    FUN=quantile,
                                                    probs=0.75)
  check$rank_ll95[(1+(116*(i-1))):(116*i)] <- apply(check2,
                                                    MARGIN=1,
                                                    FUN=quantile,
                                                    probs=0.025)
  check$rank_ul95[(1+(116*(i-1))):(116*i)] <- apply(check2,
                                                    MARGIN=1,
                                                    FUN=quantile,
                                                    probs=0.975)
  check$prob_top10[(1+(116*(i-1))):(116*i)] <- apply(check2,
                                                     MARGIN=1,
                                                     function(x) round(sum(x<=10)/6000,2))
  check$prob_bottom10[(1+(116*(i-1))):(116*i)] <- apply(check2,
                                                        MARGIN=1,
                                                        function(x) round(sum(x>=107)/6000,2))
}



# get start of week and end of week: 
check <- merge(temp, check, 
               by="weeknr")

write.csv(check, 
          paste0("weeklyranking_cismodel_swabs_cis_time_INLA_intIV_rw1",".csv"))


###########################
# do something with average value instead of ranking: 
cistime_results <- poststratI%>%group_by(cis_num, CIS.code, Region_Name, region_num, weeknr, weeknr.cis_num)%>%
  summarise(n=n())%>%dplyr::select(-n)%>%arrange(weeknr.cis_num)

cistime_results$median_prob <- apply(cisnum_time_mrp_samples_icar,
                                     MARGIN=2,
                                     FUN=quantile,
                                     probs=0.5)*100
cistime_results$ll <- apply(cisnum_time_mrp_samples_icar,
                            MARGIN=2,
                            FUN=quantile,
                            probs=0.025)*100
cistime_results$ul <- apply(cisnum_time_mrp_samples_icar,
                            MARGIN=2,
                            FUN=quantile,
                            probs=0.975)*100
cistime_results$sd <- apply(cisnum_time_mrp_samples_icar,
                            MARGIN=2,
                            FUN=sd)*100
cistime_results <- cistime_results%>%arrange(cis_num, weeknr)


cistime_results <- merge(temp, cistime_results, 
                         by="weeknr")


write.csv(cistime_results, 
          paste0("positivityestimates_cismodel_swabs_cis_time_INLA_intIV_rw1",format(min(as.Date(cistime_results$start_date,"%Y-%m-%d")),"%m%d%Y"),
                 "_",format(max(as.Date(cistime_results$start_date,"%Y-%m-%d")),"%m%d%Y"),"_",
                 run_until,".csv"))

