## continuous-binary model - compile models
rm(list=ls())

# set working directory
wdir <- file.path("~/bayes_cop_pow")

## make sim array
set.seed(82606)
nreps <- 100

# placebo group same for all
mu_1 <- -150
sigma2_1 <- 100^2
p_1 <- 0.1  
rho_1 <- 0.1

ns <- c(50, 100, 200, 400) # sample size (n per arm) 
rho_2s <- c(0.1, 0.3, 0.5) # correlation in treatment group
mu_2s <-mu_1 + c(0, 100, 150) # effectiveness in trt
sigma2_2s <- 100^2
p_2s <-p_1 + c(0, 0.3, 0.6) # prob AE in trt

a1 <- expand.grid(ns,rho_2s,mu_2s,sigma2_2s,p_2s)
names(a1) <- c("n","rho_2","mu_2","sigma2_2","p_2")

# replicate each scenario nreps times
a2 <- do.call("rbind", replicate(nreps, a1, simplify = FALSE))

# reorder and add rep id
a3<-a2[order(a2$n,a2$rho_2,a2$mu_2, a2$sigma2_2, a2$p_2),]
a3$rep_id<-rep(1:nreps,nrow(a1))

# randomize order
cb_sim_array <- a3[sample(1:nrow(a3)),]

# add seeds and simulation id
cb_sim_array$dat_seed <- sample.int(1e7,size=nrow(cb_sim_array))
cb_sim_array$samp_seed <- sample.int(1e7,size=nrow(cb_sim_array))
cb_sim_array$sim_id <- 1:nrow(cb_sim_array)

## temp
getTheta <- function(rho,p){  (rho*sqrt(p*(1-p))) / dnorm(qnorm(p)) }

# treatment
mu_2 <- sim_params$mu_2
sigma2_2 <- sim_params$sigma2_2
p_2 <- sim_params$p_2
rho_2 <- sim_params$rho_2

foo<-cb_sim_array[,c("rho_2","p_2")]

thets<-mapply(getTheta,foo[,1],foo[,2])

max(thets)

# save simulation array
saveRDS(cb_sim_array, file = file.path(wdir,"cb_sim_array.rds"))
