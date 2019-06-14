## continuous-continuous model - compile models
rm(list=ls())

# set working directory
wdir <- file.path("~/bayes_cop_pow")

## make sim array
set.seed(82606)
nreps <- 100

# placebo group same for all
mu_e_1 <- -150
sigma2_e_1 <- 100^2
mu_s_1 <- 0
sigma2_s_1 <- 100^2
rho_1 <- 0.1

ns <- c(50, 100, 200, 400) # sample size (n per arm) 
rho_2s <- c(0.1, 0.3, 0.5) # correlation in treatment group
mu_e_2s <-mu_e_1 + c(0, 100, 150) # effectiveness in trt
sigma2_e_2s <- 100^2

mu_s_2s <-mu_s_1 + c(0, 100, 150) # safety in trt
sigma2_s_2s <- 100^2

a1 <- expand.grid(ns,rho_2s,mu_e_2s,sigma2_e_2s,mu_s_2s,sigma2_s_2s)
names(a1) <- c("n","rho_2","mu_e_2","sigma2_e_2","mu_s_2","sigma2_s_2")


# replicate each scenario nreps times
a2 <- do.call("rbind", replicate(nreps, a1, simplify = FALSE))

# reorder and add rep id
a3<-a2[order(a2$n,a2$rho_2,a2$mu_e_2, a2$sigma2_e_2, a2$mu_s_2, a2$sigma2_s_2),]
a3$rep_id<-rep(1:nreps,nrow(a1))

# randomize order
cc_sim_array <- a3[sample(1:nrow(a3)),]

# add seeds and simulation id
cc_sim_array$dat_seed <- sample.int(1e7,size=nrow(cc_sim_array))
cc_sim_array$samp_seed <- sample.int(1e7,size=nrow(cc_sim_array))
cc_sim_array$sim_id <- 1:nrow(cc_sim_array)

# save simulation array
saveRDS(cc_sim_array, file = file.path(wdir,"cc_sim_array.rds"))
