## binary-binary model - sim array and compile models
rm(list=ls())
 
# set working directory
wdir <- file.path("~/bayes_cop_pow")

## make sim array
set.seed(112884)
nreps <- 100

# same for all
p_e1 <- 0.2 # prob effective
p_s1 <- 0.1 # prob AE
rho_1 <- 0.1 # tetrachoric corr

ns <- c(50, 100, 200, 400) # sample size (n per arm) 
rho_2s <- c(0.1, 0.35, 0.6) # tetrachoric correlation in treatment group
p_e2s <-p_e1 + c(0, 0.3, 0.6) # prob effective in trt
p_s2s <-p_s1 + c(0, 0.3, 0.6) # prob AE in trt

# equivalent increase in eff and safety
#p_e2s <- c(0, 0.3, 0.6)
#p_s2s <- c(0, 0.3, 0.6)

#expand.grid(rho_2s,p_e2s,p_s2s)
a1 <- expand.grid(ns,rho_2s,p_e2s,p_s2s)
names(a1) <- c("n","rho_2","p_e2","p_s2")

# replicate each scenario nreps times
a2 <- do.call("rbind", replicate(nreps, a1, simplify = FALSE))

# reorder and add rep id
a3<-a2[order(a2$n,a2$rho_2,a2$p_e2,a2$p_s2),]
a3$rep_id<-rep(1:nreps,nrow(a1))

# randomize order
bb_sim_array <- a3[sample(1:nrow(a3)),]

# add seeds and simulation id
bb_sim_array$dat_seed <- sample.int(1e7,size=nrow(bb_sim_array))
bb_sim_array$samp_seed <- sample.int(1e7,size=nrow(bb_sim_array))
bb_sim_array$sim_id <- 1:nrow(bb_sim_array)

# save simulation array
saveRDS(bb_sim_array, file = file.path(wdir,"bb_sim_array.rds"))
