## binary-binary model - simulate data and sample from posterior 

libs<-c("copula", "magrittr", "rstan", "sfsmisc")
invisible(lapply(libs, library, character.only = TRUE))

sessionInfo()
head(Sys.cpuinfo(),19)

# set working directory
wdir <- file.path("~/bayes_cop_pow")

# load sim array
bb_sim_array <- readRDS(file.path(wdir,"bb_sim_array.rds"))

#! remove for accre
#bb_sim_array <- readRDS(file.path("/home/nathan/Dropbox/njames/school/PhD/misc/conferences/ENAR2019/code/bb_sim_array.rds"))


#select parameters for given sim id using SLURM_ARRAY_TASK_ID
arg <- commandArgs(trailing=TRUE)
s_id <- as.integer(arg[1])

#! remove for accre
#s_id <- 1 

sim_params <- bb_sim_array[bb_sim_array$sim_id == s_id,]

## Simulate data
set.seed(sim_params$dat_seed)

# number of samples per arm
n <- sim_params$n

# placebo group
p_e1 <- 0.2 # prob effective
p_s1 <- 0.1 # prob AE
rho_1 <- 0.1 # tetrachoric correlation, can estimate with polycor::polychor

# normal copula
nc_p <- normalCopula(rho_1)
pbo_dist <- mvdc(nc_p, margins = c("binom", "binom"),
                 paramMargins = list(list(size = 1, prob = p_e1), 
                                     list(size = 1, prob = p_s1)) )

pbo_samps <- rMvdc(n, pbo_dist)

# treatment group
p_e2 <- sim_params$p_e2 # prob effective
p_s2 <- sim_params$p_s2 # prob AE
rho_2 <- sim_params$rho_2 # tetrachoric corr

# normal copula
nc_t <- normalCopula(rho_2)
trt_dist <- mvdc(nc_t, margins = c("binom", "binom"),
                 paramMargins = list(list(size = 1, prob = p_e2), 
                                     list(size = 1, prob = p_s2)) )

trt_samps <- rMvdc(n, trt_dist)

#combine placebo and treatment data
dat_bb <- rbind(pbo_samps,trt_samps) %>% cbind(sort(rep(c(0,1),n)),
                                               sort(rep(c(0,1),n),decreasing=TRUE),
                                               sort(rep(c(0,1),n))) %>% as.data.frame() 
names(dat_bb) <- c("efficacy","safety","treatment","trt1","trt2")

dat_bb_short <- plyr::count(dat_bb,vars=names(dat_bb))

# format data for stan
mod_data_bb <- list(N=nrow(dat_bb), 
                    x=dat_bb[,c("trt1","trt2")], 
                    y_e=dat_bb$efficacy, 
                    y_s=dat_bb$safety)

mod_data_bb_short <- list(N=nrow(dat_bb_short),
                          x=dat_bb_short[,c("trt1","trt2")],
                          y_e=dat_bb_short$efficacy,
                          y_s=dat_bb_short$safety,
                          freq=dat_bb_short$freq)

## run Stan models

# MCMC parameters
options(mc.cores = parallel::detectCores())
n_chains <- 2
n_warmup <- 3000
n_iter <- n_warmup+2500
samp_seed <- sim_params$samp_seed

# Load pre-compiled models from bbmod_init.R
mod_bb_e <- readRDS(file.path(wdir,"mod_bb_e.rds"))
mod_bb_s <- readRDS(file.path(wdir,"mod_bb_s.rds"))
mod_bb_jnt <- readRDS(file.path(wdir,"mod_bb_jnt.rds"))

#mod_bb_e <- readRDS(file.path("/home/nathan/Dropbox/njames/school/PhD/misc/conferences/ENAR2019/code/mod_bb_e.rds"))
#mod_bb_s <- readRDS(file.path("/home/nathan/Dropbox/njames/school/PhD/misc/conferences/ENAR2019/code/mod_bb_s.rds"))
#mod_bb_jnt <- readRDS(file.path("/home/nathan/Dropbox/njames/school/PhD/misc/conferences/ENAR2019/code/mod_bb_jnt.rds"))


## Sample from compiled models

# efficacy marginal model
fit_bb_e <- sampling(mod_bb_e, data=mod_data_bb, seed=samp_seed, 
                 iter=n_iter, warmup=n_warmup, chains=n_chains)

assign(paste0("summ_bb_e_",s_id), summary(fit_bb_e))
assign(paste0("samp_bb_e_",s_id), as.matrix(fit_bb_e, pars=c("p_e")))

# safety marginal model
fit_bb_s <- sampling(mod_bb_s, data=mod_data_bb, seed=samp_seed,
                 iter=n_iter, warmup=n_warmup, chains=n_chains)

assign(paste0("summ_bb_s_",s_id), summary(fit_bb_s))
assign(paste0("samp_bb_s_",s_id), as.matrix(fit_bb_s, pars=c("p_s")))

## joint copula model
# efficacy marginal model MLE for initialization
mle_bb_e <- glm(efficacy~trt1+trt2-1, data=dat_bb, family=binomial(link="probit"))

# safety marginal model MLE for initialization
mle_bb_s <- glm(safety~trt1+trt2-1, data=dat_bb, family=binomial(link="probit"))

#initalize margins at jittered MLE estimate??
init_list0 <- list(beta_e=mle_bb_e$coefficients, beta_s=mle_bb_s$coefficients)
init_list <- lapply(1:n_chains, function(x) lapply(init_list0 ,jitter, amount=1.5))

fit_bb_jnt <- sampling(mod_bb_jnt, data=mod_data_bb_short, seed=samp_seed,
                   iter=n_iter, chains=n_chains, warmup=n_warmup,
                   init=init_list, control = list(adapt_delta = 0.95))

assign(paste0("summ_bb_jnt_",s_id), summary(fit_bb_jnt))
assign(paste0("samp_bb_jnt_",s_id), as.matrix(fit_bb_jnt, pars=c("omega","p_e","p_s")))


if (0){
# store results in scratch
sdir <- file.path("/gpfs23/scratch/jamesnt")

# keep 3 fits, dataset (dat_bb), sim parameters
save(fit_bb_e, fit_bb_s, fit_bb_jnt, dat_bb, sim_params, 
     file=file.path(sdir,"bbsims", paste0("bb_sim_",s_id,".RData")))
}


# keep summaries and samples, dataset (dat_bb), sim parameters
assign(paste0("dat_bb_",s_id),dat_bb)
assign(paste0("sim_params_",s_id),sim_params)

savelist <- c("summ_bb_e", "samp_bb_e", "summ_bb_s", "samp_bb_s", 
"summ_bb_jnt", "samp_bb_jnt", "dat_bb", "sim_params")

save(list=paste0(savelist,"_",s_id), 
     file=file.path(wdir,"bbsims", paste0("bb_sim_",s_id,".RData")))

