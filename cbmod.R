## continuous-binary model - simulate data and sample from posterior 

libs<-c("copula", "magrittr", "rstan", "sfsmisc")
invisible(lapply(libs, library, character.only = TRUE))

sessionInfo()
head(Sys.cpuinfo(),19)

# set working directory
wdir <- file.path("~/bayes_cop_pow")

# load sim array
cb_sim_array <- readRDS(file.path(wdir,"cb_sim_array.rds"))

#select parameters for given sim id using SLURM_ARRAY_TASK_ID
arg <- commandArgs(trailing=TRUE)
s_id <- as.integer(arg[1])

sim_params <- cb_sim_array[cb_sim_array$sim_id == s_id,]


# Simulate data
set.seed(sim_params$dat_seed)

# function to get copula parameter given rho and p; see Costa section 3.1.2
getTheta <- function(rho,p){  (rho*sqrt(p*(1-p))) / dnorm(qnorm(p)) }

# number of samples per arm
n <- sim_params$n

# placebo group
mu_1 <- -150
sigma2_1 <- 100^2
p_1 <- 0.1  
rho_1 <- 0.1

# normal copula
nc_p <- normalCopula( getTheta(rho=rho_1, p=p_1)  )

pbo_dist <- mvdc(nc_p, margins = c("norm","binom"),
                 paramMargins = list(list(mean = mu_1, sd = sqrt(sigma2_1)), 
                                     list(size = 1, prob = p_1)) )

pbo_samps<-rMvdc(n, pbo_dist)

# treatment
mu_2 <- sim_params$mu_2
sigma2_2 <- sim_params$sigma2_2
p_2 <- sim_params$p_2
rho_2 <- sim_params$rho_2

# normal copula
nc_t <- normalCopula( getTheta(rho=rho_2, p=p_2)  )

trt_dist <- mvdc(nc_t, margins = c("norm","binom"),
                 paramMargins = list(list(mean = mu_2, sd = sqrt(sigma2_2)), 
                                     list(size = 1, prob = p_2)) )

trt_samps <- rMvdc(n, trt_dist)

#combine placebo and treatment data
dat_cb <- rbind(pbo_samps,trt_samps) %>% cbind(sort(rep(c(0,1),n)),
                                            sort(rep(c(0,1),n),decreasing=TRUE),
                                            sort(rep(c(0,1),n))) %>% as.data.frame() 
names(dat_cb) <- c("efficacy","safety","treatment","trt1","trt2")

# format data for stan
mod_data_cb <- list(N=nrow(dat_cb),
                    x=dat_cb[,c("trt1","trt2")],
                    y1=dat_cb$efficacy,
                    y2=dat_cb$safety)

## run Stan models

# MCMC parameters
options(mc.cores = parallel::detectCores())
n_chains <- 2
n_warmup <- 3000
n_iter <- n_warmup+2500
samp_seed <- sim_params$samp_seed

# Load pre-compiled models from cbmod_init.R
mod_cb_e <- readRDS(file.path(wdir,"mod_cb_e.rds"))
mod_cb_s <- readRDS(file.path(wdir,"mod_cb_s.rds"))
mod_cb_jnt <- readRDS(file.path(wdir,"mod_cb_jnt.rds"))


# Sample from compiled models

# efficacy marginal model
fit_cb_e <- sampling(mod_cb_e, data=mod_data_cb, seed=samp_seed, 
                 iter=n_iter, warmup=n_warmup, chains=n_chains)

assign(paste0("summ_cb_e_",s_id), summary(fit_cb_e))
assign(paste0("samp_cb_e_",s_id), as.matrix(fit_cb_e, pars=c("beta_e","sigma")))


# safety marginal model
fit_cb_s <- sampling(mod_cb_s, data=mod_data_cb, seed=samp_seed,
                 iter=n_iter, warmup=n_warmup, chains=n_chains)

assign(paste0("summ_cb_s_",s_id), summary(fit_cb_s))
assign(paste0("samp_cb_s_",s_id), as.matrix(fit_cb_s, pars=c("p")))

## joint copula model
# efficacy marginal model MLE for initialization
mle_cb_e<-summary(lm(efficacy~trt1+trt2-1,data=dat_cb))

# safety marginal model MLE for initialization
mle_cb_s<-glm(safety~trt1+trt2-1,data=dat_cb,family=binomial(link="probit"))

#initalize margins at jittered MLE estimate??
init_list0 <- list(beta_e=mle_cb_e$coefficients[,1], beta_s=mle_cb_s$coefficients,
                   s=rep(mle_cb_e$sigma,2))
init_list <- lapply(1:n_chains, function(x) lapply(init_list0 ,jitter))

fit_cb_jnt <- sampling(mod_cb_jnt, data=mod_data_cb, seed=samp_seed,
                   iter=n_iter, chains=n_chains, warmup=n_warmup,
                   init=init_list, control = list(adapt_delta = 0.95))

assign(paste0("summ_cb_jnt_",s_id), summary(fit_cb_jnt))
assign(paste0("samp_cb_jnt_",s_id), as.matrix(fit_cb_jnt, pars=c("mu","p","theta","rho")))

# keep summaries and samples, dataset (dat_cb), sim parameters
assign(paste0("dat_cb_",s_id),dat_cb)
assign(paste0("sim_params_",s_id),sim_params)

savelist <- c("summ_cb_e", "samp_cb_e", "summ_cb_s", "samp_cb_s", 
              "summ_cb_jnt", "samp_cb_jnt", "dat_cb", "sim_params")

save(list=paste0(savelist,"_",s_id), 
     file=file.path(wdir,"cbsims", paste0("cb_sim_",s_id,".RData")))

if (0){
  # store results in scratch
  sdir <- file.path("/gpfs23/scratch/jamesnt")
  
  save(list=paste0(savelist,"_",s_id), 
       file=file.path(sdir,"cbsims", paste0("cb_sim_",s_id,".RData")))
}
