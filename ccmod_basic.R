## continuous-continuous model - simulate data and sample from posterior 

libs<-c("copula", "magrittr", "rstan", "sfsmisc")
invisible(lapply(libs, library, character.only = TRUE))

#sessionInfo()
#head(Sys.cpuinfo(),19)

# set working directory
wdir <- file.path("~/bayes_cop_pow")

# for basic version don't call sim array with multiple scenarios, just choose one scenario
if (0){
# load sim array
cc_sim_array <- readRDS(file.path(wdir,"cc_sim_array.rds"))

#! remove for accre
#cc_sim_array <- readRDS(file.path("/home/nathan/Dropbox/njames/school/PhD/misc/conferences/ENAR2019/code/cc_sim_array.rds"))


#select parameters for given sim id using SLURM_ARRAY_TASK_ID
arg <- commandArgs(trailing=TRUE)
s_id <- as.integer(arg[1])

#! remove for accre
#s_id <- 1 

sim_params <- bb_sim_array[bb_sim_array$sim_id == s_id,]
}

## Simulate data
#set.seed(sim_params$dat_seed)
set.seed(1384)


# number of samples per arm
#n <- sim_params$n
n <- 500


# placebo group
mu_e1 <- -150 # efficacy
sigma2_e1 <- 100^2
mu_s1 <- 0 # safety
sigma2_s1 <- 100^2
rho_1 <- 0.1 # correlation

# normal copula
nc_p <- normalCopula(rho_1)

pbo_dist <- mvdc(nc_p, margins = c("norm","norm"),
                 paramMargins = list(list(mean = mu_e1, sd = sqrt(sigma2_e1)), 
                                     list(mean = mu_s1, sd = sqrt(sigma2_s1))) )

pbo_samps <- rMvdc(n, pbo_dist)

# cor(pbo_samps) #check correlation

# treatment group
mu_e2 <- -150+50 # efficacy
sigma2_e2 <- 100^2
mu_s2 <- 0+50 # safety
sigma2_s2 <- 100^2
rho_2 <- 0.1 # correlation !! update to be diff (more correlated for treatment)

# normal copula
nc_t <- normalCopula(rho_2)
trt_dist <- mvdc(nc_t, margins = c("norm","norm"),
                    paramMargins = list(list(mean = mu_e2, sd = sqrt(sigma2_e2)), 
                                        list(mean = mu_s2, sd = sqrt(sigma2_s2))) )

trt_samps <- rMvdc(n, trt_dist)

#combine placebo and treatment data
dat_cc <- rbind(pbo_samps,trt_samps) %>% cbind(sort(rep(c(0,1),n)),
                                               sort(rep(c(0,1),n),decreasing=TRUE),
                                               sort(rep(c(0,1),n))) %>% as.data.frame() 

names(dat_cc) <- c("efficacy","safety","treatment","trt1","trt2")

#dat_cc %>% ggplot(aes(y=efficacy,group=treatment))+geom_boxplot()


## Frequentist MLE

#logLik( lm(efficacy~trt1+trt2-1,data=dat_cc) )
#lm(safety~trt1+trt2-1,data=dat_cc)


Xmat<-model.matrix(~treatment, data=dat_cc)

nmLL <- function(beta.m, y, x, pobs=FALSE){
  
  np <- ncol(x)
  mu_z <- x %*% beta.m[1:np]
  
  # least squares 
  SSE <- t(y) %*% y - t(beta.m) %*% t(x) %*% y
  s2_hat <- SSE/(nrow(x)-np-1)
 
  # try using MLE
  #s2_hat <- (t(y-x%*%beta.m) %*% (y-x%*%beta.m)) / nrow(x)
  
  nLL <- -sum(dnorm(y, mean=mu_z, sd=sqrt(s2_hat), log=TRUE))

  if (!pobs) nLL else
    list(nLL = nLL, U = pnorm(y, mean=mu_z, sd=sqrt(s2_hat)) ) 
}

#nmLL(c(-150,100),dat_cc[,'efficacy'],Xmat)
logLik(lm(efficacy~treatment,data=dat_cc))
summary(lm(efficacy~treatment,data=dat_cc))
optim(c(-150,50),nmLL,y=dat_cc[,'efficacy'],x=Xmat, method="BFGS", hessian=TRUE)

#nmLL(c(0,100),dat_cc[,'safety'],Xmat)
logLik(lm(safety~treatment,data=dat_cc))
summary(lm(safety~treatment,data=dat_cc))
optim(c(0,50),nmLL,y=dat_cc[,'safety'],x=Xmat, method="BFGS")


nfLL <- function(par, y, x, copula){
  
  copnam <- substr(substitute(copula),1,11)[1]
  
  beta <- par[1] # copula param !! need to make this depend on params
  
  # try to set parameters for non-Independence copula
  if (copnam != 'indepCopula'){
    tc <-tryCatch(copula <- setTheta(copula, beta), 
                error = function(e) NULL)
    if(is.null(tc)) return(-Inf) # in case of failure, return -Inf
  }
  
  p <- ncol(x) # number of params per margin
  beta.1 <- par[1 + 1:p]
  beta.2 <- par[p+1 + 1:p]
  
  # marginal log-likelihood
  nmLL.1 <- nmLL(beta.1, y[,'efficacy'], x, pobs=TRUE)
  nmLL.2 <- nmLL(beta.2, y[,'safety'], x, pobs=TRUE)
  
  ## In case of invalid evaluation of the likelihoods, return -Inf
  if(any(is.na(c(nmLL.1$nLL, nmLL.2$nLL)))) return(-Inf)
  
  # parametric pseudo-obs
  U <- cbind(nmLL.1$U, nmLL.2$U)
  
  -sum(dCopula(u=U, copula=copula, log=TRUE)) + nmLL.1$nLL + nmLL.2$nLL
}

#nfLL(c(0.1,-150,50,0,50), dat_cc, Xmat, copula=normalCopula(dim=2))
#nfLL(c(0,-150,50,0,50), dat_cc, Xmat,copula=indepCopula(dim=2))


res.n <- optim(c(0.1,-150,50,0,50),nfLL,y=dat_cc,x=Xmat, copula=normalCopula(dim=2), 
               method="Nelder-Mead", hessian=TRUE)

cov_n_fit <-solve(res.n$hessian)
se_n_fit <-sqrt(diag(cov_n_fit))

# normal copula model
cbind(res.n$par, se_n_fit)

res.i <- optim(c(0,-150,50,0,50),nfLL,y=dat_cc,x=Xmat, copula=indepCopula(dim=2), 
               method="BFGS", hessian=TRUE, control=list(reltol=1e-16))
cov_i_fit <-solve(res.i$hessian[2:5,2:5])
se_i_fit <-sqrt(diag(cov_i_fit))

#independence copula (should match lm())
cbind(res.i$par[2:5], se_i_fit)

summary(lm(efficacy~treatment,data=dat_cc))
summary(lm(safety~treatment,data=dat_cc))


# CI too narrow --> increase type I error, increase power



# Bayesian
## run Stan models

# format data for stan
mod_data_cc <- list(N=nrow(dat_cc), 
                    x=dat_cc[,c("trt1","trt2")], 
                    y_e=dat_cc$efficacy, 
                    y_s=dat_cc$safety)


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

