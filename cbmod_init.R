## continuous-binary model - compile models
rm(list=ls())
 
library(rstan)

# set working directory
wdir <- file.path("~/bayes_cop_pow")

## Define models

# marginal efficacy model
mod_cb_e_code <- "
data {
int N;
matrix[N, 2] x;
vector[N] y1;
}
parameters {
// params for continuous (efficacy) outcome
vector[2] beta_e;
real<lower=0> sigma;  
}
model {
vector[N] mu;

// priors
beta_e ~ normal(0,1000);
sigma ~ inv_gamma(0.001,0.001); 

// marginal for continuous (efficacy) outcome
mu = x*beta_e;
y1 ~ normal(mu, sigma);
}
"

# marginal safety model
mod_cb_s_code <- "
data {
int N;
matrix[N, 2] x;
int<lower=0, upper=1> y2[N];
}
parameters {
//params for binary (safety) outcome
vector[2] beta_s;
}
model {
vector[N] p;

// priors
beta_s ~ normal(0,1000); 

// marginal for binary (safety) outcome
y2 ~ bernoulli(Phi(x*beta_s)); 
}
generated quantities {
vector[2] p;
p = Phi(beta_s);
}
"

# joint model
mod_cb_jnt_code <- "
functions {
real binorm_cop_lp(real y1, real mu, real sigma, real y2, real p, real theta) {
real targ = normal_lpdf(y1|mu,sigma);
if (y2==0) {
targ = targ + normal_lcdf((inv_Phi(1-p)-theta*
inv_Phi(normal_cdf(y1,mu,sigma)))/sqrt(1-theta^2)|0,1);
} else {
targ = targ + normal_lccdf((inv_Phi(1-p)-theta*
inv_Phi(normal_cdf(y1,mu,sigma)))/sqrt(1-theta^2)|0,1);
}
return targ;
}
}
data {
int N;
matrix[N, 2] x;
vector[N] y1;
int<lower=0, upper=1> y2[N];
}
parameters {
// params for continuous (efficacy) outcome
vector[2] beta_e;
vector<lower=0>[2] s;  

//params for binary (safety) outcome
vector[2] beta_s;

// copula dependence param
vector<lower=-1, upper=1>[2] omega;  
}
model {
vector[N] mu;
vector[N] sigma;
vector[N] p;
vector[N] theta;

// priors
beta_e ~ normal(0,1000);
beta_s ~ normal(0,1000); 
s ~ inv_gamma(0.001,0.001); 

// marginal for continuous (efficacy) outcome
mu = x*beta_e;
sigma = x*s;

// marginal for binary (safety) outcome
p = Phi(x*beta_s);

// copula dependence parameter
theta = x*omega;

// build log-likelihood
{
  vector[N] loglik;  // vectorize summation
  for (i in 1:N)
    loglik[i] = binorm_cop_lp(y1[i],mu[i],sigma[i], y2[i], p[i], theta[i]);;
  target += sum(loglik);
}


}
generated quantities {
vector[2] mu;
vector[2] p;
vector[2] theta;
vector[2] rho;

mu = beta_e;
p = Phi(beta_s);
theta = omega;

rho[1] = theta[1]*exp(normal_lpdf(inv_Phi(p[1])|0,1))/sqrt(p[1]*(1-p[1]));
rho[2] = theta[2]*exp(normal_lpdf(inv_Phi(p[2])|0,1))/sqrt(p[2]*(1-p[2]));
}
"

## Compile and save models

# efficacy marginal model
mod_cb_e <- stan_model(model_code = mod_cb_e_code)
saveRDS(mod_cb_e, file = file.path(wdir,"mod_cb_e.rds"))

# safety marginal model
mod_cb_s <- stan_model(model_code = mod_cb_s_code)
saveRDS(mod_cb_s, file = file.path(wdir,"mod_cb_s.rds"))

# joint model
mod_cb_jnt <- stan_model(model_code = mod_cb_jnt_code)
saveRDS(mod_cb_jnt, file = file.path(wdir,"mod_cb_jnt.rds"))
