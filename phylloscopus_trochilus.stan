data {
  int<lower=0> N;   
  int<lower=0> p;
  int<lower=0> J;

  int<lower=0> y[N];    
  
  int<lower=0> intercept[N];
  int<lower=0> hab_Co[N];
  int<lower=0> hab_Op[N];
  int<lower=0> hab_Urb[N];
  int<lower=0> hab_We[N];
  real apr_may[N];
  int<lower=0> year[N];
  
  int<lower=0, upper=J> route[N];
  
  real<lower=0> rec_psi;
}

parameters {
  real beta[p];
  real delta[J];
  real<lower=0> sigma;
}

transformed parameters {
  real lp[N];
  real<lower=0> mu[N];
  
  for (i in 1:N) {
    lp[i] = beta[1]*intercept[i] + beta[2]*hab_Co[i] + beta[3]*hab_Op[i]
            + beta[4]*hab_Urb[i] + beta[5]*hab_We[i] + beta[6]*apr_may[i]
            + beta[7]*year[i] + delta[route[i]];
    mu[i] = exp(lp[i]);   
  }
}

model {
  delta ~ normal(0, sigma);
  beta ~ normal(0, 10);
  target += rec_psi*poisson_lpmf(y | mu);
}
