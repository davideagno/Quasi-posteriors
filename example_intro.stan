data {
  int<lower=0> N;
  int<lower=0> p;
  
  int<lower=0> intercept[N];
  real slope[N];
  
  real y[N];
  
  real<lower=0> rec_psi;
}

parameters {
  real beta[p];
}

transformed parameters {
  real mu[N];
  
  for (i in 1:N) {
    mu[i] = beta[1]*intercept[i] + beta[2]*slope[i];
  }
}

model {
  for (i in 1:N) {
    target += rec_psi * (mu[i] - y[i] + 1) * exp(-mu[i]);
  }
}
