data {
  int<lower=0> N;
  int<lower=0> p;
  
  int<lower=0> intercept[N];
  real v1[N];
  real v2[N];
  real v3[N];
  
  real y[N];
  
  real<lower=0> rec_psi;
}

parameters {
  real beta[p];
}

transformed parameters {
  real mu[N];
  for (i in 1:N) {
    mu[i] = beta[1]*intercept[i] + beta[2]*v1[i] + beta[3]*v2[i] + beta[4]*v3[i];
  }
}

model {
  for (i in 1:N) {
    target += rec_psi * (mu[i] - y[i] + 1) * exp(-mu[i]);
  }
}
