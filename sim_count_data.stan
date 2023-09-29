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
  real eta[N];
  real mu[N];
  for (i in 1:N) {
    eta[i] = beta[1]*intercept[i] + beta[2]*v1[i] + beta[3]*v2[i] + beta[4]*v3[i];
    mu[i] = exp(eta[i]);
  }
}

model {
  for (i in 1:N) {
    target += rec_psi * (y[i]*eta[i] - mu[i]);
  }
}
