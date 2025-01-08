data {
  int<lower=0> N;
  int<lower=0> p;
  
  matrix[N, p] X;
  vector[N] y;
  
  real<lower=0> rec_psi;
}

parameters {
  vector[p] beta;
}

transformed parameters {
  vector[N] eta;
  eta = X * beta;
}

model {
  target += rec_psi * sum(y.*eta - exp(eta));
}
