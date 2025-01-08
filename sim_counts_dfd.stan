
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

model {
  vector[N] mu;
  mu = exp(X * beta);
  target += rec_psi * sum( - ((y^2)./(mu^2)) + 2*((y+1)./mu));
}
