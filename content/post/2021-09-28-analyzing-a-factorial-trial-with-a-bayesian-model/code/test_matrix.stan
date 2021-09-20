
data {
  int nrow;
  int ncol;
  
  real y[10];
  int x[nrow, ncol];
}

parameters {
  real mu;
  real<lower=0> sigma;
}

model {
  y ~ normal(mu, sigma);
}

