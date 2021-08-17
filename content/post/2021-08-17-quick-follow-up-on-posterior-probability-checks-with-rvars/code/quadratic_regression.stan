
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}

transformed data{
  vector[N] x2;
  
  for (i in 1:N) {
    x2[i] = x[i]*x[i];
  };
  
}

parameters {
  real alpha;
  real beta;
  real gamma;
  real<lower=0> sigma;
}

model {
  y ~ normal(alpha + beta*x + gamma*x2, sigma);
}