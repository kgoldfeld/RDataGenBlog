data {
  
  int<lower=0> S;
  int<lower=0> y[S]; // numerator
  int<lower=0> n[S]; // total observsations (denominator)
  
}

parameters {
  
  real<lower=0,upper=1> theta[S];
  real<lower=0> alpha;
  real<lower=0> beta;
  
}

model {
  
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
  
  theta ~ beta(alpha, beta);
  
  y ~ binomial(n, theta);
  
}

generated quantities {
  
  real mu;
  real s2;
  real disp;
  
  mu = alpha/(alpha + beta);
  s2 = alpha*beta/((alpha + beta)^2 * (alpha + beta + 1));
  disp = (s2 - mu)/(mu^2);
}


