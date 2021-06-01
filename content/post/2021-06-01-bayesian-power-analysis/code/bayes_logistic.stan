data {
  int<lower=0> N;
  int<lower=0,upper=1> y[N];
  vector[N] x;
  real s;
  real mu;
}

parameters {
  real alpha;
  real beta;
}

model {
  
  beta ~ student_t(3, mu, s);
  y ~ bernoulli_logit(alpha + beta * x);
  
}

