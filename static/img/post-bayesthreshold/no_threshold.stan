data {

    int<lower=1> N;                // number of observations
    real t[N];                     // threshold
    real y[N];                     // outcomes
}

parameters {

  real alpha;
  real beta;
  real<lower=0> sigma;

}

transformed parameters {
  vector[N] yhat;
  
  for (n in 1:N) {
    yhat[n] = alpha * (1 - t[n]) + beta * t[n];
  }
  
}


model {
  
  alpha ~ student_t(3, 0, 5);
  beta ~ student_t(3, 0, 5);
  sigma ~ exponential(1);
  
  y ~ normal(yhat, sigma);

}
