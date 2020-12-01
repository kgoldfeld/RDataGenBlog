data {

    int<lower=1> N;                // number of observations
    real x[N];                     // antibody measures
    real y[N];                     // outcomes
    
    int<lower=1> M;                // number of candidate thresholds
    real c[M];                     // candidate thresholds
  
}

transformed data {

  real lambda;
  lambda = -log(M);
  
}

parameters {

  real alpha;
  real beta;
  real<lower=0> sigma;

}

transformed parameters {
  
  vector[M] lp;
  lp = rep_vector(lambda, M);
  
  for (m in 1:M)
    for (n in 1:N)
      lp[m] = lp[m] + normal_lpdf(y[n] | x[n] < c[m] ? alpha : beta, sigma);

}

model {
  
  alpha ~ student_t(3, 0, 2.5);
  beta ~ student_t(3, 0, 2.5);
  sigma ~ exponential(1);
  
  target += log_sum_exp(lp);

}
