data {
  int<lower=1> N;  // number of observations
  int<lower=1> K;  // number of clusters
  int<lower=1> M;  // number of basis functions
  array[N] int<lower=1, upper=K> cluster;  // cluster indicator
  matrix[N, M] X_spline;  // spline basis functions
  matrix[N, M] D2_spline;  // second derivative of spline basis functions
  vector[N] y;  // response variable
  real<lower=0> lambda;
}

parameters {
  matrix[K, M] beta;  // cluster-specific spline coefficients
  real<lower=0> sigma_y;  // observation noise
  real<lower=0> sigma_beta;  // prior standard deviation for beta
}

model {
  sigma_y ~ normal(0, 1);
  sigma_beta ~ normal(0, 1);

  // Priors for beta
  for (k in 1:K) {
    beta[k] ~ normal(0, sigma_beta);
  }
  
  //Penalization
  for (k in 1:K) {
    target += -lambda * sum(square(D2_spline * beta[k]'));
  }
  
  // Likelihood
  for (n in 1:N) {
    y[n] ~ normal(X_spline[n] * beta[cluster[n]]', sigma_y);
  }
}

generated quantities {
  
  vector[N] y_pred;                    // Vector of observations.
  
  for (n in 1:N) {
    y_pred[n] = normal_rng(X_spline[n] * beta[cluster[n]]', sigma_y);
  }
}

