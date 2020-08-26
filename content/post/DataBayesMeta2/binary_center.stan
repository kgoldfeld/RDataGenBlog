data {
  
  int<lower=0> N;                // number of observations
  int<lower=0> C;                // number of control types
  int<lower=1> K;                // number of studies
  int y[N];                      // vector of categorical outcomes
  int<lower=1,upper=K> kk[N];    // site for individual
  int<lower=0,upper=1> ctrl[N];  // treatment or control
  int<lower=1,upper=C> cc[K];    // specific control for site
  
}

parameters {
  
  real alpha;               // overall intercept for treatment
  vector[K] beta_0;         // site specific intercept
  real<lower=0> sigma_b;    // sd of site intercepts

  vector[K] delta_k;        // site specific treatment effect
  real<lower=0>  eta_0;     // sd of delta_k (around delta)

  vector[C] delta_c;        // control-specific effect
  real Delta;               // overall control effect
}

transformed parameters{ 
  
  vector[N] yhat;

  for (i in 1:N)  
      yhat[i] = alpha +  beta_0[kk[i]] + (ctrl[i] * (delta_k[kk[i]]));

}

model {
  
  // priors
  
  alpha ~ student_t(3, 0, 2.5);

  beta_0 ~ normal(0, sigma_b);
  
  sigma_b ~ cauchy(0, 1);
  eta_0 ~ cauchy(0, 1);

  for (k in 1:K)
      delta_k[k] ~ normal(delta_c[cc[k]], eta_0);

  delta_c ~ normal(Delta, 0.5);
  Delta ~ normal(0, 10);
  
  
  // outcome model
  
  y ~ bernoulli_logit(yhat);
}
