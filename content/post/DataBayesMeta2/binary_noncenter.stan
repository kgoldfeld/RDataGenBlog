data {
  
  int<lower=0> N;                 // number of observations
  int<lower=0> C;                 // number of control types
  int<lower=1> K;                 // number of studies
  int y[N];                       // vector of categorical outcomes
  int<lower=1,upper=K> kk[N];     // site for individual
  int<lower=0,upper=1> ctrl[N];   // treatment or control
  int<lower=1,upper=C> cc[K];     // specific control for site

}

parameters {
  real alpha;               // overall intercept for treatment
   
  real<lower=0> sigma_b;
  real<lower=0>  eta_0;     // sd of delta_k (around delta)

  real Delta;               // overall control effect
  
   // non-centered paramerization
  
  vector[K] z_ran_rx;   // site level random effects (by period)
  vector[K] z_ran_int;  // individual level random effects 
  vector[C] z_ran_c;  // individual level random effects 

}

transformed parameters{ 
  
  vector[N] yhat;
  vector[K] beta_0;
  vector[K] delta_k;        // site specific treatment effect
  vector[C] delta_c;

  beta_0 = sigma_b * z_ran_int + alpha;
  
  for (i in 1:C)
    delta_c[i] = 0.5 * z_ran_c[i] + Delta; 
  
  for (i in 1:K)
    delta_k[i] = eta_0 * z_ran_rx[i] + delta_c[cc[i]]; 
  
  for (i in 1:N)  
      yhat[i] = beta_0[kk[i]] + ctrl[i] * delta_k[kk[i]];

}

model {
  
  // priors
  
  alpha ~ student_t(3, 0, 2.5);

  z_ran_c ~ std_normal();
  z_ran_int ~ std_normal();  
  z_ran_rx ~ std_normal();  

  sigma_b ~ cauchy(0, 1);
  eta_0 ~ cauchy(0, 1);
  
  Delta ~ normal(0, 10);
  
  // outcome model
  
  y ~ bernoulli_logit(yhat);
}

