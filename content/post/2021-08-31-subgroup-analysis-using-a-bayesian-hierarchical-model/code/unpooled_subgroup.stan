data {
  
  int<lower=0> N;
  
  int<lower=0,upper=1> rx[N];
  int<lower=1,upper=8> sub_grp[N];
  real y[N];
  
  real<lower=0> prior_sigma;
  
}

parameters {
  
  real alpha[8]; // subgroup-specific intercept
  real<lower=0> sigma;

  real t_0;
  real t_a;
  real t_b;
  real t_c;
  real t_ab;
  real t_ac;
  real t_bc;
  real t_abc;
  
}

transformed parameters {
  
  
  real theta[8]; // subgroup-specific effect
  
  real yhat[N];
  
  theta[1] = t_0;
  theta[2] = t_0 + t_a;
  theta[3] = t_0 + t_b;
  theta[4] = t_0 + t_c;
  theta[5] = t_0 + t_a + t_b + t_ab;
  theta[6] = t_0 + t_a + t_c + t_ac;
  theta[7] = t_0 + t_b + t_c + t_bc;
  theta[8] = t_0 + t_a + t_b + t_c + t_ab + t_bc + t_bc + t_abc;
  
  for (i in 1:N){
    yhat[i] = alpha[sub_grp[i]] + theta[sub_grp[i]] * rx[i];
  }

}

model {
  
  alpha ~ normal(0, prior_sigma);
  sigma ~ normal(0, prior_sigma);
  
  t_0 ~ normal(0, prior_sigma);
  t_a ~ normal(0, prior_sigma);
  t_b ~ normal(0, prior_sigma);
  t_c ~ normal(0, prior_sigma);
  t_ab ~ normal(0, prior_sigma);
  t_ac ~ normal(0, prior_sigma);
  t_bc ~ normal(0, prior_sigma);
  t_abc ~ normal(0, prior_sigma);

  y ~ normal(yhat, sigma);
  
}
