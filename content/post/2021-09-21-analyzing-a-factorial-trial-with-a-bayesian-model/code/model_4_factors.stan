data {
  
  int<lower=0> N;
  int<lower=1, upper=4> a[N];     // intervention a
  int<lower=1, upper=4> b[N];     // intervention b
  int<lower=1, upper=16> ab[N];   // interaction a & b
  vector[N] y;
  
}

parameters {
  
  real t_0;
  
  vector[3] z_a_raw;
  vector[3] z_b_raw;
  vector[15] z_ab_raw;
  
  real<lower=0> sigma;
  real<lower=0> sigma_a;
  real<lower=0> sigma_b;
  real<lower=0> sigma_ab;
  
  real<lower=0> sigma_main;

}

transformed parameters {
  
  // constrain parameters to sum to 0
  
  vector[4] z_a = append_row(z_a_raw, -sum(z_a_raw));
  vector[4] z_b = append_row(z_b_raw, -sum(z_b_raw));
  vector[16] z_ab = append_row(z_ab_raw, -sum(z_ab_raw));
  
  vector[4] t_a;
  vector[4] t_b;
  vector[16] t_ab;
  
  vector[N] yhat;
   
  for (i in 1:4) 
    t_a[i] = sigma_a * z_a[i];
    
  for (i in 1:4) 
    t_b[i] = sigma_b * z_b[i];
    
  for (i in 1:16) 
    t_ab[i] = sigma_ab * z_ab[i];
  
  // yhat
  
  for (i in 1:N){
    yhat[i] = t_0 + t_a[a[i]] + t_b[b[i]] + t_ab[ab[i]];
  }
}

model {
  
  sigma ~ normal(0, 5);
  sigma_a ~ normal(0, sigma_main);
  sigma_b ~ normal(0, sigma_main);
  sigma_ab ~ normal(0, 5);
  
  sigma_main ~ normal(0, 5);
  
  t_0 ~ normal(0, 5);

  z_a_raw ~ std_normal();
  z_b_raw ~ std_normal();
  z_ab_raw ~ std_normal();

  y ~ normal(yhat, sigma);

}
