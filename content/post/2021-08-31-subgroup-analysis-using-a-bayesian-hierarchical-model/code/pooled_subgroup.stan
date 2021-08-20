data {
  
  int<lower=0> N;
  
  int<lower=0,upper=1> rx[N];
  int<lower=1,upper=8> sub_grp[N];
  real y[N];
  
}

parameters {
  
  real delta_a;
  real delta_m;
  real delta_x;
  
  real<lower=0> sigma;
  real<lower=0> sigma_m;
  real<lower=0> sigma_x;
  
  real alpha[8]; // subgroup-specific intercept

  real t_0;
  real z_a;
  real z_b;
  real z_c;
  real z_ab;
  real z_ac;
  real z_bc;
  real t_abc;
  
}

transformed parameters {
  
  real t_a = sigma_m * z_a + delta_m;
  real t_b = sigma_m * z_b + delta_m;
  real t_c = sigma_m * z_c + delta_m;
  
  real t_ab = sigma_x * z_ab + delta_x;
  real t_ac = sigma_x * z_ac + delta_x;
  real t_bc = sigma_x * z_bc + delta_x;
  
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
  
  alpha ~ normal(delta_a, 10);
  
  sigma ~ normal(0, 10);
  sigma_m ~ normal(0, 1);
  sigma_x ~ normal(0, 1);

  delta_a ~ normal(0, 2);
  delta_m ~ normal(0, 2);
  delta_x ~ normal(0, 2);
  
  t_0 ~ normal(0, 2);
  t_abc ~ normal(0, 2);

  z_a ~ std_normal();
  z_b ~ std_normal();
  z_c ~ std_normal();

  z_ab ~ std_normal();
  z_ac ~ std_normal();
  z_bc ~ std_normal();

  y ~ normal(yhat, sigma);
  
}