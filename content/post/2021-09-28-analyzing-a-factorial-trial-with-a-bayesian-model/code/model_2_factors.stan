data {
  
  int<lower=1> N;        // number of observations
  int<lower=1> I;        // number of interventions
  int<lower=1> X2;       // number of 2-way interactions
  int main[N, I];        // interventions
  int x[N, X2];          // interactions - provide levels for each intervention?
  vector[N] y;           // outcome
  
}

parameters {
  
  real t_0;
  
  vector[3] z_raw[I];
  vector[15] z_x_raw[X2];
  
  real<lower=0> sigma;
  real<lower=0> sigma_m[I];
  real<lower=0> sigma_x[X2];
  
  real<lower=0> sigma_main;

}

transformed parameters {
  
  // constrain parameters to sum to 0
  
  vector[4] z[I]; 
  vector[16] z_x[X2]; 
  
  vector[4] t[I];
  vector[16] t_x[X2];
  
  vector[N] yhat;
  
  for (i in 1:I) {
    z[i] = append_row(z_raw[i], -sum(z_raw[i]));    
  }
  
  for (i in 1:X2) {
    z_x[i] = append_row(z_x_raw[i], -sum(z_x_raw[i]));    
  }

  for (i in 1:I) 
     for (j in 1:4) 
        t[i, j] = sigma_m[i] * z[i, j];
        
  for (i in 1:X2) 
     for (j in 1:16) 
        t_x[i, j] = sigma_x[i] * z_x[i, j];
     
  // yhat
  
  for (n in 1:N) {
    real ytemp; 
    ytemp = t_0;
    for (i in 1:I) ytemp = ytemp + t[i, main[n, i]]; // 2 sets of main effects
    for (i in 1:X2) ytemp = ytemp + t_x[i, x[n, i]]; // 1 set of interaction effects
    yhat[n] = ytemp;
  }
}

model {
  
  sigma ~ normal(0, 5);
  sigma_m ~ normal(0, sigma_main);
  sigma_x ~ normal(0, 5);
  
  sigma_main ~ normal(0, 5);
  
  t_0 ~ normal(0, 5);

  for (i in 1:I) z_raw[i] ~ std_normal();
  for (i in 1:X2) z_x_raw[i] ~ std_normal();

  y ~ normal(yhat, sigma);

}
