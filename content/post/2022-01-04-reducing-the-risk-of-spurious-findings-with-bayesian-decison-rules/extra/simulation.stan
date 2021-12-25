//---------------------------------------------//
//  rstan model for Plasma interim evaluation  //
//                                             //
//   author:                KSG                //
//   last modified date:    05/29/2020         //
//                                             //
//---------------------------------------------//

data {
  int<lower=0> N;              // number of observations
  int<lower=1> L;              // number of levels
  int<lower=0,upper=1> y[N];   // vector of categorical outcomes
  int<lower=0,upper=1> a[N];   // treatment arm for individual
  int<lower=1,upper=4> grp[N]; // grp for individual  
}

parameters {
  vector[L] alpha_g;           // group effect
  real delta_g[L];             // group treatment effects
  real delta;                  // overall treatment effect
}

transformed parameters{ 
  
  vector[N] yhat;

  for (i in 1:N)  
    yhat[i] = alpha_g[grp[i]] + a[i] * delta_g[grp[i]];
}

model {
  
  // priors
  
  alpha_g ~ normal(0, 10);
  delta_g ~ normal(delta, 0.3537);
  delta ~ normal(0, 0.3537);
  
  // outcome model
  
  for (i in 1:N)
    y[i] ~  bernoulli_logit(yhat[i]);
}

generated quantities {
  
  real<lower = 0> OR_g[L];
  real<lower = 0> Odds_g[L];

  real<lower = 0> OR;

  for (i in 1:L) {
    OR_g[i] = exp(delta_g[i]);   
    Odds_g[i] = exp(alpha_g[i]);
  }
  
  OR = exp(delta);
  
}

