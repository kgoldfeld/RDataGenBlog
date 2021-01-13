//---------------------------------------------//
//  rstan model predictive probabilities       //
//                                             //
//   author:                KSG                //
//   last modified date:    01/11/2021         //
//                                             //
//---------------------------------------------//
  
data {
  int<lower=0> N;                // number of observations
  int<lower=2> L;                // number of WHO categories
  int<lower=1,upper=L> y[N];     // vector of categorical outcomes
  int<lower=0,upper=1> rx[N];    // treatment or control
  int<lower=1> D;                // number of covariates
  row_vector[D] x[N];            // matrix of covariates  N x D matrix
}

parameters {
  
  vector[D] beta;           // covariate estimates 
  real delta;               // overall control effect
  ordered[L-1] tau;         // cut-points for cumulative odds model ([L-1] vector)
  
}

transformed parameters{ 
  
  vector[N] yhat;

  for (i in 1:N){
    yhat[i] = x[i] * beta + rx[i] * delta;
  }
}

model {
  
  // priors
  
  beta ~ student_t(3, 0, 10);
  delta ~ student_t(3, 0, 2);
  tau ~ student_t(3, 0, 8);
      
  // outcome model
  
  for (i in 1:N)
    y[i] ~ ordered_logistic(yhat[i], tau);
}

generated quantities {
  real OR = exp(delta);
}
