data {
  
  // number of observations
  int N ;                       // obs
  int nbasis ;                  // basis functions
  int nclust ;                  // clusters
  
  // response variables (binary, 0-1)
  int y_mort_pred[N] ;
  int y_mort_oth[N] ;
  int y_cens[N] ;
  
  // time (required for estimating time-variant hazard);
  // time in week of the year (1-52)
  vector[N] t ;       
  
  // spline matrix
  matrix[N, nbasis] basis ;
  
  // covariates (categorical)
  int clust[N];             // index for cluster (1-4)
  
  }
  
  parameters {
  
  // baseline hazards
  // mort_pred
  real a0_pred ;                     // spline intercept
  real w0_mean ;                     // mean of weights
  real<lower=0> sigma ;              // standard deviation of weights
  vector[nbasis] zw_pred ;           // non-centered scaling factors
  real<lower=0> l_raw ;              // raw smoothing parameter (must be exponentiated)
  
  }
  
  transformed parameters {
    
  // baseline hazard spline
  // actual spline coefficients 
  vector[nbasis] w0_pred ;   
  real l ;                           // smoothing parameter
  
  l = exp(l_raw) ;
  
  w0_pred = w0_mean + (sigma * (zw_pred / l)) ;
    
    
  }
  
  model {
  
  // priors (these are on a normal scale)
  // normal scale priors for spline parameters
  a0_pred ~ normal(0, 2);                     // intercept
  w0_mean ~ normal(0, 0.5) ;                  // mean weight
  sigma ~ exponential(1) ;                    // standard deviation
  zw_pred ~ normal(0, 0.5);                   // scaling factors
  l_raw ~ exponential(1) ;                    // smoothing parameter
  
  // model
  // mort_pred likelihood
  // linear predictor
  
  // baseline hazard spline
  vector[N] bw_pred ;
  vector[N] lambda_pred ;
  
  bw_pred = to_vector(basis * w0_pred);
  
  for (i in 1:N) {
  
  lambda_pred[i] = exp(a0_pred + bw_pred[i]) ;
  
  }
  
  y_mort_pred ~ poisson(lambda_pred) ;
  
  }