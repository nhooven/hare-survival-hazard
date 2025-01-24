data {
  
  // number of observations
  int N ;                            // obs
  int nbasis ;                       // basis functions
  int nclust ;                       // clusters
  
  // response variables (binary, 0-1)
  int y_mort_pred[N] ;
  
  // spline matrix
  matrix[N, nbasis] basis ;
  
  // covariates (categorical)
  int clust[N] ;                     // index for cluster (1-4)
  
  }
  
  parameters {
  
  // baseline hazards
  // log-normal scale
  // mort_pred
  // spline intercept (varies by cluster)
  real a0_mean ;                     // mean intercept
  real<lower=0> a0_sigma ;           // SD
  vector[nclust] a0_z ;              // non-centered scaling factor
  
  // spline coefficient mean (constant)
  real w_mean ;                      // mean coefficient, shared by all
  
  // spline coefficient SD, partially pooled with all clusters
  real<lower=0> w_sigma_mean ;       // mean SD
  real<lower=0> w_sigma_sigma ;      // SD of SD (hah!)
  vector[nclust] w_sigma_z ;         // scaling factor for SD
  
  // raw smoothing parameter, partially pooled with all clusters
  real<lower=0> l_raw_mean ;         // mean raw smoothing parameter
  real<lower=0> l_raw_sigma ;        // SD of raw smoothing parameter
  vector[nclust] l_raw_z ;           // scaling factor for smoothing parameter
  
  // spline coefficient scaling factors (vary by cluster and basis function)
  matrix[nclust, nbasis] w_z ;       // non-centered scaling factors

  }
  
  transformed parameters {
  
  // baseline hazards
  // mort_pred 
  matrix[nclust, nbasis] w_pen ;     // penalized spline coefficients
  vector<lower=0>[nclust] w_sigma ;  // coefficient SD by cluster
  vector<lower=1>[nclust] l ;        // smoothing parameter by cluster
  
  // calculate transformed SD and smoothing parameter
  for (i in 1:nclust) {
    
    w_sigma[i] = w_sigma_mean + w_sigma_sigma * w_sigma_z[i] ;
    
    l[i] = exp(l_raw_mean + l_raw_sigma * l_raw_z[i]) ;
    
  }
  
  // loop through clusters (rows) and basis functions (columns)
  for (i in 1:nclust) {
    
    for (j in 1:nbasis) {
      
      w_pen[i, j] = w_mean + (w_sigma[i] * (w_z[i, j] / l[i])) ;
      
    }
    
  }
  
  }
  
  model {
  
  // priors (log-normal scale)
  // mort_pred
  // spline intercept (varies by cluster)
  a0_mean ~ normal(0, 1) ;           // mean intercept
  a0_sigma ~ exponential(1) ;        // SD
  a0_z ~ normal(0, 1) ;              // non-centered scaling factor
  
  // spline coefficients (vary by cluster)
  w_mean ~ normal(0, 1) ;            // mean coefficient by cluster
  w_sigma_mean ~ exponential(1) ;    // mean SD 
  w_sigma_sigma ~ exponential(1) ;   // SD SD 
  w_sigma_z ~ normal(0, 1) ;         // SD scaling factor
  l_raw_mean ~ exponential(1) ;      // mean raw smoothing parameter
  l_raw_sigma ~ exponential(1) ;     // SD raw smoothing parameter
  l_raw_z ~ normal(0, 1) ;           // raw smoothing parameter scaling factor
  
  // spline coefficient scaling factors (vary by cluster and basis function)
  to_vector(w_z) ~ normal(0, 1) ;               // non-centered scaling factors
  
  // sum-to-zero constraint for identifiability
  sum(w_pen) ~ normal(0, 0.0001 * N) ;      
  
  // model
  // mort_pred likelihood
  // linear predictor
  matrix[N, nclust] bw ;             // basis functions weighted by coefficients
  matrix[N, nclust] lambda ;         // Poisson expectation
  
  // multiply basis functions by penalized coefficients
  for (j in 1:nclust) {
      
      bw[ , j] = to_vector(basis * to_vector(w_pen[j, ])) ;
      
  }
      
  // add to intercept and transform (i.e., exponentiate)
  for (i in 1:N) {
    
    for (j in 1:nclust) {
      
      lambda[i, j] = exp(a0_mean + (a0_sigma * a0_z[j]) + bw[i, j]) ;
      
    }
    
  }
      
  // likelihood
  for (i in 1:N) {
      
    target += poisson_lpmf(y_mort_pred[i] | lambda[i, clust[i]]) ;
      
  }
    
}
  