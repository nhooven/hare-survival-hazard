data {
  
  // number of observations
  int N ;                            // obs
  int nbasis ;                       // basis functions
  int nclust ;                       // clusters
  
  // response variables (binary, 0-1)
  int y_mort[N] ;
  int y_cens[N] ;
  
  // spline matrix
  matrix[N, nbasis] basis ;
  
  // cluster variable for varying effects
  int clust[N] ;                     // index for cluster (1-4)
  
  // censoring covariate
  int collar[N] ;                    // (0 = VHF?)
  
  // covariates (intrinsic)
  int sex[N] ;                       // indicator for sex (0 = female)
  real mas[N] ;                      // body mass at capture
  real hfl[N] ;                      // hind foot length at capture
  
  // covariates (treatment)
  int trt[N] ;                       // indicator for before/after treatment
  int ret[N] ;                       // indicator for retention treatment
  int pil[N] ;                       // indicator for piling treatment
  
  }
  
  parameters {
  
  // baseline hazards
  // "normal"" scale
  // mort
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
  
  // censoring 
  real h0_cens_norm ;            // baseline hazard of censoring
  real b_col ;                   // coefficient for collar
  
  // covariates (non-centered varying slopes)
  // intrinsic
  real b_sex_mean ;
  real<lower=0> b_sex_sigma ;
  vector[nclust] b_sex_z ;
  
  //real b_mas_mean ;
  //real<lower=0> b_mas_sigma ;
  //vector[nclust] b_mas_z ;
  
  real b_hfl_mean ;
  real<lower=0> b_hfl_sigma ;
  vector[nclust] b_hfl_z ;
  
  // treatment
  // pre-post interaction coefficients
  real b_trt_ret_mean ;
  real<lower=0> b_trt_ret_sigma ;
  vector[nclust] b_trt_ret_z ;
  
  real b_trt_pil_mean ;
  real<lower=0> b_trt_pil_sigma ;
  vector[nclust] b_trt_pil_z ;
  
  // treatment coefficients
  real b_ret_mean ;
  real<lower=0> b_ret_sigma ;
  vector[nclust] b_ret_z ;
  
  real b_pil_mean ;
  real<lower=0> b_pil_sigma ;
  vector[nclust] b_pil_z ;

  }
  
  transformed parameters {
  
  // baseline hazards
  // mort
  matrix[nclust, nbasis] w_pen ;     // penalized spline coefficients
  vector[nclust] w_sigma ;           // coefficient SD by cluster
  vector[nclust] l ;                 // smoothing parameter by cluster
  
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
  
  // censoring
  real h0_cens ;        // correctly scaled baseline hazard of censoring
  
  h0_cens = exp(h0_cens_norm) ;
    
  // varying slope coefficients (by cluster)
  // intrinsic
  vector[nclust] b_sex ;
  //vector[nclust] b_mas ;
  vector[nclust] b_hfl ;
  
  // treatment
  vector[nclust] b_trt_ret ;
  vector[nclust] b_trt_pil ;
  vector[nclust] b_ret ;
  vector[nclust] b_pil ;
  
  for (i in 1:nclust) {
    
    b_sex[i] = b_sex_mean + (b_sex_sigma * b_sex_z[i]) ;
    //b_mas[i] = b_mas_mean + (b_mas_sigma * b_mas_z[i]) ;
    b_hfl[i] = b_hfl_mean + (b_hfl_sigma * b_hfl_z[i]) ;
    
    b_trt_ret[i] = b_trt_ret_mean + (b_trt_ret_sigma * b_trt_ret_z[i]) ;
    b_trt_pil[i] = b_trt_pil_mean + (b_trt_pil_sigma * b_trt_pil_z[i]) ;
    b_ret[i] = b_ret_mean + (b_ret_sigma * b_ret_z[i]) ;
    b_pil[i] = b_pil_mean + (b_pil_sigma * b_pil_z[i]) ;
    
  }
  
  }
  
  model {
  
  // priors
  // mort
  // spline intercept (varies by cluster)
  a0_mean ~ normal(log(0.01), 1) ;   // mean intercept (centered on 0.01 baseline hazard)
  a0_sigma ~ exponential(1) ;        // SD
  a0_z ~ normal(0, 0.5) ;            // non-centered scaling factor
  
  // spline coefficients (vary by cluster)
  w_mean ~ normal(0, 1) ;            // mean coefficient by cluster
  w_sigma_mean ~ exponential(1) ;    // mean SD 
  w_sigma_sigma ~ exponential(1) ;   // SD SD 
  w_sigma_z ~ normal(0, 0.5) ;       // SD scaling factor
  l_raw_mean ~ exponential(1) ;      // mean raw smoothing parameter
  l_raw_sigma ~ exponential(1) ;     // SD raw smoothing parameter
  l_raw_z ~ normal(0, 0.5) ;         // raw smoothing parameter scaling factor
  
  // spline coefficient scaling factors (vary by cluster and basis function)
  to_vector(w_z) ~ normal(0, 0.5) ;  // non-centered scaling factors
  
  // sum-to-zero constraint for identifiability
  sum(w_pen) ~ normal(0, 0.0001 * N) ;    
  
  // censoring
  h0_cens_norm ~ normal(log(0.01), 1) ;  // baseline hazard of censoring (0.01 baseline hazard)
  b_col ~ cauchy(0, 2.5) ;               // coefficient for collar
  
  // coefficients (non-centered varying slopes)
  // intrinsic
  b_sex_mean ~ cauchy(0, 2.5) ;
  b_sex_sigma ~ exponential(1) ;
  b_sex_z ~ normal(0, 0.5) ;
  
  //b_mas_mean ~ normal(0, 2.5) ;
  //b_mas_sigma ~ exponential(1) ;
  //b_mas_z ~ normal(0, 0.5) ;
  
  b_hfl_mean ~ cauchy(0, 2.5) ;
  b_hfl_sigma ~ exponential(1) ;
  b_hfl_z ~ normal(0, 0.5) ;
  
  // treatment
  // pre-post interaction coefficients
  b_trt_ret_mean ~ cauchy(0, 2.5) ;
  b_trt_ret_sigma ~ exponential(1) ;
  b_trt_ret_z ~ normal(0, 0.5) ;
  
  b_trt_pil_mean ~ cauchy(0, 2.5) ;
  b_trt_pil_sigma ~ exponential(1) ;
  b_trt_pil_z ~ normal(0, 0.5) ;
  
  // treatment coefficients
  b_ret_mean ~ cauchy(0, 2.5) ;
  b_ret_sigma ~ exponential(1) ;
  b_ret_z ~ normal(0, 0.5) ;
  
  b_pil_mean ~ cauchy(0, 2.5) ;
  b_pil_sigma ~ exponential(1) ;
  b_pil_z ~ normal(0, 0.5) ;
  
  // model
  // mort
  // linear predictor
  matrix[N, nclust] bw ;             // basis functions weighted by coefficients
  matrix[N, nclust] lambda_mort ;    // Poisson expectation
  
  // multiply basis functions by penalized coefficients
  for (j in 1:nclust) {
      
      bw[ , j] = to_vector(basis * to_vector(w_pen[j, ])) ;
      
  }
      
  // add to intercept and transform (i.e., exponentiate), multiply by coefficient vector
  for (i in 1:N) {
    
    for (j in 1:nclust) {
      
      lambda_mort[i, j] = exp(a0_mean + (a0_sigma * a0_z[j]) + bw[i, j]) *  
                          exp(b_sex[j] * sex[i] +
                              //b_mas[j] * mas[i] +
                              b_hfl[j] * mas[i] +
                              b_ret[j] * ret[i] + b_trt_ret[j] * trt[i] * ret[i] +
                              b_pil[j] * pil[i] + b_trt_pil[j] * trt[i] * pil[i]) ;
      
    }
    
  }
  
  // censoring
  vector[N] lambda_cens ;
  
  for (i in 1:N) {
    
    lambda_cens[i] = h0_cens * exp(b_col * collar[i]) ;
    
  }
      
  // "joint" likelihood
  for (i in 1:N) {
      
    target += poisson_lpmf(y_mort[i] | lambda_mort[i, clust[i]]) +
              poisson_lpmf(y_cens[i] | lambda_cens[i]) ;
      
  }
  
  } 
  
  generated quantities {
    
    // hazard ratios (exp(coefs))
    real hr_sex ; 
    //real hr_mas ;
    real hr_hfl ;
    real hr_col ;
    
    // total treatment effects (also hazard ratios)
    real hr_ret_total ;
    real hr_pil_total ;
    
    hr_sex = exp(b_sex_mean) ;
    //hr_mas = exp(b_mas_mean) ;
    hr_hfl = exp(b_hfl_mean) ;
    hr_col = exp(b_col) ;
      
    hr_ret_total = exp(b_ret_mean + b_trt_ret_mean) ;
    hr_pil_total = exp(b_pil_mean + b_trt_pil_mean) ;
    
  }
  