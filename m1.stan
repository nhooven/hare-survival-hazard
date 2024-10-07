data {
  
  // number of observations
  int N ;                       // obs
  int nbasis ;                  // basis functions
  int nclust ;                  // clusters
  
  // response variable (binary, 0-1)
  int y_mort[N] ;
  
  // time (required for estimating time-variant hazard);
  // time in week of the year (1-52)
  vector[N] t ;       
  
  // spline matrix
  matrix[N, nbasis] basis ;
  
  // covariates (continuous)
  real mas[N] ;             // standardized mass (kg)
  real hfl[N] ;             // standardized hfl (cm)
  real bci[N] ;             // standardized body condition index (mass/hfl)
  
  // covariates (categorical)
  int clust[N] ;            // index for cluster (1-4)
  int sex[N] ;              // indicator for sex (0 = F, 1 = M)
  int trt[N] ;              // indicator for post-treatment
  int ret[N] ;              // indicator for retention
  int pil[N] ;              // indicator for piling
  
  }
  
  parameters {
  
  real a0 ;                          // spline intercept
  vector[nclust] c0 ;                // spline random intercepts
  real<lower=0> sigma_c0 ;           // random intercept SD
  vector[nbasis] w0 ;                // weight intercepts
  matrix[nbasis, nclust] wx ;        // random weights
  real<lower=0> sigma_w0 ;           // random weight SD
  real b_sex ;                       // slope for sex
  real b_trt_r ;                     // slope for treatment interaction with ret
  real b_trt_p ;                     // slope for treatment interaction with pil
  real b_ret ;                       // slope for retention
  real b_pil ;                       // slope for piling
  real b_mas ;                       // slope for mass
  real b_hfl ;                       // slope for hfl
  real b_bci ;                       // slope for body condition
  
  }
  
  model {
  
  // priors (these are on a normal scale)
  // normal scale priors for spline parameters
  a0 ~ normal(0, 2) ;                     // intercept
  c0[clust] ~ normal(0, 1) ;              // random intercepts
  sigma_c0 ~ exponential(1) ;             // random intercept standard deviation
  w0 ~ normal(0, 1) ;                     // weight intercepts
  wx[nbasis, clust] ~ normal(0, 1) ;      // random weights
  sigma_w0 ~ exponential(1) ;             // random weight standard deviation
    
    
  // coefficients
  b_sex ~ normal(0, 2.5) ;                // normal prior on b_sex
  b_trt_r ~ normal(0, 2.5) ;              // normal prior on b_trt_r
  b_trt_p ~ normal(0, 2.5) ;              // normal prior on b_trt_p
  b_ret ~ normal(0, 2.5) ;                // normal prior on b_ret
  b_pil ~ normal(0, 2.5) ;                // normal prior on b_pil
  b_mas ~ normal(0, 2.5) ;                // normal prior on b_mas
  b_hfl ~ normal(0, 2.5) ;                // normal prior on b_hfl
  b_bci ~ normal(0, 2.5) ;                // normal prior on b_bci
  
  // model
  // linear predictor
  matrix[nbasis, N] w ; 
  vector[N] bw ;
  vector[N] lambda ;
  vector[N] a ;
  
  for (i in 1:N) {
    
  w[ , i] = to_vector(w0 + sigma_w0 * wx[ , clust[i]]) ;      // non-centered random spline weights
  
  bw = to_vector(basis * w) ;                           // apply to basis functions
  
  a[i] = a0 + sigma_c0 * c0[clust[i]] ;                          // non-centered random spline intercept
  
  lambda[i] = exp(a[i] + bw[i]) *
  
            exp(b_sex * sex[i] + 
                b_ret * ret[i] +
                b_trt_r * trt[i] * ret[i] + 
                b_pil * pil[i] + 
                b_trt_p * trt[i] * pil[i] +
                b_mas * mas[i] +
                b_hfl * hfl[i] +
                b_bci * bci[i]) ;
  
  }
  
  y_mort ~ poisson(lambda) ;
  
  }