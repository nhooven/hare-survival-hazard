data {
  
  // number of observations
  int N ;                       // obs
  int nbasis ;               // basis functions
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
  
  // covariates (continuous)
  real mas[N] ;             // standardized mass (kg)
  real hfl[N] ;             // standardized hfl (cm)
  real bci[N] ;             // standardized body condition index (mass/hfl)
  
  // covariates (categorical)
  int collar[N] ;           // indicator for collar type (0 = VHF, 1 = GPS)
  int clust[N];             // index for cluster (1-4)
  int sex[N];               // indicator for sex (0 = F, 1 = M)
  int trt[N];               // indicator for post-treatment
  int ret[N];               // indicator for retention
  int pil[N];               // indicator for piling
  
  }
  
  parameters {
  
  // baseline hazards
  // mort_pred
  real a0_pred ;                     // spline intercept
  vector[nbasis] w0_pred ;        // weights
  
  // mort_oth
  real<lower=0> oth0 ;              // other baseline hazard
  
  // prior for baseline hazard of censoring
  real<lower=0> cens0 ;              // censoring baseline hazard
  
  // coefficients
  // mort_pred
  real bpred_sex ;                       // slope for sex
  real bpred_ret ;                       // slope for retention
  real bpred_trt_r;                      // slope for treatment interaction with ret
  real bpred_pil ;                       // slope for piling
  real bpred_trt_p;                      // slope for treatment interaction with pil
  real bpred_mas ;                       // slope for mass
  real bpred_hfl ;                       // slope for hfl
  real bpred_bci ;                       // slope for body condition
  
  // mort_oth
  real both_sex ;                       // slope for sex
  real both_ret ;                       // slope for retention
  real both_trt_r;                      // slope for treatment interaction with ret
  real both_pil ;                       // slope for piling
  real both_trt_p;                      // slope for treatment interaction with pil
  real both_mas ;                       // slope for mass
  real both_hfl ;                       // slope for hfl
  real both_bci ;                       // slope for body condition
  
  // censoring
  real bcens_col ;                       // slope for collar type
  
  }
  
  model {
  
  // priors (these are on a normal scale)
  // normal scale priors for spline parameters
  w0_pred ~ normal(0, 0.5);                   // weights
  a0_pred ~ normal(0, 2);                     // intercept
  
  oth0 ~ exponential(1) ;
  
  // prior for baseline hazard of censoring
  cens0 ~ exponential(1) ;
  
  // coefficients
  // mort_pred
  bpred_sex ~ normal(0, 2.5);                // normal prior on b_sex
  bpred_ret ~ normal(0, 2.5);                // normal prior on b_ret
  bpred_trt_r ~ normal(0, 2.5);              // normal prior on b_trt_r
  bpred_pil ~ normal(0, 2.5);                // normal prior on b_pil
  bpred_trt_p ~ normal(0, 2.5);              // normal prior on b_trt_p
  bpred_mas ~ normal(0, 2.5);                // normal prior on b_mas
  bpred_hfl ~ normal(0, 2.5);                // normal prior on b_hfl
  bpred_bci ~ normal(0, 2.5);                // normal prior on b_bci
  
  // mort_oth
  both_sex ~ normal(0, 2.5);                // normal prior on b_sex
  both_ret ~ normal(0, 2.5);                // normal prior on b_ret
  both_trt_r ~ normal(0, 2.5);              // normal prior on b_trt_r
  both_pil ~ normal(0, 2.5);                // normal prior on b_pil
  both_trt_p ~ normal(0, 2.5);              // normal prior on b_trt_p
  both_mas ~ normal(0, 2.5);                // normal prior on b_mas
  both_hfl ~ normal(0, 2.5);                // normal prior on b_hfl
  both_bci ~ normal(0, 2.5);                // normal prior on b_bci
  
  // cens
  bcens_col ~ normal(0, 2.5);                // normal prior on b_col
  
  // model
  // mort_pred likelihood
  // linear predictor
  vector[N] bw_pred ;
  vector[N] lambda_pred ;
  vector[N] a_pred ;
  
  bw_pred = to_vector(basis * w0_pred);
  
  for (i in 1:N) {
  
  a_pred[i] = a0_pred ;
  
  lambda_pred[i] = exp(a_pred[i] + bw_pred[i]) *
  
            exp(bpred_sex * sex[i] + 
                bpred_ret * ret[i] + 
                bpred_trt_r * trt[i] * ret[i] +
                bpred_pil * pil[i] +
                bpred_trt_p * trt[i] * pil[i] +
                bpred_mas * mas[i] +
                bpred_hfl * hfl[i] +
                bpred_bci * bci[i]) ;
  
  }
  
  y_mort_pred ~ poisson(lambda_pred) ;
  
  // mort_oth likelihood
  vector[N] lambda_oth ;
  
  for (i in 1:N) {
  
  lambda_oth[i] = oth0 ;
  
  }
  
  y_mort_oth ~ poisson(lambda_oth) ;
  
  // censoring likelihood
  vector[N] lambda_cens ;
  
  for (i in 1:N) {
  
  lambda_cens[i] = cens0 *
  
                exp(bcens_col * collar[i]) ;
  
  }
  
  y_cens ~ poisson(lambda_cens) ;
  
  
  }