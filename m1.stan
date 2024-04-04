data {
  
  // number of observations
  int N;                       // obs
  int num_basis;               // basis functions
  int nclust;                  // clusters
  
  // response variable (binary, 0-1)
  int y_mort[N];
  
  // time (required for estimating time-variant hazard);
  // time in week of the year (1-52)
  vector[N] t;       
  
  // spline matrix
  matrix[N, num_basis] basis;
  
  // covariates (continuous)
  real mas[N];             // standardized mass (kg)
  real hfl[N];             // standardized hfl (cm)
  real bci[N];             // standardized body condition index (mass/hfl)
  
  // covariates (categorical)
  int clust[N];            // index for cluster (1-4)
  int sex[N];              // indicator for sex (0 = F, 1 = M)
  int ret[N];              // indicator for retention
  int pil[N];              // indicator for piling
  
  }
  
  parameters {
  
  real a0;                          // spline intercept
  vector[nclust] c0;                // spline random intercepts
  real<lower=0> sigma;              // random intercept SD
  vector[num_basis] w0;             // weights
  real b_sex;                       // slope for sex
  real b_ret;                       // slope for retention
  real b_pil;                       // slope for piling
  real b_mas;                       // slope for mass
  real b_hfl;                       // slope for hfl
  real b_bci;                       // slope for body condition
  
  }
  
  model {
  
  // priors (these are on a normal scale)
  // normal scale priors for spline parameters
  w0 ~ normal(0, 1);                     // weights
  a0 ~ normal(0, 2);                     // intercept
  c0[clust] ~ normal(0, 1);              // random intercepts
  sigma ~ exponential(1);                // random intercept standard deviation
  
  // coefficients
  b_sex ~ normal(0, 2.5);                // normal prior on b_sex
  b_ret ~ normal(0, 2.5);                // normal prior on b_ret
  b_pil ~ normal(0, 2.5);                // normal prior on b_pil
  b_mas ~ normal(0, 2.5);                // normal prior on b_mas
  b_hfl ~ normal(0, 2.5);                // normal prior on b_hfl
  b_bci ~ normal(0, 2.5);                // normal prior on b_bci
  
  // model
  // linear predictor
  vector[N] bw;
  vector[N] lambda;
  vector[N] a;
  
  bw = to_vector(basis * w0);
  
  for (i in 1:N) {
  
  a[i] = a0 + sigma * c0[clust[i]];
  
  lambda[i] = exp(a[i] + bw[i]) *
  
            exp(b_sex * sex[i] + 
                b_ret * ret[i] + 
                b_pil * pil[i] +
                b_mas * mas[i] +
                b_hfl * hfl[i] +
                b_bci * bci[i]);
  
  }
  
  y_mort ~ poisson(lambda);
  
  }
  
  generated quantities {
  
  // cluster-specific intercepts
  real a_c1 = exp(a0 + sigma * c0[1]);
  real a_c2 = exp(a0 + sigma * c0[2]);
  real a_c3 = exp(a0 + sigma * c0[3]);
  real a_c4 = exp(a0 + sigma * c0[4]);
  
  // hazard ratios
  real hr_sex = exp(b_sex);
  real hr_ret = exp(b_ret);
  real hr_pil = exp(b_pil);
  real hr_mas = exp(b_mas);
  real hr_hfl = exp(b_hfl);
  real hr_bci = exp(b_bci);
  
  }