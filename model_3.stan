data {
  
  // constants
  int N_session_station ; 
  int N_session ;
  
  // observations
  int session[N_session_station] ;
  int session_station[N_session_station] ;
  
  // REM variables
  real lens ; 
  real days ;
  real day_range[N_session_station] ;
  
  // EDD data
  real det_prob ;
  real max_dist ;
  
  // correction factors
  matrix[9, N_session] ssf_pred ; 
  matrix[9, N_session] issf_pred ; 
  
  // response variables
  int count[N_session_station] ;
  
  }
  
  parameters {

  // density (one cell per camera, 9 x 12 matrix)
  matrix<lower=0>[9, N_session] D_station ;
  matrix<lower=0>[9, N_session] D_station_ssf ;
  matrix<lower=0>[9, N_session] D_station_issf ;
  
  }
  
  transformed parameters {
  
  // r* - effective detection distance
  real<lower=0, upper=3.5> r_star[N_session_station] ;
  
  for (i in 1:N_session_station) {
  
    r_star[i] = sqrt(det_prob * max_dist^2) ;
  
  }
  
  }
  
  model {
  
  // priors
    
    // density (prior scale should be hares/ha)
    to_row_vector(D_station) ~ gamma(1.25, 1.50) ;
    to_row_vector(D_station_ssf) ~ gamma(1.25, 1.50) ;
    to_row_vector(D_station_issf) ~ gamma(1.25, 1.50) ;
    
  // likelihoods
  // naive density
  for (i in 1:N_session_station) {
    
    // Poisson - counts
    target += poisson_lpmf(count[i] | log((to_row_vector(D_station)[i] * 100) * // should be in hares/km2
                                          days *
                                          (day_range[i]) *
                                          (2.0 + lens) *
                                          (r_star[i] * 1000) /
                                          3.14159265)) ;
  
   }
   
   // density corrected with ssf weights
   for (i in 1:N_session_station) {
    
    // Poisson - counts
    target += poisson_lpmf(count[i] | log((to_row_vector(D_station_ssf)[i] * 100) * // should be in hares/km2
                                          days *
                                          (day_range[i]) *
                                          (2.0 + lens) *
                                          (r_star[i] * 1000) /
                                          3.14159265) *
                                          (1 / to_row_vector(ssf_pred)[i])) ;       // correction factor
  
   }
   
   // density corrected with issf UD weights
   for (i in 1:N_session_station) {
    
    // Poisson - counts
    target += poisson_lpmf(count[i] | log((to_row_vector(D_station_issf)[i] * 100) * // should be in hares/km2
                                          days *
                                          (day_range[i]) *
                                          (2.0 + lens) *
                                          (r_star[i] * 1000) /
                                          3.14159265) *
                                          (1 / to_row_vector(issf_pred)[i])) ;           // correction factor in likelihood) ;
  
   }
  
  }
  
  generated quantities {
    
  // calculate site density/CV as a derived parameter - row-wise from the matrix
  real D_site[N_session] ;
  real D_site_ssf[N_session] ;
  real D_site_issf[N_session] ;
  real CV_site[N_session] ;
  real CV_site_ssf[N_session] ;
  real CV_site_issf[N_session] ;
  
  for (i in 1:N_session) {
    
    D_site[i] = mean(D_station[ , i]) ;
    D_site_ssf[i] = mean(D_station_ssf[ , i]) ;
    D_site_issf[i] = mean(D_station_issf[ , i]) ;
    
    CV_site[i] = sd(D_station[ , i]) / D_site[i] ;
    CV_site_ssf[i] = sd(D_station_ssf[ , i]) / D_site_ssf[i] ;
    CV_site_issf[i] = sd(D_station_issf[ , i]) / D_site_issf[i] ;
    
  }
  
  }
  