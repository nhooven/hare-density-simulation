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
  real ssf_pred[N_session_station] ; 
  real issf_pred[N_session_station] ; 
  
  // response variables
  int count[N_session_station] ;
  
  }
  
  parameters {

  // density (one cell per camera, 9 x 21 matrix)
  matrix<lower=0>[9, N_session] D_station ;
  
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
    
    // density (prior scale should be hares/ha) - mean is shape / rate
    to_row_vector(D_station) ~ gamma(5, 1) ;
    
  // likelihoods
  for (i in 1:N_session_station) {
    
    // Poisson - counts
    target += poisson_lpmf(count[i] | log((to_row_vector(D_station)[i] * 100) * // should be in hares/km2
                                          days *
                                          (day_range[i] / 1000) *
                                          (2.0 + lens) *
                                          (r_star[i] * 1000) /
                                          3.14159265)) ;
  
   }
  
  }
  
  generated quantities {
    
  // calculate site density/CV as a derived parameter - row-wise from the matrix
  real D_site[N_session] ;
  real CV_site[N_session] ;
  
  for (i in 1:N_session) {
    
    D_site[i] = mean(D_station[ , i]) ;
    
    CV_site[i] = sd(D_station[ , i]) / D_site[i] ;
    
  }
  
  }
  