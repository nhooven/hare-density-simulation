// corrected densities
real D_site_ssf[N_session] ;
real D_site_issf[N_session] ;
real CV_site_ssf[N_session] ;
real CV_site_issf[N_session] ;

for (i in 1:N_session) {
  
  for (j in 1:)
    
    
    D_site_ssf[i] = mean(mean(ssf_pred) / ssf_pred[i] * D_station[ , i])
  
  
  
  ;
  
  CV_site[i] = sd(D_station[ , i]) / D_site[i] ;
  
}