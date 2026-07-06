# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 09 - Random encounter model calculations
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 28 Apr 2025
# Date completed: 29 Apr 2025
# Date last modified: 06 Jul 2026
# R version: 4.5.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation

#_______________________________________________________________________
# 2. Purpose ----
#_______________________________________________________________________

# In this script, we will use both the aggregated detections and speed sampling distributions
# to calculate the random encounter model density estimate.

# Our approach will be to use non-parametric bootstrapping coupled with semi-parametric
# Monte Carlo sampling to fully incorporate uncertainty in speed estimates while
# also constructing 95% confidence envelopes on our mean REM estimates

# First, we'll loop through each replicate, resampling 9 cameras (with replacement)
# for bootstrapping (how many iterations??).

# Then, once the cameras are sampled, we'll take a random draw for the speed estimate 
# and then calculate the mean D for that iteration of the bootstrap. Once all bootstrap
# iterations are completed, we can then summarize with a mean D, CV, and CIs

# Recall that we have 1,000 replicates with 4 abundances and 12 scenarios in each
# so we will need to conduct bootstrapping within each of these 48,000 combinations (per question)

#_______________________________________________________________________
# 3. Read in data ----
#_______________________________________________________________________

# which speeds and contacts
which.ones <- "1_TV1"

# aggregated contacts
contacts <- readRDS(paste0(getwd(), "/data_derived/aggregated_contacts/agg_contacts_", which.ones, ".rds"))

# speeds
pooled.speeds.mean.day <- readRDS(paste0(getwd(), "/data_derived/speed_distributions/speeds_", which.ones, "_mean.rds"))
pooled.speeds.rms.day <- readRDS(paste0(getwd(), "/data_derived/speed_distributions/speeds_", which.ones, "_rms.rds"))

#_______________________________________________________________________
# 4. Define function ----

# for reference:
# i = replicate
# j = scenario (the speed distributions are the same across abundances)
# k = abundance
# l = bootstrap iteration

# YES I know that this is a disgusting nested loop. Deal with it. 

# Benchmarking: with 1,000 bootstrapped iterations, this takes ~ 1 min per replicate
# giving 1000 minutes = ~ 17 hours

# we'll try one with just 100 replicates to start

#_______________________________________________________________________

boot_rem <- function (contacts,
                      speeds,
                      n.iter = 1000) {
  
  # loop through replicates
  start.time <- Sys.time()
  
  all.summary.i <- data.frame()
  
  for (i in 1:1000) {
    
    # subset
    contacts.i <- contacts %>% filter(rep == i)
    speeds.i <- speeds[ , i, ]
    
    # loop through scenarios
    all.summary.j <- data.frame()
    
    for (j in 1:12) {
      
      # ensure that speeds.j has no NAs so we don't need a while loop
      speeds.j <- speeds.i[ , j]
      
      speeds.j <- speeds.j[is.na(speeds.j) == FALSE]
      
      # speeds have to exist
      if (length(speeds.j) == 0) { 
        
        all.summary.j <- data.frame(rep = i,
                                    scenario = j,
                                    abund = c(32, 16, 8, 4),
                                    true.D = c(32, 16, 8, 4) / 10,
                                    mean.REM.D = NA,
                                    cv.REM.D = NA,
                                    l95.REM.D = NA,
                                    u95.REM.D = NA)
        
      } else {
        
        # loop through abundances
        all.summary.k <- data.frame()
        
        for (k in c(32, 16, 8, 4)) {
          
          contacts.k <- contacts.i %>% filter(abund == k)
          
          # bootstrapping
          df.boot.k <- data.frame()
          
          for (l in 1:n.iter) {
            
            # sample cameras with replacement
            samp.cams <- contacts.k[sample(1:9, size = 9, replace = T), ]
            
            # sample speeds
            samp.cams$speed <- sample(speeds.j, size = 9, replace = F)
            
            # calculate mean REM
            samp.cams <- samp.cams %>%
              
              mutate(REM = (points / 28) *  # we'll try points first
                       (pi / (speed * (3.5 / 1000) * (2.0 + ((57.3 * pi) / 180)))) /
                       100)   # per hectare
            
            # iteration df
            df.iter <- data.frame(iter = l,
                                  true.D = k / 10,
                                  REM.D = mean(samp.cams$REM))
            
            # bind in to all iterations
            df.boot.k <- rbind(df.boot.k, df.iter)
            
          }
          
          # calculate mean, CV, and 95% CI from bootstraps
          summary.k <- data.frame(rep = i,
                                  scenario = j,
                                  abund = k,
                                  true.D = k / 10,
                                  mean.REM.D = mean(df.boot.k$REM.D, na.rm = T),
                                  cv.REM.D = sd(df.boot.k$REM.D, na.rm = T) / mean(df.boot.k$REM.D, na.rm = T),
                                  l95.REM.D = quantile(df.boot.k$REM.D, probs = 0.025, na.rm = T),
                                  u95.REM.D = quantile(df.boot.k$REM.D, probs = 0.975, na.rm = T))
          
          # bind in
          all.summary.k <- rbind(all.summary.k, summary.k)
          
        } # k
        
        # bind in
        all.summary.j <- rbind(all.summary.j, all.summary.k)
        
      } # else
      
    } # j
    
    # bind in
    all.summary.i <- rbind(all.summary.i, all.summary.j)
    
    if (i %% 25 == 0) {
      
      elapsed.time <- round(as.numeric(difftime(Sys.time(), 
                                                start.time, 
                                                units = "mins")), 
                            digits = 1)
      
      print(paste0("Completed replicate ", 
                   i, 
                   " of ", 
                   1000, 
                   " - ", 
                   elapsed.time, 
                   " mins"))
      
    }
    
  } # i
  
  # return
  return(all.summary.i)
  
}

#_______________________________________________________________________
# 5. Use function ----
#_______________________________________________________________________

boot.rem.mean <- boot_rem(contacts, pooled.speeds.mean.day, n.iter = 1000)
boot.rem.rms <- boot_rem(contacts, pooled.speeds.rms.day, n.iter = 1000)

#_______________________________________________________________________
# 6. Write files ----
#_______________________________________________________________________

saveRDS(boot.rem.mean, paste0(getwd(), "/data_derived/REM_results/", which.ones, ".rds"))
