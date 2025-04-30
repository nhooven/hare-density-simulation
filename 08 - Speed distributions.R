# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 08 - Speed distributions
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 02 Apr 2025
# Date completed: 
# Date last modified: 28 Apr 2025
# R version: 4.4.3

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation

#_______________________________________________________________________
# 2. Purpose ----
#_______________________________________________________________________

# Here to save computation time during the REM sampling/calculation procedure,
# we'll pre-load samples from the Gaussian sampling distributions of speeds
# and save them to 4 arrays that we can easily sample from within the
# Monte Carlo sampling step

# We can also use these data for visualization, showing what the models imply

# Here we'll sample 500 draws, giving a 500 x 500 x 12 array (cols are indivs, rows are draws, dim 3 is scenario)

#_______________________________________________________________________
# 3. Read in data ----
#_______________________________________________________________________

speeds.1T <- read.csv(paste0(getwd(), "/Derived data/Sampled - CTMM speeds/ctmm_speeds_1T.csv"))
speeds.1NT <- read.csv(paste0(getwd(), "/Derived data/Sampled - CTMM speeds/ctmm_speeds_1NT.csv"))
speeds.2T.1 <- read.csv(paste0(getwd(), "/Derived data/Sampled - CTMM speeds/ctmm_speeds_2T.csv"))
speeds.2T.2 <- read.csv(paste0(getwd(), "/Derived data/Sampled - CTMM speeds/ctmm_speeds_2T_2.csv"))
speeds.2NT <- read.csv(paste0(getwd(), "/Derived data/Sampled - CTMM speeds/ctmm_speeds_2NT.csv"))

# bind 2T together
speeds.2T <- rbind(speeds.2T.1, speeds.2T.2)

rm(speeds.2T.1, speeds.2T.2)

#_______________________________________________________________________
# 4. Parametric SEs ----

# here we'll assume a normal sampling distribution (in the future we're do this non-parametrically)

#_______________________________________________________________________

# define function
parametric_se <- function (df) {
  
  df.1 <- df %>%
  
  # back-calculate the SD of the sampling dist. (i.e., the SE)
  # we'll use the mean absolute difference from the mean to the 95% CIs
  mutate(parametric.se = (((speed.mean.hi - speed.mean.est) +
                          (speed.mean.est - speed.mean.lo)) / 2) / 1.96)
  
  return(df.1)
  
}

# use function
speeds.1T <- parametric_se(speeds.1T)
speeds.1NT <- parametric_se(speeds.1NT)
speeds.2T <- parametric_se(speeds.2T)
speeds.2NT <- parametric_se(speeds.2NT)

#_______________________________________________________________________
# 5. Create arrays ----

# x = 500 (draws)
# y = 500 (indivs)
# z = 12 (scenarios)

#_______________________________________________________________________

# define function
speed_array <- function (df) {
  
  speeds.all <- array(data = NA,
                      dim = c(500, 500, 12))
  
  for (z in 1:12) {
    
    # subset by scenario
    df.z <- df %>% filter(scenario == z)
    
    speeds.z <- matrix(data = NA,
                       nrow = 500,
                       ncol = 500)
    
    for (y in 1:max(df.z$indiv)) {
      
      df.y <- df.z %>% filter(indiv == y)
      
      # sample using mean and parametric SE
      # ensure that non-velocity models get NAs here
      if (is.na(df.y$speed.mean.est) == F) {
        
        speeds.y <- rnorm(n = 500, mean = df.y$speed.mean.est, sd = df.y$parametric.se)
        
      } else {
        
        speeds.y <- rep(NA, times = 500)
        
      }
      
      # bind in by column
      speeds.z[ , y] <- speeds.y
      
    }
    
    # bind in by depth of array
    speeds.all[ , , z] <- speeds.z
    
  }
  
  # return array
  return(speeds.all)
  
}

# use function
speed.array.1T <- speed_array(speeds.1T)
speed.array.1NT <- speed_array(speeds.1NT)
speed.array.2T <- speed_array(speeds.2T)
speed.array.2NT <- speed_array(speeds.2NT)

#_______________________________________________________________________
# 5. Create pooled speed distributions ----

# Here we will aggregate individuals for each replicate/scenario combination
# based upon the 12 individuals collared in each approach
# Thus, we will need to include both T and NT individuals

# These pooled distributions are essentially mixtures of normals

#_______________________________________________________________________
# 5a. Read in lookup tables ----
#_______________________________________________________________________

collared.1 <- read.csv(paste0(getwd(), "/Derived data/Sampled reps/collared_1.csv"))
collared.2 <- read.csv(paste0(getwd(), "/Derived data/Sampled reps/collared_2.csv"))
collared.3 <- read.csv(paste0(getwd(), "/Derived data/Sampled reps/collared_3.csv"))

#_______________________________________________________________________
# 5b. Define function ----

# this will create a 6000 (12 x 500) x 1000 x 12 array
# where x = sample draws
# y = replicates
# z = scenarios

#_______________________________________________________________________

pooled_speeds <- function (speed.array.T,
                           speed.array.NT,
                           q = 1) {
  
  # use correct "collared" file
  if (q == 1) {
    
    collared <- collared.1
    
  } 
  
  if (q == 2) {
    
    collared <- collared.2
    
  } 
  
  if (q == 3) {
    
    collared <- collared.3
    
  } 
  
  # define blank array
  all.pooled <- array(data = NA,
                      dim = c(6000, 1000, 12))
  
  # loop through reps
  for (y in 1:1000) {
    
    # subset collared individuals
    focal.collared.T <- collared %>% filter(rep == y, source == "T")
    focal.collared.NT <- collared %>% filter(rep == y, source == "NT")
    
    # extract correct individuals (columns) from the speed arrays
    # target
    indivs.T <- speed.array.T[ , focal.collared.T$indiv, ]
    indivs.NT <- speed.array.NT[ , focal.collared.NT$indiv, ]
  
    # loop through scenarios
    all.speeds.z <- array(data = NA,
                          dim = c(6000, 1, 12))
    
    for (z in 1:12) {
      
      # extract as 1D vector and bind together
      speeds.z <- c(as.vector(indivs.T[ , , z]),
                    as.vector(indivs.NT[ , , z]))
      
      # bind into array
      all.speeds.z[ , , z] <- speeds.z
      
    }
    
    # bind columns
    all.pooled[ , y, ] <- all.speeds.z
  
  }
  
  # return
  return(all.pooled)
  
}

#_______________________________________________________________________
# 5c. Use function ----
#_______________________________________________________________________

pooled.speeds.1 <- pooled_speeds(speed.array.1T, speed.array.1NT, q = 1)
pooled.speeds.2 <- pooled_speeds(speed.array.2T, speed.array.2NT, q = 2)
pooled.speeds.3 <- pooled_speeds(speed.array.2T, speed.array.2NT, q = 3)

#_______________________________________________________________________
# 6. Convert from m/s to km/day ----
#_______________________________________________________________________

# conversion factor is 1000 / 86400 ~ 0.0115574
pooled.speeds.1.day <- pooled.speeds.1 / (1000 / 86400)
pooled.speeds.2.day <- pooled.speeds.2 / (1000 / 86400)
pooled.speeds.3.day <- pooled.speeds.3 / (1000 / 86400)

#_______________________________________________________________________
# 7. Write to file ----
#_______________________________________________________________________

save(pooled.speeds.1.day, file = paste0(getwd(), "/Derived data/Speed distributions/speeds_1.RData"))
save(pooled.speeds.2.day, file = paste0(getwd(), "/Derived data/Speed distributions/speeds_2.RData"))
save(pooled.speeds.3.day, file = paste0(getwd(), "/Derived data/Speed distributions/speeds_3.RData"))
