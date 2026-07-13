# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 08 - Speed distributions
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 02 Apr 2025
# Date completed: 
# Date last modified: 13 Jul 2026
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

speeds.1T.TV1 <- readRDS(paste0(getwd(), "/data_derived/sampled_ctmm_speeds/1T_TV1.rds"))
speeds.1T.TV2 <- readRDS(paste0(getwd(), "/data_derived/sampled_ctmm_speeds/1T_TV2.rds"))
speeds.1T.TV3 <- readRDS(paste0(getwd(), "/data_derived/sampled_ctmm_speeds/1T_TV3.rds"))

speeds.1NT.TV1 <- readRDS(paste0(getwd(), "/data_derived/sampled_ctmm_speeds/1NT_TV1.rds"))
speeds.1NT.TV2 <- readRDS(paste0(getwd(), "/data_derived/sampled_ctmm_speeds/1NT_TV2.rds"))
speeds.1NT.TV3 <- readRDS(paste0(getwd(), "/data_derived/sampled_ctmm_speeds/1NT_TV3.rds"))

#_______________________________________________________________________
# 4. Parametric SEs ----

# here we'll assume a normal sampling distribution (in the future we're do this non-parametrically)

#_______________________________________________________________________

# define function
parametric_se <- function (df) {
  
  df.1 <- df %>%
  
  # back-calculate the SD of the sampling dist. (i.e., the SE)
  # we'll use the mean absolute difference from the mean to the 95% CIs
  mutate(parametric.se.mean = (((speed.mean.hi - speed.mean.est) +
                              (speed.mean.est - speed.mean.lo)) / 2) / 1.96,
         parametric.se.rms = (((rms.speed.hi - rms.speed.est) +
                              (rms.speed.est - rms.speed.lo)) / 2) / 1.96)
  
  return(df.1)
  
}

# use function
speeds.1T.TV1.1 <- parametric_se(speeds.1T.TV1)
speeds.1T.TV2.1 <- parametric_se(speeds.1T.TV2)
speeds.1T.TV3.1 <- parametric_se(speeds.1T.TV3)

speeds.1NT.TV1.1 <- parametric_se(speeds.1NT.TV1)
speeds.1NT.TV2.1 <- parametric_se(speeds.1NT.TV2)
speeds.1NT.TV3.1 <- parametric_se(speeds.1NT.TV3)

#_______________________________________________________________________
# 5. Create arrays ----

# x = 500 (draws)
# y = 500 (indivs)
# z = 12 (scenarios)

#_______________________________________________________________________

# define function
speed_array <- function (df, .which = "mean") {
  
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
        
        if (.which == "mean") {
          
          speeds.y <- rnorm(n = 500, mean = df.y$speed.mean.est, sd = df.y$parametric.se.mean)
          
        } else {
          
          speeds.y <- rnorm(n = 500, mean = df.y$rms.speed.est, sd = df.y$parametric.se.rms)
          
        }
        
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
speed.array.1T.TV1.mean <- speed_array(speeds.1T.TV1.1, "mean")
speed.array.1T.TV2.mean <- speed_array(speeds.1T.TV2.1, "mean")
speed.array.1T.TV3.mean <- speed_array(speeds.1T.TV3.1, "mean")

speed.array.1NT.TV1.mean <- speed_array(speeds.1NT.TV1.1, "mean")
speed.array.1NT.TV2.mean <- speed_array(speeds.1NT.TV2.1, "mean")
speed.array.1NT.TV3.mean <- speed_array(speeds.1NT.TV3.1, "mean")

speed.array.1T.TV1.rms <- speed_array(speeds.1T.TV1.1, "rms")
speed.array.1T.TV2.rms <- speed_array(speeds.1T.TV2.1, "rms")
speed.array.1T.TV3.rms <- speed_array(speeds.1T.TV3.1, "rms")

speed.array.1NT.TV1.rms <- speed_array(speeds.1NT.TV1.1, "rms")
speed.array.1NT.TV2.rms <- speed_array(speeds.1NT.TV2.1, "rms")
speed.array.1NT.TV3.rms <- speed_array(speeds.1NT.TV3.1, "rms")

#_______________________________________________________________________
# 5. Create pooled speed distributions ----

# Here we will aggregate individuals for each replicate/scenario combination
# based upon the 12 individuals collared in each approach
# Thus, we will need to include both T and NT individuals

# These pooled distributions are essentially mixtures of normals

#_______________________________________________________________________
# 5a. Read in lookup tables ----
#_______________________________________________________________________

collared.1 <- readRDS(paste0(getwd(), "/data_derived/sampled_reps/collared_1.rds"))
collared.2 <- readRDS(paste0(getwd(), "/data_derived/sampled_reps/collared_2.rds"))
collared.3 <- readRDS(paste0(getwd(), "/data_derived/sampled_reps/collared_3.rds"))

#_______________________________________________________________________
# 5b. Define function ----

# this will create a 6000 (12 x 500) x 1000 x 12 array
# where x = sample draws
# y = replicates
# z = scenarios

#_______________________________________________________________________

pooled_speeds <- function (.collared,
                           .speed.array.T,
                           .speed.array.NT) {
  
  # define blank array
  all.pooled <- array(data = NA,
                      dim = c(6000, 1000, 12))
  
  # loop through reps
  for (y in 1:1000) {
    
    # subset collared individuals
    focal.collared.T <- .collared %>% filter(rep == y, source == "T")
    focal.collared.NT <- .collared %>% filter(rep == y, source == "NT")
    
    # extract correct individuals (columns) from the speed arrays
    # target
    indivs.T <- .speed.array.T[ , focal.collared.T$indiv, ]
    indivs.NT <- .speed.array.NT[ , focal.collared.NT$indiv, ]
  
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

pooled.speeds.1.TV1.mean <- pooled_speeds(collared.1, speed.array.1T.TV1.mean, speed.array.1NT.TV1.mean)
pooled.speeds.1.TV2.mean <- pooled_speeds(collared.1, speed.array.1T.TV2.mean, speed.array.1NT.TV2.mean)
pooled.speeds.1.TV3.mean <- pooled_speeds(collared.1, speed.array.1T.TV3.mean, speed.array.1NT.TV3.mean)

pooled.speeds.1.TV1.rms <- pooled_speeds(collared.1, speed.array.1T.TV1.rms, speed.array.1NT.TV1.rms)
pooled.speeds.1.TV2.rms <- pooled_speeds(collared.1, speed.array.1T.TV2.rms, speed.array.1NT.TV2.rms)
pooled.speeds.1.TV3.rms <- pooled_speeds(collared.1, speed.array.1T.TV3.rms, speed.array.1NT.TV3.rms)

#_______________________________________________________________________
# 6. Convert from m/s to km/day ----
#_______________________________________________________________________

# conversion factor is 1000 / 86400 ~ 0.0115574
pooled.speeds.1.TV1.mean.day <- pooled.speeds.1.TV1.mean / (1000 / 86400)
pooled.speeds.1.TV2.mean.day <- pooled.speeds.1.TV2.mean / (1000 / 86400)
pooled.speeds.1.TV3.mean.day <- pooled.speeds.1.TV3.mean / (1000 / 86400)

pooled.speeds.1.TV1.rms.day <- pooled.speeds.1.TV1.rms / (1000 / 86400)
pooled.speeds.1.TV2.rms.day <- pooled.speeds.1.TV2.rms / (1000 / 86400)
pooled.speeds.1.TV3.rms.day <- pooled.speeds.1.TV3.rms / (1000 / 86400)

#_______________________________________________________________________
# 7. Write to file ----
#_______________________________________________________________________

saveRDS(pooled.speeds.1.TV1.mean.day, paste0(getwd(), "/data_derived/speed_distributions/speeds_1_TV1_mean.rds"))
saveRDS(pooled.speeds.1.TV2.mean.day, paste0(getwd(), "/data_derived/speed_distributions/speeds_1_TV2_mean.rds"))
saveRDS(pooled.speeds.1.TV3.mean.day, paste0(getwd(), "/data_derived/speed_distributions/speeds_1_TV3_mean.rds"))

saveRDS(pooled.speeds.1.TV1.rms.day, paste0(getwd(), "/data_derived/speed_distributions/speeds_1_TV1_rms.rds"))
saveRDS(pooled.speeds.1.TV2.rms.day, paste0(getwd(), "/data_derived/speed_distributions/speeds_1_TV2_rms.rds"))
saveRDS(pooled.speeds.1.TV3.rms.day, paste0(getwd(), "/data_derived/speed_distributions/speeds_1_TV3_rms.rds"))

#_______________________________________________________________________
# 8. Summaries ----
#_______________________________________________________________________
# 8a. Table S3 - percent bias by scenario
#_______________________________________________________________________

# true speeds
true.speeds.1T.TV1 <- readRDS(paste0(getwd(), "/data_derived/sampled_speeds/speeds_1T_TV1.rds"))
true.speeds.1T.TV2 <- readRDS(paste0(getwd(), "/data_derived/sampled_speeds/speeds_1T_TV2.rds"))
true.speeds.1T.TV3 <- readRDS(paste0(getwd(), "/data_derived/sampled_speeds/speeds_1T_TV3.rds"))

true.speeds.1NT.TV1 <- readRDS(paste0(getwd(), "/data_derived/sampled_speeds/speeds_1NT_TV1.rds"))
true.speeds.1NT.TV2 <- readRDS(paste0(getwd(), "/data_derived/sampled_speeds/speeds_1NT_TV2.rds"))
true.speeds.1NT.TV3 <- readRDS(paste0(getwd(), "/data_derived/sampled_speeds/speeds_1NT_TV3.rds"))

# function
speed_bias <- function (.trueT, .trueNT, .estT, .estNT) {
  
  # start with target individuals
  all.speeds <- .trueT |>
    
    # drop unneeded columns
    dplyr::select(-c(question, target)) |>
    
    # join in estimates
    left_join(
      
      .estT |>
        
        # keep only mean estimates for simplicity
        dplyr::select(indiv, scenario, speed.mean.est, sld.speed.mean.est, rms.speed.est)
      
    ) |>
    
    # drop NAs (we can't use them anyway)
    drop_na() |>
    
    # bind in NT individuals
    bind_rows(
      
      .trueNT |>
        
        # drop unneeded columns
        dplyr::select(-c(question, target)) |>
        
        # join in estimates
        left_join(
          
          .estNT |>
            
            # keep only mean estimates for simplicity
            dplyr::select(indiv, scenario, speed.mean.est, sld.speed.mean.est, rms.speed.est)
          
        ) |>
        
        # drop NAs (we can't use them anyway)
        drop_na()
      
    ) |>
    
    # calculate percent bias
    mutate(
      
      bias.estimated = ((speed.mean.est - true.speed) / true.speed * 100),
      bias.rms = ((rms.speed.est - true.speed) / true.speed * 100),
      bias.sld = ((sld.speed.mean.est - true.speed) / true.speed * 100)
      
    ) |>
    
    # summarize
    group_by(scenario) |>
    
    summarize(
      
      bias.estimated.mean = mean(bias.estimated),
      bias.estimated.sd = sd(bias.estimated),
      
      bias.rms.mean = mean(bias.rms),
      bias.rms.sd = sd(bias.rms),
      
      bias.sld.mean = mean(bias.sld),
      bias.sld.sd = sd(bias.sld),
      
    )
  
  return(all.speeds)
  
}

# use
speed.bias.TV1 <- speed_bias(true.speeds.1T.TV1,
                             true.speeds.1NT.TV1,
                             speeds.1T.TV1,
                             speeds.1NT.TV1)

speed.bias.TV2 <- speed_bias(true.speeds.1T.TV2,
                             true.speeds.1NT.TV2,
                             speeds.1T.TV2,
                             speeds.1NT.TV2)

speed.bias.TV3 <- speed_bias(true.speeds.1T.TV3,
                             true.speeds.1NT.TV3,
                             speeds.1T.TV3,
                             speeds.1NT.TV3)

# write to clipboard
write.table(speed.bias.TV1, "clipboard", sep = "\t")
write.table(speed.bias.TV2, "clipboard", sep = "\t")
write.table(speed.bias.TV3, "clipboard", sep = "\t")

# 07-13-2026
# we'll still need to do the Q2 individuals

#_______________________________________________________________________
# 8b. Table S4 - Share of top-performing CTMM types
#_______________________________________________________________________

top_ctmms <- function (.estT, .estNT) {
  
  top.ctmms <- .estT |>
    
    bind_rows(.estNT) |>
    
    # keep relevant columns
    dplyr::select(scenario, top.model) |>
    
    group_by(scenario, top.model) |>
    
    summarize(total = n()) |>
    
    # pivot wider
    pivot_wider(names_from = top.model,
                values_from = total)
  
  return(top.ctmms)
  
}

# use
top.TV1 <- top_ctmms(speeds.1T.TV1, speeds.1NT.TV1)
top.TV2 <- top_ctmms(speeds.1T.TV2, speeds.1NT.TV2)
top.TV3 <- top_ctmms(speeds.1T.TV3, speeds.1NT.TV3)

# write to clipboard
write.table(top.TV1, "clipboard", sep = "\t")
write.table(top.TV2, "clipboard", sep = "\t")
write.table(top.TV3, "clipboard", sep = "\t")

# 07-13-2026
# we'll still need to do the Q2 individuals