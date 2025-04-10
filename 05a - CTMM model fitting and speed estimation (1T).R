# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 05a - CTMM model fitting and speed estimation (1T)
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 10 Apr 2025
# Date completed: 10 Apr 2025
# Date last modified: 10 Apr 2025
# R version: 4.4.3

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation
library(lubridate)            # work with dates
library(ctmm)                 # CTSP movement modeling
library(amt)                  # work with animal movement tracks
library(sf)                   # spatial operations

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

# "GPS" tracks
all.tracks <- read.csv("Derived data/Sampled - GPS tracks/tracks_1T.csv")

#_______________________________________________________________________
# 3. Define scenarios ----
#_______________________________________________________________________

# S1: 0.5 hr fix rate, 100% fix success, 4 week duration (no need to modify track)
# S2: 0.5 hr fix rate, 100% fix success, 2 week duration 
# S3: 0.5 hr fix rate, 60% fix success, 4 week duration 
# S4: 0.5 hr fix rate, 60% fix success, 2 week duration 
# S5: 1 hr fix rate, 100% fix success, 4 week duration 
# S6: 1 hr fix rate, 100% fix success, 2 week duration
# S7: 1 hr fix rate, 60% fix success, 4 week duration 
# S8: 1 hr fix rate, 60% fix success, 2 week duration 
# S9: 4 hr fix rate, 100% fix success, 4 week duration 
# S10: 4 hr fix rate, 100% fix success, 2 week duration 
# S11: 4 hr fix rate, 60% fix success, 4 week duration 
# S12: 4 hr fix rate, 60% fix success, 2 week duration 

scenarios <- expand.grid("duration" = (c(4, 2)),
                         "success" = c(100, 60),
                         "rate" = c(0.5, 1, 4))

scenarios$scenario <- 1:12

# we'll subset each track to these scenarios within the function to save workspace memory

#_______________________________________________________________________
# 4. Write functions ----
#_______________________________________________________________________
# 4a. create_track ----

# this one will take a xy df and return a track for subsetting

#_______________________________________________________________________

create_track <- function (df) {
  
  df.track <- df %>% 
    
    mutate(t_ = ymd_hms(t_)) %>%
    
    make_track(.x = x_,
               .y = y_,
               .t = t_,
               crs = "EPSG:32611")
  
  return(df.track)
  
}

#_______________________________________________________________________
# 4b. prep_telem ----

# this one will prepare a correctly-formatted telemetry object for ctmm
# AFTER the track is processed by scenario (i.e., reduced in quality)

#_______________________________________________________________________

prep_telem <- function (track) {
  
  # coerce to telemetry 
  # convert to lat long so as.telemetry doesn't mess it up
  coords <- track %>%
    
    st_as_sf(coords = c("x_", 
                        "y_"),
             crs = "EPSG:32611") %>%
    
    st_transform(crs = "EPSG:4326") %>%
    
    st_coordinates()
  
  # add columns
  track$location.long <- coords[ , 1]
  track$location.lat <- coords[ , 2]
  
  # create telemetry object
  indiv.telem <- xpectr::suppress_mw(
    
    as.telemetry(object = data.frame("timestamp" = track$t_,
                                     "location.lat" = track$location.lat,
                                     "location.long" = track$location.long),
                timeformat = "auto",
                keep = FALSE)
    
    )
  
  # return
  return(indiv.telem)
  
}

#_______________________________________________________________________
# 4c. ctmm_speeds ----
#_______________________________________________________________________

# this function will:
# (1) iterate over all 500 individuals
# (2) iterate over each scenario, within individual
# (3) subset the track to the correct scenario
# (4) select the best performing CTMM
# (5) extract the stationary speed estimate (closely follows the simulated speeds and is much faster)

#_______________________________________________________________________

ctmm_speeds <- function (df) {
  
  # start time
  start.time <- Sys.time()
  
  # iterate over individual
  all.speeds <- data.frame()
  
  for (i in 1:length(unique(df$indiv))) {
    
    focal.indiv <- df %>% filter(indiv == unique(df$indiv)[i])
    
    # iterate over scenarios
    focal.i <- data.frame()
    
    for (j in 1:nrow(scenarios)) {
      
      focal.scenario <- scenarios %>% slice(j)
      
      # create a track
      focal.track <- create_track(focal.indiv)
      
      # TRACK DURATION
      if (focal.scenario$duration == 2) {
        
        focal.track <- focal.track %>% slice(1:(n() / 2))
        
      }
      
      # FIX RATE
      if (focal.scenario$rate %in% c(1, 4)) {
        
      focal.track <- focal.track %>%
        
        track_resample(period = hours(focal.scenario$rate),
                       tolerance = minutes(0))  
        
      }
      
      # FIX SUCCESS
      if (focal.scenario$success == 60) {
        
        focal.track <- focal.track %>%
          
          slice(sample(1:nrow(focal.track), size = round(nrow(focal.track) * 0.6))) %>%
          
          arrange(t_)
        
      }
      
      # create telemetry object
      focal.telem <- prep_telem(focal.track)
      
      # guesstimate parameters
      guess.param <- variogram.fit(variogram(focal.telem), 
                                   name = "guess.param", 
                                   interactive = FALSE)
      
      # model selection
      fitted.mods <- ctmm.select(focal.telem,
                                 CTMM = guess.param,
                                 verbose = TRUE)
      
      # speed estimation (stationary from a correlated velo model)
      # this is only valid if the top model contains Tau2
      if (names(fitted.mods)[1] %in% c("OUF",
                                       "OUF anisotropic",
                                       "OUf",
                                       "OUf anisotropic",
                                       "IOU")) {
        
        # mean speed of Gaussian movement process
        speed.mean <- ctmm::speed(fitted.mods[[1]], 
                                  units = FALSE)
        
        # store in df
        focal.i.scenario <- data.frame(indiv = i,
                                       scenario = j,
                                       top.model = names(fitted.mods)[1],
                                       speed.mean.lo = speed.mean$CI[1],
                                       speed.mean.est = speed.mean$CI[2],
                                       speed.mean.hi = speed.mean$CI[3]) 
        
      } else {
        
        # store in df
        focal.i.scenario <- data.frame(indiv = i,
                                       scenario = j,
                                       top.model = names(fitted.mods)[1],
                                       speed.mean.lo = NA,
                                       speed.mean.est = NA,
                                       speed.mean.hi = NA) 
        
      }
      
      # SLD speed estimation
      focal.i.scenario$sld.speed.mean.lo <- mean(speed(focal.track), na.rm = T) - (sd(speed(focal.track), na.rm = T) / sqrt(nrow(focal.track)))
      focal.i.scenario$sld.speed.mean.est <- mean(speed(focal.track), na.rm = T)
      focal.i.scenario$sld.speed.mean.hi <- mean(speed(focal.track), na.rm = T) + (sd(speed(focal.track), na.rm = T) / sqrt(nrow(focal.track)))
      
      # bind into df
      focal.i <- rbind(focal.i, focal.i.scenario)
      
    }
    
    # bind into df
    all.speeds <- rbind(all.speeds, focal.i)
    
    # print status message after every 10 iterationd is completed
    # AND write to disk to save what's left of my sanity
    if (i %% 1 == 0) {
      
      elapsed.time <- round(as.numeric(difftime(Sys.time(), 
                                                start.time, 
                                                units = "mins")), 
                            digits = 1)
      
      print(paste0("Completed individual ", 
                   i, 
                   " of ", 
                   length(unique(df$indiv)), 
                   " - ", 
                   elapsed.time, 
                   " mins"))
      
      # write .csv
      write.csv(all.speeds, file = paste0(getwd(), "/Derived data/Sampled - CTMM speeds/ctmm_speeds_1T.csv"))
      
  }
  
  }
  
  return(all.speeds)
  
}

#_______________________________________________________________________
# 5. Use function ----

# benchmarking
#ctmm.speeds.test <- ctmm_speeds(df = all.tracks %>% filter(indiv == 1))
# takes 11.2 minutes for one individual
# this equates to ~ 93.3 hours (~ 3.9 days)

#_______________________________________________________________________

ctmm.speeds <- ctmm_speeds(df = all.tracks)

#_______________________________________________________________________
# 6. Write to .csv ----

# (even though it should be written already)

#_______________________________________________________________________

write.csv(ctmm.speeds, file = paste0(getwd(), "/Derived data/Sampled - CTMM speeds/ctmm_speeds_1T.csv"))
