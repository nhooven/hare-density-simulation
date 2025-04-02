# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 06a - Speed estimation (fold 1)
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 26 Mar 2025
# Date completed: 
# Date last modified: 01 Apr 2025
# R version: 4.4.3

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation
library(ctmm)
library(amt)                  # work with tracks

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

tracks.lo1 <- read.csv(paste0(getwd(), "/Derived data/Simulated tracks/tracks_lo1.csv"))
tracks.col.lo1 <- read.csv(paste0(getwd(), "/Derived data/Simulated tracks/tracks_col_lo1.csv"))

tracks.hi1 <- read.csv(paste0(getwd(), "/Derived data/Simulated tracks/tracks_hi1.csv"))
tracks.col.hi1 <- read.csv(paste0(getwd(), "/Derived data/Simulated tracks/tracks_col_hi1.csv"))

#_______________________________________________________________________
# 3. Extract collared tracks only ----
#_______________________________________________________________________
# 3a. Q1 ----
#_______________________________________________________________________

tracks.lo1.1 <- tracks.lo1 %>% filter(indiv %in% 1:5)

#_______________________________________________________________________
# 3b. Q2 ----
#_______________________________________________________________________

tracks.hi1.1 <- tracks.hi1 %>% filter(indiv %in% 1:5)

#_______________________________________________________________________
# 4. Rarefy tracks for each scenario ----

# define scenarios
# these range from ideal (hourly fixes, 100% success, full month tracking)
# to poor (4-hour fixes, 60% success, 2 weeks tracking)

scenarios <- expand.grid(fix.rate = c(1, 2, 4),
                         fix.success = c(100, 80),      # assume that 60% fix success is too bad
                         track.dur = c(4, 2))           # assume that 1 weekk is far too few data

# in this step, we will loop through:
  # iteration (j)
  # scenario (k)
  # individual (l)

scenario_rarefy <- function (
    
  cam.df,   # 5 individuals from the camera contact dataset
  col.df    # 10 individuals from "out of sample"
  
  ) {
  
  # bind together
  both.df <- rbind(cam.df, col.df)
  
  # create unique identifier
  both.df$ID <- paste0(both.df$use, "_", both.df$indiv)
  
  # df to bind to
  all.tracks <- data.frame()
  
  # NESTED LOOPS
  # loop 1 - iteration (j)
  for (j in 1:max(both.df$iter)) {
    
    df.j <- both.df %>% filter(iter == j)
    
    # loop 2 - scenario (k)
    for (k in 1:nrow(scenarios)) {
      
      row.k <- scenarios %>% slice(k)
      
      # loop 3 - individual (l)
      for (l in 1:length(unique(df.j$ID))) {
        
        df.l <- df.j %>% filter(ID == unique(df.j$ID)[l])
        
        # NESTED IF-ELSE
        # track duration
        if (row.k$track.dur == 2) {
          
          df.l.1 <- df.l %>% slice(1:(n() / 2))
          
        } else {
          
          df.l.1 <- df.l
          
          }
        
        # fix success
        if (row.k$fix.success == 80) {
          
          df.l.2 <- df.l.1 %>% slice(sample(1:nrow(df.l.1), size = nrow(df.l.1) * 0.80))
          
        } else {
            
            df.l.2 <- df.l.1
            
          }
          
        # fix rate
        # create a track
        track.l <- df.l.2 %>%
          
          mutate(t_ = ymd_hms(t_)) %>%
          
          make_track(.x = x_,
                     .y = y_,
                     .t = t_)
        
        if (row.k$fix.rate == 2) {
          
          track.l.1 <- track.l %>% track_resample(rate = hours(2), tolerance = minutes(0))
          
        } else {
          
          if (row.k$fix.rate == 4) {
            
            track.l.1 <- track.l %>% track_resample(rate = hours(4), tolerance = minutes(0))
            
          } else {
            
            track.l.1 <- track.l
            
          }
          
        }
        
        # add in identifiers and bind in
        track.l.2 <- track.l.1 %>%
          
          # drop "burst"
          dplyr::select(x_, y_, t_) %>%
          
          mutate(iter = df.l$iter[1],
                 ID = df.l$ID[1],
                 scenario = k)
      
        all.tracks <- rbind(all.tracks, track.l.2)
        
      }
      
    }
    
  }
  
  # return
  return(all.tracks)
  
}

#_______________________________________________________________________
# 4b. Use function ----
#_______________________________________________________________________

q1.rare <- scenario_rarefy(tracks.lo1.1, tracks.col.lo1)

#_______________________________________________________________________
# 5. Select movement models and simulate speeds ----

# this will take FOREVER
# the hope is that as the data quality get worse, this will go quicker
# I'll make sure to tally which model was the top one 

#_______________________________________________________________________
# 5a. Q1 ----

# here we will be conducting model selection for 100 x 15 x 12 things
# AND extracting speeds for a subset of those

# all combinations from q1.rare
q1.rare.combo <- expand.grid(iter = 1:max(q1.rare$iter),
                             ID = unique(q1.rare$ID),
                             scenario = 1:max(q1.rare$scenario))

#_______________________________________________________________________

# list of top models (we can always extract parameters from these later)
top.models <- list()

# df of speeds
all.speeds <- data.frame()

# start.time
start.time <- Sys.time()

for (i in 1:nrow(q1.rare.combo)) {
  
  # focal row
  focal.i <- q1.rare.combo %>% slice(i)
  
  # subset rarefied tracks
  focal.track <- q1.rare %>%
    
    filter(iter == focal.i$iter,
           ID == focal.i$ID,
           scenario == focal.i$scenario)
  
  # coerce to telemetry 
  # convert to lat long so as.telemetry doesn't mess it up
  focal.coords <- focal.track %>%
    
    st_as_sf(coords = c("x_", "y_"),
             crs = "EPSG:32611") %>%
    
    st_transform(crs = "EPSG:4326") %>%
    
    st_coordinates()
  
  # add columns
  focal.track$location.long <- focal.coords[ , 1]
  focal.track$location.lat <- focal.coords[ , 2]
  
  # create telemetry object
  indiv.telem <- as.telemetry(object = data.frame("timestamp" = focal.track$t_,
                                                  "location.lat" = focal.track$location.lat,
                                                  "location.long" = focal.track$location.long),
                              timeformat = "auto",
                              keep = FALSE)
  
  # guesstimate parameters
  guess.param <- variogram.fit(variogram(indiv.telem), 
                               name = "guess.param", 
                               interactive = FALSE)
  
  # model selection
  fitted.mods <- ctmm.select(indiv.telem,
                             CTMM = guess.param,
                             verbose = TRUE)
  
  # bind top model into list
  top.models[[i]] <- fitted.mods[[1]]
  
  # Speed estimation
  # this is only valid if the top model contains Tau2
  if (names(fitted.mods)[1] %in% c("OUF",
                                   "OUF anisotropic",
                                   "OUf",
                                   "OUf anisotropic",
                                   "IOU")) {
    
    # mean speed of Gaussian movement process
    speed.mean <- ctmm::speed(fitted.mods[[1]], 
                              units = FALSE)
    
    # speed conditional on relocations (takes ~ 25 secs on scenario 1)
    #speed.cond <- ctmm::speed(fitted.mods[[1]],        
    #                          data = indiv.telem,    # condition on data
    #                          units = FALSE,         # keep in m/s
    #                          fast = TRUE,           # central limit theorem, prohibitively slow otherwise
    #                          robust = TRUE,         
    #                          trace = FALSE)          # progress bar
    
    # store in focal.i
    focal.i.1 <- focal.i %>%
      
      mutate(top.model = names(fitted.mods)[1],
             speed.mean.lo = speed.mean$CI[1],
             speed.mean.est = speed.mean$CI[2],
             speed.mean.hi = speed.mean$CI[3])
    #,
    #         speed.cond.lo = speed.cond$CI[1],
    #         speed.cond.est = speed.cond$CI[2],
    #         speed.cond.hi = speed.cond$CI[3])
    
  } else {
    
    # store in focal.i
    focal.i.1 <- focal.i %>%
      
      mutate(top.model = names(fitted.mods)[1],
             speed.mean.lo = NA,
             speed.mean.est = NA,
             speed.mean.hi = NA)
    
    #,
    #         speed.cond.lo = NA,
    #         speed.cond.est = NA,
    #         speed.cond.hi = NA)
    
  }
  
  # bind into df
  all.speeds <- rbind(all.speeds, focal.i.1)
  
  # print status message after every iteration is completed
  # AND write to disk to save what's left of my sanity
  if (i %% 180 == 0) {
    
    elapsed.time <- round(as.numeric(difftime(Sys.time(), 
                                              start.time, 
                                              units = "mins")), 
                          digits = 1)
    
    print(paste0("Completed iter ", 
                 i / 180, 
                 " of ", 
                 max(q1.rare.combo$iter), 
                 " - ", 
                 elapsed.time, 
                 " mins"))
    
    # write .csv
    write.csv(all.speeds, file = paste0(getwd(), "/Derived data/CTMM speeds/speeds_lo1.csv"))
    
    # save models
    save(top.models, file = paste0(getwd(), "/Derived data/CTMMs/top_models_lo1.RData"))
    
  }

}

#_______________________________________________________________________
# 5b. Q2 ----

# here we will be conducting model selection for 100 x 15 things
# AND extracting speeds for a subset of those

# bind df together
all.tracks.hi1 <- rbind(tracks.hi1.1, tracks.col.hi1)

# add unique identifier
all.tracks.hi1$ID <- paste0(all.tracks.hi1$indiv, "_", all.tracks.hi1$use)

#_______________________________________________________________________

# list of top models (we can always extract parameters from these later)
top.models.q2 <- list()

# df of speeds
all.speeds.q2 <- data.frame()

# start.time
start.time <- Sys.time()

for (i in 1:length(unique(all.tracks.hi1$ID))) {
  
  # focal track
  focal.track <- all.tracks.hi1 %>%
    
    filter(ID == unique(all.tracks.hi1$ID)[i])
  
  # focal row
  focal.i <- data.frame(iter = focal.track$iter[1],
                        ID = focal.track$ID[1])
  
  # coerce to telemetry 
  # convert to lat long so as.telemetry doesn't mess it up
  focal.coords <- focal.track %>%
    
    st_as_sf(coords = c("x_", "y_"),
             crs = "EPSG:32611") %>%
    
    st_transform(crs = "EPSG:4326") %>%
    
    st_coordinates()
  
  # add columns
  focal.track$location.long <- focal.coords[ , 1]
  focal.track$location.lat <- focal.coords[ , 2]
  
  # create telemetry object
  indiv.telem <- as.telemetry(object = data.frame("timestamp" = focal.track$t_,
                                                  "location.lat" = focal.track$location.lat,
                                                  "location.long" = focal.track$location.long),
                              timeformat = "auto",
                              keep = FALSE)
  
  # guesstimate parameters
  guess.param <- variogram.fit(variogram(indiv.telem), 
                               name = "guess.param", 
                               interactive = FALSE)
  
  # model selection
  fitted.mods <- ctmm.select(indiv.telem,
                             CTMM = guess.param,
                             verbose = TRUE)
  
  # bind top model into list
  top.models.q2[[i]] <- fitted.mods[[1]]
  
  # Speed estimation
  # this is only valid if the top model contains Tau2
  if (names(fitted.mods)[1] %in% c("OUF",
                                   "OUF anisotropic",
                                   "OUf",
                                   "OUf anisotropic",
                                   "IOU")) {
    
    # mean speed of Gaussian movement process
    speed.mean <- ctmm::speed(fitted.mods[[1]], 
                              units = FALSE)
    
    # speed conditional on relocations (takes ~ 25 secs on scenario 1)
    #speed.cond <- ctmm::speed(fitted.mods[[1]],        
    #                          data = indiv.telem,    # condition on data
    #                          units = FALSE,         # keep in m/s
    #                          fast = TRUE,           # central limit theorem, prohibitively slow otherwise
    #                          robust = TRUE,         
    #                          trace = FALSE)          # progress bar
    
    # store in focal.i
    focal.i.1 <- focal.i %>%
      
      mutate(top.model = names(fitted.mods)[1],
             speed.mean.lo = speed.mean$CI[1],
             speed.mean.est = speed.mean$CI[2],
             speed.mean.hi = speed.mean$CI[3])
    
    #,
    #         speed.cond.lo = speed.cond$CI[1],
    #         speed.cond.est = speed.cond$CI[2],
    #         speed.cond.hi = speed.cond$CI[3])
    
  } else {
    
    # store in focal.i
    focal.i.1 <- focal.i %>%
      
      mutate(top.model = names(fitted.mods)[1],
             speed.mean.lo = NA,
             speed.mean.est = NA,
             speed.mean.hi = NA)
    
    #,
    #         speed.cond.lo = NA,
    #         speed.cond.est = NA,
    #         speed.cond.hi = NA)
    
  }
  
  # bind into df
  all.speeds.q2 <- rbind(all.speeds.q2, focal.i.1)
  
  # print status message after every iteration is completed
  # AND write to disk to save what's left of my sanity
  if (i %% 15 == 0) {
    
    elapsed.time <- round(as.numeric(difftime(Sys.time(), 
                                              start.time, 
                                              units = "mins")), 
                          digits = 1)
    
    print(paste0("Completed iter ", 
                 i / 15, 
                 " of ", 
                 max(q1.rare.combo$iter), 
                 " - ", 
                 elapsed.time, 
                 " mins"))
    
    # write .csv
    write.csv(all.speeds.q2, file = paste0(getwd(), "/Derived data/CTMM speeds/speeds_hi1.csv"))
    
    # save models
    save(top.models.q2, file = paste0(getwd(), "/Derived data/CTMMs/top_models_hi1.RData"))
    
  }
  
}

#_______________________________________________________________________
# 6. Take a breath  ----
#_______________________________________________________________________

# it's over!