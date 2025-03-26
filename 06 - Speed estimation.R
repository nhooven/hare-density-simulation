# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 06 - Speed estimation
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 26 Mar 2025
# Date completed: 
# Date last modified: 26 Mar 2025
# R version: 4.4.3

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation
library(ctmm)

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

tracks.lo <- read.csv(paste0(getwd(), "/Derived data/Simulated tracks/tracks_lo.csv"))
tracks.collar.lo <- read.csv(paste0(getwd(), "/Derived data/Simulated tracks/tracks_collar_lo.csv"))

#_______________________________________________________________________
# 3. Extract collared tracks only ----
#_______________________________________________________________________
# 3a. Q1 ----
#_______________________________________________________________________

tracks.lo.1 <- tracks.lo %>%
  
  filter(collared == 1)

#_______________________________________________________________________
# 4. Fit movement models ----
#_______________________________________________________________________
# 4. Q1 ----

# because our data here are all "high quality", we'll assume that the 
# OUF model will be the best one

# each replicate gets 5 "camera" individuals and 10 additional ones

#_______________________________________________________________________

# define function
q1_ctmm <- function () {
  
  # loop through replicates
  start.time <- Sys.time()
  
  all.params.i <- data.frame()
  all.ctmm.i <- list()
  
  for (i in 1:6) {
    
    # subset "camera" tracks
    tracks.camera.focal <- tracks.lo.1 %>%
      
      filter(rep == i) %>%
      
      # remove extraneous columns
      dplyr::select(-c(abund.group,
                       collared))
    
    # subset "collared" tracks
    tracks.collar.focal <- tracks.collar.lo %>%
      
      filter(rep == i)
    
    # bind together
    tracks.focal <- rbind(tracks.camera.focal, tracks.collar.focal)
    
    # add unique ID
    tracks.focal$ID <- paste0(tracks.focal$indiv, "_", tracks.focal$use)
    
    # loop through indivs
    all.params.j <- data.frame()
    all.ctmm.j <- list()
    
    for (j in 1:length(unique(tracks.focal$ID))) {
      
      focal.indiv <- tracks.focal %>%
        
        filter(ID == unique(tracks.focal$ID)[j]) %>%
      
      # only keep 4 weeks' worth of data for Q1
      slice(1:(4 * 7 * 24))
      
      # convert to lat long so as.telemetry doesn't mess it up
      focal.coords <- focal.indiv %>%
        
        st_as_sf(coords = c("x_", "y_"),
                 crs = "EPSG:32611") %>%
        
        st_transform(crs = "EPSG:4326") %>%
        
        st_coordinates()
      
      # add columns
      focal.indiv$location.long <- focal.coords[ , 1]
      focal.indiv$location.lat <- focal.coords[ , 2]
      
      # create telemetry object
      indiv.telem <- as.telemetry(object = data.frame("timestamp" = focal.indiv$t_,
                                                      "location.lat" = focal.indiv$location.lat,
                                                      "location.long" = focal.indiv$location.long),
                                  timeformat = "auto",
                                  keep = FALSE)
      
      # fit OUF model
      # guesstimated model parameters from the variogram
      guess.param <- variogram.fit(variogram(indiv.telem), 
                                   name = "guess.param", 
                                   interactive = FALSE)
      
      # define model
      indiv.ctmm <- ctmm(tau = guess.param$tau,
                         sigma = guess.param$sigma,
                         isotropic = FALSE,
                         range = TRUE)
      
      # fit model with pHREML
      indiv.ctmm.fit <- ctmm.fit(data = indiv.telem,
                                 CTMM = indiv.ctmm)
      
      # save parameters to df (for comparison later)
      all.params.j <- rbind(all.params.j,
                            data.frame(indiv = focal.indiv$indiv[1],
                                       rep = focal.indiv$rep[1],
                                       use = focal.indiv$use[1],
                                       ID = focal.indiv$ID[1],
                                       mean1 = indiv.ctmm.fit$mu[1],
                                       mean2 = indiv.ctmm.fit$mu[2],
                                       tau1 = indiv.ctmm.fit$tau[1],
                                       tau2 = indiv.ctmm.fit$tau[2],
                                       sigma.maj = indiv.ctmm.fit$sigma@par[1],
                                       sigma.min = indiv.ctmm.fit$sigma@par[2],
                                       angle = indiv.ctmm.fit$sigma@par[3]))
      
      # save fitted CTMM
      all.ctmm.j[[j]] <- indiv.ctmm.fit
      
    }
    
    # bind again
    all.params.i <- rbind(all.params.i, all.params.j)
    
    # and save to list
    all.ctmm.i[[i]] <- all.ctmm.j
    
    # print status message
    elapsed.time <- round(as.numeric(difftime(Sys.time(), 
                                              start.time, 
                                              units = "mins")), 
                          digits = 1)
    
    print(paste0("Completed replicate ", i, " of ", 6, " - ", elapsed.time))
  
  }
  
  # return
  return(list(all.params.i,
              all.ctmm.i))
  
}

# use function
q1.ctmm <- q1_ctmm()


