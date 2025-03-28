# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 04a - CTMM simulations (camera - low)
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 25 Mar 2025
# Date completed: 26 Mar 2025
# Date last modified: 28 Mar 2025
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

indivs <- read.csv(paste0(getwd(), "/Derived data/Individuals and parameters/camera_lo.csv"))

# unit boundary
unit.bound <- st_read(paste0(getwd(), "/Derived data/Shapefiles/unit_bound.shp"))

# camera viewsheds
vs <- st_read(paste0(getwd(), "/Derived data/Shapefiles/cams_vs.shp"))

vs$cam.id <- 1:nrow(vs)

st_crs(vs) <- "epsg:32611"

#_______________________________________________________________________
# 3. Simulation parameters ----
#_______________________________________________________________________

# simulation duration
sim.duration.wk <- 8

sim.duration.sec <- sim.duration.wk * 7 * 24 * 60 * 60 

# "sampling rate" for fundamental (i.e., straight line) steps
sim.samp.rate <- 20

# time steps to simulate on
sim.timestep <- seq(1, sim.duration.sec, sim.samp.rate)

#_______________________________________________________________________
# 4. Function ----
#_______________________________________________________________________

ctmm_sim <- function(data = indivs) {
  
  # loop through all rows
  all.passes <- data.frame()
  all.relocs <- data.frame()
  all.speeds <- data.frame()
  
  # start time
  start.time <- Sys.time()
  
  for (i in 1:nrow(data)) {
    
    # subset focal row
    focal.row <- data %>% slice(i)
    
    # define OUF model
    ouf.mod <- ctmm(tau = c(focal.row$tau1,             # positional autocorr time
                            focal.row$tau2),            # velocity autocorr time
                    isotropic = FALSE,                  # anisotropic
                    sigma = c(focal.row$sigma.maj,      # variance along major axis
                              focal.row$sigma.min,      # variance along minor axis
                              focal.row$angle),         # angle of axis
                    mu = c(focal.row$mean1,             # x coord
                           focal.row$mean2))            # y coord
    
    # run simulation
    ouf.sim <- simulate(object = ouf.mod, 
                        t = sim.timestep,
                        complete = T)
    
    # coerce telemetry object to sf
    telem.df <- data.frame(t = ouf.sim$t,
                           x = ouf.sim$x,
                           y = ouf.sim$y,
                           timestamp = as.POSIXct(ouf.sim$timestamp))
    
    telem.sf <- st_as_sf(telem.df,
                         coords = c("x", "y"),
                         crs = "epsg:32611")
    
    # tally camera contacts ("passes")
    # in this case, 20-second locations that fall within each viewshed
    # these will be independent detections (> 20 minutes apart)
    passes <- st_intersection(vs, 
                              telem.sf) %>%
      
      # drop the geometry
      st_drop_geometry() %>% 
      
      # group by camera
      group_by(cam.id) %>% 
      
      # add a time difference column
      mutate(time.diff = as.numeric(t - lag(t))) %>%
      
      # drop anything < 20 minutes
      filter(time.diff > 20 %#% "minutes") %>%
      
      # tally all intersections
      tally() %>%
      
      # and add identifiers
      mutate(indiv = focal.row$indiv,
             rep = focal.row$rep,
             use = focal.row$use,
             abund.group = focal.row$abund.group,
             collared = focal.row$collared)
    
    # bind in to passes df
    all.passes <- rbind(all.passes, passes)
    
    # sample to a "GPS track"
    # this will be a regular 1-hr track, which can be subsampled for Q2
    telem.track <- telem.df %>%
      
      make_track(.x = x,
                 .y = y,
                 .t = timestamp)
    
    # calculate "true" speed (m/s)
    telem.track.speed <- telem.track %>%
      
      steps() 
    
    # 8 week
    speed.8wk <- sum(telem.track.speed$sl_) / sim.duration.sec
    
    # 4 week
    telem.track.speed.4 <- telem.track.speed %>% slice(1:(n() / 2))
    
    speed.4wk <- sum(telem.track.speed.4$sl_) / (sim.duration.sec / 2)
    
    # resample track for telemetering
    telem.track.1 <- telem.track %>%
      
      track_resample(rate = hours(1),
                     tolerance = minutes(0)) %>%
      
      # remove "burst"
      dplyr::select(-burst_) %>%
      
      # and add identifiers
      mutate(indiv = focal.row$indiv,
             rep = focal.row$rep,
             use = focal.row$use,
             abund.group = focal.row$abund.group,
             collared = focal.row$collared)
    
    # bind in to relocations df
    all.relocs <- rbind(all.relocs, telem.track.1)
    
    # bind into speeds df
    all.speeds <- rbind(all.speeds,
                        data.frame(indiv = focal.row$indiv,
                                   rep = focal.row$rep,
                                   use = focal.row$use,
                                   abund.group = focal.row$abund.group,
                                   collared = focal.row$collared,
                                   speed.8wk = speed.8wk,
                                   speed.4wk = speed.4wk))
    
    # status message (every 10 iterations)
    if (i %% 10 == 0) {
      
      elapsed.time <- round(as.numeric(difftime(Sys.time(), 
                                                start.time, 
                                                units = "mins")), 
                            digits = 1)
      
      print(paste0("Completed sim ", i, " of ", nrow(data), " - ", elapsed.time, " mins"))
      
    }
    
  }
  
  # pack into list and return
  passes.relocs.speeds <- list(all.passes,
                               all.relocs,
                               all.speeds)
  
  return(passes.relocs.speeds)
  
}

#_______________________________________________________________________
# 5. Use function ----
#_______________________________________________________________________

sim.all <- ctmm_sim(data = indivs)

#_______________________________________________________________________
# 6. Write to .csv ----
#_______________________________________________________________________

write.csv(sim.all[[1]], paste0(getwd(), "/Derived data/Camera detections/detections_lo.csv"))

write.csv(sim.all[[2]], paste0(getwd(), "/Derived data/Simulated tracks/tracks_lo.csv"))

write.csv(sim.all[[3]], paste0(getwd(), "/Derived data/Speeds/speeds_lo.csv"))
