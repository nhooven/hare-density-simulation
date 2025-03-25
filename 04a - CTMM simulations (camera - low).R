# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 04a - CTMM simulations (camera - low)
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 25 Mar 2025
# Date completed: 
# Date last modified: 
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
  for (i in 1:nrow(data)) {
    
    # subset focal row
    focal.row <- data %>% slice(i)
    
    # define OUF model
    ouf.mod <- ctmm(tau = c(focal.row$tau1,
                            focal.row$tau2), 
                    isotropic = FALSE, 
                    sigma = c(focal.row$sigma.maj,
                              focal.row$sigma.min,
                              focal.row$angle), 
                    mu = c(focal.row$mean1, 
                           focal.row$mean2))
    
    # run simulation
    ouf.sim <- simulate(object = ouf.mod, 
                        t = sim.timestep)
    
    # tally camera contacts
    # coerce telemetry object to sf
    telem.df <- data.frame(t = ouf.sim$t,
                           x = ouf.sim$x,
                           y = ouf.sim$y)
    
    telem.sf <- st_as_sf(telem.df,
                         coords = c("x", "y"),
                         crs = "epsg:32611")
    
    # tally "passes"
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
      tally()
  
    
  }
  
}