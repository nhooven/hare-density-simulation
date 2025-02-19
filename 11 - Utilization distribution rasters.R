# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 11 - Utilization distribution rasters
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 12 Feb 2025
# Date completed: 12 Feb 2025
# Date last modified: 18 Feb 2025
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)          
library(lubridate)          # work with dates
library(terra)              # rasters
library(sf)                 # polygons
library(amt)                # tracks
library(ctmm)               # ctmms/aKDEs

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

# read in a dataset
sims.B1 <- read.csv(paste0(getwd(), "/UD sims/sims_B1.csv"))
sims.B2 <- read.csv(paste0(getwd(), "/UD sims/sims_B2.csv"))
sims.B3 <- read.csv(paste0(getwd(), "/UD sims/sims_B3.csv"))
sims.A1 <- read.csv(paste0(getwd(), "/UD sims/sims_A1.csv"))
sims.A2 <- read.csv(paste0(getwd(), "/UD sims/sims_A2.csv"))
sims.A3 <- read.csv(paste0(getwd(), "/UD sims/sims_A3.csv"))

# unit boundary
unit.bound <- st_read(paste0(getwd(), "/Derived_data/Shapefiles/unit_bound.shp"))

# define a raster grid on which to calculate the aKDE
landscape.covs.B1 <- rast("Rasters/Scaled covariates/B1.tif")

base.grid <- raster::raster(landscape.covs.B1$forage)

#_______________________________________________________________________
# 2. Define function ----
#_______________________________________________________________________

# function
akde_count <- function(sims) {
  
  # here we want to calculate expected use at any one time step, per pixel
  # we will take one replicate from EACH individual to do this calculation
  # this will give us a raster of expected counts for each replicate 
  # we can then calculate a pixel-wise mean and SD, invert to create a correction factor, etc.
  
  # loop through replicates m (sim.rep from the sims df)
  # define base raster
  all.rast.m <- crop(rast(nrows = nrow(base.grid),
                          ncols = ncol(base.grid),
                          resolution = res(base.grid),
                          extent = raster::extent(base.grid),
                          vals = NA,
                          crs = "EPSG:32611"),
                     unit.bound)
  
  # loop
  for (m in 1:max(sims$sim.rep)) {
    
    # subset
    focal.sims <- sims %>%
      
      filter(sim.rep == m) %>%
      
      # format time correctly (anything sampled at 00:00 doesn't have a time component)
      mutate(t_ = case_when(nchar(t_) == 19 ~  t_,
                            nchar(t_) == 10 ~ paste0(t_, " 00:00:00"))) %>%
      
      mutate(t_ = ymd_hms(t_))
    
    # now we have 100 individuals n, with a focal m
    # we'll need to loop through these individuals to fit aKDEs
    
    # define base raster
    all.rast.n <- rast(nrows = nrow(base.grid),
                       ncols = ncol(base.grid),
                       resolution = res(base.grid),
                       extent = extent(base.grid),
                       vals = NA,
                       crs = "EPSG:32611")
    
    # loop
    for (n in 1:max(focal.sims$indiv)) {
      
      # subset individual and coerce to track
      indiv.sim <- focal.sims %>%
        
        filter(indiv == n) %>%
        
        make_track(.x = x_, 
                   .y = y_, 
                   .t = t_, 
                   crs = 32611)
      
      # convert to telemetry object
      indiv.telem <- as_telemetry(indiv.sim,
                                  timeformat = "auto",
                                  timezone = "America/Los_Angeles",
                                  keep = TRUE)
      
      # change projection
      projection(indiv.telem) <- "EPSG:32611"
      
      # fit CTSP models
      # guesstimated model parameters from the variogram
      guess.param <- variogram.fit(variogram(indiv.telem), 
                                   name = "guess.param", 
                                   interactive = FALSE)
      
      # baseline model (we'll use the Ornstein-Uhlenbeck process)
      ctmm.model.1 <- ctmm(tau = guess.param$tau[1],
                           omega = FALSE,
                           range = TRUE,
                           error = FALSE,
                           isotropic = TRUE,
                           data = indiv.telem)
      
      # fit aKDE
      focal.akde <- akde(data = indiv.telem,
                         CTMM = ctmm.fit(data = indiv.telem,
                                         CTMM = ctmm.model.1),
                         grid = base.grid)
      
      # convert to SpatRaster and normalize
      focal.akde.pdf <- rast(base.grid)
      values(focal.akde.pdf) <- focal.akde$PDF / sum(focal.akde$PDF)  # must sum to 1
      
      # name layer (unique individual)
      names(focal.akde.pdf) <- n
      
      # bind into raster stack
      all.rast.n <- c(all.rast.n, focal.akde.pdf)
      
    }
    
    # sum all rasters
    all.rast.n.sum <- sum(all.rast.n, na.rm = TRUE)
    
    # crop
    all.rast.n.sum.crop <- crop(all.rast.n.sum, unit.bound)
    
    # add name of replicate
    names(all.rast.n.sum.crop) <- m
    
    # and add into a stack of length m
    all.rast.m <- c(all.rast.m, all.rast.n.sum.crop)
    
  }
  
  # return
  return(all.rast.m)
  
}

#_______________________________________________________________________
# 3. Use function ----
#_______________________________________________________________________

 <- akde_count()

# completed these in multiple sessions - completed 02-18-2025
