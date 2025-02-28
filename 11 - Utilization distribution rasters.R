# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 11 - Utilization distribution rasters
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 12 Feb 2025
# Date completed: 12 Feb 2025
# Date last modified: 28 Feb 2025
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)          
library(lubridate)          # work with dates
library(terra)              # rasters
library(sf)                 # polygons
library(amt)                # tracks
library(adehabitatHR)               

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

# define a raster grid on which to calculate the KDE (should be a SpatialPixels)
templ.rast <- rast("Rasters/Scaled covariates/B1.tif")
templ.rast.sp <- as(raster::raster(templ.rast), "SpatialPixels")

#_______________________________________________________________________
# 2. Define function ----

# here we want to:
# (1) iterate through all 100 simulation replicates
# (2) fit simple KDEs for all individuals
# (3) sum KDEs within each replicate
# (2) tally the number of locations per pixel 
#         - of possible (n.indiv * 337)
# (3) calculate mean and SD over all replicates

#_______________________________________________________________________

sim_count <- function(sims) {
  
  # extract indices
  id.trt <- sims$trt[1]
  id.rep <- sims$rep[1]
  
  # loop through all replicates
  # base raster to bind into
  base.rast <- rast(crop(templ.rast$fora, 
                         unit.bound),
                    vals = NA)
  
  start.time <- Sys.time()
  
  for (i in 1:max(sims$sim.rep)) {
    
    # subset 
    sims.focal <- sims %>%
      
      filter(sim.rep == i)
    
    # loop through individuals
    base.rast.rep <- rast(crop(templ.rast$fora, 
                               unit.bound),
                          vals = NA)
    
    for (j in 1:max(sims.focal$indiv)) {
      
      # subset 
      sims.focal.indiv <- sims.focal %>%
        
        filter(indiv == j)
      
      # create sf points
      sims.focal.sf <- st_as_sf(sims.focal.indiv,
                                coords = c("x_", "y_"),
                                crs = "EPSG:32611")
      
     # convert sf to SpatialPoints objects
     sims.sp <- SpatialPoints(as(sims.focal.sf, "Spatial"))
     sims.sp@proj4string <- CRS("EPSG:32611") 
     
     # fit kernel with reference h parameter (way faster than LSCV)
     mean.kernel <- kernelUD(xy = sims.sp, h = "href", grid = templ.rast.sp)
     
     # convert to SpatRaster and crop to unit boundary
     mean.kernel.rast <- crop(rast(mean.kernel), unit.bound)
     
     # change name to individual
     names(mean.kernel.rast) <- j
     
     # bind in
     base.rast.rep <- c(base.rast.rep, mean.kernel.rast)
      
    }
    
    # calculate the sum
    sum.rast.rep <- sum(base.rast.rep, na.rm = TRUE)
    
    # "normalize" to calculate the RSS
    rss.rast.rep <- sum.rast.rep / mean(values(sum.rast.rep))
    
    # and calculate the correction factor
    cf.rast.rep <- 1 / rss.rast.rep
    
    # change layer name to replicate
    names(cf.rast.rep) <- i
    
    # bind into base raster
    base.rast <- c(base.rast, cf.rast.rep)
    
    # status message (every 10 replicates)
    if (i %% 10 == 0) {
      
      elapsed.time <- round(as.numeric(difftime(Sys.time(), 
                                                start.time, 
                                                units = "mins")), 
                            digits = 1)
      
      print(paste0("Completed rep ", i, " of ", max(sims$sim.rep), " - ", elapsed.time, " mins"))
      
    }
    
  }
  
  # remove first layer
  total.rast <- base.rast[[-1]]
  
  # calculate mean and SD
  mean.rast <- c(mean(total.rast, na.rm = T), stdev(total.rast, na.rm = T))
  
  # return
  return(mean.rast)
  
}

#_______________________________________________________________________
# 3. Use function ----
#_______________________________________________________________________

sim.count.B1 <- sim_count(sims.B1)
sim.count.B2 <- sim_count(sims.B2)
sim.count.B3 <- sim_count(sims.B3)
sim.count.A1 <- sim_count(sims.A1)
sim.count.A2 <- sim_count(sims.A2)
sim.count.A3 <- sim_count(sims.A3)

#_______________________________________________________________________
# 4. Examine rasters ----
#_______________________________________________________________________

plot(sim.count.B1)
plot(sim.count.B2)
plot(sim.count.B3)
plot(sim.count.A1)
plot(sim.count.A2)
plot(sim.count.A3)

#_______________________________________________________________________
# 5. Write rasters ----
#_______________________________________________________________________

writeRaster(sim.count.B1, paste0(getwd(), "/Rasters/UD predictions/B1.tif"))
