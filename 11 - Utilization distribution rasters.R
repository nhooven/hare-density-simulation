# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 11 - Utilization distribution rasters
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 12 Feb 2025
# Date completed: 12 Feb 2025
# Date last modified: 26 Feb 2025
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)          
library(lubridate)          # work with dates
library(terra)              # rasters
library(sf)                 # polygons

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
templ.rast<- rast("Rasters/Scaled covariates/B1.tif")

#_______________________________________________________________________
# 2. Define function ----

# here we want to:
# (1) iterate through all 100 simulation replicates
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
  base.rast <- rast(crop(templ.rast$fora, unit.bound),
                    vals = NA)

  for (i in 1:max(sims$sim.rep)) {
    
    # subset 
    sims.focal <- sims %>%
      
      filter(sim.rep == i)
    
    # create sf points
    sims.focal.sf <- st_as_sf(sims.focal,
                              coords = c("x_", "y_"),
                              crs = "EPSG:32611")
    
    # rasterize the counts
    focal.rast <- rasterize(sims.focal.sf,
                            templ.rast,
                            fun = "count")
    
    # crop
    focal.rast.crop <- crop(focal.rast, unit.bound)
    
    # change layer name
    names(focal.rast.crop) <- i
    
    # "normalize" to calculate the RSS
    focal.rast.rss <- focal.rast.crop / mean(values(focal.rast.crop))
    
    # and calculate the correction factor
    focal.rast.cf <- 1 / focal.rast.rss
    
    # bind into base raster
    base.rast <- c(base.rast, focal.rast.cf)
    
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
