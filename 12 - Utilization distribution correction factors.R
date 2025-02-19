# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 12 - Utilization distribution correction factors
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 18 Feb 2025
# Date completed: 18 Feb 2025
# Date last modified: 18 Feb 2025
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)          
library(lubridate)          # work with dates
library(terra)              # rasters
library(sf)                 # polygons

#_______________________________________________________________________
# 2. Read in rasters ----
#_______________________________________________________________________

rast.B1 <- rast(paste0(getwd(), "/Rasters/UD rasters/B1.tif"))
rast.B2 <- rast(paste0(getwd(), "/Rasters/UD rasters/B2.tif"))
rast.B3 <- rast(paste0(getwd(), "/Rasters/UD rasters/B3.tif"))
rast.A1 <- rast(paste0(getwd(), "/Rasters/UD rasters/A1.tif"))
rast.A2 <- rast(paste0(getwd(), "/Rasters/UD rasters/A2.tif"))
rast.A3 <- rast(paste0(getwd(), "/Rasters/UD rasters/A3.tif"))

#_______________________________________________________________________
# 3. Calculate mean and SD of correction factor ----
#_______________________________________________________________________
# 3a. Define function ----
#_______________________________________________________________________

calc_cf <- function(x) {
  
  # initialize a new raster
  new.rast <- x$lyr.1
  
  # loop through each replicate, calculating the correction factor
  for (i in 2:nlyr(x)) {
    
    # subset out focal raster
    focal.rast <- x[[i]]
    
    # divide by the mean value
    focal.rast.rss <- focal.rast / mean(values(focal.rast))
    
    # and take the reciprocal of the RSS for a correction factor
    focal.rast.cf <- 1 / focal.rast.rss
    
    # and bind into the new raster
    new.rast <- c(new.rast, focal.rast.cf)
    
  }
  
  # calculate the mean and SD of the correction factor
  cf.rast <- c(mean(new.rast, na.rm = TRUE),
               stdev(new.rast, na.rm = TRUE))
  
  # return
  return(cf.rast)
  
}

#_______________________________________________________________________
# 3b. Use function ----
#_______________________________________________________________________

cf.B1 <- calc_cf(rast.B1)
cf.B2 <- calc_cf(rast.B2)
cf.B3 <- calc_cf(rast.B3)
cf.A1 <- calc_cf(rast.A1)
cf.A2 <- calc_cf(rast.A2)
cf.A3 <- calc_cf(rast.A3)

#_______________________________________________________________________
# 4. Write rasters ----
#_______________________________________________________________________

writeRaster(cf.B1, paste0(getwd(), "/Rasters/UD predictions/B1.tif"), overwrite = TRUE)
writeRaster(cf.B2, paste0(getwd(), "/Rasters/UD predictions/B2.tif"), overwrite = TRUE)
writeRaster(cf.B3, paste0(getwd(), "/Rasters/UD predictions/B3.tif"), overwrite = TRUE)
writeRaster(cf.A1, paste0(getwd(), "/Rasters/UD predictions/A1.tif"), overwrite = TRUE)
writeRaster(cf.A2, paste0(getwd(), "/Rasters/UD predictions/A2.tif"), overwrite = TRUE)
writeRaster(cf.A3, paste0(getwd(), "/Rasters/UD predictions/A3.tif"), overwrite = TRUE)
