# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 11 - Utilization distribution rasters
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 08 Jan 2025
# Date completed: 08 Jan 2025
# Date last modified: 08 Jan 2025
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation
library(sf)              # read in shapefiles
library(terra)           # rasters
library(amt)             # simulate tracks
library(lubridate)       # work with time
library(sp)              # spatial points
library(raster)          # rasters
library(adehabitatHR)    # fit kernels
library(matrixStats)     # rowwise SDs

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

# simulated start and endpoints
sims.S1L <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_UD_S1L.csv"))
sims.S2L <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_UD_S2L.csv"))
sims.S3L <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_UD_S3L.csv"))
sims.S1H <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_UD_S1H.csv"))
sims.S2H <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_UD_S2H.csv"))
sims.S3H <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_UD_S3H.csv"))

sims.C1L <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_UD_C1L.csv"))
sims.C2L <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_UD_C2L.csv"))
sims.C3L <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_UD_C3L.csv"))
sims.C1H <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_UD_C1H.csv"))
sims.C2H <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_UD_C2H.csv"))
sims.C3H <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_UD_C3H.csv"))

# unit boundary
unit.bound <- st_read(paste0(getwd(), "/Derived_data/Shapefiles/unit_bound.shp"))

# read one raster in as a template
landscape.template <- rast("Rasters/Scaled covariates/S1L.tif")

#_______________________________________________________________________
# 3. Keep only the endpoints ----
#_______________________________________________________________________

sims.S1L <- sims.S1L %>% filter(which.point == "end")
sims.S2L <- sims.S2L %>% filter(which.point == "end")
sims.S3L <- sims.S3L %>% filter(which.point == "end")
sims.S1H <- sims.S1H %>% filter(which.point == "end")
sims.S2H <- sims.S2H %>% filter(which.point == "end")
sims.S3H <- sims.S3H %>% filter(which.point == "end")

sims.C1L <- sims.C1L %>% filter(which.point == "end")
sims.C2L <- sims.C2L %>% filter(which.point == "end")
sims.C3L <- sims.C3L %>% filter(which.point == "end")
sims.C1H <- sims.C1H %>% filter(which.point == "end")
sims.C2H <- sims.C2H %>% filter(which.point == "end")
sims.C3H <- sims.C3H %>% filter(which.point == "end")

#_______________________________________________________________________
# 4. Determine template size ----
#_______________________________________________________________________
# 4a. Bind all locations together ----
#_______________________________________________________________________

sims.all <- rbind(sims.S1L, sims.S2L, sims.S3L,
                  sims.S1H, sims.S2H, sims.S3H,
                  sims.C1L, sims.C2L, sims.C3L,
                  sims.C1H, sims.C2H, sims.C3H)

# promote to sf
sims.all.sf <- st_as_sf(sims.all,
                        coords = c("x_", 
                                   "y_"),
                        crs = "epsg:32611")

#_______________________________________________________________________
# 4b. Create bounding box polygon ----
#_______________________________________________________________________

# extract bounding box
bbox <- st_bbox(sims.all.sf)

# add corners to df
bbox.df <- data.frame("x" = c(bbox[1], 
                              bbox[1],
                              bbox[3],
                              bbox[3]),
                      "y" = c(bbox[2], 
                              bbox[4],
                              bbox[4],
                              bbox[2]))

# make a polygon
bbox.poly <- bbox.df %>%
  
  st_as_sf(coords = c("x", "y"),
           crs = "epsg:32611") %>%
  
  summarize(geomety = st_combine(geometry)) %>%
  
  st_cast("POLYGON")

#_______________________________________________________________________
# 4c. Crop template raster ----
#_______________________________________________________________________

# crop
template <- crop(landscape.template$stem, bbox.poly)

# initialize a SpatialPixels object equivalent to the raster
template.sp <- as(as(template, "Raster"), "SpatialPixels")

#_______________________________________________________________________
# 5. Fit kernel UDs ----
#_______________________________________________________________________
# 5a. Define function ----
#_______________________________________________________________________

sim_kernel <- function (id.landscape,
                        id.variability,
                        id.rep,
                        template.sp)     {
  
  # subset data
  sims.focal <- sims.all.sf %>%
    
    filter(landscape == id.landscape,
           variability == id.variability,
           rep == id.rep)
  
  # convert sf to SpatialPoints objects
  sims.sp <- as(sims.focal, "Spatial") 
  
  # fit kernel with reference h parameter (way faster than LSCV)
  mean.kernel <- kernelUD(xy = sims.sp, h = "href", grid = template.sp)
  
  # convert to SpatRaster and crop to unit boundary
  mean.kernel.rast <- crop(rast(mean.kernel), unit.bound)
  
  # convert to "RSS" by scaling the mean to be 1
  mean.kernel.rss <- mean.kernel.rast / mean(values(mean.kernel.rast))
  
  # bootstrap to calculate SE
  # we'll use 500 bootstrapped (replaced) samples
  pred.samples <- matrix(data = NA, nrow = 5000, ncol = 500)
  
  # fill matrix
  for (i in 1:ncol(pred.samples)) {
    
    pred.samples[ , i] <- sample(1:5000, size = nrow(pred.samples), replace = TRUE)
    
  }
  
  # fit kernels for each, calculate RSS, then cast the SD to a raster
  pred.values <- matrix(data = NA, nrow = nrow(values(mean.kernel.rss)), ncol = 500)
  
  for (j in 1:ncol(pred.values)) {
    
    # subset data
    sims.sampled <- sims.focal %>%
      
      slice(pred.samples[ , j])
    
    # convert sf to SpatialPoints object
    sims.sampled.sp <- SpatialPoints(as(sims.sampled, "Spatial"))
    sims.sampled.sp@proj4string <- CRS("+proj=utm +zone=11 +datum=WGS84 +units=m +no_defs +type=crs")
    
    # fit kernel with reference h parameter (way faster than LSCV)
    mean.kernel.sampled <- kernelUD(xy = sims.sampled.sp, h = "href", grid = template.sp)
    
    # convert to SpatRaster and crop to unit boundary
    mean.kernel.sampled.rast <- crop(rast(mean.kernel.sampled), unit.bound)
    
    # convert to "RSS" by scaling the mean to be 1
    mean.kernel.sampled.rss <- mean.kernel.sampled.rast / mean(values(mean.kernel.sampled.rast))
    
    # add the values to the matrix
    pred.values[ , j] <- values(mean.kernel.sampled.rss)
    
  }
  
  # add rowwise SDs to a new raster
  se.rast <- rast(mean.kernel.rss, vals = rowSds(pred.values))
  
  # bind together
  pred.all <- c(mean.kernel.rss, se.rast)
  names(pred.all) <- c("mean", "se")
  
  # return
  return(pred.all)
  
}

#_______________________________________________________________________
# 5b. Use function ----
#_______________________________________________________________________

kernel.S1L <- sim_kernel("simple", "low", 1, template.sp)
kernel.S2L <- sim_kernel("simple", "low", 2, template.sp)
kernel.S3L <- sim_kernel("simple", "low", 3, template.sp)
kernel.S1H <- sim_kernel("simple", "high", 1, template.sp)
kernel.S2H <- sim_kernel("simple", "high", 2, template.sp)
kernel.S3H <- sim_kernel("simple", "high", 3, template.sp)

kernel.C1L <- sim_kernel("complex", "low", 1, template.sp)
kernel.C2L <- sim_kernel("complex", "low", 2, template.sp)
kernel.C3L <- sim_kernel("complex", "low", 3, template.sp)
kernel.C1H <- sim_kernel("complex", "high", 1, template.sp)
kernel.C2H <- sim_kernel("complex", "high", 2, template.sp)
kernel.C3H <- sim_kernel("complex", "high", 3, template.sp)

#_______________________________________________________________________
# 6. Write rasters ----
#_______________________________________________________________________

writeRaster(kernel.S1L, paste0(getwd(), "/Rasters/UD rasters/S1L.tif"), overwrite = TRUE)
writeRaster(kernel.S2L, paste0(getwd(), "/Rasters/UD rasters/S2L.tif"), overwrite = TRUE)
writeRaster(kernel.S3L, paste0(getwd(), "/Rasters/UD rasters/S3L.tif"), overwrite = TRUE)
writeRaster(kernel.S1H, paste0(getwd(), "/Rasters/UD rasters/S1H.tif"), overwrite = TRUE)
writeRaster(kernel.S2H, paste0(getwd(), "/Rasters/UD rasters/S2H.tif"), overwrite = TRUE)
writeRaster(kernel.S3H, paste0(getwd(), "/Rasters/UD rasters/S3H.tif"), overwrite = TRUE)

writeRaster(kernel.C1L, paste0(getwd(), "/Rasters/UD rasters/C1L.tif"), overwrite = TRUE)
writeRaster(kernel.C2L, paste0(getwd(), "/Rasters/UD rasters/C2L.tif"), overwrite = TRUE)
writeRaster(kernel.C3L, paste0(getwd(), "/Rasters/UD rasters/C3L.tif"), overwrite = TRUE)
writeRaster(kernel.C1H, paste0(getwd(), "/Rasters/UD rasters/C1H.tif"), overwrite = TRUE)
writeRaster(kernel.C2H, paste0(getwd(), "/Rasters/UD rasters/C2H.tif"), overwrite = TRUE)
writeRaster(kernel.C3H, paste0(getwd(), "/Rasters/UD rasters/C3H.tif"), overwrite = TRUE)
