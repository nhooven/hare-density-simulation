# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 09 - Generate predictive rasters
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 24 Dec 2024
# Date completed: 27 Dec 2024
# Date last modified: 27 Dec 2024
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation
library(terra)           # rasters
library(sf)

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

# SSF/iSSF coefficients
all.params <- read.csv(paste0(getwd(), "/Derived_data/Model parameters/all_params.csv"))

# scaling factors (mean and sd)
all.scale <- read.csv(paste0(getwd(), "/Derived_data/Model parameters/mean_sd_covs.csv"))

# unit boundary
unit.bound <- st_read(paste0(getwd(), "/Derived_data/Shapefiles/unit_bound.shp"))

# unprocessed covariate rasters
landscape.covs.S1 <- rast("Rasters/S1.tif")
landscape.covs.S2 <- rast("Rasters/S2.tif")
landscape.covs.S3 <- rast("Rasters/S3.tif")
landscape.covs.C1 <- rast("Rasters/C1.tif")
landscape.covs.C2 <- rast("Rasters/C2.tif")
landscape.covs.C3 <- rast("Rasters/C3.tif")

# sl/ta distributions
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_S1L.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_S2L.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_S3L.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_S1H.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_S2H.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_S3H.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_C1L.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_C2L.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_C3L.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_C1H.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_C2H.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_C3H.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/ta_dist.RData"))

#_______________________________________________________________________
# 3. Movement-naive SSF predictions ----
#_______________________________________________________________________
# 3a. Define function ----
#_______________________________________________________________________

SSF_pred <- function (landscape.covs,     # raster
                      id.landscape,
                      id.variability,
                      id.rep) {
  
  # subset all.params
  focal.params <- all.params %>%
    
    filter(landscape == id.landscape,
           variability == id.variability,
           rep == id.rep,
           type == "SSF")
  
  # subset scaling factors
  focal.scale <- all.scale %>%
    
    filter(landscape == id.landscape,
           variability == id.variability,
           rep == id.rep)
  
  # scale raster correctly
  landscape.covs.1 <- c((landscape.covs$stem - focal.scale$mean.stem) / focal.scale$sd.stem,
                        (landscape.covs$edge - focal.scale$mean.edge) / focal.scale$sd.edge,
                        landscape.covs$mature)
  
  # calculate mean, lower 95, and upper 95 predictions
  # NOTE: These will not be exponentiated yet
  # mean
  pred.mean <- landscape.covs.1$stem * focal.params$estimate[focal.params$term == "stem.s"] +
               landscape.covs.1$edge * focal.params$estimate[focal.params$term == "edge.s"] +
               landscape.covs.1$mature * focal.params$estimate[focal.params$term == "mature"]
  
  # lower and upper 95
  # parameterize distributions
  stem.low <- qnorm(p = 0.025, 
                    mean = focal.params$estimate[focal.params$term == "stem.s"], 
                    sd = focal.params$estimate[focal.params$term == "sd__stem.s"])
  
  stem.high <- qnorm(p = 0.975, 
                     mean = focal.params$estimate[focal.params$term == "stem.s"], 
                     sd = focal.params$estimate[focal.params$term == "sd__stem.s"])
  
  edge.low <- qnorm(p = 0.025, 
                    mean = focal.params$estimate[focal.params$term == "edge.s"], 
                    sd = focal.params$estimate[focal.params$term == "sd__edge.s"])
  
  edge.high <- qnorm(p = 0.975, 
                     mean = focal.params$estimate[focal.params$term == "edge.s"], 
                     sd = focal.params$estimate[focal.params$term == "sd__edge.s"])
  
  mature.low <- qnorm(p = 0.025, 
                      mean = focal.params$estimate[focal.params$term == "mature"], 
                      sd = focal.params$estimate[focal.params$term == "sd__mature"])
  
  mature.high <- qnorm(p = 0.975, 
                       mean = focal.params$estimate[focal.params$term == "mature"], 
                       sd = focal.params$estimate[focal.params$term == "sd__mature"])
  
  # calculate
  pred.low <- landscape.covs.1$stem * stem.low +
              landscape.covs.1$edge * edge.low +
              landscape.covs.1$mature * mature.low
  
  pred.high <- landscape.covs.1$stem * stem.high +
               landscape.covs.1$edge * edge.high +
               landscape.covs.1$mature * mature.high
  
  # crop and exponentiate
  pred.mean.crop <- exp(crop(pred.mean, unit.bound))
  pred.low.crop <- exp(crop(pred.low, unit.bound))
  pred.high.crop <- exp(crop(pred.high, unit.bound))
  
  # calculate RSS for correction factor
  pred.mean.rss <- pred.mean.crop / mean(values(pred.mean.crop))
  pred.low.rss <- pred.low.crop / mean(values(pred.low.crop))
  pred.high.rss <- pred.high.crop / mean(values(pred.high.crop))
  
  # bind together
  pred.all <- c(pred.low.rss, pred.mean.rss, pred.high.rss)
  names(pred.all) <- c("low", "mean", "high")
  
  # return predictive raster
  return(pred.all)
  
}

#_______________________________________________________________________
# 3b. Use function ----
#_______________________________________________________________________

SSF.pred.S1L <- SSF_pred(landscape.covs.S1, "simple", "low", 1)
SSF.pred.S2L <- SSF_pred(landscape.covs.S2, "simple", "low", 2)
SSF.pred.S3L <- SSF_pred(landscape.covs.S3, "simple", "low", 3)
SSF.pred.S1H <- SSF_pred(landscape.covs.S1, "simple", "high", 1)
SSF.pred.S2H <- SSF_pred(landscape.covs.S2, "simple", "high", 2)
SSF.pred.S3H <- SSF_pred(landscape.covs.S3, "simple", "high", 3)

SSF.pred.C1L <- SSF_pred(landscape.covs.C1, "complex", "low", 1)
SSF.pred.C2L <- SSF_pred(landscape.covs.C2, "complex", "low", 2)
SSF.pred.C3L <- SSF_pred(landscape.covs.C3, "complex", "low", 3)
SSF.pred.C1H <- SSF_pred(landscape.covs.C1, "complex", "high", 1)
SSF.pred.C2H <- SSF_pred(landscape.covs.C2, "complex", "high", 2)
SSF.pred.C3H <- SSF_pred(landscape.covs.C3, "complex", "high", 3)

#_______________________________________________________________________
# 3c. Write rasters ----
#_______________________________________________________________________

writeRaster(SSF.pred.S1L, paste0(getwd(), "/Rasters/SSF predictions/S1L.tif"), overwrite = TRUE)
writeRaster(SSF.pred.S2L, paste0(getwd(), "/Rasters/SSF predictions/S2L.tif"), overwrite = TRUE)
writeRaster(SSF.pred.S3L, paste0(getwd(), "/Rasters/SSF predictions/S3L.tif"), overwrite = TRUE)
writeRaster(SSF.pred.S1H, paste0(getwd(), "/Rasters/SSF predictions/S1H.tif"), overwrite = TRUE)
writeRaster(SSF.pred.S2H, paste0(getwd(), "/Rasters/SSF predictions/S2H.tif"), overwrite = TRUE)
writeRaster(SSF.pred.S3H, paste0(getwd(), "/Rasters/SSF predictions/S3H.tif"), overwrite = TRUE)

writeRaster(SSF.pred.C1L, paste0(getwd(), "/Rasters/SSF predictions/C1L.tif"), overwrite = TRUE)
writeRaster(SSF.pred.C2L, paste0(getwd(), "/Rasters/SSF predictions/C2L.tif"), overwrite = TRUE)
writeRaster(SSF.pred.C3L, paste0(getwd(), "/Rasters/SSF predictions/C3L.tif"), overwrite = TRUE)
writeRaster(SSF.pred.C1H, paste0(getwd(), "/Rasters/SSF predictions/C1H.tif"), overwrite = TRUE)
writeRaster(SSF.pred.C2H, paste0(getwd(), "/Rasters/SSF predictions/C2H.tif"), overwrite = TRUE)
writeRaster(SSF.pred.C3H, paste0(getwd(), "/Rasters/SSF predictions/C3H.tif"), overwrite = TRUE)

#_______________________________________________________________________
# 4. Day range "adjustment" rasters ----

# we'll generate these using the habitat-specific step length
# distribution adjustments from the fitted iSSFs

# the main assumption in the REM is that day range is equal across
# the entire camera array. Each day range corresponds to a certain
# gamma SL distribution with E(x) = mean step length * 12

#_______________________________________________________________________
# 4a. Define function ----
#_______________________________________________________________________

sl_raster <- function(landscape.covs,     # raster
                      sl.dist,
                      id.landscape,
                      id.variability,
                      id.rep) {
  
  # calculate expectation (shape * scale)
  mean.sl <- sl.dist$params$shape * sl.dist$params$scale
  
  # create raster
  base.sl.rast <- rast(landscape.covs$stem)
  
  # assign value
  values(base.sl.rast) <- mean.sl
  
  # extract iSSF parameters
  focal.params <- all.params %>%
    
    filter(landscape == id.landscape,
           variability == id.variability,
           rep == id.rep,
           type == "iSSF")
  
  # subset scaling factors
  focal.scale <- all.scale %>%
    
    filter(landscape == id.landscape,
           variability == id.variability,
           rep == id.rep)
  
  # scale raster correctly
  landscape.covs.1 <- c((landscape.covs$stem - focal.scale$mean.stem) / focal.scale$sd.stem,
                        (landscape.covs$edge - focal.scale$mean.edge) / focal.scale$sd.edge,
                        landscape.covs$mature)
  
  # calculate spatially-explicit SL shape parameter adjustments
  pred.mean <- focal.params$estimate[focal.params$term == "log(sl_)"] +
               landscape.covs.1$stem * focal.params$estimate[focal.params$term == "stem.s:log(sl_)"] +
               landscape.covs.1$mature * focal.params$estimate[focal.params$term == "mature:log(sl_)"]
  
  # adjust the gamma distributions with the new shape parameters (in a raster)
  shape.adjust <- sl.dist$params$shape + pred.mean
  
  # and calculate the new expectation
  adjusted.mean.sl <- shape.adjust * sl.dist$params$scale 
  
  # crop
  adjusted.mean.sl.crop <- crop(adjusted.mean.sl, unit.bound)
  
  # transform to day range (in km)
  adjusted.day.range <- (adjusted.mean.sl.crop * 12) / 1000
  
  # return
  return(adjusted.day.range)
  
}

#_______________________________________________________________________
# 4b. Use function ----
#_______________________________________________________________________

adj.dr.S1L <- sl_raster(landscape.covs.S1, sl.dist.S1L, "simple", "low", 1)
adj.dr.S2L <- sl_raster(landscape.covs.S2, sl.dist.S2L, "simple", "low", 2)
adj.dr.S3L <- sl_raster(landscape.covs.S3, sl.dist.S3L, "simple", "low", 3)
adj.dr.S1H <- sl_raster(landscape.covs.S1, sl.dist.S1H, "simple", "high", 1)
adj.dr.S2H <- sl_raster(landscape.covs.S2, sl.dist.S2H, "simple", "high", 2)
adj.dr.S3H <- sl_raster(landscape.covs.S3, sl.dist.S3H, "simple", "high", 3)

adj.dr.C1L <- sl_raster(landscape.covs.C1, sl.dist.C1L, "complex", "low", 1)
adj.dr.C2L <- sl_raster(landscape.covs.C2, sl.dist.C2L, "complex", "low", 2)
adj.dr.C3L <- sl_raster(landscape.covs.C3, sl.dist.C3L, "complex", "low", 3)
adj.dr.C1H <- sl_raster(landscape.covs.C1, sl.dist.C1H, "complex", "high", 1)
adj.dr.C2H <- sl_raster(landscape.covs.C2, sl.dist.C2H, "complex", "high", 2)
adj.dr.C3H <- sl_raster(landscape.covs.C3, sl.dist.C3H, "complex", "high", 3)

#_______________________________________________________________________
# 4c. Write rasters ----
#_______________________________________________________________________

writeRaster(adj.dr.S1L, paste0(getwd(), "/Rasters/Adjusted day range/S1L.tif"), overwrite = TRUE)
writeRaster(adj.dr.S2L, paste0(getwd(), "/Rasters/Adjusted day range/S2L.tif"), overwrite = TRUE)
writeRaster(adj.dr.S3L, paste0(getwd(), "/Rasters/Adjusted day range/S3L.tif"), overwrite = TRUE)
writeRaster(adj.dr.S1H, paste0(getwd(), "/Rasters/Adjusted day range/S1H.tif"), overwrite = TRUE)
writeRaster(adj.dr.S2H, paste0(getwd(), "/Rasters/Adjusted day range/S2H.tif"), overwrite = TRUE)
writeRaster(adj.dr.S3H, paste0(getwd(), "/Rasters/Adjusted day range/S3H.tif"), overwrite = TRUE)

writeRaster(adj.dr.C1L, paste0(getwd(), "/Rasters/Adjusted day range/C1L.tif"), overwrite = TRUE)
writeRaster(adj.dr.C2L, paste0(getwd(), "/Rasters/Adjusted day range/C2L.tif"), overwrite = TRUE)
writeRaster(adj.dr.C3L, paste0(getwd(), "/Rasters/Adjusted day range/C3L.tif"), overwrite = TRUE)
writeRaster(adj.dr.C1H, paste0(getwd(), "/Rasters/Adjusted day range/C1H.tif"), overwrite = TRUE)
writeRaster(adj.dr.C2H, paste0(getwd(), "/Rasters/Adjusted day range/C2H.tif"), overwrite = TRUE)
writeRaster(adj.dr.C3H, paste0(getwd(), "/Rasters/Adjusted day range/C3H.tif"), overwrite = TRUE)

#_______________________________________________________________________
# 5. Scaled covariate rasters for UD simulation ----
#_______________________________________________________________________
# 5a. Define function ----
#_______________________________________________________________________

scaled_covs <- function(landscape.covs,     # raster
                        id.landscape,
                        id.variability,
                        id.rep) {
  
  # subset scaling factors
  focal.scale <- all.scale %>%
    
    filter(landscape == id.landscape,
           variability == id.variability,
           rep == id.rep)
  
  # scale raster correctly
  landscape.covs.1 <- c((landscape.covs$stem - focal.scale$mean.stem) / focal.scale$sd.stem,
                        (landscape.covs$edge - focal.scale$mean.edge) / focal.scale$sd.edge,
                        landscape.covs$mature)
  
  # return
  return(landscape.covs.1)
  
}

#_______________________________________________________________________
# 5b. Use function ----
#_______________________________________________________________________

scaled.covs.S1L <- scaled_covs(landscape.covs.S1, "simple", "low", 1)
scaled.covs.S2L <- scaled_covs(landscape.covs.S2, "simple", "low", 2)
scaled.covs.S3L <- scaled_covs(landscape.covs.S3, "simple", "low", 3)
scaled.covs.S1H <- scaled_covs(landscape.covs.S1, "simple", "high", 1)
scaled.covs.S2H <- scaled_covs(landscape.covs.S2, "simple", "high", 2)
scaled.covs.S3H <- scaled_covs(landscape.covs.S3, "simple", "high", 3)

scaled.covs.C1L <- scaled_covs(landscape.covs.C1, "complex", "low", 1)
scaled.covs.C2L <- scaled_covs(landscape.covs.C2, "complex", "low", 2)
scaled.covs.C3L <- scaled_covs(landscape.covs.C3, "complex", "low", 3)
scaled.covs.C1H <- scaled_covs(landscape.covs.C1, "complex", "high", 1)
scaled.covs.C2H <- scaled_covs(landscape.covs.C2, "complex", "high", 2)
scaled.covs.C3H <- scaled_covs(landscape.covs.C3, "complex", "high", 3)

#_______________________________________________________________________
# 5c. Write rasters ----
#_______________________________________________________________________

writeRaster(scaled.covs.S1L, paste0(getwd(), "/Rasters/Scaled covariates/S1L.tif"), overwrite = TRUE)
writeRaster(scaled.covs.S2L, paste0(getwd(), "/Rasters/Scaled covariates/S2L.tif"), overwrite = TRUE)
writeRaster(scaled.covs.S3L, paste0(getwd(), "/Rasters/Scaled covariates/S3L.tif"), overwrite = TRUE)
writeRaster(scaled.covs.S1H, paste0(getwd(), "/Rasters/Scaled covariates/S1H.tif"), overwrite = TRUE)
writeRaster(scaled.covs.S2H, paste0(getwd(), "/Rasters/Scaled covariates/S2H.tif"), overwrite = TRUE)
writeRaster(scaled.covs.S3H, paste0(getwd(), "/Rasters/Scaled covariates/S3H.tif"), overwrite = TRUE)

writeRaster(scaled.covs.C1L, paste0(getwd(), "/Rasters/Scaled covariates/C1L.tif"), overwrite = TRUE)
writeRaster(scaled.covs.C2L, paste0(getwd(), "/Rasters/Scaled covariates/C2L.tif"), overwrite = TRUE)
writeRaster(scaled.covs.C3L, paste0(getwd(), "/Rasters/Scaled covariates/C3L.tif"), overwrite = TRUE)
writeRaster(scaled.covs.C1H, paste0(getwd(), "/Rasters/Scaled covariates/C1H.tif"), overwrite = TRUE)
writeRaster(scaled.covs.C2H, paste0(getwd(), "/Rasters/Scaled covariates/C2H.tif"), overwrite = TRUE)
writeRaster(scaled.covs.C3H, paste0(getwd(), "/Rasters/Scaled covariates/C3H.tif"), overwrite = TRUE)
