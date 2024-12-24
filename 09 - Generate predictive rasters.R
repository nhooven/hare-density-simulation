# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 09 - Generate predictive rasters
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 24 Dec 2024
# Date completed: 
# Date last modified: 
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