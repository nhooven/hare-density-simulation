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
  # mean
  pred.mean <- exp(landscape.covs.1$stem * focal.params$estimate[focal.params$term == "stem.s"] +
                   landscape.covs.1$edge * focal.params$estimate[focal.params$term == "edge.s"] +
                   landscape.covs.1$mature * focal.params$estimate[focal.params$term == "mature"])
  
  # lower 95
  # parameterize distribution
  qnorm()
  
  
}