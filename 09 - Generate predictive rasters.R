# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 09 - Generate predictive rasters
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 24 Dec 2024
# Date completed: 27 Dec 2024
# Date last modified: 10 Feb 2025
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation
library(terra)           # rasters
library(sf)
library(mvtnorm)         # multivariate normal distribution
library(matrixStats)     # rowwise SDs

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

# SSF/iSSF coefficients
all.params <- read.csv(paste0(getwd(), "/Derived_data/Model parameters/all_params.csv"))

# scaling factors (mean and sd)
all.scale <- read.csv(paste0(getwd(), "/Derived_data/Model parameters/mean_sd_covs.csv"))

# add in trt because I forgot it
all.scale$trt <- rep(c("before", "after", "before", "after"), each = 3)

# variance-covariance matrices
load(paste0(getwd(), "/Derived_data/Model parameters/vcov.RData"))

# vcov lookup table
vcov.lookup <- read.csv(paste0(getwd(), "/Derived_data/Lookup/vcov.csv"))

# unit boundary
unit.bound <- st_read(paste0(getwd(), "/Derived_data/Shapefiles/unit_bound.shp"))

# unprocessed covariate rasters
landscape.covs.B1 <- rast("Rasters/B1.tif")
landscape.covs.B2 <- rast("Rasters/B2.tif")
landscape.covs.B3 <- rast("Rasters/B3.tif")
landscape.covs.A1 <- rast("Rasters/A1.tif")
landscape.covs.A2 <- rast("Rasters/A2.tif")
landscape.covs.A3 <- rast("Rasters/A3.tif")

# sl/ta distributions
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_B1.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_B2.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_B3.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_A1.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_A2.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_A3.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/ta_dist.RData"))

#_______________________________________________________________________
# 3. Sample draws from the multivariate normal ----
#_______________________________________________________________________
# 3a. Define function ----
#_______________________________________________________________________

sample_mvn <- function (n.samples,
                        id.trt,
                        id.rep,
                        id.type) {
  
  # subset all.params
  focal.params <- all.params %>%
    
    filter(trt == id.trt,
           rep == id.rep,
           type == id.type)
  
  # subset variance-covariance matrix
  vcov.lookup.1 <- vcov.lookup %>% 
    
    filter(trt == id.trt,
           rep == id.rep,
           type == id.type)
  
  # subset list
  vcov.focal <- all.vcov[[vcov.lookup.1$index]]
  
  # sample from the multivariate normal
  if (id.type == "RSF") {
    
    mvn.samples <- as.data.frame(rmvnorm(n = n.samples, 
                                         mean = focal.params$estimate[c(1:3)], 
                                         sigma = vcov.focal))
    
  } else {
    
    mvn.samples <- as.data.frame(rmvnorm(n = n.samples, 
                                         mean = focal.params$estimate[c(1:6)], 
                                         sigma = vcov.focal))
    
  }
  
  # return
  return(mvn.samples)
  
}

#_______________________________________________________________________
# 3b. Use function ----
#_______________________________________________________________________

# RSF
mvn.sampled.RSF <- list(sample_mvn(100, "before", 1, "RSF"),
                        sample_mvn(100, "before", 2, "RSF"),
                        sample_mvn(100, "before", 3, "RSF"),
                        sample_mvn(100, "after", 1, "RSF"),
                        sample_mvn(100, "after", 2, "RSF"),
                        sample_mvn(100, "after", 3, "RSF"))

# iSSF
mvn.sampled.iSSF <- list(sample_mvn(100, "before", 1, "iSSF"),
                         sample_mvn(100, "before", 2, "iSSF"),
                         sample_mvn(100, "before", 3, "iSSF"),
                         sample_mvn(100, "after", 1, "iSSF"),
                         sample_mvn(100, "after", 2, "iSSF"),
                         sample_mvn(100, "after", 3, "iSSF"))

#_______________________________________________________________________
# 4. RSF predictions ----
#_______________________________________________________________________
# 4a. Define function ----
#_______________________________________________________________________

RSF_pred <- function (landscape.covs,     # raster
                      id.trt,
                      id.rep) {
  
  # subset all.params
  focal.params <- all.params %>%
    
    filter(trt == id.trt,
           rep == id.rep,
           type == "RSF")
  
  # subset scaling factors
  focal.scale <- all.scale %>%
    
    filter(trt == id.trt,
           rep == id.rep,
           model == "RSF")
  
  # subset variance-covariance matrix
  vcov.lookup.1 <- vcov.lookup %>% 
    
    filter(trt == id.trt,
           rep == id.rep,
           type == "RSF")
  
  # subset MVN samples
  mvn.samples <- mvn.sampled.RSF[[vcov.lookup.1$index]]
  
  # scale raster correctly
  landscape.covs.1 <- c((landscape.covs$forage - focal.scale$mean.forage) / focal.scale$sd.forage,
                        (landscape.covs$edge - focal.scale$mean.edge) / focal.scale$sd.edge,
                        landscape.covs$open)
  
  # calculate mean prediction
  # NOTE: These will not be exponentiated yet
  # mean
  pred.mean <- landscape.covs.1$forage * focal.params$estimate[focal.params$term == "forage.s"] +
               landscape.covs.1$edge * focal.params$estimate[focal.params$term == "edge.s"] +
               landscape.covs.1$open * focal.params$estimate[focal.params$term == "open"]
  
  # crop and exponentiate
  pred.mean.crop <- exp(crop(pred.mean, unit.bound))

  # calculate RSS for correction factor
  pred.mean.rss <- pred.mean.crop / mean(values(pred.mean.crop))
  
  # convert to correction factor
  pred.mean.cf <- 1 / pred.mean.rss
  
  # bootstrap for SE raster
  pred.boot <- matrix(data = NA, nrow = nrow(values(pred.mean.cf)), ncol = 100)
  
  for (i in 1:100) {
    
    # extract samples
    mvn.samples.focal <- mvn.samples[i, ]
    
    # calculate raster
    pred.sample <- landscape.covs.1$forage * mvn.samples.focal$V1 +
                   landscape.covs.1$edge * mvn.samples.focal$V2 +
                   landscape.covs.1$open * mvn.samples.focal$V3
    
    # crop and exponentiate
    pred.sample.crop <- exp(crop(pred.sample, unit.bound))
    
    # calculate RSS for correction factor
    pred.sample.rss <- pred.sample.crop / mean(values(pred.sample.crop))
    
    # convert to correction factor
    pred.sample.cf <- 1 / pred.sample.rss
    
    # extract values and bind into matrix
    pred.values <- values(pred.sample.cf)
    
    # bind into matrix
    pred.boot[ , i] <- pred.values
    
  }
  
  # add rowwise SDs to a new raster
  se.rast <- rast(pred.mean.cf, vals = rowSds(pred.boot))
  
  # bind together
  pred.all <- c(pred.mean.rss, pred.mean.cf, se.rast)
  names(pred.all) <- c("rss", "cf", "se")
  
  # return predictive raster
  return(pred.all)
  
}

#_______________________________________________________________________
# 4b. Use function ----
#_______________________________________________________________________

RSF.pred.B1 <- RSF_pred(landscape.covs.B1, "before", 1)
RSF.pred.B2 <- RSF_pred(landscape.covs.B2, "before", 2)
RSF.pred.B3 <- RSF_pred(landscape.covs.B3, "before", 3)
RSF.pred.A1 <- RSF_pred(landscape.covs.A1, "after", 1)
RSF.pred.A2 <- RSF_pred(landscape.covs.A2, "after", 2)
RSF.pred.A3 <- RSF_pred(landscape.covs.A3, "after", 3)

#_______________________________________________________________________
# 4c. Write rasters ----
#_______________________________________________________________________

writeRaster(RSF.pred.B1, paste0(getwd(), "/Rasters/RSF predictions/B1.tif"), overwrite = TRUE)
writeRaster(RSF.pred.B2, paste0(getwd(), "/Rasters/RSF predictions/B2.tif"), overwrite = TRUE)
writeRaster(RSF.pred.B3, paste0(getwd(), "/Rasters/RSF predictions/B3.tif"), overwrite = TRUE)
writeRaster(RSF.pred.A1, paste0(getwd(), "/Rasters/RSF predictions/A1.tif"), overwrite = TRUE)
writeRaster(RSF.pred.A2, paste0(getwd(), "/Rasters/RSF predictions/A2.tif"), overwrite = TRUE)
writeRaster(RSF.pred.A3, paste0(getwd(), "/Rasters/RSF predictions/A3.tif"), overwrite = TRUE)

#_______________________________________________________________________
# 5. Day range "adjustment" rasters ----

# we'll generate these using the habitat-specific step length
# distribution adjustments from the fitted iSSFs

# the main assumption in the REM is that day range is equal across
# the entire camera array. Each day range corresponds to a certain
# gamma SL distribution with E(x) = mean step length * 12

# first, we'll calculate distributions of day range for each scenario,
# then generate spatially-explicit expectation rasters

#_______________________________________________________________________
# 5a. Hyperdistribution of step lengths ----
#_______________________________________________________________________

# simulated data
sim.data.all <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/init_sims.csv"))

# collared lookup tables
load(paste0(getwd(), "/Derived_data/Lookup/collared_lookup_1.RData"))
load(paste0(getwd(), "/Derived_data/Lookup/collared_lookup_2.RData"))
load(paste0(getwd(), "/Derived_data/Lookup/collared_lookup_3.RData"))

# define function
static_dr <- function (id.trt,
                       id.rep) {
  
  # subset lookup table
  if (id.rep == 1) {
    
    lookup.focal <- collared.1
    
  } else {
    
    if (id.rep == 2) {
      
      lookup.focal <- collared.2
      
    } else {
      
      lookup.focal <- collared.3
      
    }
    
  }

  # subset all data
  sim.data.focal <- sim.data.all %>%
    
    filter(trt == id.trt &
           rep == id.rep &
           indiv %in% lookup.focal)
  
  # loop through all individuals
  all.static.dr <- data.frame()
  
  for (i in unique(sim.data.focal$indiv)) {
    
    # filter and keep only columns we need 
    sim.data.focal.1 <- sim.data.focal %>% 
      
      filter(indiv == i) %>%
      
      dplyr::select(x_, y_, t_, indiv) %>%
      
      # coerce t_ to POSIXct
      mutate(t_ = as.POSIXct(t_))
    
    # make track
    focal.track <- sim.data.focal.1  %>% make_track(.x = x_,
                                                    .y = y_,
                                                    .t = t_,
                                                    all_cols = TRUE,
                                                    check_duplicates = TRUE,
                                                    crs = 32611)
    
    # steps
    focal.steps <- steps(focal.track, keep_cols = "start")
    
    # fit SL distribution
    sl.dist.focal <- fit_distr(focal.steps$sl_, dist_name = "gamma")
    
    # calculate expected SL and multiply by 12
    mean.sl.focal <- sl.dist.focal$params$shape * sl.dist.focal$params$scale
    
    dr.focal <- mean.sl.focal * 12
    
    # create df and bind together
    all.static.dr <- rbind(all.static.dr, 
                           data.frame(trt = id.trt,
                                      rep = id.rep,
                                      indiv = i,
                                      mean.dr = dr.focal))
    
  }
  
  # return
  return(all.static.dr)
  
}

# apply function
static.dr.all <- rbind(static_dr("before", 1),
                       static_dr("before", 2),
                       static_dr("before", 3),
                       static_dr("after", 1),
                       static_dr("after", 2),
                       static_dr("after", 3))

# calculate means and SDs (and divide by 1000 since we really want this to be in km)
static.dr.summary <- static.dr.all %>%
  
  group_by(trt, rep) %>%
  
  summarize(dr.mean = mean(mean.dr) / 1000,
            dr.sd = sd(mean.dr) / 1000)

static.dr.summary

# write to .csv
write.csv(static.dr.summary, "Derived_data/Model parameters/static_dr.csv")

#_______________________________________________________________________
# 5b. Define function ----
#_______________________________________________________________________

sl_raster <- function(landscape.covs,     # raster
                      sl.dist,
                      id.trt,
                      id.rep) {
  
  # calculate expectation (shape * scale)
  mean.sl <- sl.dist$params$shape * sl.dist$params$scale
  
  # create raster
  base.sl.rast <- rast(landscape.covs$forage)
  
  # assign value
  values(base.sl.rast) <- mean.sl
  
  # extract iSSF parameters
  focal.params <- all.params %>%
    
    filter(trt == id.trt,
           rep == id.rep,
           type == "iSSF")
  
  # subset scaling factors
  focal.scale <- all.scale %>%
    
    filter(trt == id.trt,
           rep == id.rep,
           model == "SSF")
  
  # subset variance-covariance matrix
  vcov.lookup.1 <- vcov.lookup %>% 
    
    filter(trt == id.trt,
           rep == id.rep,
           type == "iSSF")
  
  # subset MVN samples
  mvn.samples <- mvn.sampled.iSSF[[vcov.lookup.1$index - 6]]
  
  # scale raster correctly
  landscape.covs.1 <- c((landscape.covs$forage - focal.scale$mean.forage) / focal.scale$sd.forage,
                        (landscape.covs$edge - focal.scale$mean.edge) / focal.scale$sd.edge,
                        landscape.covs$open)
  
  # calculate spatially-explicit SL shape parameter adjustments
  # this is a bit weird since one adjustment is on the start location and one is on the end location
  # we'll just roll with it
  pred.mean <- focal.params$estimate[focal.params$term == "log(sl_)"] +
               landscape.covs.1$forage * focal.params$estimate[focal.params$term == "forage.s:log(sl_)"] +
               landscape.covs.1$open * focal.params$estimate[focal.params$term == "open_start:log(sl_)"]
  
  # adjust the gamma distributions with the new shape parameters (in a raster)
  shape.adjust <- sl.dist$params$shape + pred.mean
  
  # and calculate the new expectation
  adjusted.mean.sl <- shape.adjust * sl.dist$params$scale 
  
  # crop
  adjusted.mean.sl.crop <- crop(adjusted.mean.sl, unit.bound)
  
  # transform to day range (in km)
  adjusted.day.range <- (adjusted.mean.sl.crop * 12) / 1000
  
  # bootstrap for SE raster
  pred.boot <- matrix(data = NA, nrow = nrow(values(adjusted.day.range)), ncol = 100)
  
  for (i in 1:100) {
    
    # extract samples
    mvn.samples.focal <- mvn.samples[i, ]
    
    # calculate spatially-explicit SL shape parameter adjustments
    pred.sample <- mvn.samples.focal$V4 +
                   landscape.covs.1$forage * mvn.samples.focal$V5 +
                   landscape.covs.1$open * mvn.samples.focal$V6
    
    # adjust the gamma distributions with the new shape parameters (in a raster)
    shape.adjust.sample <- sl.dist$params$shape + pred.sample
    
    # and calculate the new expectation
    adjusted.mean.sl.sample <- shape.adjust.sample * sl.dist$params$scale 
    
    # crop
    adjusted.mean.sl.crop.sample <- crop(adjusted.mean.sl.sample, unit.bound)
    
    # transform to day range (in km)
    adjusted.day.range.sample <- (adjusted.mean.sl.crop.sample * 12) / 1000
    
    # extract values
    pred.values <- values(adjusted.day.range.sample)
    
    # bind into matrix
    pred.boot[ , i] <- pred.values
    
  }
  
  # add rowwise SDs to a new raster
  se.rast <- rast(adjusted.day.range, vals = rowSds(pred.boot))
  
  # bind together
  pred.all <- c(adjusted.day.range, se.rast)
  names(pred.all) <- c("mean", "se")
  
  # return predictive raster
  return(pred.all)

}

#_______________________________________________________________________
# 5b. Use function ----
#_______________________________________________________________________

adj.dr.B1 <- sl_raster(landscape.covs.B1, sl.dist.B1, "before", 1)
adj.dr.B2 <- sl_raster(landscape.covs.B2, sl.dist.B2, "before", 2)
adj.dr.B3 <- sl_raster(landscape.covs.B3, sl.dist.B3, "before", 3)
adj.dr.A1 <- sl_raster(landscape.covs.A1, sl.dist.A1, "after", 1)
adj.dr.A2 <- sl_raster(landscape.covs.A2, sl.dist.A2, "after", 2)
adj.dr.A3 <- sl_raster(landscape.covs.A3, sl.dist.A3, "after", 3)

#_______________________________________________________________________
# 5c. Write rasters ----
#_______________________________________________________________________

writeRaster(adj.dr.B1, paste0(getwd(), "/Rasters/Adjusted day range/B1.tif"), overwrite = TRUE)
writeRaster(adj.dr.B2, paste0(getwd(), "/Rasters/Adjusted day range/B2.tif"), overwrite = TRUE)
writeRaster(adj.dr.B3, paste0(getwd(), "/Rasters/Adjusted day range/B3.tif"), overwrite = TRUE)
writeRaster(adj.dr.A1, paste0(getwd(), "/Rasters/Adjusted day range/A1.tif"), overwrite = TRUE)
writeRaster(adj.dr.A2, paste0(getwd(), "/Rasters/Adjusted day range/A2.tif"), overwrite = TRUE)
writeRaster(adj.dr.A3, paste0(getwd(), "/Rasters/Adjusted day range/A3.tif"), overwrite = TRUE)

#_______________________________________________________________________
# 6. Scaled covariate rasters for UD simulation ----
#_______________________________________________________________________
# 6a. Define function ----
#_______________________________________________________________________

scaled_covs <- function(landscape.covs,     # raster
                        id.trt,
                        id.rep) {
  
  # subset scaling factors
  focal.scale <- all.scale %>%
    
    filter(trt == id.trt,
           rep == id.rep,
           model == "SSF")
  
  # scale raster correctly
  landscape.covs.1 <- c((landscape.covs$forage - focal.scale$mean.forage) / focal.scale$sd.forage,
                        (landscape.covs$edge - focal.scale$mean.edge) / focal.scale$sd.edge,
                        landscape.covs$open)
  
  # return
  return(landscape.covs.1)
  
}

#_______________________________________________________________________
# 6b. Use function ----
#_______________________________________________________________________

scaled.covs.B1 <- scaled_covs(landscape.covs.B1, "before", 1)
scaled.covs.B2 <- scaled_covs(landscape.covs.B2, "before", 2)
scaled.covs.B3 <- scaled_covs(landscape.covs.B3, "before", 3)
scaled.covs.A1 <- scaled_covs(landscape.covs.A1, "after", 1)
scaled.covs.A2 <- scaled_covs(landscape.covs.A2, "after", 2)
scaled.covs.A3 <- scaled_covs(landscape.covs.A3, "after", 3)

#_______________________________________________________________________
# 6c. Write rasters ----
#_______________________________________________________________________

writeRaster(scaled.covs.B1, paste0(getwd(), "/Rasters/Scaled covariates/B1.tif"), overwrite = TRUE)
writeRaster(scaled.covs.B2, paste0(getwd(), "/Rasters/Scaled covariates/B2.tif"), overwrite = TRUE)
writeRaster(scaled.covs.B3, paste0(getwd(), "/Rasters/Scaled covariates/B3.tif"), overwrite = TRUE)
writeRaster(scaled.covs.A1, paste0(getwd(), "/Rasters/Scaled covariates/A1.tif"), overwrite = TRUE)
writeRaster(scaled.covs.A2, paste0(getwd(), "/Rasters/Scaled covariates/A2.tif"), overwrite = TRUE)
writeRaster(scaled.covs.A3, paste0(getwd(), "/Rasters/Scaled covariates/A3.tif"), overwrite = TRUE)

