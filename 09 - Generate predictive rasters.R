# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 09 - Generate predictive rasters
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 24 Dec 2024
# Date completed: 27 Dec 2024
# Date last modified: 30 Jan 2025
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

# variance-covariance matrices
load(paste0(getwd(), "/Derived_data/Model parameters/vcov.RData"))

# vcov lookup table
vcov.lookup <- read.csv(paste0(getwd(), "/Derived_data/Lookup/vcov.csv"))

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
# 3. Sample draws from the multivariate normal ----
#_______________________________________________________________________
# 3a. Define function ----
#_______________________________________________________________________

sample_mvn <- function (n.samples,
                        id.landscape,
                        id.variability,
                        id.rep,
                        id.type) {
  
  # subset all.params
  focal.params <- all.params %>%
    
    filter(landscape == id.landscape,
           variability == id.variability,
           rep == id.rep,
           type == id.type)
  
  # subset variance-covariance matrix
  vcov.lookup.1 <- vcov.lookup %>% 
    
    filter(landscape == id.landscape,
           variability == id.variability,
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
mvn.sampled.RSF <- list(sample_mvn(100, "simple", "low", 1, "RSF"),
                        sample_mvn(100, "simple", "low", 2, "RSF"),
                        sample_mvn(100, "simple", "low", 3, "RSF"),
                        sample_mvn(100, "simple", "high", 1, "RSF"),
                        sample_mvn(100, "simple", "high", 2, "RSF"),
                        sample_mvn(100, "simple", "high", 3, "RSF"),
                        sample_mvn(100, "complex", "low", 1, "RSF"),
                        sample_mvn(100, "complex", "low", 2, "RSF"),
                        sample_mvn(100, "complex", "low", 3, "RSF"),
                        sample_mvn(100, "complex", "high", 1, "RSF"),
                        sample_mvn(100, "complex", "high", 2, "RSF"),
                        sample_mvn(100, "complex", "high", 3, "RSF"))

# iSSF
mvn.sampled.iSSF <- list(sample_mvn(100, "simple", "low", 1, "iSSF"),
                         sample_mvn(100, "simple", "low", 2, "iSSF"),
                         sample_mvn(100, "simple", "low", 3, "iSSF"),
                         sample_mvn(100, "simple", "high", 1, "iSSF"),
                         sample_mvn(100, "simple", "high", 2, "iSSF"),
                         sample_mvn(100, "simple", "high", 3, "iSSF"),
                         sample_mvn(100, "complex", "low", 1, "iSSF"),
                         sample_mvn(100, "complex", "low", 2, "iSSF"),
                         sample_mvn(100, "complex", "low", 3, "iSSF"),
                         sample_mvn(100, "complex", "high", 1, "iSSF"),
                         sample_mvn(100, "complex", "high", 2, "iSSF"),
                         sample_mvn(100, "complex", "high", 3, "iSSF"))

#_______________________________________________________________________
# 4. RSF predictions ----
#_______________________________________________________________________
# 4a. Define function ----
#_______________________________________________________________________

RSF_pred <- function (landscape.covs,     # raster
                      id.landscape,
                      id.variability,
                      id.rep) {
  
  # subset all.params
  focal.params <- all.params %>%
    
    filter(landscape == id.landscape,
           variability == id.variability,
           rep == id.rep,
           type == "RSF")
  
  # subset scaling factors
  focal.scale <- all.scale %>%
    
    filter(landscape == id.landscape,
           variability == id.variability,
           rep == id.rep,
           model == "RSF")
  
  # subset variance-covariance matrix
  vcov.lookup.1 <- vcov.lookup %>% 
    
    filter(landscape == id.landscape,
           variability == id.variability,
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

RSF.pred.S1L <- RSF_pred(landscape.covs.S1, "simple", "low", 1)
RSF.pred.S2L <- RSF_pred(landscape.covs.S2, "simple", "low", 2)
RSF.pred.S3L <- RSF_pred(landscape.covs.S3, "simple", "low", 3)
RSF.pred.S1H <- RSF_pred(landscape.covs.S1, "simple", "high", 1)
RSF.pred.S2H <- RSF_pred(landscape.covs.S2, "simple", "high", 2)
RSF.pred.S3H <- RSF_pred(landscape.covs.S3, "simple", "high", 3)

RSF.pred.C1L <- RSF_pred(landscape.covs.C1, "complex", "low", 1)
RSF.pred.C2L <- RSF_pred(landscape.covs.C2, "complex", "low", 2)
RSF.pred.C3L <- RSF_pred(landscape.covs.C3, "complex", "low", 3)
RSF.pred.C1H <- RSF_pred(landscape.covs.C1, "complex", "high", 1)
RSF.pred.C2H <- RSF_pred(landscape.covs.C2, "complex", "high", 2)
RSF.pred.C3H <- RSF_pred(landscape.covs.C3, "complex", "high", 3)

#_______________________________________________________________________
# 4c. Write rasters ----
#_______________________________________________________________________

writeRaster(RSF.pred.S1L, paste0(getwd(), "/Rasters/RSF predictions/S1L.tif"), overwrite = TRUE)
writeRaster(RSF.pred.S2L, paste0(getwd(), "/Rasters/RSF predictions/S2L.tif"), overwrite = TRUE)
writeRaster(RSF.pred.S3L, paste0(getwd(), "/Rasters/RSF predictions/S3L.tif"), overwrite = TRUE)
writeRaster(RSF.pred.S1H, paste0(getwd(), "/Rasters/RSF predictions/S1H.tif"), overwrite = TRUE)
writeRaster(RSF.pred.S2H, paste0(getwd(), "/Rasters/RSF predictions/S2H.tif"), overwrite = TRUE)
writeRaster(RSF.pred.S3H, paste0(getwd(), "/Rasters/RSF predictions/S3H.tif"), overwrite = TRUE)

writeRaster(RSF.pred.C1L, paste0(getwd(), "/Rasters/RSF predictions/C1L.tif"), overwrite = TRUE)
writeRaster(RSF.pred.C2L, paste0(getwd(), "/Rasters/RSF predictions/C2L.tif"), overwrite = TRUE)
writeRaster(RSF.pred.C3L, paste0(getwd(), "/Rasters/RSF predictions/C3L.tif"), overwrite = TRUE)
writeRaster(RSF.pred.C1H, paste0(getwd(), "/Rasters/RSF predictions/C1H.tif"), overwrite = TRUE)
writeRaster(RSF.pred.C2H, paste0(getwd(), "/Rasters/RSF predictions/C2H.tif"), overwrite = TRUE)
writeRaster(RSF.pred.C3H, paste0(getwd(), "/Rasters/RSF predictions/C3H.tif"), overwrite = TRUE)

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
sim.data.all <- rbind(read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_S1L.csv")),
                      read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_S2L.csv")),
                      read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_S3L.csv")),
                      read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_S1H.csv")),
                      read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_S2H.csv")),
                      read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_S3H.csv")),
                      read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_C1L.csv")),
                      read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_C2L.csv")),
                      read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_C3L.csv")),
                      read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_C1H.csv")),
                      read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_C2H.csv")),
                      read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_C3H.csv")))

# lookup table
lookup <- read.csv(paste0(getwd(), "/Derived_data/Lookup/collared_lookup.csv"))

# define function
static_dr <- function (id.landscape,
                       id.variability,
                       id.rep) {
  
  # subset lookup table
  lookup.focal <- lookup %>%
    
    filter(landscape == id.landscape &
           variability == id.variability &
           rep == id.rep)
  
  # subset all data
  sim.data.focal <- sim.data.all %>%
    
    filter(landscape == id.landscape &
             variability == id.variability &
             rep == id.rep &
             indiv %in% lookup.focal$indiv)
  
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
                           data.frame(landscape = id.landscape,
                                      variability = id.variability,
                                      rep = id.rep,
                                      indiv = i,
                                      mean.dr = dr.focal))
    
  }
  
  # return
  return(all.static.dr)
  
}

# apply function
static.dr.all <- rbind(static_dr("simple", "low", 1),
                       static_dr("simple", "low", 2),
                       static_dr("simple", "low", 3),
                       static_dr("simple", "high", 1),
                       static_dr("simple", "high", 2),
                       static_dr("simple", "high", 3),
                       static_dr("complex", "low", 1),
                       static_dr("complex", "low", 2),
                       static_dr("complex", "low", 3),
                       static_dr("complex", "high", 1),
                       static_dr("complex", "high", 2),
                       static_dr("complex", "high", 3))

# calculate means and SDs
static.dr.summary <- static.dr.all %>%
  
  group_by(landscape, variability, rep) %>%
  
  summarize(dr.mean = mean(mean.dr),
            dr.sd = sd(mean.dr))

static.dr.summary

# write to .csv
write.csv(static.dr.summary, "Derived_data/Model parameters/static_dr.csv")

#_______________________________________________________________________
# 5b. Define function ----
#_______________________________________________________________________

sl_raster <- function(landscape.covs,     # raster
                      sl.dist,
                      id.landscape,
                      id.variability,
                      id.rep) {
  
  # calculate expectation (shape * scale)
  mean.sl <- sl.dist$params$shape * sl.dist$params$scale
  
  # create raster
  base.sl.rast <- rast(landscape.covs$forage)
  
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
           rep == id.rep,
           model == "SSF")
  
  # subset variance-covariance matrix
  vcov.lookup.1 <- vcov.lookup %>% 
    
    filter(landscape == id.landscape,
           variability == id.variability,
           rep == id.rep,
           type == "iSSF")
  
  # subset MVN samples
  mvn.samples <- mvn.sampled.iSSF[[vcov.lookup.1$index - 12]]
  
  # scale raster correctly
  landscape.covs.1 <- c((landscape.covs$forage - focal.scale$mean.forage) / focal.scale$sd.forage,
                        (landscape.covs$edge - focal.scale$mean.edge) / focal.scale$sd.edge,
                        landscape.covs$open)
  
  # calculate spatially-explicit SL shape parameter adjustments
  pred.mean <- focal.params$estimate[focal.params$term == "log(sl_)"] +
               landscape.covs.1$forage * focal.params$estimate[focal.params$term == "forage.s:log(sl_)"] +
               landscape.covs.1$open * focal.params$estimate[focal.params$term == "open:log(sl_)"]
  
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
# 5c. Write rasters ----
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
# 6. Scaled covariate rasters for UD simulation ----
#_______________________________________________________________________
# 6a. Define function ----
#_______________________________________________________________________

scaled_covs <- function(landscape.covs,     # raster
                        id.landscape,
                        id.variability,
                        id.rep) {
  
  # subset scaling factors
  focal.scale <- all.scale %>%
    
    filter(landscape == id.landscape,
           variability == id.variability,
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
# 6c. Write rasters ----
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

