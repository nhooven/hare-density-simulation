# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 08 - Fit habitat models
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 09 Dec 2024
# Date completed: 09 Dec 2024
# Date last modified: 30 Jan 2025
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation
library(terra)           # rasters
library(amt)
library(sp)              # mcps
library(glmmTMB)         # modeling
library(broom.mixed)     # tidy model outputs

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

# lookup table
lookup <- read.csv(paste0(getwd(), "/Derived_data/Lookup/collared_lookup.csv"))

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
                      
# rasters
landscape.covs.S1 <- rast("Rasters/S1.tif")
landscape.covs.S2 <- rast("Rasters/S2.tif")
landscape.covs.S3 <- rast("Rasters/S3.tif")
landscape.covs.C1 <- rast("Rasters/C1.tif")
landscape.covs.C2 <- rast("Rasters/C2.tif")
landscape.covs.C3 <- rast("Rasters/C3.tif")

#_______________________________________________________________________
# 3. Clean data ----
#_______________________________________________________________________
# 3a. Define function to create tracks and steps ----
#_______________________________________________________________________

loop_steps <- function (id.landscape,
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
  all.focal.steps <- data.frame()
  
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
    
    # bind together
    all.focal.steps <- rbind(all.focal.steps, focal.steps)
    
  }
  
  # return
  return(all.focal.steps)
  
}

#_______________________________________________________________________
# 3b. Create steps ----
#_______________________________________________________________________

steps.S1L <- loop_steps("simple", "low", 1)
steps.S2L <- loop_steps("simple", "low", 2)
steps.S3L <- loop_steps("simple", "low", 3)
steps.S1H <- loop_steps("simple", "high", 1)
steps.S2H <- loop_steps("simple", "high", 2)
steps.S3H <- loop_steps("simple", "high", 3)

steps.C1L <- loop_steps("complex", "low", 1)
steps.C2L <- loop_steps("complex", "low", 2)
steps.C3L <- loop_steps("complex", "low", 3)
steps.C1H <- loop_steps("complex", "high", 1)
steps.C2H <- loop_steps("complex", "high", 2)
steps.C3H <- loop_steps("complex", "high", 3)

#_______________________________________________________________________
# 4. Step length and turning angle distributions ----
#_______________________________________________________________________
# 4a. Fit distributions ----
#_______________________________________________________________________

# sl
sl.dist.S1L <- fit_distr(steps.S1L$sl_, dist_name = "gamma")
sl.dist.S2L <- fit_distr(steps.S2L$sl_, dist_name = "gamma")
sl.dist.S3L <- fit_distr(steps.S3L$sl_, dist_name = "gamma")
sl.dist.S1H <- fit_distr(steps.S1H$sl_, dist_name = "gamma")
sl.dist.S2H <- fit_distr(steps.S2H$sl_, dist_name = "gamma")
sl.dist.S3H <- fit_distr(steps.S3H$sl_, dist_name = "gamma")

sl.dist.C1L <- fit_distr(steps.C1L$sl_, dist_name = "gamma")
sl.dist.C2L <- fit_distr(steps.C2L$sl_, dist_name = "gamma")
sl.dist.C3L <- fit_distr(steps.C3L$sl_, dist_name = "gamma")
sl.dist.C1H <- fit_distr(steps.C1H$sl_, dist_name = "gamma")
sl.dist.C2H <- fit_distr(steps.C2H$sl_, dist_name = "gamma")
sl.dist.C3H <- fit_distr(steps.C3H$sl_, dist_name = "gamma")

# ta
ta.dist <- make_unif_distr(min = -pi, max = pi)

#_______________________________________________________________________
# 4b. Save distributions ----
#_______________________________________________________________________

save(sl.dist.S1L, file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_S1L.RData"))
save(sl.dist.S2L, file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_S2L.RData"))
save(sl.dist.S3L, file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_S3L.RData"))
save(sl.dist.S1H, file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_S1H.RData"))
save(sl.dist.S2H, file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_S2H.RData"))
save(sl.dist.S3H, file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_S3H.RData"))
save(sl.dist.C1L, file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_C1L.RData"))
save(sl.dist.C2L, file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_C2L.RData"))
save(sl.dist.C3L, file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_C3L.RData"))
save(sl.dist.C1H, file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_C1H.RData"))
save(sl.dist.C2H, file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_C2H.RData"))
save(sl.dist.C3H, file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_C3H.RData"))

save(ta.dist, file = paste0(getwd(), "/Derived_data/Model parameters/ta_dist.RData"))

#_______________________________________________________________________
# 5. RSF - Sample random locations and extract covariates ----
#_______________________________________________________________________
# 5a. Define function ----
#_______________________________________________________________________

sample_extract_RSF <- function (steps.df,
                                landscape) {
  
  # loop through all individuals
  RSF.samples <- data.frame()
  
  for (i in unique(steps.df$indiv)) {
    
    # subset for MCP fitting
    steps.df.1 <- steps.df %>% 
      
      filter(indiv == i)
    
    # promote to SP
    steps.sp <- sp::SpatialPoints(coords = data.frame(x = steps.df.1$x1_,
                                                      y = steps.df.1$y1_),
                                  proj4string = sp::CRS("EPSG:32611"))
    
    # fit a 100% MCP
    focal.mcp <- adehabitatHR::mcp(steps.sp,
                                   percent = 100)
    
    # sample 30 random locations for each used location within MCP
    random.sp <- spsample(focal.mcp, 
                          n = length(steps.sp) * 30,
                          type = "regular")
    
    # promote both to SPDF and add case ID
    used.df <- data.frame(x = steps.sp@coords[ , 1],
                          y = steps.sp@coords[ , 2],
                          case = 1)
    
    used.spdf <- SpatialPointsDataFrame(coords = steps.sp@coords,
                                        data = used.df,
                                        proj4string = sp::CRS("EPSG:32611"))
    
    random.df <- data.frame(x = random.sp@coords[ , 1],
                            y = random.sp@coords[ , 2],
                            case = 0)
    
    random.spdf <- SpatialPointsDataFrame(coords = random.sp@coords,
                                          data = random.df,
                                          proj4string = sp::CRS("EPSG:32611"))
    
    # bind together and promote to sf
    all.sf <- st_as_sf(rbind(used.spdf, random.spdf))
    
    # extract values from rasters
    focal.forage <- terra::extract(landscape$forage, all.sf, ID = FALSE)
    focal.edge <- terra::extract(landscape$edge, all.sf, ID = FALSE)
    focal.open <- terra::extract(landscape$open, all.sf, ID = FALSE)
    
    # and bind in
    all.sf.1 <- cbind(all.sf, focal.forage, focal.edge, focal.open)

    # bind together
    RSF.samples <- as.data.frame(rbind(RSF.samples, all.sf.1))
    
  }
  
  # scale all continuous covariates and add weights for RSF
  RSF.samples <- RSF.samples %>%
    
    mutate(forage.s = as.numeric(scale(forage)),
           edge.s = as.numeric(scale(edge)),
           weight = ifelse(case == 0,
                           5000,
                           1),
           indiv = i)
    
  # return
  return(RSF.samples)
  
}

#_______________________________________________________________________
# 5b. Use function ----
#_______________________________________________________________________

RSF.S1L <- sample_extract_RSF(steps.S1L, landscape.covs.S1)
RSF.S2L <- sample_extract_RSF(steps.S2L, landscape.covs.S2)
RSF.S3L <- sample_extract_RSF(steps.S3L, landscape.covs.S3)
RSF.S1H <- sample_extract_RSF(steps.S1H, landscape.covs.S1)
RSF.S2H <- sample_extract_RSF(steps.S2H, landscape.covs.S2)
RSF.S3H <- sample_extract_RSF(steps.S3H, landscape.covs.S3)

RSF.C1L <- sample_extract_RSF(steps.C1L, landscape.covs.C1)
RSF.C2L <- sample_extract_RSF(steps.C2L, landscape.covs.C2)
RSF.C3L <- sample_extract_RSF(steps.C3L, landscape.covs.C3)
RSF.C1H <- sample_extract_RSF(steps.C1H, landscape.covs.C1)
RSF.C2H <- sample_extract_RSF(steps.C2H, landscape.covs.C2)
RSF.C3H <- sample_extract_RSF(steps.C3H, landscape.covs.C3)

#_______________________________________________________________________
# 6. SSF - Sample random steps and extract covariates ----
#_______________________________________________________________________
# 6a. Define function ----
#_______________________________________________________________________

sample_extract_SSF <- function (steps.df,
                                sl.dist,
                                landscape) {
  
  # loop through all individuals
  all.steps <- data.frame()
  
  for (i in unique(steps.df$indiv)) {

    # subset and extract
    steps.df.1 <- steps.df %>% 
      
      filter(indiv == i) %>%
      
      # sample random steps
      random_steps(n_control = 30,
                   sl_distr = sl.dist,
                   ta_distr = ta.dist) %>%
      
      # extract covariates
      extract_covariates(landscape)
    
    # bind together
    all.steps <- rbind(all.steps, steps.df.1)
    
  }
  
  # scale all continuous covariates
  all.steps <- all.steps %>%
    
    mutate(forage.s = as.numeric(scale(forage)),
           edge.s = as.numeric(scale(edge)))
  
  # return
  return(all.steps)
  
}

#_______________________________________________________________________
# 6b. Use function ----
#_______________________________________________________________________

SSF.S1L <- sample_extract_SSF(steps.S1L, sl.dist.S1L, landscape.covs.S1)
SSF.S2L <- sample_extract_SSF(steps.S2L, sl.dist.S2L, landscape.covs.S2)
SSF.S3L <- sample_extract_SSF(steps.S3L, sl.dist.S3L, landscape.covs.S3)
SSF.S1H <- sample_extract_SSF(steps.S1H, sl.dist.S1H, landscape.covs.S1)
SSF.S2H <- sample_extract_SSF(steps.S2H, sl.dist.S2H, landscape.covs.S2)
SSF.S3H <- sample_extract_SSF(steps.S3H, sl.dist.S3H, landscape.covs.S3)

SSF.C1L <- sample_extract_SSF(steps.C1L, sl.dist.C1L, landscape.covs.C1)
SSF.C2L <- sample_extract_SSF(steps.C2L, sl.dist.C2L, landscape.covs.C2)
SSF.C3L <- sample_extract_SSF(steps.C3L, sl.dist.C3L, landscape.covs.C3)
SSF.C1H <- sample_extract_SSF(steps.C1H, sl.dist.C1H, landscape.covs.C1)
SSF.C2H <- sample_extract_SSF(steps.C2H, sl.dist.C2H, landscape.covs.C2)
SSF.C3H <- sample_extract_SSF(steps.C3H, sl.dist.C3H, landscape.covs.C3)

#_______________________________________________________________________
# 7. Means and SDs of continuous variables ----
#_______________________________________________________________________
# 7a. Define function ----
#_______________________________________________________________________

extract_mean_sd <- function(sampled.steps,
                            id.landscape,
                            id.variability,
                            id.rep,
                            id.model) {
  
  return(data.frame(mean.forage = mean(sampled.steps$forage),
                    sd.forage = sd(sampled.steps$forage),
                    mean.edge = mean(sampled.steps$edge),
                    sd.edge = sd(sampled.steps$edge),
                    landscape = id.landscape,
                    variability = id.variability,
                    rep = id.rep,
                    model = id.model))
  
}

#_______________________________________________________________________
# 7b. Use function ----
#_______________________________________________________________________

mean.sd.all <- rbind(#RSF
                     extract_mean_sd(RSF.S1L, "simple", "low", 1, "RSF"),
                     extract_mean_sd(RSF.S2L, "simple", "low", 2, "RSF"),
                     extract_mean_sd(RSF.S3L, "simple", "low", 3, "RSF"),
                     extract_mean_sd(RSF.S1H, "simple", "high", 1, "RSF"),
                     extract_mean_sd(RSF.S2H, "simple", "high", 2, "RSF"),
                     extract_mean_sd(RSF.S3H, "simple", "high", 3, "RSF"),
                     extract_mean_sd(RSF.C1L, "complex", "low", 1, "RSF"),
                     extract_mean_sd(RSF.C2L, "complex", "low", 2, "RSF"),
                     extract_mean_sd(RSF.C3L, "complex", "low", 3, "RSF"),
                     extract_mean_sd(RSF.C1H, "complex", "high", 1, "RSF"),
                     extract_mean_sd(RSF.C2H, "complex", "high", 2, "RSF"),
                     extract_mean_sd(RSF.C3H, "complex", "high", 3, "RSF"),
                     
                     # SSF
                     extract_mean_sd(SSF.S1L, "simple", "low", 1, "SSF"),
                     extract_mean_sd(SSF.S2L, "simple", "low", 2, "SSF"),
                     extract_mean_sd(SSF.S3L, "simple", "low", 3, "SSF"),
                     extract_mean_sd(SSF.S1H, "simple", "high", 1, "SSF"),
                     extract_mean_sd(SSF.S2H, "simple", "high", 2, "SSF"),
                     extract_mean_sd(SSF.S3H, "simple", "high", 3, "SSF"),
                     extract_mean_sd(SSF.C1L, "complex", "low", 1, "SSF"),
                     extract_mean_sd(SSF.C2L, "complex", "low", 2, "SSF"),
                     extract_mean_sd(SSF.C3L, "complex", "low", 3, "SSF"),
                     extract_mean_sd(SSF.C1H, "complex", "high", 1, "SSF"),
                     extract_mean_sd(SSF.C2H, "complex", "high", 2, "SSF"),
                     extract_mean_sd(SSF.C3H, "complex", "high", 3, "SSF"))

#_______________________________________________________________________
# 8. Fit RSF models ----
#_______________________________________________________________________
# 8a. Define functions ----
#_______________________________________________________________________

# full random slopes
fit_RSF <- function (sampled.points) {
  
  RSF.struc <- glmmTMB(case ~ forage.s +
                              edge.s +
                              open +
                              (0 + forage.s | indiv) +
                              (0 + edge.s | indiv) +
                              (0 + open | indiv) +
                              (1 | indiv),
                             weights = weight,
                             family = binomial,
                             data = sampled.points,
                             doFit = FALSE) 
  
  RSF.struc$parameters$theta[4] <- log(1e3)
  RSF.struc$mapArg <- list(theta = factor(c(1:3, NA)))
  
  RSF.model <- fitTMB(RSF.struc)
  
  # return
  return(RSF.model)
  
}

# no random slopes to fix convergence issues
fit_naive_ssf_1 <- function(sampled.steps) {
  
  SSF.model <- fit_issf(data = sampled.steps,
                        formula = case_ ~ stem.s +
                                          edge.s +
                                          mature +
                                          log(sl_) +
                                          strata(step_id_))
  
  # return
  return(SSF.model)
  
}

#_______________________________________________________________________
# 8b. Fit models ----
#_______________________________________________________________________

RSF.S1L <- fit_RSF(RSF.S1L)
RSF.S2L <- fit_RSF(RSF.S2L)
RSF.S3L <- fit_RSF(RSF.S3L)
RSF.S1H <- fit_RSF(RSF.S1H)
RSF.S2H <- fit_RSF(RSF.S2H)
RSF.S3H <- fit_RSF(RSF.S3H)
RSF.C1L <- fit_RSF(RSF.C1L)
RSF.C2L <- fit_RSF(RSF.C2L)
RSF.C3L <- fit_RSF(RSF.C3L)
RSF.C1H <- fit_RSF(RSF.C1H)
RSF.C2H <- fit_RSF(RSF.C2H)
RSF.C3H <- fit_RSF(RSF.C3H)

#_______________________________________________________________________
# 9. Fit iSSF models ----
#_______________________________________________________________________
# 9a. Define function ----
#_______________________________________________________________________

# full random slopes
fit_issf <- function(sampled.steps) {
  
  iSSF.struc <- glmmTMB(case_ ~ forage.s +
                                forage.s:log(sl_) +
                                edge.s +
                                open +
                                open:log(sl_) +
                                log(sl_) +
                                (0 + forage.s | indiv) +
                                (0 + edge.s | indiv) +
                                (0 + open | indiv) +
                                (1 | step_id_),
                              family = poisson,
                              data = sampled.steps,
                              doFit = FALSE) 
  
  iSSF.struc$parameters$theta[4] <- log(1e3)
  iSSF.struc$mapArg <- list(theta = factor(c(1:3, NA)))
  
  iSSF.model <- fitTMB(iSSF.struc)
  
  # return
  return(iSSF.model)
  
}

# no random slopes to fix convergence issues
fit_issf_1 <- function(sampled.steps) {
  
  iSSF.model <- amt::fit_issf(data = sampled.steps,
                              formula = case_ ~ forage.s +
                                                forage.s:log(sl_) +
                                                edge.s +
                                                open +
                                                open:log(sl_) +
                                                log(sl_) +
                                strata(step_id_))
  
  # return
  return(iSSF.model)
  
}

#_______________________________________________________________________
# 9b. Fit models ----
#_______________________________________________________________________

iSSF.S1L <- fit_issf_1(SSF.S1L)
iSSF.S2L <- fit_issf(SSF.S2L)
iSSF.S3L <- fit_issf(SSF.S3L)
iSSF.S1H <- fit_issf(SSF.S1H)
iSSF.S2H <- fit_issf(SSF.S2H)
iSSF.S3H <- fit_issf_1(SSF.S3H)
iSSF.C1L <- fit_issf(SSF.C1L)
iSSF.C2L <- fit_issf(SSF.C2L)
iSSF.C3L <- fit_issf(SSF.C3L)
iSSF.C1H <- fit_issf(SSF.C1H)
iSSF.C2H <- fit_issf(SSF.C2H)
iSSF.C3H <- fit_issf(SSF.C3H)

#_______________________________________________________________________
# 10. Extract and bind HS coefficients together ----
#_______________________________________________________________________
# 10a. Define function ----
#_______________________________________________________________________

extract_coefs <- function(model,
                          id.type,
                          id.landscape,
                          id.variability,
                          id.rep) {
  
  # ask which class it is
  if (class(model)[1] == "glmmTMB") {
    
    # use tidy function
    tidy.model <- broom.mixed::tidy(model) %>%
      
      # add id columns
      mutate(landscape = id.landscape,
             variability = id.variability,
             rep = id.rep,
             type = id.type) %>%
      
      # remove first three columns
      dplyr::select(-c(effect, component, group)) %>%
      
      # drop sd__(Intercept) term
      filter(term != "sd__(Intercept)")
    
  } else {
    
    tidy.model <- tidy(model$model)
    
    # bind in 0 SD terms
    tidy.model <- rbind(tidy.model,
                        data.frame(term = "sd__forage.s",
                                   estimate = 0.0,
                                   std.error = NA,
                                   statistic = NA,
                                   p.value = NA),
                        data.frame(term = "sd__edge.s",
                                   estimate = 0.0,
                                   std.error = NA,
                                   statistic = NA,
                                   p.value = NA),
                        data.frame(term = "sd__open",
                                   estimate = 0.0,
                                   std.error = NA,
                                   statistic = NA,
                                   p.value = NA))
    
    # add in id columns
    tidy.model <- tidy.model %>%
      
      # add id columns
      mutate(landscape = id.landscape,
             variability = id.variability,
             rep = id.rep,
             type = id.type)
    
  }
  
    # return
    return(tidy.model)
  
}

#_______________________________________________________________________
# 10b. Use function and bind ----
#_______________________________________________________________________

all.tidy <- rbind(extract_coefs(RSF.S1L, "RSF", "simple", "low", 1),
                  extract_coefs(RSF.S2L, "RSF", "simple", "low", 2),
                  extract_coefs(RSF.S3L, "RSF", "simple", "low", 3),
                  extract_coefs(RSF.S1H, "RSF", "simple", "high", 1),
                  extract_coefs(RSF.S2H, "RSF", "simple", "high", 2),
                  extract_coefs(RSF.S3H, "RSF", "simple", "high", 3),
                  extract_coefs(RSF.C1L, "RSF", "complex", "low", 1),
                  extract_coefs(RSF.C2L, "RSF", "complex", "low", 2),
                  extract_coefs(RSF.C3L, "RSF", "complex", "low", 3),
                  extract_coefs(RSF.C1H, "RSF", "complex", "high", 1),
                  extract_coefs(RSF.C2H, "RSF", "complex", "high", 2),
                  extract_coefs(RSF.C3H, "RSF", "complex", "high", 3),
                  extract_coefs(iSSF.S1L, "iSSF", "simple", "low", 1),
                  extract_coefs(iSSF.S2L, "iSSF", "simple", "low", 2),
                  extract_coefs(iSSF.S3L, "iSSF", "simple", "low", 3),
                  extract_coefs(iSSF.S1H, "iSSF", "simple", "high", 1),
                  extract_coefs(iSSF.S2H, "iSSF", "simple", "high", 2),
                  extract_coefs(iSSF.S3H, "iSSF", "simple", "high", 3),
                  extract_coefs(iSSF.C1L, "iSSF", "complex", "low", 1),
                  extract_coefs(iSSF.C2L, "iSSF", "complex", "low", 2),
                  extract_coefs(iSSF.C3L, "iSSF", "complex", "low", 3),
                  extract_coefs(iSSF.C1H, "iSSF", "complex", "high", 1),
                  extract_coefs(iSSF.C2H, "iSSF", "complex", "high", 2),
                  extract_coefs(iSSF.C3H, "iSSF", "complex", "high", 3))

# drop intercept term
all.tidy <- all.tidy %>% filter(term != "(Intercept)")

#_______________________________________________________________________
# 10c. Plot ----
#_______________________________________________________________________

# subset for plotting
all.tidy.beta <- all.tidy %>% filter(term %in% unique(all.tidy$term)[c(1:3, 7:9)])

# plot
ggplot(all.tidy.beta,
       aes(y = term,
           x = estimate,
           shape = landscape,
           color = variability)) +
  
  facet_grid(landscape ~ variability) +
  
  geom_vline(xintercept = 0) +
  
  theme_bw() +
  
  geom_point() +
  
  theme(panel.grid = element_blank(),
        legend.position = "none")

#_______________________________________________________________________
# 11. Extract and store variance-covariance matrices ----
#_______________________________________________________________________

# define function
extract_vcov <- function(model,
                         id.type,
                         id.landscape,
                         id.variability,
                         id.rep) {
  
  # ask which class it is
  if (class(model)[1] == "glmmTMB") {
    
    # extract the matrix component and remove the intercept
    vcov.focal <- vcov(model)[[1]][-1, -1]
    
  } else {
    
    vcov.focal <- model$model$var
    
  }
  
  # return
  return(vcov.focal)
  
}

# extract and bind together
all.vcov <- list(extract_vcov(RSF.S1L, "RSF", "simple", "low", 1),
                 extract_vcov(RSF.S2L, "RSF", "simple", "low", 2),
                 extract_vcov(RSF.S3L, "RSF", "simple", "low", 3),
                 extract_vcov(RSF.S1H, "RSF", "simple", "high", 1),
                 extract_vcov(RSF.S2H, "RSF", "simple", "high", 2),
                 extract_vcov(RSF.S3H, "RSF", "simple", "high", 3),
                 extract_vcov(RSF.C1L, "RSF", "complex", "low", 1),
                 extract_vcov(RSF.C2L, "RSF", "complex", "low", 2),
                 extract_vcov(RSF.C3L, "RSF", "complex", "low", 3),
                 extract_vcov(RSF.C1H, "RSF", "complex", "high", 1),
                 extract_vcov(RSF.C2H, "RSF", "complex", "high", 2),
                 extract_vcov(RSF.C3H, "RSF", "complex", "high", 3),
                 extract_vcov(iSSF.S1L, "iSSF", "simple", "low", 1),
                 extract_vcov(iSSF.S2L, "iSSF", "simple", "low", 2),
                 extract_vcov(iSSF.S3L, "iSSF", "simple", "low", 3),
                 extract_vcov(iSSF.S1H, "iSSF", "simple", "high", 1),
                 extract_vcov(iSSF.S2H, "iSSF", "simple", "high", 2),
                 extract_vcov(iSSF.S3H, "iSSF", "simple", "high", 3),
                 extract_vcov(iSSF.C1L, "iSSF", "complex", "low", 1),
                 extract_vcov(iSSF.C2L, "iSSF", "complex", "low", 2),
                 extract_vcov(iSSF.C3L, "iSSF", "complex", "low", 3),
                 extract_vcov(iSSF.C1H, "iSSF", "complex", "high", 1),
                 extract_vcov(iSSF.C2H, "iSSF", "complex", "high", 2),
                 extract_vcov(iSSF.C3H, "iSSF", "complex", "high", 3)) 

# save to file
save(all.vcov, file = paste0(getwd(), "/Derived_data/Model parameters/vcov.RData"))

# write lookup table to csv
vcov.lookup <- cbind(expand.grid(rep = 1:3,
                                 variability = c("low", "high"),
                                 landscape = c("simple", "complex"),
                                 type = c("RSF", "iSSF")),
                     index = 1:24)

write.csv(vcov.lookup, paste0(getwd(), "/Derived_data/Lookup/vcov.csv"))

#_______________________________________________________________________
# 12. Write to .csv ----
#_______________________________________________________________________

write.csv(all.tidy, paste0(getwd(), "/Derived_data/Model parameters/all_params.csv"))
write.csv(mean.sd.all, paste0(getwd(), "/Derived_data/Model parameters/mean_sd_covs.csv"))
