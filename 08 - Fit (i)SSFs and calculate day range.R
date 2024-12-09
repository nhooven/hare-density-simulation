# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 08 - Fit (i)SSFs and calculate day range
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 09 Dec 2024
# Date completed: 09 Dec 2024
# Date last modified: 09 Dec 2024
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation
library(amt)
library(glmmTMB)         # modeling

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

# lookup table
lookup <- read.csv(paste0(getwd(), "/Derived_data/Lookup/lookup_collared.csv"))

# simulated data
sim.data.SL <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_simple_low.csv"))
sim.data.CL <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_complex_low.csv"))
sim.data.SH <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_simple_high.csv"))
sim.data.CH <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_complex_high.csv"))

# rasters
landscape.covs.simple <- rast("Rasters/simple.tif")
landscape.covs.complex <- rast("Rasters/complex.tif")

#_______________________________________________________________________
# 3. Clean data ----
#_______________________________________________________________________
# 3a. Split lookup tables ----
#_______________________________________________________________________

lookup.SL <- lookup %>% filter(landscape == "simple" & variability == "low")
lookup.CL <- lookup %>% filter(landscape == "complex" & variability == "low")
lookup.SH <- lookup %>% filter(landscape == "simple" & variability == "high")
lookup.CH <- lookup %>% filter(landscape == "complex" & variability == "high")

#_______________________________________________________________________
# 3b. Subset ----
#_______________________________________________________________________

sim.data.SL.1 <- sim.data.SL %>% filter(indiv %in% lookup.SL$indiv)
sim.data.CL.1 <- sim.data.CL %>% filter(indiv %in% lookup.CL$indiv)
sim.data.SH.1 <- sim.data.SH %>% filter(indiv %in% lookup.SH$indiv)
sim.data.CH.1 <- sim.data.CH %>% filter(indiv %in% lookup.CH$indiv)

#_______________________________________________________________________
# 3c. Define function to create tracks and steps ----
#_______________________________________________________________________

loop_steps <- function (x) {
  
  # loop through all individuals
  all.steps <- data.frame()
  
  for (i in unique(x$indiv)) {
    
    # filter and keep only columns we need 
    x.1 <- x %>% 
      
      filter(indiv == i) %>%
      
      dplyr::select(x_, y_, t_, indiv) %>%
      
      # coerce t_ to POSIXct
      mutate(t_ = as.POSIXct(t_))
    
    # make track
    x.track <- x.1 %>% make_track(.x = x_,
                                  .y = y_,
                                  .t = t_,
                                  all_cols = TRUE,
                                  check_duplicates = TRUE,
                                  crs = 32611)
    
    # steps
    x.steps <- steps(x.track, keep_cols = "start")
    
    # bind together
    all.steps <- rbind(x.steps, all.steps)
    
  }
  
  # return
  return(all.steps)
  
}

#_______________________________________________________________________
# 3d. Create steps ----
#_______________________________________________________________________

steps.SL <- loop_steps(sim.data.SL.1)
steps.CL <- loop_steps(sim.data.CL.1)
steps.SH <- loop_steps(sim.data.SH.1)
steps.CH <- loop_steps(sim.data.CH.1)

#_______________________________________________________________________
# 4. Calculate mean day range from 2-hr step lengths ----
#_______________________________________________________________________

day.range <- data.frame(landscape = c("simple", "complex", "simple", "complex"),
                        variability = c("low", "low", "high", "high"),
                        day.range = c(mean(steps.SL$sl_) * 12 / 1000,
                                      mean(steps.CL$sl_) * 12 / 1000,
                                      mean(steps.SH$sl_) * 12 / 1000,
                                      mean(steps.CH$sl_) * 12 / 1000))

# write to .csv
write.csv(day.range, paste0(getwd(), "/Derived_data/Model parameters/day_range.csv"))

#_______________________________________________________________________
# 5. Step length and turning angle distributions ----
#_______________________________________________________________________

# sl
sl.dist.SL <- fit_distr(steps.SL$sl_, dist_name = "gamma")
sl.dist.CL <- fit_distr(steps.CL$sl_, dist_name = "gamma")
sl.dist.SH <- fit_distr(steps.SH$sl_, dist_name = "gamma")
sl.dist.CH <- fit_distr(steps.CH$sl_, dist_name = "gamma")

# ta
ta.dist <- make_unif_distr(min = -pi, max = pi)

# save distributions
save(sl.dist.SL, file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_SL.RData"))
save(sl.dist.CL, file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_CL.RData"))
save(sl.dist.SH, file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_SH.RData"))
save(sl.dist.CH, file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_CH.RData"))
save(ta.dist, file = paste0(getwd(), "/Derived_data/Model parameters/ta_dist.RData"))

#_______________________________________________________________________
# 6. Sample random steps and extract covariates ----
#_______________________________________________________________________
# 6a. Define function ----
#_______________________________________________________________________

sample_extract <- function (x,
                            sl.dist,
                            landscape) {
  
  # loop through all individuals
  all.steps <- data.frame()
  
  for (i in unique(x$indiv)) {
    
    x.1 <- x %>% 
      
      filter(indiv == i) %>%
      
      # sample random steps
      random_steps(n_control = 20,
                   sl_distr = sl.dist,
                   ta_distr = ta.dist) %>%
      
      # extract covariates
      extract_covariates(landscape)
    
    # bind together
    all.steps <- rbind(x.1, all.steps)
    
  }
  
  # return
  return(all.steps)
  
}

#_______________________________________________________________________
# 6b. Use function ----
#_______________________________________________________________________

sampled.SL <- sample_extract(steps.SL, sl.dist.SL, landscape.covs.simple)
sampled.CL <- sample_extract(steps.CL, sl.dist.CL, landscape.covs.complex)
sampled.SH <- sample_extract(steps.SH, sl.dist.SH, landscape.covs.simple)
sampled.CH <- sample_extract(steps.CH, sl.dist.CH, landscape.covs.complex)

#_______________________________________________________________________
# 7. Fit movement-naive models ----
#_______________________________________________________________________
# 7a. Simple - low ----
#_______________________________________________________________________

SSF.SL.struc <- glmmTMB(case_ ~ stem +
                                edge +
                                mature +
                                log(sl_) +
                                (0 + stem | indiv) +
                                (0 + edge | indiv) +
                                (0 + mature | indiv) +
                                (1 | step_id_),
                        family = poisson,
                        data = sampled.SL,
                        doFit = FALSE) 

SSF.SL.struc$parameters$theta[4] <- log(1e3)
SSF.SL.struc$mapArg <- list(theta = factor(c(1:3, NA)))

SSF.SL.model <- fitTMB(SSF.SL.struc)

summary(SSF.SL.model)

#_______________________________________________________________________
# 7b. Complex - low ----
#_______________________________________________________________________

SSF.CL.struc <- glmmTMB(case_ ~ stem +
                                edge +
                                mature +
                                log(sl_) +
                                (0 + stem | indiv) +
                                (0 + edge | indiv) +
                                (0 + mature | indiv) +
                                (1 | step_id_),
                        family = poisson,
                        data = sampled.CL,
                        doFit = FALSE) 

SSF.CL.struc$parameters$theta[4] <- log(1e3)
SSF.CL.struc$mapArg <- list(theta = factor(c(1:3, NA)))

SSF.CL.model <- fitTMB(SSF.CL.struc)

summary(SSF.CL.model)

#_______________________________________________________________________
# 7c. Simple - high ----
#_______________________________________________________________________

SSF.SH.struc <- glmmTMB(case_ ~ stem +
                                edge +
                                mature +
                                log(sl_) +
                                (0 + stem | indiv) +
                                (0 + edge | indiv) +
                                (0 + mature | indiv) +
                                (1 | step_id_),
                        family = poisson,
                        data = sampled.SH,
                        doFit = FALSE) 

SSF.SH.struc$parameters$theta[4] <- log(1e3)
SSF.SH.struc$mapArg <- list(theta = factor(c(1:3, NA)))

SSF.SH.model <- fitTMB(SSF.SH.struc)

summary(SSF.SH.model)

#_______________________________________________________________________
# 7d. Complex - high ----
#_______________________________________________________________________

SSF.CH.struc <- glmmTMB(case_ ~ stem +
                                edge +
                                mature +
                                log(sl_) +
                                (0 + stem | indiv) +
                                (0 + edge | indiv) +
                                (0 + mature | indiv) +
                                (1 | step_id_),
                        family = poisson,
                        data = sampled.CH,
                        doFit = FALSE) 

SSF.CH.struc$parameters$theta[4] <- log(1e3)
SSF.CH.struc$mapArg <- list(theta = factor(c(1:3, NA)))

SSF.CH.model <- fitTMB(SSF.CH.struc)

summary(SSF.CH.model)

#_______________________________________________________________________
# 7e. Extract and bind HS coefficients together ----
#_______________________________________________________________________

SSF.params <- rbind(data.frame(names = c("stem", "edge", "mature"),
                               betas = SSF.SL.model$fit$par[2:4],
                               landscape = "simple",
                               variability = "low"),
                    data.frame(names = c("stem", "edge", "mature"),
                               betas = SSF.CL.model$fit$par[2:4],
                               landscape = "complex",
                               variability = "low"),
                    data.frame(names = c("stem", "edge", "mature"),
                               betas = SSF.SH.model$fit$par[2:4],
                               landscape = "simple",
                               variability = "high"),
                    data.frame(names = c("stem", "edge", "mature"),
                               betas = SSF.CH.model$fit$par[2:4],
                               landscape = "complex",
                               variability = "high"))

# plot
ggplot(SSF.params,
       aes(x = names,
           y = betas,
           shape = landscape,
           color = variability)) +
  
  theme_bw() +
  
  geom_point()

#_______________________________________________________________________
# 8. Movement-naive raster predictions ----

# here we'll calculate log-RSS (the movement-free selection kernel) from these
# movement-naive SSFs across the study area and then extract these values at 
# each camera as the "correction factor"

#_______________________________________________________________________
# 8a. Simple - low ----
#_______________________________________________________________________

SSF.params.1 <- SSF.params %>% 
  
  filter(landscape == "simple",
         variability == "low")

pred.raster.SL <- landscape.covs.simple$stem * SSF.params.1$betas[1] +
                  landscape.covs.simple$edge * SSF.params.1$betas[2] +
                  landscape.covs.simple$mature * SSF.params.1$betas[3]

terra::plot(pred.raster.SL)

#_______________________________________________________________________
# 8b. Complex - low ----
#_______________________________________________________________________

SSF.params.2 <- SSF.params %>% 
  
  filter(landscape == "complex",
         variability == "low")

pred.raster.CL <- landscape.covs.complex$stem * SSF.params.2$betas[1] +
                  landscape.covs.complex$edge * SSF.params.2$betas[2] +
                  landscape.covs.complex$mature * SSF.params.2$betas[3]

terra::plot(pred.raster.CL)

#_______________________________________________________________________
# 8c. Simple - high ----
#_______________________________________________________________________

SSF.params.3 <- SSF.params %>% 
  
  filter(landscape == "simple",
         variability == "high")

pred.raster.SH <- landscape.covs.simple$stem * SSF.params.3$betas[1] +
                  landscape.covs.simple$edge * SSF.params.3$betas[2] +
                  landscape.covs.simple$mature * SSF.params.3$betas[3]

terra::plot(pred.raster.SH)

#_______________________________________________________________________
# 8d. Complex - high ----
#_______________________________________________________________________

SSF.params.4 <- SSF.params %>% 
  
  filter(landscape == "complex",
         variability == "high")

pred.raster.CH <- landscape.covs.complex$stem * SSF.params.4$betas[1] +
                  landscape.covs.complex$edge * SSF.params.4$betas[2] +
                  landscape.covs.complex$mature * SSF.params.4$betas[3]

terra::plot(pred.raster.CH)

#_______________________________________________________________________
# 8e. Bind rasters together, rename, and save ----
#_______________________________________________________________________

SSF.pred.rasters <- c(pred.raster.SL, 
                      pred.raster.SH,
                      pred.raster.CL,
                      pred.raster.CH)

names(SSF.pred.rasters) <- c("SL", "SH", "CL", "CH")

writeRaster(SSF.pred.rasters, filename = "Rasters/SSF_pred.tif", overwrite = T)

#_______________________________________________________________________
# 9. Fit correctly-specified iSSFs ----
#_______________________________________________________________________
# 9a. Simple - low ----
#_______________________________________________________________________

iSSF.SL.struc <- glmmTMB(case_ ~ stem +
                                 stem:log(sl_) +
                                 edge +
                                 mature +
                                 mature:log(sl_) +
                                 log(sl_) +
                                 (0 + stem | indiv) +
                                 (0 + edge | indiv) +
                                 (0 + mature | indiv) +
                                 (1 | step_id_),
                        family = poisson,
                        data = sampled.SL,
                        doFit = FALSE) 

iSSF.SL.struc$parameters$theta[4] <- log(1e3)
iSSF.SL.struc$mapArg <- list(theta = factor(c(1:3, NA)))

iSSF.SL.model <- fitTMB(iSSF.SL.struc)

summary(iSSF.SL.model)

#_______________________________________________________________________
# 9b. Complex - low ----
#_______________________________________________________________________

iSSF.CL.struc <- glmmTMB(case_ ~ stem +
                                 stem:log(sl_) +
                                 edge +
                                 mature +
                                 mature:log(sl_) +
                                 log(sl_) +
                                 (0 + stem | indiv) +
                                 (0 + edge | indiv) +
                                 (0 + mature | indiv) +
                                 (1 | step_id_),
                        family = poisson,
                        data = sampled.CL,
                        doFit = FALSE) 

iSSF.CL.struc$parameters$theta[4] <- log(1e3)
iSSF.CL.struc$mapArg <- list(theta = factor(c(1:3, NA)))

iSSF.CL.model <- fitTMB(iSSF.CL.struc)

summary(iSSF.CL.model)

#_______________________________________________________________________
# 9c. Simple - high ----
#_______________________________________________________________________

iSSF.SH.struc <- glmmTMB(case_ ~ stem +
                                 stem:log(sl_) +
                                 edge +
                                 mature +
                                 mature:log(sl_) +
                                 log(sl_) +
                                 (0 + stem | indiv) +
                                 (0 + edge | indiv) +
                                 (0 + mature | indiv) +
                                 (1 | step_id_),
                        family = poisson,
                        data = sampled.SH,
                        doFit = FALSE) 

iSSF.SH.struc$parameters$theta[4] <- log(1e3)
iSSF.SH.struc$mapArg <- list(theta = factor(c(1:3, NA)))

iSSF.SH.model <- fitTMB(iSSF.SH.struc)

summary(iSSF.SH.model)

#_______________________________________________________________________
# 9d. Complex - high ----
#_______________________________________________________________________

iSSF.CH.struc <- glmmTMB(case_ ~ stem +
                                 stem:log(sl_) +
                                 edge +
                                 mature +
                                 mature:log(sl_) +
                                 log(sl_) +
                                 (0 + stem | indiv) +
                                 (0 + edge | indiv) +
                                 (0 + mature | indiv) +
                                 (1 | step_id_),
                        family = poisson,
                        data = sampled.CH,
                        doFit = FALSE) 

iSSF.CH.struc$parameters$theta[4] <- log(1e3)
iSSF.CH.struc$mapArg <- list(theta = factor(c(1:3, NA)))

iSSF.CH.model <- fitTMB(iSSF.CH.struc)

summary(iSSF.CH.model)

#_______________________________________________________________________
# 9e. Extract and bind coefficients together ----
#_______________________________________________________________________

iSSF.params <- rbind(data.frame(names = c("stem", "edge", "mature", "log_sl", "stem_sl", "mature_sl"),
                                betas = iSSF.SL.model$fit$par[2:7],
                                landscape = "simple",
                                variability = "low"),
                     data.frame(names = c("stem", "edge", "mature", "log_sl", "stem_sl", "mature_sl"),
                                betas = iSSF.CL.model$fit$par[2:7],
                                landscape = "complex",
                                variability = "low"),
                     data.frame(names = c("stem", "edge", "mature", "log_sl", "stem_sl", "mature_sl"),
                                betas = iSSF.SH.model$fit$par[2:7],
                                landscape = "simple",
                                variability = "high"),
                     data.frame(names = c("stem", "edge", "mature", "log_sl", "stem_sl", "mature_sl"),
                                betas = iSSF.CH.model$fit$par[2:7],
                                landscape = "complex",
                                variability = "high"))

# bind with movement-naive SSF parameters
SSF.params$type <- "SSF"
iSSF.params$type <- "iSSF"

all.params <- rbind(SSF.params, iSSF.params)

#_______________________________________________________________________
# 10. Write to .csv ----
#_______________________________________________________________________

write.csv(all.params, paste0(getwd(), "/Derived_data/Model parameters/all_params.csv"))
