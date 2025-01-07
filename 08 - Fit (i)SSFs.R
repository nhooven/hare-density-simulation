# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 08 - Fit (i)SSFs
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 09 Dec 2024
# Date completed: 09 Dec 2024
# Date last modified: 07 Jan 2025
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation
library(terra)           # rasters
library(amt)
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
# 5. Sample random steps and extract covariates ----
#_______________________________________________________________________
# 5a. Define function ----
#_______________________________________________________________________

sample_extract <- function (steps.df,
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
  
  # scale all continuous covariates and bind mean/SD into lookup table
  all.steps <- all.steps %>%
    
    mutate(stem.s = as.numeric(scale(stem)),
           edge.s = as.numeric(scale(edge)))
  
  # return
  return(all.steps)
  
}

#_______________________________________________________________________
# 5b. Use function ----
#_______________________________________________________________________

sampled.S1L <- sample_extract(steps.S1L, sl.dist.S1L, landscape.covs.S1)
sampled.S2L <- sample_extract(steps.S2L, sl.dist.S2L, landscape.covs.S2)
sampled.S3L <- sample_extract(steps.S3L, sl.dist.S3L, landscape.covs.S3)
sampled.S1H <- sample_extract(steps.S1H, sl.dist.S1H, landscape.covs.S1)
sampled.S2H <- sample_extract(steps.S2H, sl.dist.S2H, landscape.covs.S2)
sampled.S3H <- sample_extract(steps.S3H, sl.dist.S3H, landscape.covs.S3)

sampled.C1L <- sample_extract(steps.C1L, sl.dist.C1L, landscape.covs.C1)
sampled.C2L <- sample_extract(steps.C2L, sl.dist.C2L, landscape.covs.C2)
sampled.C3L <- sample_extract(steps.C3L, sl.dist.C3L, landscape.covs.C3)
sampled.C1H <- sample_extract(steps.C1H, sl.dist.C1H, landscape.covs.C1)
sampled.C2H <- sample_extract(steps.C2H, sl.dist.C2H, landscape.covs.C2)
sampled.C3H <- sample_extract(steps.C3H, sl.dist.C3H, landscape.covs.C3)

#_______________________________________________________________________
# 6. Means and SDs of continuous variables ----
#_______________________________________________________________________
# 6a. Define function ----
#_______________________________________________________________________

extract_mean_sd <- function(sampled.steps,
                            id.landscape,
                            id.variability,
                            id.rep) {
  
  return(data.frame(mean.stem = mean(sampled.steps$stem),
                    sd.stem = sd(sampled.steps$stem),
                    mean.edge = mean(sampled.steps$edge),
                    sd.edge = sd(sampled.steps$edge),
                    landscape = id.landscape,
                    variability = id.variability,
                    rep = id.rep))
  
}

#_______________________________________________________________________
# 6b. Use function ----
#_______________________________________________________________________

mean.sd.all <- rbind(extract_mean_sd(sampled.S1L, "simple", "low", 1),
                     extract_mean_sd(sampled.S2L, "simple", "low", 2),
                     extract_mean_sd(sampled.S3L, "simple", "low", 3),
                     extract_mean_sd(sampled.S1H, "simple", "high", 1),
                     extract_mean_sd(sampled.S2H, "simple", "high", 2),
                     extract_mean_sd(sampled.S3H, "simple", "high", 3),
                     extract_mean_sd(sampled.C1L, "complex", "low", 1),
                     extract_mean_sd(sampled.C2L, "complex", "low", 2),
                     extract_mean_sd(sampled.C3L, "complex", "low", 3),
                     extract_mean_sd(sampled.C1H, "complex", "high", 1),
                     extract_mean_sd(sampled.C2H, "complex", "high", 2),
                     extract_mean_sd(sampled.C3H, "complex", "high", 3))

#_______________________________________________________________________
# 7. Fit movement-naive models ----
#_______________________________________________________________________
# 7a. Define functions ----
#_______________________________________________________________________

# full random slopes
fit_naive_ssf <- function (sampled.steps) {
  
  SSF.struc <- glmmTMB(case_ ~ stem.s +
                               edge.s +
                               mature +
                               log(sl_) +
                               (0 + stem.s | indiv) +
                               (0 + edge.s | indiv) +
                               (0 + mature | indiv) +
                               (1 | step_id_),
                             family = poisson,
                             data = sampled.steps,
                             doFit = FALSE) 
  
  SSF.struc$parameters$theta[4] <- log(1e3)
  SSF.struc$mapArg <- list(theta = factor(c(1:3, NA)))
  
  SSF.model <- fitTMB(SSF.struc)
  
  # return
  return(SSF.model)
  
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
# 7b. Fit models ----
#_______________________________________________________________________

SSF.S1L <- fit_naive_ssf(sampled.S1L)
SSF.S2L <- fit_naive_ssf(sampled.S2L)
SSF.S3L <- fit_naive_ssf(sampled.S3L)
SSF.S1H <- fit_naive_ssf(sampled.S1H)
SSF.S2H <- fit_naive_ssf(sampled.S2H)
SSF.S3H <- fit_naive_ssf(sampled.S3H)
SSF.C1L <- fit_naive_ssf(sampled.C1L)
SSF.C2L <- fit_naive_ssf_1(sampled.C2L)
SSF.C3L <- fit_naive_ssf(sampled.C3L)
SSF.C1H <- fit_naive_ssf(sampled.C1H)
SSF.C2H <- fit_naive_ssf(sampled.C2H)
SSF.C3H <- fit_naive_ssf(sampled.C3H)

#_______________________________________________________________________
# 8. Fit correctly-specified iSSFs ----
#_______________________________________________________________________
# 8a. Define function ----
#_______________________________________________________________________

# full random slopes
fit_mv_issf <- function(sampled.steps) {
  
  iSSF.struc <- glmmTMB(case_ ~ stem.s +
                                stem.s:log(sl_) +
                                edge.s +
                                mature +
                                mature:log(sl_) +
                                log(sl_) +
                                (0 + stem.s | indiv) +
                                (0 + edge.s | indiv) +
                                (0 + mature | indiv) +
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
fit_mv_issf_1 <- function(sampled.steps) {
  
  iSSF.model <- fit_issf(data = sampled.steps,
                         formula = case_ ~ stem.s +
                                           stem.s:log(sl_) +
                                           edge.s +
                                           mature +
                                           mature:log(sl_) +
                                           log(sl_) +
                           strata(step_id_))
  
  # return
  return(iSSF.model)
  
}

#_______________________________________________________________________
# 7b. Fit models ----
#_______________________________________________________________________

iSSF.S1L <- fit_mv_issf(sampled.S1L)
iSSF.S2L <- fit_mv_issf(sampled.S2L)
iSSF.S3L <- fit_mv_issf_1(sampled.S3L)
iSSF.S1H <- fit_mv_issf(sampled.S1H)
iSSF.S2H <- fit_mv_issf_1(sampled.S2H)
iSSF.S3H <- fit_mv_issf(sampled.S3H)
iSSF.C1L <- fit_mv_issf(sampled.C1L)
iSSF.C2L <- fit_mv_issf_1(sampled.C2L)
iSSF.C3L <- fit_mv_issf(sampled.C3L)
iSSF.C1H <- fit_mv_issf(sampled.C1H)
iSSF.C2H <- fit_mv_issf(sampled.C2H)
iSSF.C3H <- fit_mv_issf(sampled.C3H)

#_______________________________________________________________________
# 8. Extract and bind HS coefficients together ----
#_______________________________________________________________________

# define function
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
                        data.frame(term = "sd__stem.s",
                                   estimate = 0.0,
                                   std.error = NA,
                                   statistic = NA,
                                   p.value = NA),
                        data.frame(term = "sd__edge.s",
                                   estimate = 0.0,
                                   std.error = NA,
                                   statistic = NA,
                                   p.value = NA),
                        data.frame(term = "sd__mature",
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

# extract and bind together
all.tidy <- rbind(extract_coefs(SSF.S1L, "SSF", "simple", "low", 1),
                  extract_coefs(SSF.S2L, "SSF", "simple", "low", 2),
                  extract_coefs(SSF.S3L, "SSF", "simple", "low", 3),
                  extract_coefs(SSF.S1H, "SSF", "simple", "high", 1),
                  extract_coefs(SSF.S2H, "SSF", "simple", "high", 2),
                  extract_coefs(SSF.S3H, "SSF", "simple", "high", 3),
                  extract_coefs(SSF.C1L, "SSF", "complex", "low", 1),
                  extract_coefs(SSF.C2L, "SSF", "complex", "low", 2),
                  extract_coefs(SSF.C3L, "SSF", "complex", "low", 3),
                  extract_coefs(SSF.C1H, "SSF", "complex", "high", 1),
                  extract_coefs(SSF.C2H, "SSF", "complex", "high", 2),
                  extract_coefs(SSF.C3H, "SSF", "complex", "high", 3),
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

# subset for plotting
all.tidy.beta <- all.tidy %>% filter(term %in% unique(all.tidy$term)[c(1:4, 8:9)])

# plot
ggplot(all.tidy.beta,
       aes(x = term,
           y = estimate,
           shape = landscape,
           color = variability)) +
  
  theme_bw() +
  
  geom_point()

#_______________________________________________________________________
# 9. Extract and store variance-covariance matrices ----
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
all.vcov <- list(extract_vcov(SSF.S1L, "SSF", "simple", "low", 1),
                 extract_vcov(SSF.S2L, "SSF", "simple", "low", 2),
                 extract_vcov(SSF.S3L, "SSF", "simple", "low", 3),
                 extract_vcov(SSF.S1H, "SSF", "simple", "high", 1),
                 extract_vcov(SSF.S2H, "SSF", "simple", "high", 2),
                 extract_vcov(SSF.S3H, "SSF", "simple", "high", 3),
                 extract_vcov(SSF.C1L, "SSF", "complex", "low", 1),
                 extract_vcov(SSF.C2L, "SSF", "complex", "low", 2),
                 extract_vcov(SSF.C3L, "SSF", "complex", "low", 3),
                 extract_vcov(SSF.C1H, "SSF", "complex", "high", 1),
                 extract_vcov(SSF.C2H, "SSF", "complex", "high", 2),
                 extract_vcov(SSF.C3H, "SSF", "complex", "high", 3),
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
                                 type = c("SSF", "iSSF")),
                     index = 1:24)

write.csv(vcov.lookup, paste0(getwd(), "/Derived_data/Lookup/vcov.csv"))

#_______________________________________________________________________
# 10. Write to .csv ----
#_______________________________________________________________________

write.csv(all.tidy, paste0(getwd(), "/Derived_data/Model parameters/all_params.csv"))
write.csv(mean.sd.all, paste0(getwd(), "/Derived_data/Model parameters/mean_sd_covs.csv"))
