# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 08 - Fit habitat models
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 09 Dec 2024
# Date completed: 09 Dec 2024
# Date last modified: 10 Feb 2025
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
load(paste0(getwd(), "/Derived_data/Lookup/collared_lookup_1.RData"))
load(paste0(getwd(), "/Derived_data/Lookup/collared_lookup_2.RData"))
load(paste0(getwd(), "/Derived_data/Lookup/collared_lookup_3.RData"))

# simulated data
sim.data.all <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/init_sims.csv"))
                      
# rasters
landscape.covs.B1 <- rast("Rasters/B1.tif")
landscape.covs.B2 <- rast("Rasters/B2.tif")
landscape.covs.B3 <- rast("Rasters/B3.tif")
landscape.covs.A1 <- rast("Rasters/A1.tif")
landscape.covs.A2 <- rast("Rasters/A2.tif")
landscape.covs.A3 <- rast("Rasters/A3.tif")

#_______________________________________________________________________
# 3. Clean data ----
#_______________________________________________________________________
# 3a. Define function to create tracks and steps ----
#_______________________________________________________________________

loop_steps <- function (id.trt,
                        id.rep) {
  
  # subset lookup 
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

steps.B1 <- loop_steps("before", 1)
steps.B2 <- loop_steps("before", 2)
steps.B3 <- loop_steps("before", 3)

steps.A1 <- loop_steps("after", 1)
steps.A2 <- loop_steps("after", 2)
steps.A3 <- loop_steps("after", 3)

#_______________________________________________________________________
# 4. Step length and turning angle distributions ----
#_______________________________________________________________________
# 4a. Fit distributions ----
#_______________________________________________________________________

# sl
sl.dist.B1 <- fit_distr(steps.B1$sl_, dist_name = "gamma")
sl.dist.B2 <- fit_distr(steps.B2$sl_, dist_name = "gamma")
sl.dist.B3 <- fit_distr(steps.B3$sl_, dist_name = "gamma")

sl.dist.A1 <- fit_distr(steps.A1$sl_, dist_name = "gamma")
sl.dist.A2 <- fit_distr(steps.A2$sl_, dist_name = "gamma")
sl.dist.A3 <- fit_distr(steps.A3$sl_, dist_name = "gamma")

# ta
ta.dist <- make_unif_distr(min = -pi, max = pi)

#_______________________________________________________________________
# 4b. Save distributions ----
#_______________________________________________________________________

save(sl.dist.B1, file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_B1.RData"))
save(sl.dist.B2, file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_B2.RData"))
save(sl.dist.B3, file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_B3.RData"))
save(sl.dist.A1, file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_A1.RData"))
save(sl.dist.A2, file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_A2.RData"))
save(sl.dist.A3, file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_A3.RData"))

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

RSF.B1 <- sample_extract_RSF(steps.B1, landscape.covs.B1)
RSF.B2 <- sample_extract_RSF(steps.B2, landscape.covs.B2)
RSF.B3 <- sample_extract_RSF(steps.B3, landscape.covs.B3)
RSF.A1 <- sample_extract_RSF(steps.A1, landscape.covs.A1)
RSF.A2 <- sample_extract_RSF(steps.A2, landscape.covs.A2)
RSF.A3 <- sample_extract_RSF(steps.A3, landscape.covs.A3)

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
      
      # extract covariates (end and start of step)
      extract_covariates(landscape, where = "both")
    
    # bind together
    all.steps <- rbind(all.steps, steps.df.1)
    
  }
  
  # scale all continuous covariates
  all.steps <- all.steps %>%
    
    mutate(forage.s = as.numeric(scale(forage_end)),
           edge.s = as.numeric(scale(edge_end)))
  
  # return
  return(all.steps)
  
}

#_______________________________________________________________________
# 6b. Use function ----
#_______________________________________________________________________

SSF.B1 <- sample_extract_SSF(steps.B1, sl.dist.B1, landscape.covs.B1)
SSF.B2 <- sample_extract_SSF(steps.B2, sl.dist.B2, landscape.covs.B2)
SSF.B3 <- sample_extract_SSF(steps.B3, sl.dist.B3, landscape.covs.B3)
SSF.A1 <- sample_extract_SSF(steps.A1, sl.dist.A1, landscape.covs.A1)
SSF.A2 <- sample_extract_SSF(steps.A2, sl.dist.A2, landscape.covs.A2)
SSF.A3 <- sample_extract_SSF(steps.A3, sl.dist.A3, landscape.covs.A3)

#_______________________________________________________________________
# 7. Means and SDs of continuous variables ----
#_______________________________________________________________________
# 7a. Define function ----
#_______________________________________________________________________

extract_mean_sd <- function(sampled.steps,
                            id.trt,
                            id.rep,
                            id.model) {
  
  if (id.model == "RSF") {
    
    return(data.frame(mean.forage = mean(sampled.steps$forage),
                    sd.forage = sd(sampled.steps$forage),
                    mean.edge = mean(sampled.steps$edge),
                    sd.edge = sd(sampled.steps$edge),
                    rep = id.rep,
                    model = id.model))
    
  } else {
    
    return(data.frame(mean.forage = mean(sampled.steps$forage_end),
                      sd.forage = sd(sampled.steps$forage_end),
                      mean.edge = mean(sampled.steps$edge_end),
                      sd.edge = sd(sampled.steps$edge_end),
                      rep = id.rep,
                      model = id.model))
    
  }
  
}

#_______________________________________________________________________
# 7b. Use function ----
#_______________________________________________________________________

mean.sd.all <- rbind(#RSF
                     extract_mean_sd(RSF.B1, "before", 1, "RSF"),
                     extract_mean_sd(RSF.B2, "before", 2, "RSF"),
                     extract_mean_sd(RSF.B3, "before", 3, "RSF"),
                     extract_mean_sd(RSF.A1, "after", 1, "RSF"),
                     extract_mean_sd(RSF.A2, "after", 2, "RSF"),
                     extract_mean_sd(RSF.A3, "after", 3, "RSF"),
                     
                     # SSF
                     extract_mean_sd(SSF.B1, "before", 1, "SSF"),
                     extract_mean_sd(SSF.B2, "before", 2, "SSF"),
                     extract_mean_sd(SSF.B3, "before", 3, "SSF"),
                     extract_mean_sd(SSF.A1, "after", 1, "SSF"),
                     extract_mean_sd(SSF.A2, "after", 2, "SSF"),
                     extract_mean_sd(SSF.A3, "after", 3, "SSF"))

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

#_______________________________________________________________________
# 8b. Fit models ----
#_______________________________________________________________________

RSF.model.B1 <- fit_RSF(RSF.B1)
RSF.model.B2 <- fit_RSF(RSF.B2)
RSF.model.B3 <- fit_RSF(RSF.B3)
RSF.model.A1 <- fit_RSF(RSF.A1)
RSF.model.A2 <- fit_RSF(RSF.A2)
RSF.model.A3 <- fit_RSF(RSF.A3)

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
                                open_start +
                                open_start:log(sl_) +
                                log(sl_) +
                                (0 + forage.s | indiv) +
                                (0 + forage.s:log(sl_) | indiv) +
                                (0 + edge.s | indiv) +
                                (0 + open_start | indiv) +
                                (0 + open_start:log(sl_) | indiv) +
                                (0 + log(sl_) | indiv) +
                                (1 | step_id_),
                              family = poisson,
                              data = sampled.steps,
                              doFit = FALSE) 
  
  iSSF.struc$parameters$theta[4] <- log(1e3)
  iSSF.struc$mapArg <- list(theta = factor(c(1:6, NA)))
  
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

iSSF.model.B1 <- fit_issf(SSF.B1)
iSSF.model.B2 <- fit_issf(SSF.B2)
iSSF.model.B3 <- fit_issf(SSF.B3)
iSSF.model.A1 <- fit_issf(SSF.A1)
iSSF.model.A2 <- fit_issf(SSF.A2)
iSSF.model.A3 <- fit_issf(SSF.A3)

#_______________________________________________________________________
# 10. Extract and bind HS coefficients together ----
#_______________________________________________________________________
# 10a. Define function ----
#_______________________________________________________________________

extract_coefs <- function(model,
                          id.type,
                          id.trt,
                          id.rep) {
  
  # ask which type 
  # RSF
  if (id.type == "RSF") {
    
    # use tidy function
    tidy.model <- broom.mixed::tidy(model) %>%
      
      # add id columns
      mutate(trt = id.trt,
             rep = id.rep,
             type = id.type) %>%
      
      # remove first three columns
      dplyr::select(-c(effect, component, group)) %>%
      
      # drop sd__(Intercept) term
      filter(term != "sd__(Intercept)")
    
  # iSSF
  } else {
    
    # use tidy function
    tidy.model <- broom.mixed::tidy(model) %>%
      
      # add id columns
      mutate(trt = id.trt,
             rep = id.rep,
             type = id.type) %>%
      
      # remove first three columns
      dplyr::select(-c(effect, component, group)) %>%
      
      # drop sd__(Intercept) term
      filter(term != "sd__(Intercept)")
  
  }
  
    # return
    return(tidy.model)
  
}

#_______________________________________________________________________
# 10b. Use function and bind ----
#_______________________________________________________________________

all.tidy <- rbind(extract_coefs(RSF.model.B1, "RSF", "before", 1),
                  extract_coefs(RSF.model.B2, "RSF", "before", 2),
                  extract_coefs(RSF.model.B3, "RSF", "before", 3),
                  extract_coefs(RSF.model.A1, "RSF", "after", 1),
                  extract_coefs(RSF.model.A2, "RSF", "after", 2),
                  extract_coefs(RSF.model.A3, "RSF", "after", 3),
                  extract_coefs(iSSF.model.B1, "iSSF", "before", 1),
                  extract_coefs(iSSF.model.B2, "iSSF", "before", 2),
                  extract_coefs(iSSF.model.B3, "iSSF", "before", 3),
                  extract_coefs(iSSF.model.A1, "iSSF", "after", 1),
                  extract_coefs(iSSF.model.A2, "iSSF", "after", 2),
                  extract_coefs(iSSF.model.A3, "iSSF", "after", 3))

# drop intercept term
all.tidy <- all.tidy %>% filter(term != "(Intercept)")

# change open_start to open
all.tidy$term[all.tidy$term == "open_start"] <- "open"
all.tidy$term[all.tidy$term == "sd__open_start"] <- "sd__open"

#_______________________________________________________________________
# 10c. Plot betas ----
#_______________________________________________________________________

# subset for plotting
all.tidy.beta <- all.tidy %>% 
  
  filter(term %in% unique(all.tidy$term)[c(1:3, 7:9)]) %>%
  
  # reorder and label factor
  mutate(term = factor(term,
                       levels = c("forage.s", "forage.s:log(sl_)", "edge.s",
                                  "open", "open_start:log(sl_)", "log(sl_)"),
                       labels = c("forage", "forage:log(sl)", "edge",
                                  "open", "open:log(sl)", "log(sl)")))

# plot
ggplot(all.tidy.beta,
       aes(y = term,
           x = estimate,
           color = type,
           group = rep)) +
  
  facet_grid(type ~ trt,
             scales = "free_y") +
  
  geom_vline(xintercept = 0) +
  
  theme_bw() +
  
  # error
  geom_errorbarh(aes(y = term,
                     xmin = estimate - (std.error * 1.96),
                     xmax = estimate + (std.error * 1.96),
                     color = type,
                     group = rep),
                 height = 0,
                 position = position_dodge(width = 0.5)) +
  
  geom_point(size = 2,
             position = position_dodge(width = 0.5)) +
  
  theme(panel.grid = element_blank())

#_______________________________________________________________________
# 10d. Plot RE SDs ----
#_______________________________________________________________________

# subset for plotting
all.tidy.sd <- all.tidy %>% 
  
  filter(term %in% unique(all.tidy$term)[c(4:6, 10:12)]) %>%
  
  # reorder and label factor
  mutate(term = factor(term,
                       levels = c("sd__forage.s", "sd__forage.s:log(sl_)", "sd__edge.s",
                                  "sd__open", "sd__open_start:log(sl_)", "sd__log(sl_)"),
                       labels = c("forage", "forage:log(sl)", "edge",
                                  "open", "open:log(sl)", "log(sl)")))

# plot
ggplot(all.tidy.sd,
       aes(y = trt,
           x = estimate,
           fill = type,
           group = rep)) +
  
  facet_grid(type ~ term,
             scales = "free_x") +
  
  geom_vline(xintercept = 0) +
  
  theme_bw() +
  
  # error
  geom_col(position = position_dodge(width = 1.0),
           color = "black") +
  
  theme(panel.grid = element_blank())

# seems like the models struggled to recover variance on all parameters
# makes sense since we simulated without individual variability!

#_______________________________________________________________________
# 11. Extract and store variance-covariance matrices ----
#_______________________________________________________________________

# define function
extract_vcov <- function(model,
                         id.type,
                         id.trt,
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
all.vcov <- list(extract_vcov(RSF.model.B1, "RSF", "before", 1),
                 extract_vcov(RSF.model.B2, "RSF", "before", 2),
                 extract_vcov(RSF.model.B3, "RSF", "before", 3),
                 extract_vcov(RSF.model.A1, "RSF", "after", 1),
                 extract_vcov(RSF.model.A2, "RSF", "after", 2),
                 extract_vcov(RSF.model.A3, "RSF", "after", 3),
                 extract_vcov(iSSF.model.B1, "iSSF", "before", 1),
                 extract_vcov(iSSF.model.B2, "iSSF", "before", 2),
                 extract_vcov(iSSF.model.B3, "iSSF", "before", 3),
                 extract_vcov(iSSF.model.A1, "iSSF", "after", 1),
                 extract_vcov(iSSF.model.A2, "iSSF", "after", 2),
                 extract_vcov(iSSF.model.A3, "iSSF", "after", 3)) 

# save to file
save(all.vcov, file = paste0(getwd(), "/Derived_data/Model parameters/vcov.RData"))

# write lookup table to csv
vcov.lookup <- cbind(expand.grid(rep = 1:3,
                                 trt = c("before", "after"),
                                 type = c("RSF", "iSSF")),
                     index = 1:12)

write.csv(vcov.lookup, paste0(getwd(), "/Derived_data/Lookup/vcov.csv"))

#_______________________________________________________________________
# 12. Write to .csv ----
#_______________________________________________________________________

write.csv(all.tidy, paste0(getwd(), "/Derived_data/Model parameters/all_params.csv"))
write.csv(mean.sd.all, paste0(getwd(), "/Derived_data/Model parameters/mean_sd_covs.csv"))
