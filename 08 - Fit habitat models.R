# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 08 - Fit habitat models
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 09 Dec 2024
# Date completed: 09 Dec 2024
# Date last modified: 21 Feb 2025
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation
library(terra)           # rasters
library(amt)
library(sf)
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
ta.dist.B1 <- fit_distr(steps.B1$ta_, dist_name = "vonmises")
ta.dist.B2 <- fit_distr(steps.B2$ta_, dist_name = "vonmises")
ta.dist.B3 <- fit_distr(steps.B3$ta_, dist_name = "vonmises")

ta.dist.A1 <- fit_distr(steps.A1$ta_, dist_name = "vonmises")
ta.dist.A2 <- fit_distr(steps.A2$ta_, dist_name = "vonmises")
ta.dist.A3 <- fit_distr(steps.A3$ta_, dist_name = "vonmises")

#_______________________________________________________________________
# 4b. Extract parameters ----
#_______________________________________________________________________

sl.dist.params <- data.frame(trt = c("before", "before", "before", "before", "before", "before",
                                     "after", "after", "after", "after", "after", "after"),
                             rep = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3),
                             shape = c(sl.dist.B1$params$shape,
                                       sl.dist.B2$params$shape,
                                       sl.dist.B3$params$shape,
                                       sl.dist.A1$params$shape,
                                       sl.dist.A2$params$shape,
                                       sl.dist.A3$params$shape),
                             scale = c(sl.dist.B1$params$scale,
                                       sl.dist.B2$params$scale,
                                       sl.dist.B3$params$scale,
                                       sl.dist.A1$params$scale,
                                       sl.dist.A2$params$scale,
                                       sl.dist.A3$params$scale))

ta.dist.params <- data.frame(trt = c("before", "before", "before", "before", "before", "before",
                                     "after", "after", "after", "after", "after", "after"),
                             rep = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3),
                             kappa = c(ta.dist.B1$params$kappa,
                                       ta.dist.B2$params$kappa,
                                       ta.dist.B3$params$kappa,
                                       ta.dist.A1$params$kappa,
                                       ta.dist.A2$params$kappa,
                                       ta.dist.A3$params$kappa))

#_______________________________________________________________________
# 4c. Save to lookup ----
#_______________________________________________________________________

write.csv(sl.dist.params, paste0(getwd(), "/Derived_data/Lookup/sl_params.csv"))
write.csv(ta.dist.params, paste0(getwd(), "/Derived_data/Lookup/ta_params.csv"))

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
    focal.fora <- terra::extract(landscape$fora, all.sf, ID = FALSE)
    focal.elev <- terra::extract(landscape$elev, all.sf, ID = FALSE)
    focal.open <- terra::extract(landscape$open, all.sf, ID = FALSE)
    
    # and bind in
    all.sf.1 <- cbind(all.sf, focal.fora, focal.elev, focal.open)

    # bind together
    RSF.samples <- as.data.frame(rbind(RSF.samples, all.sf.1))
    
  }
  
  # scale all continuous covariates and add weights for RSF
  RSF.samples <- RSF.samples %>%
    
    mutate(fora.s = as.numeric(scale(fora)),
           elev.s = as.numeric(scale(elev)),
           open.s = as.numeric(scale(open)),
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
# 6. iSSF - Sample random steps and extract covariates ----
#_______________________________________________________________________
# 6a. Define function ----
#_______________________________________________________________________

sample_extract_SSF <- function (steps.df,
                                sl.dist,
                                ta.dist,
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
    
    mutate(fora.end.s = as.numeric(scale(fora_end)),
           fora.start.s = as.numeric(scale(fora_start)),
           elev.s = as.numeric(scale(elev_end)),
           open.start.s = as.numeric(scale(open_start)),
           open.end.s = as.numeric(scale(open_end))) %>%
    
    # create unique stratum
    mutate(stratum = paste0(indiv, step_id_))
  
  # return
  return(all.steps)
  
}

#_______________________________________________________________________
# 6b. Use function ----
#_______________________________________________________________________

SSF.B1 <- sample_extract_SSF(steps.B1, sl.dist.B1, ta.dist.B1, landscape.covs.B1)
SSF.B2 <- sample_extract_SSF(steps.B2, sl.dist.B2, ta.dist.B2, landscape.covs.B2)
SSF.B3 <- sample_extract_SSF(steps.B3, sl.dist.B3, ta.dist.B3, landscape.covs.B3)
SSF.A1 <- sample_extract_SSF(steps.A1, sl.dist.A1, ta.dist.A1, landscape.covs.A1)
SSF.A2 <- sample_extract_SSF(steps.A2, sl.dist.A2, ta.dist.A2, landscape.covs.A2)
SSF.A3 <- sample_extract_SSF(steps.A3, sl.dist.A3, ta.dist.A3, landscape.covs.A3)

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
    
    return(data.frame(mean.fora = mean(sampled.steps$fora),
                      sd.fora = sd(sampled.steps$fora),
                      mean.elev = mean(sampled.steps$elev),
                      sd.elev = sd(sampled.steps$elev),
                      mean.open = mean(sampled.steps$open),
                      sd.open = sd(sampled.steps$open),
                      rep = id.rep,
                      model = id.model))
    
  } else {
    
    return(data.frame(mean.fora.start = mean(sampled.steps$fora_start),
                      sd.fora.start = sd(sampled.steps$fora_start),
                      mean.fora.end = mean(sampled.steps$fora_end),
                      sd.fora.end = sd(sampled.steps$fora_end),
                      mean.elev = mean(sampled.steps$elev_end),
                      sd.elev = sd(sampled.steps$elev_end),
                      mean.open.start = mean(sampled.steps$open_start),
                      sd.open.start = sd(sampled.steps$open_start),
                      mean.open.end = mean(sampled.steps$open_end),
                      sd.open.end = sd(sampled.steps$open_end),
                      rep = id.rep,
                      model = id.model))
    
  }
  
}

#_______________________________________________________________________
# 7b. Use function ----
#_______________________________________________________________________

mean.sd.rsf <- rbind(extract_mean_sd(RSF.B1, "before", 1, "RSF"),
                     extract_mean_sd(RSF.B2, "before", 2, "RSF"),
                     extract_mean_sd(RSF.B3, "before", 3, "RSF"),
                     extract_mean_sd(RSF.A1, "after", 1, "RSF"),
                     extract_mean_sd(RSF.A2, "after", 2, "RSF"),
                     extract_mean_sd(RSF.A3, "after", 3, "RSF"))
                     
mean.sd.ssf <- rbind(extract_mean_sd(SSF.B1, "before", 1, "SSF"),
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
  
  RSF.struc <- glmmTMB(case ~ fora.s +
                              elev.s +
                              I(elev.s^2) +
                              open.s +
                              (1 | indiv) +
                              (0 + fora.s | indiv) +
                              (0 + elev.s | indiv) +
                              (0 + I(elev.s^2) | indiv) +
                              (0 + open.s | indiv),
                             weights = weight,
                             family = binomial,
                             data = sampled.points,
                             doFit = FALSE) 
  
  RSF.struc$parameters$theta[1] <- log(1e3)
  RSF.struc$mapArg <- list(theta = factor(c(NA, 1:4)))     # account for two slopes for elevation
  
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
issf_fit <- function(sampled.steps) {
  
  iSSF.struc <- glmmTMB(case_ ~ fora.start.s:log(sl_) +
                                fora.end.s +
                                fora.end.s:cos(ta_) +
                                elev.s +
                                I(elev.s^2) +
                                open.start.s:log(sl_)  +
                                open.end.s +
                                log(sl_) +
                                cos(ta_) +
                                (1 | stratum) +                 # in the interest of convergence
                                (0 + fora.end.s | indiv) +      # we'll only include RS for base coefs
                                (0 + elev.s | indiv) +
                                (0 + I(elev.s^2) | indiv) +
                                (0 + open.end.s | indiv),
                              family = poisson,
                              data = sampled.steps,
                              doFit = FALSE) 
  
  iSSF.struc$parameters$theta[1] <- log(1e3)
  iSSF.struc$mapArg <- list(theta = factor(c(NA, 1:4)))
  
  iSSF.model <- fitTMB(iSSF.struc)
  
  # return
  return(iSSF.model)
  
}

# no random slopes to fix convergence issues
issf_fit_1 <- function(sampled.steps) {
  
  iSSF.model <- amt::fit_issf(data = sampled.steps,
                              formula = case_ ~ fora.start.s:log(sl_) +
                                                fora.end.s +
                                                fora.end.s:cos(ta_) +
                                                elev.s +
                                                I(elev.s^2) +
                                                open.start.s:log(sl_)  +
                                                open.end.s +
                                                log(sl_) +
                                                cos(ta_) +
                                                strata(stratum))
  
  # return
  return(iSSF.model)
  
}

#_______________________________________________________________________
# 9b. Fit models ----
#_______________________________________________________________________

iSSF.model.B1 <- issf_fit_1(SSF.B1)
iSSF.model.B2 <- issf_fit(SSF.B2)
iSSF.model.B3 <- issf_fit_1(SSF.B3)
iSSF.model.A1 <- issf_fit_1(SSF.A1)
iSSF.model.A2 <- issf_fit_1(SSF.A2)
iSSF.model.A3 <- issf_fit_1(SSF.A3)

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
    
    if ("glmmTMB" %in% class(model)) {
      
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
      
    } else {
      
      # use tidy function
      tidy.model <- broom.mixed::tidy(model$model) %>%
        
        # add id columns
        mutate(trt = id.trt,
               rep = id.rep,
               type = id.type)
      
    }
  
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

#_______________________________________________________________________
# 10c. Plot betas ----
#_______________________________________________________________________

# I SHOULD ADD IN THE SIMULATION VALUES

# subset for plotting
all.tidy.beta <- all.tidy %>% 
  
  filter(term %in% unique(all.tidy$term)[c(1:4, 9:15)]) %>%
  
  # replace names
  mutate(term = recode(term, 
                       "fora.end.s" = "fora.s",
                       "open.end.s" = "open.s")) %>%
  
  # add original betas from simulation
  add_row(term = "fora.s",
          estimate = 1.5,
          std.error = 0,
          statistic = 0,
          p.value = 0,
          trt = c("before"),
          rep = NA,
          type = "sim") %>%
  
  add_row(term = "fora.end.s:cos(ta_)",
          estimate = -0.5,
          std.error = 0,
          statistic = 0,
          p.value = 0,
          trt = c("before"),
          rep = NA,
          type = "sim") %>%
  
  add_row(term = "log(sl_):fora.start.s",
          estimate = -0.05,
          std.error = 0,
          statistic = 0,
          p.value = 0,
          trt = c("before"),
          rep = NA,
          type = "sim") %>%
  
  add_row(term = "elev.s",
          estimate = 0.75,
          std.error = 0,
          statistic = 0,
          p.value = 0,
          trt = c("before"),
          rep = NA,
          type = "sim") %>%
  
  add_row(term = "I(elev.s^2)",
          estimate = -0.75,
          std.error = 0,
          statistic = 0,
          p.value = 0,
          trt = c("before"),
          rep = NA,
          type = "sim") %>%
  
  add_row(term = "log(sl_):open.start.s",
          estimate = 0.25,
          std.error = 0,
          statistic = 0,
          p.value = 0,
          trt = c("before"),
          rep = NA,
          type = "sim") %>%
  
  add_row(term = "open.s",
          estimate = -1.0,
          std.error = 0,
          statistic = 0,
          p.value = 0,
          trt = c("before"),
          rep = NA,
          type = "sim") %>%
  
  # reorder and label factor
  mutate(term = factor(term,
                       levels = rev(c("fora.s", 
                                  "fora.end.s:cos(ta_)",
                                  "log(sl_):fora.start.s", 
                                  "elev.s",
                                  "I(elev.s^2)", 
                                  "open.s", 
                                  "log(sl_):open.start.s", 
                                  "log(sl_)",
                                  "cos(ta_)")),
         
                        labels = rev(c("fora", 
                                       "fora:cos(ta)",
                                       "fora:log(sl)", 
                                       "elev",
                                       "elev^2", 
                                       "open", 
                                       "open:log(sl)", 
                                       "log(sl)",
                                       "cos(ta)"))))

# plot
ggplot(all.tidy.beta,
       aes(y = type,
           x = estimate,
           fill = type,
           group = rep,
           shape = trt)) +
  
  facet_wrap(~ term,
             scales = "free") +
  
  geom_vline(xintercept = 0) +
  
  theme_bw() +
  
  # error
  geom_errorbarh(aes(y = type,
                     xmin = estimate - (std.error * 1.96),
                     xmax = estimate + (std.error * 1.96),
                     color = type,
                     group = rep),
                 height = 0,
                 position = position_dodge(width = 1.0)) +
  
  geom_point(size = 2,
             position = position_dodge(width = 1.0)) +
  
  scale_shape_manual(values = c(21, 23)) +
  
  theme(panel.grid = element_blank())

# this looks terrible but it helps me figure things out

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
write.csv(mean.sd.rsf, paste0(getwd(), "/Derived_data/Model parameters/mean_sd_rsf.csv"))
write.csv(mean.sd.ssf, paste0(getwd(), "/Derived_data/Model parameters/mean_sd_ssf.csv"))
