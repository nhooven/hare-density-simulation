# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 09a - Testing UD simulations
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 10 Dec 2024
# Date completed: 
# Date last modified: 10 Dec 2024
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation
library(sf)              # read in shapefiles
library(terra)           # rasters
library(amt)             # simulate tracks
library(lubridate)       # work with time

#_______________________________________________________________________
# 2. Read in rasters and unit boundary ----
#_______________________________________________________________________

landscape.covs.simple <- rast("Rasters/simple.tif")
landscape.covs.complex <- rast("Rasters/complex.tif")

unit.bound <- st_read(paste0(getwd(), "/Derived_data/Shapefiles/unit_bound.shp"))

#_______________________________________________________________________
# 3. Define simulation parameters ----
#_______________________________________________________________________
# 3a. Movement parameter distributions ----
#_______________________________________________________________________

# we'll just compare the simple-low and complex-low scenarios here

# step lengths (gamma)
load(file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_SL.RData"))
load(file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_CL.RData"))

# turning angle
load(file = paste0(getwd(), "/Derived_data/Model parameters/ta_dist.RData"))

#_______________________________________________________________________
# 3b. Habitat selection coefficients and identifiers ----
#_______________________________________________________________________

model.params <- read.csv(paste0(getwd(), "/Derived_data/Model parameters/all_params.csv"))

# keep only iSSF parameters
model.params <- model.params %>% filter(type == "iSSF")

#_______________________________________________________________________
# 3c. Initialization space ----
#_______________________________________________________________________

# extract bounding box
unit.bbox <- st_bbox(unit.bound)

# we'll take a random draw from this area to initialize each path

#_______________________________________________________________________
# 3d. Redistribution kernel parameters ----
#_______________________________________________________________________

# control steps
rk.control <- 10000

# tolerance outside the landscape (hopefully we won't have to deal with this much)
rk.tolerance <- 0.01    # 1%

#_______________________________________________________________________
# 4. Define functions needed for simulation ----
#_______________________________________________________________________
# 4a. Start step ----
#_______________________________________________________________________

# we'll borrow the home range centroid function from original simulations

make_hrc <- function() {
  
  hrc.x <- runif(n = 1,
                 min = unit.bbox[1],
                 max = unit.bbox[3])
  
  hrc.y <- runif(n = 1,
                 min = unit.bbox[2],
                 max = unit.bbox[4])
  
  # concatenate
  hrc <- c(hrc.x, hrc.y)
  
  # return
  return(hrc)
  
}

#_______________________________________________________________________
# 5. Run simulations iteratively ----

# n of replicates
n.reps <- 300

# n of time steps (keep this at a month)
n.steps <- 336

# 11 Dec 2024
# Following Potts and Borger (2023), let's run a few hundred sims from the same
# starting location, take the end points, and smooth with a tradition KDE

#_______________________________________________________________________
# 5a. Define functions ----
#_______________________________________________________________________

# transient UDs
sim_issf_tud <- function (ls,
                         var) {
  
  # extract fit parameters
  focal.params <- model.params %>%
    
    filter(landscape == ls &
             variability == var) %>%
    
    dplyr::select(names, betas)
  
  # extract correct sl distribution
  if (ls == "simple" & var == "low") {
    
    focal.sl.dist <- sl.dist.SL
    
  } else {
    
    if (ls == "complex" & var == "low") {
      
      focal.sl.dist <- sl.dist.CL
      
    } else {
      
      if (ls == "simple" & var == "high") {
        
        focal.sl.dist <- sl.dist.SH
        
      } else {
        
        focal.sl.dist <- sl.dist.CH
        
      }
      
      
    }
    
  }
  
  # use correct raster stack
  if (ls == "simple") {
    
    landscape.covs <- landscape.covs.simple
    
  } else {
    
    landscape.covs <- landscape.covs.complex
    
  }
  
  # create df to hold all sims
  sims.df <- data.frame()
  
  # run simulations
  start.time <- Sys.time()
  
  for (i in 1:n.reps) {
    
    # define hrc
    hrc <- make_hrc()
    
    # define start step
    start.step <- make_start(x = c(hrc[1],
                                   hrc[2]),
                             ta_ = 0,
                             time = ymd_hm("2024-09-01 18:00", 
                                           tz = "America/Los_Angeles"),
                             dt = hours(2),
                             crs = crs("EPSG:32611"))
    
    # make iSSF model
    # here the terms are important to get right so redistribution_kernel() works okay
    issf.model <- make_issf_model(coefs = c("stem_end" = focal.params$betas[focal.params$names == "stem"],  
                                            "stem_end:log(sl_)" = focal.params$betas[focal.params$names == "stem_sl"],
                                            "edge_end" = focal.params$betas[focal.params$names == "edge"],
                                            "mature_start" = focal.params$betas[focal.params$names == "mature"],
                                            "log(sl_):mature_start" = focal.params$betas[focal.params$names == "mature_sl"],
                                            "log(sl_)" = focal.params$betas[focal.params$names == "log_sl"]),              
                                  sl = focal.sl.dist,
                                  ta = ta.dist)
    
    # initialize redistribution kernel
    rk <- redistribution_kernel(x = issf.model,
                                start = start.step,
                                map = landscape.covs,
                                n.control = rk.control,   
                                max.dist = get_max_dist(issf.model),
                                tolerance.outside = rk.tolerance)
    
    # run simulation
    sim.path <- simulate_path(rk,
                              n.steps = n.steps,
                              start = start.step,
                              verbose = TRUE)
    
    # extract final location for fitting uncorrelated UD
    sim.final <- sim.path %>%
      
      # add in individual identifiers (for step so we can extract them later)
      mutate(step = rownames(sim.path),
             indiv = i)
    
    # bind to df
    sims.df <- bind_rows(sims.df, sim.final)
    
    # status message
    elapsed.time <- round(as.numeric(difftime(Sys.time(), 
                                              start.time, 
                                              units = "mins")), 
                          digits = 1)
    
    print(paste0("Completed path ", i, " of ", n.reps, " - ", elapsed.time, " mins"))
    
  }
  
  # return
  return(sims.df)
  
}

# steady-state UD
sim_issf_ssud <- function (ls,
                           var,
                           n.steps.ssud) {
  
  # start time
  start.time <- Sys.time()
  
  # extract fit parameters
  focal.params <- model.params %>%
    
    filter(landscape == ls &
           variability == var) %>%
    
    dplyr::select(names, betas)
  
  # extract correct sl distribution
  if (ls == "simple" & var == "low") {
    
    focal.sl.dist <- sl.dist.SL
    
  } else {
    
    if (ls == "complex" & var == "low") {
      
      focal.sl.dist <- sl.dist.CL
      
    } else {
      
      if (ls == "simple" & var == "high") {
        
        focal.sl.dist <- sl.dist.SH
        
      } else {
        
        focal.sl.dist <- sl.dist.CH
        
      }
      
      
    }
    
  }
  
  # use correct raster stack
  if (ls == "simple") {
    
    landscape.covs <- landscape.covs.simple
    
  } else {
    
    landscape.covs <- landscape.covs.complex
    
  }
  
  # run simulations
  start.time <- Sys.time()
    
  # define start step (center of unit, SSUD should not be sensitive to initial conditions)
  start.step <- make_start(x = c(unit.bbox[1] + (unit.bbox[3] - unit.bbox[1]) / 2,
                                 unit.bbox[2] + (unit.bbox[4] - unit.bbox[2]) / 2),
                           ta_ = 0,
                           time = ymd_hm("2024-09-01 18:00", 
                                         tz = "America/Los_Angeles"),
                           dt = hours(2),
                           crs = crs("EPSG:32611"))
  
  # make iSSF model
  # here the terms are important to get right so redistribution_kernel() works okay
  issf.model <- make_issf_model(coefs = c("stem_end" = focal.params$betas[focal.params$names == "stem"],  
                                          "stem_end:log(sl_)" = focal.params$betas[focal.params$names == "stem_sl"],
                                          "edge_end" = focal.params$betas[focal.params$names == "edge"],
                                          "mature_start" = focal.params$betas[focal.params$names == "mature"],
                                          "log(sl_):mature_start" = focal.params$betas[focal.params$names == "mature_sl"],
                                          "log(sl_)" = focal.params$betas[focal.params$names == "log_sl"]),              
                                sl = focal.sl.dist,
                                ta = ta.dist)
  
  # initialize redistribution kernel
  rk <- redistribution_kernel(x = issf.model,
                              start = start.step,
                              map = landscape.covs,
                              n.control = rk.control,   
                              max.dist = get_max_dist(issf.model),
                              tolerance.outside = 0.05)         # 5%
  
  # run simulation
  sim.path <- simulate_path(rk,
                            n.steps = n.steps.ssud,
                            start = start.step,
                            verbose = TRUE)
  
  # extract final location for fitting uncorrelated UD
  sim.final <- sim.path %>%
    
    # add in individual identifiers (for step so we can extract them later)
    mutate(step = rownames(sim.path))
  
  # return
  return(sim.final)
  
  # print completed time
  elapsed.time <- round(as.numeric(difftime(Sys.time(), 
                                            start.time, 
                                            units = "mins")), 
                        digits = 1)
  
  print(paste0("Completed after ", elapsed.time, " mins"))
  
}

#_______________________________________________________________________
# 5b. Run simulations ----
#_______________________________________________________________________

sims.SL.100 <- sim_issf_ud("simple", "low")
sims.CL.100 <- sim_issf_ud("complex", "low")

sims.CL.200 <- sim_issf_ud("complex", "low")

# SSUD
sims.SL.ssud <- sim_issf_ssud("simple", "low", 100000)
sims.CL.ssud <- sim_issf_ssud("complex", "low", 100000)

# had issues with animal stepping outside the landscape
# maybe a large n with the transient UDs makes more sense? 
# a large n of uncorrelated points, smooth with a kernel, call it good?

#_______________________________________________________________________
# 6. Plot 2D density maps at - track length ----
#_______________________________________________________________________
# 6a. Simple - TUD ----
#_______________________________________________________________________

# subset df
sims.SL.subset <- sims.SL %>% 
  
  filter(step %in% c(50, 100, 150, 200, 250, 300, 337)) %>%
  
  mutate(step = factor(step,
                       levels = c(50, 100, 150, 200, 250, 300, 337)))

# subset df (200 - 500)
sims.CL.subset.200 <- sims.CL.200 %>% 
  
  filter(step %in% c(50, 100, 150, 200, 250, 300, 350, 400, 450, 500)) %>%
  
  mutate(step = factor(step,
                       levels = c(50, 100, 150, 200, 250, 300, 350, 400, 450, 500)))

# facetted density plot
ggplot(sims.CL.subset.200,
       aes(x = x_,
           y = y_)) +
  
  facet_wrap(~ step, nrow = 2) +
  
  theme_bw() +
  
  stat_density2d_filled(aes(fill = after_stat(level))) +
  
  scale_fill_viridis_d(option = "magma") +
  
  geom_point(shape = 21,
             color = "white",
             alpha = 0.15) +
  
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank())

# 10 Dec 2024
# This is showing me that the spatial pattern stabilizes somewhat quickly (over a hundred steps)
# but the shape of the UD is pretty sensitive to individual points.
# Thus, it is probably more important to simulate more tracks than run them for longer. 

#_______________________________________________________________________
# 6b. Complex - TUD ----
#_______________________________________________________________________

# subset df
sims.CL.subset <- sims.CL.100 %>% 
  
  filter(step %in% c(50, 100, 150, 200, 250, 300)) %>%
  
  mutate(step = factor(step,
                       levels = c(50, 100, 150, 200, 250, 300)))

# facetted density plot
ggplot(sims.CL.subset,
       aes(x = x_,
           y = y_)) +
  
  facet_wrap(~ step, nrow = 1) +
  
  theme_bw() +
  
  stat_density2d_filled(aes(fill = after_stat(level))) +
  
  scale_fill_viridis_d(option = "magma") +
  
  geom_point(shape = 21,
             color = "white") +
  
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank())

# this is interesting - we're not seeing clear stabilization even up to 300 steps
# should we crank up the total steps?

# maybe it's worth doing a longer simulation just to make sure we're characterizing
# this UD well. 

#_______________________________________________________________________
# 7. Plot 2D density maps at - sample size ----

# define function to create plottable df
sim_UD_sample_size <- function (df, n) {
  
  x.1 <- df %>% 
    
    filter(step == n.steps) %>%
    
    slice(1:n) %>%
    
    mutate(samp.size = n)
  
  return (x.1)
  
}

#_______________________________________________________________________
# 7a. Simple ----
#_______________________________________________________________________

# subset df
sims.SL.subset.1 <- rbind(sim_UD_sample_size(sims.SL.100, 10),
                          sim_UD_sample_size(sims.SL.100, 20),
                          sim_UD_sample_size(sims.SL.100, 30),
                          sim_UD_sample_size(sims.SL.100, 40),
                          sim_UD_sample_size(sims.SL.100, 50),
                          sim_UD_sample_size(sims.SL.100, 60),
                          sim_UD_sample_size(sims.SL.100, 70),
                          sim_UD_sample_size(sims.SL.100, 80),
                          sim_UD_sample_size(sims.SL.100, 90),
                          sim_UD_sample_size(sims.SL.100, 100))

# factor levels
sims.SL.subset.1$samp.size <- factor(sims.SL.subset.1$samp.size,
                                     levels = seq(10, 100, by = 10))


# facetted density plot
ggplot(sims.SL.subset.1,
       aes(x = x_,
           y = y_)) +
  
  facet_wrap(~ samp.size, nrow = 2) +
  
  theme_bw() +
  
  stat_density2d_filled(aes(fill = after_stat(level))) +
  
  scale_fill_viridis_d(option = "magma") +
  
  geom_point(shape = 21,
             color = "white") +
  
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank())

# yes! exactly what I wanted. Here we're not getting much more after 50 individual tracks. 

#_______________________________________________________________________
# 7b. Complex ----
#_______________________________________________________________________

# subset df
sims.CL.subset.1 <- rbind(sim_UD_sample_size(sims.CL.100, 10),
                          sim_UD_sample_size(sims.CL.100, 20),
                          sim_UD_sample_size(sims.CL.100, 30),
                          sim_UD_sample_size(sims.CL.100, 40),
                          sim_UD_sample_size(sims.CL.100, 50),
                          sim_UD_sample_size(sims.CL.100, 60),
                          sim_UD_sample_size(sims.CL.100, 70),
                          sim_UD_sample_size(sims.CL.100, 80),
                          sim_UD_sample_size(sims.CL.100, 90),
                          sim_UD_sample_size(sims.CL.100, 100))

# factor levels
sims.CL.subset.1$samp.size <- factor(sims.CL.subset.1$samp.size,
                                     levels = seq(10, 100, by = 10))


# facetted density plot
ggplot(sims.CL.subset.1,
       aes(x = x_,
           y = y_)) +
  
  facet_wrap(~ samp.size, nrow = 2) +
  
  theme_bw() +
  
  stat_density2d_filled(aes(fill = after_stat(level))) +
  
  scale_fill_viridis_d(option = "magma") +
  
  geom_point(shape = 21,
             color = "white") +
  
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank())

# it looks like the sample size required to recover a sensible UD in the complex landscape
# might be higher - this is no shock!

#_______________________________________________________________________
# 8. SSUD convergence ----

# define function to create plottable df
sim_UD_ssud_sample_size <- function (df, n) {
  
  x.1 <- df %>% 
    
    slice(1:n) %>%
    
    mutate(samp.size = n)
  
  return (x.1)
  
}

#_______________________________________________________________________
# 8a. Simple - SSUD ----
#_______________________________________________________________________

# subset df
sims.SL.ssud.subset.1 <- rbind(sim_UD_ssud_sample_size(sims.SL.ssud, 1000),
                               sim_UD_ssud_sample_size(sims.SL.ssud, 2000),
                               sim_UD_ssud_sample_size(sims.SL.ssud, 3000),
                               sim_UD_ssud_sample_size(sims.SL.ssud, 4000),
                               sim_UD_ssud_sample_size(sims.SL.ssud, 5000),
                               sim_UD_ssud_sample_size(sims.SL.ssud, 6000),
                               sim_UD_ssud_sample_size(sims.SL.ssud, 7000),
                               sim_UD_ssud_sample_size(sims.SL.ssud, 8000),
                               sim_UD_ssud_sample_size(sims.SL.ssud, 9000),
                               sim_UD_ssud_sample_size(sims.SL.ssud, 10000))

# factor levels
sims.SL.ssud.subset.1$samp.size <- factor(sims.SL.ssud.subset.1$samp.size,
                                          levels = seq(1000, 10000, by = 1000))


# facetted density plot
ggplot(sims.SL.ssud.subset.1,
       aes(x = x_,
           y = y_)) +
  
  facet_wrap(~ samp.size, nrow = 2) +
  
  theme_bw() +
  
  geom_hex() +
  
  scale_fill_viridis_c(option = "magma") +
  
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank())

# points within the camera grid
ggplot(sims.SL.ssud %>% slice(1001:10000),
       aes(x = x_,
           y = y_)) +
  
  theme_bw() +
  
  geom_point(size = 0.5) +
  
  geom_rect(aes(xmin = unit.bbox[1],
                ymin = unit.bbox[2],
                xmax = unit.bbox[3],
                ymax = unit.bbox[4]),
            color = "green",
            fill = NA,
            linewidth = 1.1) +
  
  coord_cartesian(xlim = c(900, 1300),
                  ylim = c(900, 1300))

# I'm thinking we should run this for 100,000 iterations

# let's look at the raster
library(terra)

blank.rast <- rast(nrows = nrow(landscape.covs.simple),
                   ncols = ncol(landscape.covs.simple),
                   xmin = xmin(landscape.covs.simple),
                   ymin = ymin(landscape.covs.simple),
                   xmax = xmax(landscape.covs.simple),
                   ymax = ymax(landscape.covs.simple),
                   res = 10,
                   crs = crs(landscape.covs.simple))

point.rast <- rasterize(matrix(c(sims.SL.ssud$x_,
                                 sims.SL.ssud$y_),
                               ncol = 2),
                        y = blank.rast,
                        fun = "count")

point.rast

plot(point.rast)

# we obviously can't use the native resolution for this UD raster
# but 10 looks decent, it just needs more samples

# what if we normalize it?
plot(point.rast / sum(values(point.rast), na.rm = T))

# now this UD should sum to 1, which is what we want