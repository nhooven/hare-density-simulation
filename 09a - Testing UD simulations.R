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
n.reps <- 50

# n of time steps
n.steps <- 336

# 10 Dec 2024 
# 100 steps does not seem to lead to stability... hmm

#_______________________________________________________________________
# 5a. Define function ----
#_______________________________________________________________________

sim_issf_ud <- function (ls,
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

#_______________________________________________________________________
# 5b. Run simulations ----
#_______________________________________________________________________

sims.SL <- sim_issf_ud("simple", "low")
sims.CL <- sim_issf_ud("complex", "low")

#_______________________________________________________________________
# 6. Plot 2D density maps at certain cut points ----
#_______________________________________________________________________
# 6a. Simple ----
#_______________________________________________________________________

# subset df
sims.SL.subset <- sims.SL %>% 
  
  filter(step %in% c(50, 100, 150, 200, 250, 300, 337)) %>%
  
  mutate(step = factor(step,
                       levels = c(50, 100, 150, 200, 250, 300, 337)))

# facetted density plot
ggplot(sims.SL.subset,
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

# 10 Dec 2024
# This is showing me that the spatial pattern stabilizes somewhat quickly (over a hundred steps)
# but the shape of the UD is pretty sensitive to individual points.
# Thus, it is probably more important to simulate more tracks than run them for longer. 

