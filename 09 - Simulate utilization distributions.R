# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 09 - Simulate utilization distributions
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 09 Dec 2024
# Date completed: 
# Date last modified: 
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

# step lengths (gamma)
load(file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_SL.RData"))
load(file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_CL.RData"))
load(file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_SH.RData"))
load(file = paste0(getwd(), "/Derived_data/Model parameters/sl_dist_CH.RData"))

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

# we'll take a random draw from this distribution to initialize each path

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

# 09 Dec 2024
# let's run 50 tracks to see what that gives us

# n of replicates
n.reps <- 50

# calculate approximate time in hours (assuming 336 locations)
(42 * n.reps) / (60 * 60)

# WHERE I STOPPED - 09 DEC 2024

#_______________________________________________________________________
# 5a. Create data.frame to hold all sims ----
#_______________________________________________________________________

sims.df <- data.frame()

#_______________________________________________________________________
# 5b. Run simulations ----
#_______________________________________________________________________

start.time <- Sys.time()

for (i in 1:n.reps) {
  
  # define hrc
  hrc <- make_hrc()
  
  # define start step
  start.step <- make_start(x = c(hrc[1] + rnorm(n = 1, mean = 0, sd = 50),
                                 hrc[2] + rnorm(n = 1, mean = 0, sd = 50)),
                           ta_ = 0,
                           time = ymd_hm("2024-09-01 18:00", 
                                         tz = "America/Los_Angeles"),
                           dt = hours(2),
                           crs = crs("EPSG:32611"))
  
  # calculate home ranging parameters
  hr.params <- hr_params(e.var = e.var,
                         e.var.sd = e.var.sd,
                         hrc = hrc)
  
  # make iSSF model
  # here the terms are important to get right so redistribution_kernel() works okay
  issf.model <- make_issf_model(coefs = c("stem_end" = coef.stem,  
                                          "stem_end:log(sl_)" = coef.stem.sl,
                                          "edge_end" = coef.edge,
                                          "mature_start" = coef.mature,
                                          "log(sl_):mature_start" = coef.mature.sl,
                                          x2_ = hr.params[1],
                                          y2_ = hr.params[2], 
                                          "I(x2_^2 + y2_^2)" = hr.params[3]),              
                                sl = sl.dist,
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
                            n.steps = 336,
                            start = start.step,
                            verbose = TRUE)
  
  # add in identifiers
  sim.path.1 <- sim.path %>% 
    
    mutate(indiv = i,
           landscape = id.landscape,
           variability = id.variability)
  
  # bind to df
  sims.df <- bind_rows(sims.df, sim.path.1)
  
  # status message
  elapsed.time <- round(as.numeric(difftime(Sys.time(), 
                                            start.time, 
                                            units = "mins")), 
                        digits = 1)
  
  print(paste0("Completed path ", i, " of ", n.reps, " - ", elapsed.time, " mins"))
  
}

#_______________________________________________________________________
# 6. Plot tracks ----
#_______________________________________________________________________

sims.sf <- st_as_sf(sims.df,
                    coords = c("x_", 
                               "y_"),
                    crs = "epsg:32611") 

# plot paths
ggplot() +
  
  theme_bw() +
  
  # unit boundary as a black square
  geom_sf(data = unit.bound,
          fill = NA) +
  
  # paths
  geom_path(data = sims.df,
            aes(x = x_,
                y = y_,
                color = as.factor(indiv)),
            alpha = 0.5) +
  
  # remove legend
  theme(legend.position = "none") +
  
  # coordinate system to make nice axis labels
  coord_sf(datum = st_crs(32611))

#_______________________________________________________________________
# 7. Write to .csv ----
#_______________________________________________________________________

write.csv(sims.df, paste0(getwd(), "/Derived_data/Simulated data/sims_simple_low.csv"))
