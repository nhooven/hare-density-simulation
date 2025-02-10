# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 03 - Initial movement simulations
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 15 Nov 2024
# Date completed: 21 Nov 2024
# Date last modified: 10 Feb 2025
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

B1 <- rast(paste0(getwd(), "/Rasters/B1.tif"))
B2 <- rast(paste0(getwd(), "/Rasters/B2.tif"))
B3 <- rast(paste0(getwd(), "/Rasters/B3.tif"))

A1 <- rast(paste0(getwd(), "/Rasters/A1.tif"))
A2 <- rast(paste0(getwd(), "/Rasters/A2.tif"))
A3 <- rast(paste0(getwd(), "/Rasters/A3.tif"))

unit.bound <- st_read(paste0(getwd(), "/Derived_data/Shapefiles/unit_bound.shp"))

# scale rasters
landscape.covs.B <- list(c(scale(B1$forage),
                           scale(B1$edge),
                           B1$open),
                         c(scale(B2$forage),
                           scale(B2$edge),
                           B2$open),
                         c(scale(B3$forage),
                           scale(B3$edge),
                           B3$open))

landscape.covs.A <- list(c(scale(A1$forage),
                           scale(A1$edge),
                           A1$open),
                         c(scale(A2$forage),
                           scale(A2$edge),
                           A2$open),
                         c(scale(A3$forage),
                           scale(A3$edge),
                           A3$open))

#_______________________________________________________________________
# 3. Define simulation parameters ----
#_______________________________________________________________________
# 3a. Movement parameter distributions ----
#_______________________________________________________________________

# step lengths (gamma)
# make distribution
sl.dist <- make_gamma_distr(shape = 1.2,
                            scale = 70)

# turning angles (uniform)
# make distribution
ta.dist <- make_unif_distr(min = -pi,
                           max = pi)

#_______________________________________________________________________
# 3b. Habitat selection coefficients ----
#_______________________________________________________________________

# mean coefficients
coef.forage <- 1.5          # selection for forage
coef.forage.sl <- -0.3      # shorter sl with higher forage       
coef.edge <- -0.5           # avoidance of edge distance
coef.open <- -2.5           # base avoidance of open (start of step)
coef.open.sl <- 0.5         # interaction with log(sl) (longer movements when starting in open)

#_______________________________________________________________________
# 3c. Home ranging parameters ----
#_______________________________________________________________________

# potential "home range" centroids
# allow this to be a random uniform draw from x and y within the unit boundary

# extract bounding box
unit.bbox <- st_bbox(unit.bound)

# we'll take a random draw from both distributions as each path is initialized

# bivariate normal variance
# during initial runs it appears that 5,000 gives reasonable HRs (~ 2 ha)
e.var <- 5000

# the Bx and By coefficients are then calculated from the "home range" centroid and the variance

#_______________________________________________________________________
# 3d. Redistribution kernel parameters ----
#_______________________________________________________________________

# control steps
rk.control <- 100

# tolerance outside the landscape (hopefully we won't have to deal with this much)
rk.tolerance <- 0.01    # 1%

#_______________________________________________________________________
# 4. Define functions needed for simulation ----
#_______________________________________________________________________
# 4a. Home range centroid and start step ----
#_______________________________________________________________________

# here we'll ensure that the HRC is within the unit (true target of density estimation)
# while the start step is nearby

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
# 4b. Calculate home ranging parameters ----
#_______________________________________________________________________

hr_params <- function(e.var = e.var,          # expected variance of the bivariate normal
                      hrc = hrc)              # home range centroid as previously drawn
  
  {
  
  # x2 + y2 coefficient
  b.x2y2 <- -1 / e.var
  
  # solve for x and y coefficients
  b.x <- hrc[1] * 2 * -b.x2y2 
  b.y <- hrc[2] * 2 * -b.x2y2
  
  hr.params <- c(b.x, b.y, b.x2y2)
  
  return(hr.params)
  
}

#_______________________________________________________________________
# 5. Run simulations iteratively ----

# n of replicates for each landscape replicate
n.reps <- 100

# each one takes about 3 seconds
(3 * 100) / 60      # 5 minutes??

#_______________________________________________________________________
# 5a. Create dfs to hold all sims ----
#_______________________________________________________________________

sims.df <- data.frame()

#_______________________________________________________________________
# 5b. Define function ----
#_______________________________________________________________________

init_sim <- function(id.rep) {
  
  # extract landscapes
  focal.landscape.B <- landscape.covs.B[[id.rep]]
  focal.landscape.A <- landscape.covs.A[[id.rep]]
  
  # loop through each replicate
  sims.df <- data.frame()
  
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
                           hrc = hrc)
    
    # make iSSF model
    # here the terms are important to get right so redistribution_kernel() works okay
    issf.model <- make_issf_model(coefs = c("forage_end" = coef.forage,  
                                            "forage_end:log(sl_)" = coef.forage.sl,
                                            "edge_end" = coef.edge,
                                            "open_start" = coef.open,
                                            "log(sl_):open_start" = coef.open.sl,
                                            x2_ = hr.params[1],
                                            y2_ = hr.params[2], 
                                            "I(x2_^2 + y2_^2)" = hr.params[3]),              
                                  sl = sl.dist,
                                  ta = ta.dist)
    
    # BEFORE simulations
    # initialize redistribution kernel
    rk.B <- redistribution_kernel(x = issf.model,
                                  start = start.step,
                                  map = focal.landscape.B,     # use correct landscape here
                                  n.control = rk.control,   
                                  max.dist = get_max_dist(issf.model),
                                  tolerance.outside = rk.tolerance)
    
    # run simulation
    sim.path.B <- simulate_path(rk.B,
                                n.steps = 336,
                                start = start.step,
                                verbose = TRUE)
    
    # extract endpoint for starting point of AFTER simulation
    end.step <- make_start(x = c(sim.path.B$x_[337],
                                 sim.path.B$y_[337]),
                           ta_ = 0,
                           time = sim.path.B$t_[337],
                           dt = hours(2),
                           crs = crs("EPSG:32611"))
    
    # AFTER simulations
    # initialize redistribution kernel
    rk.A <- redistribution_kernel(x = issf.model,  # same issf parameters
                                  start = end.step,
                                  map = focal.landscape.A,     # use correct landscape here
                                  n.control = rk.control,   
                                  max.dist = get_max_dist(issf.model),
                                  tolerance.outside = rk.tolerance)
    
    # run simulation
    sim.path.A <- simulate_path(rk.A,
                                n.steps = 336,
                                start = end.step,
                                verbose = TRUE)
    
    # bind together
    sim.path.B$trt <- "before"
    sim.path.A$trt <- "after"
    
    sim.path <- rbind(sim.path.B,
                      sim.path.A)
    
    sim.path.1 <- sim.path %>% 
      
      mutate(indiv = i,
             rep = id.rep)
    
    # bind to df
    sims.df <- rbind(sims.df, sim.path.1)
    
    # status message (every 10 iterations)
    if (i %% 10 == 0) {
    
    elapsed.time <- round(as.numeric(difftime(Sys.time(), 
                                              start.time, 
                                              units = "mins")), 
                          digits = 1)
    
    print(paste0("Completed path ", i, " of ", n.reps, " - ", elapsed.time, " mins"))
    
    }
  
  }
  
  # return
  return(sims.df)

}

#_______________________________________________________________________
# 5c. Use function ----
#_______________________________________________________________________

init.sim.1 <- init_sim(1)
init.sim.2 <- init_sim(2)
init.sim.3 <- init_sim(3)

#_______________________________________________________________________
# 6. Plot tracks ----
#_______________________________________________________________________

# plot paths
ggplot() +
  
  theme_bw() +
  
  # unit boundary as a black square
  geom_sf(data = unit.bound,
          fill = NA) +
  
  facet_wrap(~ trt) +
  
  # paths
  geom_point(data = init.sim.3,
            aes(x = x_,
                y = y_),
            alpha = 0.05,
            size = 0.5,
            color = "purple") +
  
  # coordinate system to make nice axis labels
  coord_sf(datum = st_crs(32611))

#_______________________________________________________________________
# 7. Write to .csvs ----
#_______________________________________________________________________

# bind together
init.sim.all <- rbind(init.sim.1, init.sim.2, init.sim.3)

write.csv(init.sim.all, paste0(getwd(), "/Derived_data/Simulated data/init_sims.csv"))
