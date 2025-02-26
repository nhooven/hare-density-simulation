# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 10 - Simulate utilization distributions
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 09 Dec 2024
# Date completed: 12 Dec 2024
# Date last modified: 26 Feb 2025
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
# 2. Read in data ----
#_______________________________________________________________________

# iSSF coefficients
all.params <- read.csv(paste0(getwd(), "/Derived_data/Model parameters/all_params.csv"))

# unit boundary
unit.bound <- st_read(paste0(getwd(), "/Derived_data/Shapefiles/unit_bound.shp"))

# scaled covariate rasters
landscape.covs.B1 <- rast("Rasters/Scaled covariates/B1.tif")
landscape.covs.B2 <- rast("Rasters/Scaled covariates/B2.tif")
landscape.covs.B3 <- rast("Rasters/Scaled covariates/B3.tif")
landscape.covs.A1 <- rast("Rasters/Scaled covariates/A1.tif")
landscape.covs.A2 <- rast("Rasters/Scaled covariates/A2.tif")
landscape.covs.A3 <- rast("Rasters/Scaled covariates/A3.tif")

# sl/ta distributions
sl.params <- read.csv(paste0(getwd(), "/Derived_data/Lookup/sl_params.csv"))
ta.params <- read.csv(paste0(getwd(), "/Derived_data/Lookup/ta_params.csv"))

#_______________________________________________________________________
# 3. Define simulation parameters ----
#_______________________________________________________________________
# 3a. Habitat selection coefficients and identifiers ----
#_______________________________________________________________________

# keep only iSSF parameters
model.params <- all.params %>% filter(type == "iSSF")

#_______________________________________________________________________
# 3b. Initialization space ----
#_______________________________________________________________________

# extract bounding box
unit.bbox <- st_bbox(unit.bound)

# define function
hr_params <- function(e.var = 2000,          # variance of the bivariate normal
                      hrc = hrc)             # home range centroid as previously drawn
  
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
# 4. Define grid of starting locations ----

# instead of picking random starting locations or starting them from
# all the same location, let's create a regular grid that encompasses the
# unit and some area around it, and simulate many tracks for each individual

# we really just want predicted probability of use over our focal area
# each pixel of the unit will be equidistant to the same number of 
# starting locations

#_______________________________________________________________________
# 4a. Define sampling space ----
#_______________________________________________________________________

# buffer unit boundary (double the area)
# our current unit is sqrt(100,000) wide and tall
# we want a square that is 2.5x that
side.length <- sqrt(250000)

# extract coordinates - find centroid of raster
rast.centroid <- c(ext(landscape.covs.B1)[2] / 2,
                   ext(landscape.covs.B1)[2] / 2)

# how many meters per side? (let's make our unit 25 ha)
m.side <- sqrt(25 * 10000)

# how many meters to add and subtract?
m.side.half <- m.side / 2

# create polygon feature
unit.buff <- st_polygon(list(cbind(c(rast.centroid[1] - m.side.half, 
                                     rast.centroid[1] + m.side.half,
                                     rast.centroid[1] + m.side.half,
                                     rast.centroid[1] - m.side.half,
                                     rast.centroid[1] - m.side.half), 
                                   c(rast.centroid[2] - m.side.half, 
                                     rast.centroid[2] - m.side.half,
                                     rast.centroid[2] + m.side.half,
                                     rast.centroid[2] + m.side.half,
                                     rast.centroid[2] - m.side.half))))

# assign to sf object
unit.buff.sf <- st_as_sf(st_sfc(unit.buff))

st_crs(unit.buff.sf) <- crs("EPSG:32611")

# plot
plot(st_geometry(unit.buff.sf))
plot(st_geometry(unit.bound), add = T)

#_______________________________________________________________________
# 4b. Sample 100 locations ----
#_______________________________________________________________________

sampled.starts <- st_sample(unit.buff.sf, 
                            size = 100,
                            type = "regular",
                            exact = TRUE)

# just get these as centered as possible

plot(st_geometry(unit.buff.sf))
plot(st_geometry(unit.bound), add = T)
plot(st_geometry(sampled.starts), add = T)

sampled.starts.coords <- st_coordinates(sampled.starts)

#_______________________________________________________________________
# 5. Run transient UD simulations ----
#_______________________________________________________________________
# 5a. Define function ----
#_______________________________________________________________________

tud_sim <- function (landscape.covs,
                     id.trt,
                     id.rep,
                     n.reps = 100,
                     n.steps = 336) {
  
  # extract fit parameters
  focal.params <- model.params %>%
    
    filter(trt == id.trt,
           rep == id.rep)
  
  # extract sl/ta parameters
  focal.sl <- sl.params %>%
    
    filter(trt == id.trt,
           rep == id.rep)
  
  focal.ta <- ta.params %>%
    
    filter(trt == id.trt,
           rep == id.rep)
  
  # make distributions
  sl.dist <- make_gamma_distr(shape = focal.sl$shape,
                              scale = focal.sl$scale)
  
  ta.dist <- make_vonmises_distr(kappa = focal.ta$kappa)
  
  # loop through individuals
  sim.all <- data.frame()
  
  start.time <- Sys.time()
  
  for (i in 1:n.reps) {
    
    # define start step
    sampled.start <- sampled.starts.coords[i, ]
    
    # hrc - common to all reps within this individual
    hrc <- c(sampled.start[1], sampled.start[2])
    
    # start steps begin with a little noise
    start.step <- make_start(x = c(rnorm(1, hrc[1], 15),
                                   rnorm(1, hrc[2], 15)),
                             ta_ = 0,
                             time = ymd_hm("2024-09-01 18:00", 
                                           tz = "America/Los_Angeles"),
                             dt = hours(2),
                             crs = crs("EPSG:32611"))
    
    # home ranging parameters
    hr.params <- hr_params(e.var = 2000,      
                           hrc = hrc)
    
    hr.params <- unname(hr.params)
    
    # make iSSF model
    # here the terms are important to get right so redistribution_kernel() works okay
    issf.model <- make_issf_model(coefs = c("sl_:fora_start" = focal.params$estimate[focal.params$term == "sl_:fora_start"],
                                            "fora_end" = focal.params$estimate[focal.params$term == "fora_end"],   
                                            "fora_end:cos(ta_)" = focal.params$estimate[focal.params$term == "fora_end:cos(ta_)"],
                                            "elev_end" = focal.params$estimate[focal.params$term == "elev.s"],
                                            "I(elev_end^2)" = focal.params$estimate[focal.params$term == "I(elev.s^2)"],
                                            "sl_:open_start" = focal.params$estimate[focal.params$term == "sl_:open_start"],
                                            "sl_" = focal.params$estimate[focal.params$term == "sl_"],
                                            "log(sl_)" = focal.params$estimate[focal.params$term == "log(sl_)"],
                                            "cos(ta_)" = focal.params$estimate[focal.params$term == "cos(ta_)"],
                                            x2_ = hr.params[1],
                                            y2_ = hr.params[2], 
                                            "I(x2_^2 + y2_^2)" = hr.params[3]),              
                                  sl = sl.dist,
                                  ta = ta.dist)
    
    # initialize redistribution kernel
    rk <- redistribution_kernel(x = issf.model,
                                start = start.step,
                                map = landscape.covs,
                                n.control = 100,   
                                max.dist = get_max_dist(issf.model),
                                tolerance.outside = 0.01)
    
    # create df to hold all sims
    sim.indiv.all <- data.frame()
    
    # run simulations - loop through n.reps
    for (j in 1:n.reps) {
      
      # run simulation
      sim.path <- simulate_path(rk,
                                n.steps = n.steps,
                                start = start.step,
                                verbose = TRUE)
      
      # add identifiers
      sim.path <- sim.path %>%
        
        mutate(trt = id.trt,
               rep = id.rep,
               indiv = i,
               sim.rep = j)
      
      # bind to df
      sim.indiv.all <- rbind(sim.indiv.all, sim.path)
      
    }
    
    # bind to df
    sim.all <- rbind(sim.all, sim.indiv.all)
    
    # status message (every 10 individuals)
    if (i %% 10 == 0) {
      
      elapsed.time <- round(as.numeric(difftime(Sys.time(), 
                                                start.time, 
                                                units = "mins")), 
                            digits = 1)
      
      print(paste0("Completed individual ", i, " of ", nrow(n.reps), " - ", elapsed.time, " mins"))
      
    }
    
  }
  
  # return
  return(sim.all)
  
}

#_______________________________________________________________________
# 5b. Run simulations ----
#_______________________________________________________________________

sims.B1 <- tud_sim(landscape.covs.B1, "before", 1)
sims.B2 <- tud_sim(landscape.covs.B2, "before", 2)
sims.B3 <- tud_sim(landscape.covs.B3, "before", 3)
sims.A1 <- tud_sim(landscape.covs.A1, "after", 1)
sims.A2 <- tud_sim(landscape.covs.A2, "after", 2)
sims.A3 <- tud_sim(landscape.covs.A3, "after", 3)
