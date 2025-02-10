# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 10 - Simulate utilization distributions
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 09 Dec 2024
# Date completed: 12 Dec 2024
# Date last modified: 31 Jan 2025
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
landscape.covs.S1L <- rast("Rasters/Scaled covariates/S1L.tif")
landscape.covs.S2L <- rast("Rasters/Scaled covariates/S2L.tif")
landscape.covs.S3L <- rast("Rasters/Scaled covariates/S3L.tif")
landscape.covs.S1H <- rast("Rasters/Scaled covariates/S1H.tif")
landscape.covs.S2H <- rast("Rasters/Scaled covariates/S2H.tif")
landscape.covs.S3H <- rast("Rasters/Scaled covariates/S3H.tif")

landscape.covs.C1L <- rast("Rasters/Scaled covariates/C1L.tif")
landscape.covs.C2L <- rast("Rasters/Scaled covariates/C2L.tif")
landscape.covs.C3L <- rast("Rasters/Scaled covariates/C3L.tif")
landscape.covs.C1H <- rast("Rasters/Scaled covariates/C1H.tif")
landscape.covs.C2H <- rast("Rasters/Scaled covariates/C2H.tif")
landscape.covs.C3H <- rast("Rasters/Scaled covariates/C3H.tif")

# sl/ta distributions
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_S1L.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_S2L.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_S3L.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_S1H.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_S2H.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_S3H.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_C1L.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_C2L.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_C3L.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_C1H.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_C2H.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/sl_dist_C3H.RData"))
load(paste0(getwd(), "/Derived_data/Model parameters/ta_dist.RData"))

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
hr_params <- function(e.var = 5000,          # variance of the bivariate normal
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
# 4. Burn-in simulations (half a "month") ----

# here we'll start 100 simulations from the landscape centroid,
# use the point estimates from the iSSF, and extract the 
# endpoints as starting steps for the subsequent sims

#_______________________________________________________________________
# 4a. Define function ----
#_______________________________________________________________________

tud_sim_burnin <- function (landscape.covs,
                            sl.dist, 
                            id.landscape,
                            id.variability,
                            id.rep,
                            n.reps = 100,
                            n.steps = 168) {
  
  # extract fit parameters
  focal.params <- model.params %>%
    
    filter(landscape == id.landscape &
           variability == id.variability,
           rep == id.rep)
  
  # create df to hold all sims
  sim.endpoint.all <- data.frame()
  
  # run simulations
  start.time <- Sys.time()
  
  # loop through n.reps
  for (i in 1:n.reps) {
    
    # home range center (centroid)
    hrc <- c(ext(landscape.covs)[2] / 2, ext(landscape.covs)[2] / 2)
    
    # define start step
    start.step <- make_start(x = c(hrc[1],
                                   hrc[2]),
                             ta_ = 0,
                             time = ymd_hm("2024-09-01 18:00", 
                                           tz = "America/Los_Angeles"),
                             dt = hours(2),
                             crs = crs("EPSG:32611"))
    
    # calculate home ranging parameters and unname
    hr.params <- hr_params(hrc = hrc)
    
    hr.params <- unname(hr.params)
    
    # make iSSF model
    # here the terms are important to get right so redistribution_kernel() works okay
    issf.model <- make_issf_model(coefs = c("forage_end" = focal.params$estimate[focal.params$term == "forage.s"],  
                                            "forage_end:log(sl_)" = focal.params$estimate[focal.params$term == "forage.s:log(sl_)"],
                                            "edge_end" = focal.params$estimate[focal.params$term == "edge.s"],
                                            "open_start" = focal.params$estimate[focal.params$term == "open"],
                                            "log(sl_):open_start" = focal.params$estimate[focal.params$term == "open:log(sl_)"],
                                            "log(sl_)" = focal.params$estimate[focal.params$term == "log(sl_)"],
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
    
    # run simulation
    sim.path <- simulate_path(rk,
                              n.steps = n.steps,
                              start = start.step,
                              verbose = TRUE)
    
    # add identifiers and extract endpoint
    sim.endpoint <- sim.path %>%
      
      slice(n()) %>%
      
      mutate(landscape = id.landscape,
             variability = id.variability,
             rep = id.rep,
             sim.rep = i)
    
    # bind to df
    sim.endpoint.all <- rbind(sim.endpoint.all, sim.endpoint)
    
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
  return(sim.endpoint.all)
  
}

test.run <- tud_sim_burnin(landscape.covs.S1L, sl.dist.S1L, "simple", "low", 1)


test.run.sf <- st_as_sf(test.run,
                        coords = c("x_", 
                                   "y_"),
                        crs = "epsg:32611")

# plot
ggplot() +
  
  theme_bw() +
  
  # points
  geom_sf(data = test.run.sf,
          size = 1.5) +
  
  # unit boundary as a black square
  geom_sf(data = unit.bound,
          fill = NA,
          color = "black",
          linewidth = 1.25) +
  
  # remove legend
  theme(legend.position = "none") +
  
  # coordinate system to make nice axis labels
  coord_sf(datum = st_crs(32611))

# looks reasonable - we can use this procedure to generate starting locations no problem

#_______________________________________________________________________
# 4b. Use function ----
#_______________________________________________________________________

burnin.S1L <- tud_sim_burnin(landscape.covs.S1L, sl.dist.S1L, "simple", "low", 1)   # completed 01-30-2025
burnin.S2L <- tud_sim_burnin(landscape.covs.S2L, sl.dist.S2L, "simple", "low", 2)   # completed 01-30-2025
burnin.S3L <- tud_sim_burnin(landscape.covs.S3L, sl.dist.S3L, "simple", "low", 3)   # completed 01-30-2025
burnin.S1H <- tud_sim_burnin(landscape.covs.S1H, sl.dist.S1H, "simple", "high", 1)  # completed 01-31-2025
burnin.S2H <- tud_sim_burnin(landscape.covs.S2H, sl.dist.S2H, "simple", "high", 2)  # completed 01-31-2025
burnin.S3h <- tud_sim_burnin(landscape.covs.S3H, sl.dist.S3H, "simple", "high", 3)  # completed 01-31-2025

burnin.C1L <- tud_sim_burnin(landscape.covs.C1L, sl.dist.C1L, "complex", "low", 1)  # completed 01-31-2025
burnin.C2L <- tud_sim_burnin(landscape.covs.C2L, sl.dist.C2L, "complex", "low", 2)  # completed 01-31-2025
burnin.C3L <- tud_sim_burnin(landscape.covs.C3L, sl.dist.C3L, "complex", "low", 3)  # completed 01-31-2025
burnin.C1H <- tud_sim_burnin(landscape.covs.C1H, sl.dist.C1H, "complex", "high", 1) # completed 01-31-2025
burnin.C2H <- tud_sim_burnin(landscape.covs.C2H, sl.dist.C2H, "complex", "high", 2) # completed 01-31-2025
burnin.C3H <- tud_sim_burnin(landscape.covs.C3H, sl.dist.C3H, "complex", "high", 3) # completed 01-31-2025

#_______________________________________________________________________
# 5. Run transient UD simulations ----
#_______________________________________________________________________
# 5a. Define function ----
#_______________________________________________________________________

tud_sim <- function (burnin,
                     landscape.covs,
                     sl.dist, 
                     id.landscape,
                     id.variability,
                     id.rep,
                     n.reps = 100,
                     n.steps = 336) {
  
  # extract fit parameters
  focal.params <- model.params %>%
    
    filter(landscape == id.landscape &
             variability == id.variability,
           rep == id.rep)
  
  # loop through individuals (sim.rep from burnin)
  sim.all <- data.frame()
  
  start.time <- Sys.time()
  
  for (i in unique(burnin$sim.rep)) {
    
    # define focal burnin row
    burnin.focal <- burnin[i, ]
    
    # draw from random slope distributions (each individual gets its own)
    focal.params.draw <- data.frame(forage = rnorm(1,
                                                   focal.params$estimate[focal.params$term == "forage.s"],
                                                   focal.params$estimate[focal.params$term == "sd__forage.s"]),
                                    edge = rnorm(1,
                                                 focal.params$estimate[focal.params$term == "edge.s"],
                                                 focal.params$estimate[focal.params$term == "sd__edge.s"]), 
                                    open = rnorm(1,
                                                 focal.params$estimate[focal.params$term == "open"],
                                                 focal.params$estimate[focal.params$term == "sd__open"]),
                                    log.sl = focal.params$estimate[focal.params$term == "log(sl_)"],
                                    forage.log.sl = focal.params$estimate[focal.params$term == "forage.s:log(sl_)"],
                                    open.log.sl = focal.params$estimate[focal.params$term == "open:log(sl_)"])
    
    # start step - common to all reps within this individual
    hrc <- c(burnin.focal$x_, burnin.focal$y_)
    
    start.step <- make_start(x = c(hrc[1],
                                   hrc[2]),
                             ta_ = 0,
                             time = ymd_hm("2024-09-01 18:00", 
                                           tz = "America/Los_Angeles"),
                             dt = hours(2),
                             crs = crs("EPSG:32611"))
    
    # home ranging parameters
    hr.params <- hr_params(e.var = 5000 * 50,      # 20 times just to ensure tracks stay close but not too close
                           hrc = hrc)
    
    hr.params <- unname(hr.params)
    
    # make iSSF model
    # here the terms are important to get right so redistribution_kernel() works okay
    issf.model <- make_issf_model(coefs = c("forage_end" = focal.params.draw$forage,  
                                            "forage_end:log(sl_)" = focal.params.draw$forage.log.sl,
                                            "edge_end" = focal.params.draw$edge,
                                            "open_start" = focal.params.draw$open,
                                            "log(sl_):open_start" = focal.params.draw$open.log.sl,
                                            "log(sl_)" = focal.params.draw$log.sl,
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
    
    # run simulations
    
    
    # run simulations - loop through n.reps
    for (j in 1:n.reps) {
      
      # run simulation
      sim.path <- simulate_path(rk,
                                n.steps = n.steps,
                                start = start.step,
                                verbose = TRUE)
      
      # add identifiers
      sim.path <- sim.path %>%
        
        mutate(landscape = id.landscape,
               variability = id.variability,
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
      
      print(paste0("Completed individual ", i, " of ", nrow(burnin), " - ", elapsed.time, " mins"))
      
    }
    
  }
  
  # return
  return(sim.all)
  
}

# 100 runs of this takes 5.3 minutes, so the full 100 x 100 would take 8.83 hours

#_______________________________________________________________________
# 5b. Run simulations ----
#_______________________________________________________________________

sims.S1L <- tud_sim(burnin.S1L, landscape.covs.S1L, sl.dist.S1L, "simple", "low", 1)
sims.S2L <- tud_sim(burnin.S2L, landscape.covs.S2L, sl.dist.S2L, "simple", "low", 2)
sims.S3L <- tud_sim(burnin.S3L, landscape.covs.S3L, sl.dist.S3L, "simple", "low", 3)
sims.S1H <- tud_sim(burnin.S1H, landscape.covs.S1H, sl.dist.S1H, "simple", "high", 1)
sims.S2H <- tud_sim(burnin.S2H, landscape.covs.S2H, sl.dist.S2H, "simple", "high", 2)
sims.S3H <- tud_sim(burnin.S3H, landscape.covs.S3H, sl.dist.S3H, "simple", "high", 3)

sims.C1L <- tud_sim(burnin.C1L, landscape.covs.C1L, sl.dist.C1L, "complex", "low", 1)
sims.C2L <- tud_sim(burnin.C2L, landscape.covs.C2L, sl.dist.C2L, "complex", "low", 2)
sims.C3L <- tud_sim(burnin.C3L, landscape.covs.C3L, sl.dist.C3L, "complex", "low", 3)
sims.C1H <- tud_sim(burnin.C1H, landscape.covs.C1H, sl.dist.C1H, "complex", "high", 1)
sims.C2H <- tud_sim(burnin.C2H, landscape.covs.C2H, sl.dist.C2H, "complex", "high", 2)
sims.C3H <- tud_sim(burnin.C3H, landscape.covs.C3H, sl.dist.C3H, "complex", "high", 3)



# 01-31-2025
# let's write a separate script for fitting aKDEs, just keeps things cleaner









#_______________________________________________________________________
# 6. Plot points ----
#_______________________________________________________________________
# 6a. Convert to sf ----
#_______________________________________________________________________

sims.S1L.sf <- st_as_sf(sims.S1L,
                        coords = c("x_", 
                                   "y_"),
                        crs = "epsg:32611") 

sims.S2L.sf <- st_as_sf(sims.S2L %>% filter(which.point == "end"),
                        coords = c("x_", 
                                   "y_"),
                        crs = "epsg:32611")

sims.S3L.sf <- st_as_sf(sims.S3L %>% filter(which.point == "end"),
                        coords = c("x_", 
                                   "y_"),
                        crs = "epsg:32611")

sims.S1H.sf <- st_as_sf(sims.S1H %>% filter(which.point == "end"),
                        coords = c("x_", 
                                   "y_"),
                        crs = "epsg:32611")

sims.S2H.sf <- st_as_sf(sims.S2H %>% filter(which.point == "end"),
                        coords = c("x_", 
                                   "y_"),
                        crs = "epsg:32611")

sims.S3H.sf <- st_as_sf(sims.S3H %>% filter(which.point == "end"),
                        coords = c("x_", 
                                   "y_"),
                        crs = "epsg:32611")

sims.C1L.sf <- st_as_sf(sims.C1L %>% filter(which.point == "end"),
                        coords = c("x_", 
                                   "y_"),
                        crs = "epsg:32611")

sims.C2L.sf <- st_as_sf(sims.C2L %>% filter(which.point == "end"),
                        coords = c("x_", 
                                   "y_"),
                        crs = "epsg:32611")

sims.C3L.sf <- st_as_sf(sims.C3L %>% filter(which.point == "end"),
                        coords = c("x_", 
                                   "y_"),
                        crs = "epsg:32611")

sims.C1H.sf <- st_as_sf(sims.C1H %>% filter(which.point == "end"),
                        coords = c("x_", 
                                   "y_"),
                        crs = "epsg:32611")

sims.C2H.sf <- st_as_sf(sims.C2H %>% filter(which.point == "end"),
                        coords = c("x_", 
                                   "y_"),
                        crs = "epsg:32611")

sims.C3H.sf <- st_as_sf(sims.C3H %>% filter(which.point == "end"),
                        coords = c("x_", 
                                   "y_"),
                        crs = "epsg:32611")

#_______________________________________________________________________
# 6b. Plots ----
#_______________________________________________________________________

# SL
ggplot() +
  
  theme_bw() +
  
  geom_sf(data = unit.buff) +
  
  # points
  geom_sf(data = sims.S1L.sf,
            alpha = 0.15,
          size = 0.50) +
  
  # unit boundary as a black square
  geom_sf(data = unit.bound,
          fill = NA,
          color = "yellow",
          linewidth = 1.25) +
  
  # remove legend
  theme(legend.position = "none") +
  
  # coordinate system to make nice axis labels
  coord_sf(datum = st_crs(32611))

#_______________________________________________________________________
# 7. Fit kernel UDs ----
#_______________________________________________________________________
# 7a. Initialize a SpatialPixels object equivalent to the main raster ----
#_______________________________________________________________________

landscape.sp <- as(as(landscape.covs.S1L$stem, "Raster"), "SpatialPixels")

#_______________________________________________________________________
# 7b. Define function ----
#_______________________________________________________________________

sim_kernel <- function (template,     # SpatialPixels object
                        sims.sf) {
  
  # convert sf to SpatialPoints objects
  sims.sp <- as(sims.sf, "Spatial") 
  
  # fit kernel with LSCV h parameter
  kernel <- kernelUD(xy = sims.sp, h = "LSCV", grid = template)
  
  # convert to SpatRaster
  kernel.rast <- rast(kernel)
  
  # return
  return(kernel.rast)
  
}

#_______________________________________________________________________
# 7c. Use function ----
#_______________________________________________________________________

kernel.S1L <- sim_kernel(landscape.sp, sims.S1L.sf)
kernel.S2L <- sim_kernel(landscape.sp, sims.S2L.sf)
kernel.S3L <- sim_kernel(landscape.sp, sims.S3L.sf)
kernel.S1H <- sim_kernel(landscape.sp, sims.S1H.sf)
kernel.S2H <- sim_kernel(landscape.sp, sims.S2H.sf)
kernel.S3H <- sim_kernel(landscape.sp, sims.S3H.sf)

kernel.C1L <- sim_kernel(landscape.sp, sims.C1L.sf)
kernel.C2L <- sim_kernel(landscape.sp, sims.C2L.sf)
kernel.C3L <- sim_kernel(landscape.sp, sims.C3L.sf)
kernel.C1H <- sim_kernel(landscape.sp, sims.C1H.sf)
kernel.C2H <- sim_kernel(landscape.sp, sims.C2H.sf)
kernel.C3H <- sim_kernel(landscape.sp, sims.C3H.sf)

# crop to unit boundary
kernel.S1L.crop <- crop(kernel.S1L, unit.bound)
kernel.S2L.crop <- crop(kernel.S2L, unit.bound)
kernel.S3L.crop <- crop(kernel.S3L, unit.bound)
kernel.S1H.crop <- crop(kernel.S1H, unit.bound)
kernel.S2H.crop <- crop(kernel.S2H, unit.bound)
kernel.S3H.crop <- crop(kernel.S3H, unit.bound)

kernel.C1L.crop <- crop(kernel.C1L, unit.bound)
kernel.C2L.crop <- crop(kernel.C2L, unit.bound)
kernel.C3L.crop <- crop(kernel.C3L, unit.bound)
kernel.C1H.crop <- crop(kernel.C1H, unit.bound)
kernel.C2H.crop <- crop(kernel.C2H, unit.bound)
kernel.C3H.crop <- crop(kernel.C3H, unit.bound)

#_______________________________________________________________________
# 7d. Plot for a sanity check ----
#_______________________________________________________________________

ggplot() +
  
  theme_bw() +
  
  tidyterra::geom_spatraster(data = kernel.C3H.crop,
                             aes(fill = ud)) +
  
  coord_sf(datum = sf::st_crs(32611)) +
  
  scale_fill_viridis_c(option = "magma")

#_______________________________________________________________________
# 8. Write to .csvs ----
#_______________________________________________________________________

write.csv(sims.S1L, paste0(getwd(), "/Derived_data/Simulated data/sims_UD_S1L.csv"))
write.csv(sims.S2L, paste0(getwd(), "/Derived_data/Simulated data/sims_UD_S2L.csv"))
write.csv(sims.S3L, paste0(getwd(), "/Derived_data/Simulated data/sims_UD_S3L.csv"))
write.csv(sims.S1H, paste0(getwd(), "/Derived_data/Simulated data/sims_UD_S1H.csv"))
write.csv(sims.S2H, paste0(getwd(), "/Derived_data/Simulated data/sims_UD_S2H.csv"))
write.csv(sims.S3H, paste0(getwd(), "/Derived_data/Simulated data/sims_UD_S3H.csv"))

write.csv(sims.C1L, paste0(getwd(), "/Derived_data/Simulated data/sims_UD_C1L.csv"))
write.csv(sims.C2L, paste0(getwd(), "/Derived_data/Simulated data/sims_UD_C2L.csv"))
write.csv(sims.C3L, paste0(getwd(), "/Derived_data/Simulated data/sims_UD_C3L.csv"))
write.csv(sims.C1H, paste0(getwd(), "/Derived_data/Simulated data/sims_UD_C1H.csv"))
write.csv(sims.C2H, paste0(getwd(), "/Derived_data/Simulated data/sims_UD_C2H.csv"))
write.csv(sims.C3H, paste0(getwd(), "/Derived_data/Simulated data/sims_UD_C3H.csv"))

#_______________________________________________________________________
# 9. Write rasters ----
#_______________________________________________________________________

writeRaster(kernel.S1L.crop, filename = paste0(getwd(), "/Rasters/kernel_S1L.tif"), overwrite = TRUE)
writeRaster(kernel.S2L.crop, filename = paste0(getwd(), "/Rasters/kernel_S2L.tif"), overwrite = TRUE)
writeRaster(kernel.S3L.crop, filename = paste0(getwd(), "/Rasters/kernel_S3L.tif"), overwrite = TRUE)
writeRaster(kernel.S1H.crop, filename = paste0(getwd(), "/Rasters/kernel_S1H.tif"), overwrite = TRUE)
writeRaster(kernel.S2H.crop, filename = paste0(getwd(), "/Rasters/kernel_S2H.tif"), overwrite = TRUE)
writeRaster(kernel.S3H.crop, filename = paste0(getwd(), "/Rasters/kernel_S3H.tif"), overwrite = TRUE)

writeRaster(kernel.C1L.crop, filename = paste0(getwd(), "/Rasters/kernel_C1L.tif"), overwrite = TRUE)
writeRaster(kernel.C2L.crop, filename = paste0(getwd(), "/Rasters/kernel_C2L.tif"), overwrite = TRUE)
writeRaster(kernel.C3L.crop, filename = paste0(getwd(), "/Rasters/kernel_C3L.tif"), overwrite = TRUE)
writeRaster(kernel.C1H.crop, filename = paste0(getwd(), "/Rasters/kernel_C1H.tif"), overwrite = TRUE)
writeRaster(kernel.C2H.crop, filename = paste0(getwd(), "/Rasters/kernel_C2H.tif"), overwrite = TRUE)
writeRaster(kernel.C3H.crop, filename = paste0(getwd(), "/Rasters/kernel_C3H.tif"), overwrite = TRUE)
