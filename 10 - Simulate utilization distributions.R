# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 10 - Simulate utilization distributions
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 09 Dec 2024
# Date completed: 12 Dec 2024
# Date last modified: 27 Dec 2024
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation
library(sf)              # read in shapefiles
library(terra)           # rasters
library(amt)             # simulate tracks
library(lubridate)       # work with time
library(sp)              # spatial points
library(raster)          
library(adehabitatHR)    # fit kernels

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

# SSF/iSSF coefficients
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
# 3c. Initialization space ----
#_______________________________________________________________________

# extract bounding box
unit.bbox <- st_bbox(unit.bound)

# we'll take the centroid from this box to initialize each path

#_______________________________________________________________________
# 3d. Calculate home ranging parameters ----
#_______________________________________________________________________

hr_params <- function(e.var = e.var,          # expected variance of the bivariate normal
                      e.var.sd = e.var.sd,    # sd of the draws for the bivariate normal variance
                      hrc = start.step)       # home range centroid as previously drawn
  
{
  
  # variance
  var.focal <- rnorm(n = 1, mean = e.var, sd = e.var.sd)
  
  # x2 + y2 coefficient
  b.x2y2 <- -1 / var.focal
  
  # solve for x and y coefficients
  b.x <- hrc$x_ * 2 * -b.x2y2 
  b.y <- hrc$y_ * 2 * -b.x2y2
  
  hr.params <- c(b.x, b.y, b.x2y2)
  
  return(hr.params)
  
}

#_______________________________________________________________________
# 3d. Redistribution kernel parameters ----
#_______________________________________________________________________

# control steps
rk.control <- 10000

# tolerance outside the landscape (hopefully we won't have to deal with this much)
rk.tolerance <- 0.05    # 5%

#_______________________________________________________________________
# 4. Run simulations iteratively ----
#_______________________________________________________________________
# 4a. Define function ----
#_______________________________________________________________________

sim_issf_ud <- function (landscape.covs,
                         sl.dist, 
                         n.reps = 300,
                         n.steps = 336,
                         id.landscape,
                         id.variability,
                         id.rep) {
  
  # extract fit parameters
  focal.params <- model.params %>%
    
    filter(landscape == id.landscape &
           variability == id.variability,
           rep == id.rep)
  
  # create df to hold all sims
  sims.df <- data.frame()
  
  # run simulations
  start.time <- Sys.time()
  
  for (i in 1:n.reps) {
    
    # define start step (centroid of unit to keep initial conditions constant)
    start.step <- make_start(x = c(unit.bbox[1] + (unit.bbox[3] - unit.bbox[1]) / 2,
                                   unit.bbox[2] + (unit.bbox[4] - unit.bbox[2]) / 2),
                                   ta_ = 0,
                                   time = ymd_hm("2024-09-01 18:00", 
                                                 tz = "America/Los_Angeles"),
                                   dt = hours(2),
                                   crs = crs("EPSG:32611"))
    
    # calculate home ranging parameters
    hr.params <- hr_params(e.var = 5000,
                           e.var.sd = 100,
                           hrc = start.step)
    
    # change names
    names(hr.params)[1:2] <- ""
                             
    # make iSSF model
    # here the terms are important to get right so redistribution_kernel() works okay
    issf.model <- make_issf_model(coefs = c("stem_end" = focal.params$estimate[focal.params$term == "stem.s"],  
                                            "stem_end:log(sl_)" = focal.params$estimate[focal.params$term == "stem.s:log(sl_)"],
                                            "edge_end" = focal.params$estimate[focal.params$term == "edge.s"],
                                            "mature_start" = focal.params$estimate[focal.params$term == "mature"],
                                            "log(sl_):mature_start" = focal.params$estimate[focal.params$term == "mature:log(sl_)"],
                                            "log(sl_)" = focal.params$estimate[focal.params$term == "log(sl_)"],
                                            "x2_" = hr.params[1],
                                            "y2_" = hr.params[2], 
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
                              n.steps = n.steps,
                              start = start.step,
                              verbose = TRUE)
    
    # 27 Dec 2024
    # Still haven't decided whether I want to use the TUD or the SSUD...
    
    # extract final location for fitting uncorrelated UD
    sim.final <- sim.path %>%
      
      slice(nrow(sim.path)) %>%
      
      # add in individual identifier
      mutate(indiv = i)
    
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
# 4b. Run simulations ----
#_______________________________________________________________________

sims.SL <- sim_issf_ud("simple", "low")
sims.CL <- sim_issf_ud("complex", "low")
sims.SH <- sim_issf_ud("simple", "high")
sims.CH <- sim_issf_ud("complex", "high")

#_______________________________________________________________________
# 6. Plot points ----
#_______________________________________________________________________
# 6a. Convert to sf ----
#_______________________________________________________________________

sims.SL.sf <- st_as_sf(sims.SL %>% drop_na(x_),
                       coords = c("x_", 
                                  "y_"),
                       crs = "epsg:32611") 

sims.CL.sf <- st_as_sf(sims.CL %>% drop_na(x_),
                       coords = c("x_", 
                                  "y_"),
                       crs = "epsg:32611") 

sims.SH.sf <- st_as_sf(sims.SH %>% drop_na(x_),
                       coords = c("x_", 
                                  "y_"),
                       crs = "epsg:32611") 

sims.CH.sf <- st_as_sf(sims.CH %>% drop_na(x_),
                       coords = c("x_", 
                                  "y_"),
                       crs = "epsg:32611") 

#_______________________________________________________________________
# 6b. Plots ----
#_______________________________________________________________________

# SL
ggplot() +
  
  theme_bw() +
  
  # unit boundary as a black square
  geom_sf(data = unit.bound,
          fill = NA) +
  
  # paths
  geom_sf(data = sims.SL.sf,
            alpha = 0.5) +
  
  # remove legend
  theme(legend.position = "none") +
  
  # coordinate system to make nice axis labels
  coord_sf(datum = st_crs(32611))

# CL
ggplot() +
  
  theme_bw() +
  
  # unit boundary as a black square
  geom_sf(data = unit.bound,
          fill = NA) +
  
  # paths
  geom_sf(data = sims.CL.sf,
          alpha = 0.5) +
  
  # remove legend
  theme(legend.position = "none") +
  
  # coordinate system to make nice axis labels
  coord_sf(datum = st_crs(32611))

# SH
ggplot() +
  
  theme_bw() +
  
  # unit boundary as a black square
  geom_sf(data = unit.bound,
          fill = NA) +
  
  # paths
  geom_sf(data = sims.SH.sf,
          alpha = 0.5) +
  
  # remove legend
  theme(legend.position = "none") +
  
  # coordinate system to make nice axis labels
  coord_sf(datum = st_crs(32611))

# CH
ggplot() +
  
  theme_bw() +
  
  # unit boundary as a black square
  geom_sf(data = unit.bound,
          fill = NA) +
  
  # paths
  geom_sf(data = sims.CH.sf,
          alpha = 0.5) +
  
  # remove legend
  theme(legend.position = "none") +
  
  # coordinate system to make nice axis labels
  coord_sf(datum = st_crs(32611))

#_______________________________________________________________________
# 7. Density maps ----
#_______________________________________________________________________
# 7a. Add identifier ----
#_______________________________________________________________________

sims.SL$scenario <- "SL"
sims.CL$scenario <- "CL"
sims.SH$scenario <- "SH"
sims.CH$scenario <- "CH"

sims.all <- rbind(sims.SL, sims.CL, sims.SH, sims.CH)

# order factor levels
sims.all$scenario <- factor(sims.all$scenario, levels = c("SL", "CL", "SH", "CH"))

#_______________________________________________________________________
# 7b. Plot ----
#_______________________________________________________________________

ggplot(sims.all,
       aes(x = x_,
           y = y_)) +
  
  facet_wrap(~ scenario) +
  
  theme_bw() +
  
  stat_density2d_filled(aes(fill = after_stat(level))) +
  
  scale_fill_viridis_d(option = "magma") +
  
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank())

# SH is weird and has many NAs because of the positive edge coefficient - 
# let's look at re-doing that model later

#_______________________________________________________________________
# 8. Fit kernel UDs ----
#_______________________________________________________________________
# 8a. Initialize a SpatialPixels object equivalent to the main raster ----
#_______________________________________________________________________

landscape.sp <- as(as(landscape.covs.simple$stem, "Raster"), "SpatialPixels")

#_______________________________________________________________________
# 8b. Define function ----
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
# 8c. Use function ----
#_______________________________________________________________________

kernel.SL <- sim_kernel(landscape.sp, sims.SL.sf)
kernel.CL <- sim_kernel(landscape.sp, sims.CL.sf)
kernel.SH <- sim_kernel(landscape.sp, sims.SH.sf)
kernel.CH <- sim_kernel(landscape.sp, sims.CH.sf)

#_______________________________________________________________________
# 8d. Plot ----
#_______________________________________________________________________

# SL
ggplot() +
  
  theme_bw() +
  
  tidyterra::geom_spatraster(data = kernel.SL,
                            aes(fill = ud)) +
  
  coord_sf(datum = sf::st_crs(32611)) +
  
  scale_fill_viridis_c(option = "magma")

# CL
ggplot() +
  
  theme_bw() +
  
  tidyterra::geom_spatraster(data = kernel.CL,
                             aes(fill = ud)) +
  
  coord_sf(datum = sf::st_crs(32611)) +
  
  scale_fill_viridis_c(option = "magma")

# SH
ggplot() +
  
  theme_bw() +
  
  tidyterra::geom_spatraster(data = kernel.SH,
                             aes(fill = ud)) +
  
  coord_sf(datum = sf::st_crs(32611)) +
  
  scale_fill_viridis_c(option = "magma")

# CL
ggplot() +
  
  theme_bw() +
  
  tidyterra::geom_spatraster(data = kernel.CL,
                             aes(fill = ud)) +
  
  coord_sf(datum = sf::st_crs(32611)) +
  
  scale_fill_viridis_c(option = "magma")

#_______________________________________________________________________
# 8e. Stack ----
#_______________________________________________________________________

kernel.all <- c(kernel.SL, kernel.CL, kernel.SH, kernel.CH)

names(kernel.all) <- c("SL", "CL", "SH", "CH")

#_______________________________________________________________________
# 9. Write to .csvs ----
#_______________________________________________________________________

write.csv(sims.SL, paste0(getwd(), "/Derived_data/Simulated data/sims_UD_SL.csv"))
write.csv(sims.CL, paste0(getwd(), "/Derived_data/Simulated data/sims_UD_CL.csv"))
write.csv(sims.SH, paste0(getwd(), "/Derived_data/Simulated data/sims_UD_SH.csv"))
write.csv(sims.CH, paste0(getwd(), "/Derived_data/Simulated data/sims_UD_CH.csv"))

#_______________________________________________________________________
# 10. Write rasters ----
#_______________________________________________________________________

writeRaster(kernel.all, filename = paste0(getwd(), "/Rasters/kernel_all_sim.tif"))
