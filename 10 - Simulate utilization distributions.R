# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 10 - Simulate utilization distributions
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 09 Dec 2024
# Date completed: 12 Dec 2024
# Date last modified: 10 Jan 2025
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
library(raster)          # rasters
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

model.params$term[model.params$term == "sd__mature.s"] <- "sd__mature"

#_______________________________________________________________________
# 3b. Initialization space ----
#_______________________________________________________________________

# extract bounding box
unit.bbox <- st_bbox(unit.bound)

# simulation "stage"
unit.buff <- st_buffer(unit.bound, dist = 1000)

plot(landscape.covs.S1L$stem)
plot(st_geometry(unit.buff), add = T)
plot(st_geometry(unit.bound), add = T)

#_______________________________________________________________________
# 3c. Convert everything outside the "stage" to the mean covariate value ----

# we're only interested in simulated movement within a small area of the 
# unit boundary. We want to avoid any large gradients pulling tracks away

#_______________________________________________________________________

# define function
define_stage <- function (landscape.covs) {
  
  # initialize three rasters
  out.stage.1 <- rast(landscape.covs$stem, vals = 0)
  out.stage.2 <- rast(landscape.covs$edge, vals = 0)
  out.stage.3 <- rast(landscape.covs$mature, vals = 0)
  
  # mask each raster
  mask.1 <- mask(landscape.covs$stem, unit.buff)
  mask.2 <- mask(landscape.covs$edge, unit.buff)
  mask.3 <- mask(landscape.covs$mature, unit.buff)
  
  # mosaic
  mosaic.1 <- mosaic(mask.1, out.stage.1)
  mosaic.2 <- mosaic(mask.2, out.stage.2)
  mosaic.3 <- mosaic(mask.3, out.stage.3)
  
  # ensure that what's "mature" is 1
  mosaic.3 <- ifel(mosaic.3 == 0.5, 1, 0)
  
  # bind together and return
  stage.rast <- c(mosaic.1, mosaic.2, mosaic.3)
  
  return(stage.rast)
  
}

# use function
landscape.covs.S1L.1 <- define_stage(landscape.covs.S1L)
landscape.covs.S2L.1 <- define_stage(landscape.covs.S2L)
landscape.covs.S3L.1 <- define_stage(landscape.covs.S3L)
landscape.covs.S1H.1 <- define_stage(landscape.covs.S1H)
landscape.covs.S2H.1 <- define_stage(landscape.covs.S2H)
landscape.covs.S3H.1 <- define_stage(landscape.covs.S3H)

landscape.covs.C1L.1 <- define_stage(landscape.covs.C1L)
landscape.covs.C2L.1 <- define_stage(landscape.covs.C2L)
landscape.covs.C3L.1 <- define_stage(landscape.covs.C3L)
landscape.covs.C1H.1 <- define_stage(landscape.covs.C1H)
landscape.covs.C2H.1 <- define_stage(landscape.covs.C2H)
landscape.covs.C3H.1 <- define_stage(landscape.covs.C3H)

#_______________________________________________________________________
# 3d. Calculate home ranging parameters ----
#_______________________________________________________________________

# determine a reasonable bivariate normal size
# we'll need to come up with a variance, simulate draws, and examine an
# outer quantile (e.g., 95%)

# load package
library(MASS)

# define means
bvn.mean <- c(ext(landscape.covs.C1H.1)[2] / 2, ext(landscape.covs.C1H.1)[2] / 2)

# variance-covariance matrix (assume they don't covary at all)
bvn.vcov <- matrix(c(300000, 0, 
                     0, 300000),
                   ncol = 2)

# sample
bvn.samples <- as.data.frame(mvrnorm(n = 10000,
                                     bvn.mean,
                                     bvn.vcov))

names(bvn.samples) <- c("x", "y")

# plot
ggplot() +
  
  geom_sf(data = st_as_sf(bvn.samples, 
                   coords = c("x", "y"),
                   crs = "epsg:32611"),
          alpha = 0.15) +
  
  geom_sf(data = unit.buff,
          fill = NA) +
  
  geom_sf(data = unit.bound) +
  
  coord_sf(datum = st_crs(32611))

# define function
hr_params <- function(e.var = 300000,        # variance of the bivariate normal
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
# 3f. Redistribution kernel parameters ----
#_______________________________________________________________________

# control steps
rk.control <- 10000

# tolerance outside the landscape (hopefully we won't have to deal with this much)
rk.tolerance <- 0.10    # 10%

#_______________________________________________________________________
# 4. Run steady-state utilization distribution (SSUD) simulations ----
#_______________________________________________________________________
# 4a. Define function ----
#_______________________________________________________________________

sim_issf_ssud <- function (landscape.covs,
                           sl.dist, 
                           id.landscape,
                           id.variability,
                           id.rep,
                           n.reps = 1,
                           n.steps = 20000) {
  
  # extract fit parameters
  focal.params <- model.params %>%
    
    filter(landscape == id.landscape &
           variability == id.variability,
           rep == id.rep)
  
  # create df to hold all sims
  sims.df <- data.frame()
  
  # run simulations
  start.time <- Sys.time()
  
  # loop through n.reps
  for (i in 1:n.reps) {
    
    # home range center (centroid)
    hrc <- c(ext(landscape.covs.C1H.1)[2] / 2, ext(landscape.covs.C1H.1)[2] / 2)
    
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
    issf.model <- make_issf_model(coefs = c("stem_end" = focal.params$estimate[focal.params$term == "stem.s"],  
                                            "stem_end:log(sl_)" = focal.params$estimate[focal.params$term == "stem.s:log(sl_)"],
                                            "edge_end" = focal.params$estimate[focal.params$term == "edge.s"],
                                            "mature_start" = focal.params$estimate[focal.params$term == "mature"],
                                            "log(sl_):mature_start" = focal.params$estimate[focal.params$term == "mature:log(sl_)"],
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
                                n.control = rk.control,   
                                max.dist = get_max_dist(issf.model),
                                tolerance.outside = rk.tolerance)
    
    # run simulation
    sim.path <- simulate_path(rk,
                              n.steps = n.steps,
                              start = start.step,
                              verbose = TRUE)
    
    # add identifiers
    sim.path.1 <- sim.path %>%
      
      mutate(landscape = id.landscape,
             variability = id.variability,
             rep = id.rep,
             sim.rep = i)
    
    # bind to df
    sims.df <- rbind(sims.df, sim.path.1)
    
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

sims.S1L <- sim_issf_ssud(landscape.covs.S1L.1, sl.dist.S1L, "simple", "low", 1)
sims.S2L <- sim_issf_ssud(landscape.covs.S2L.1, sl.dist.S2L, "simple", "low", 2)
sims.S3L <- sim_issf_ssud(landscape.covs.S3L.1, sl.dist.S3L, "simple", "low", 3)
sims.S1H <- sim_issf_ssud(landscape.covs.S1H.1, sl.dist.S1H, "simple", "high", 1)
sims.S2H <- sim_issf_ssud(landscape.covs.S2H.1, sl.dist.S2H, "simple", "high", 2)
sims.S3H <- sim_issf_ssud(landscape.covs.S3H.1, sl.dist.S3H, "simple", "high", 3)

sims.C1L <- sim_issf_ssud(landscape.covs.C1L.1, sl.dist.C1L, "complex", "low", 1)
sims.C2L <- sim_issf_ssud(landscape.covs.C2L.1, sl.dist.C2L, "complex", "low", 2)
sims.C3L <- sim_issf_ssud(landscape.covs.C3L.1, sl.dist.C3L, "complex", "low", 3)
sims.C1H <- sim_issf_ssud(landscape.covs.C1H.1, sl.dist.C1H, "complex", "high", 1)
sims.C2H <- sim_issf_ssud(landscape.covs.C2H.1, sl.dist.C2H, "complex", "high", 2)
sims.C3H <- sim_issf_ssud(landscape.covs.C3H.1, sl.dist.C3H, "complex", "high", 3)

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
