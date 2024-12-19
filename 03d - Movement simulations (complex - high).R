# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 03d - Movement simulations (complex - high)
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 15 Nov 2024
# Date completed: 21 Nov 2024
# Date last modified: 19 Dec 2024
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

landscape.covs.1 <- rast("Rasters/C1.tif")
landscape.covs.2 <- rast("Rasters/C2.tif")
landscape.covs.3 <- rast("Rasters/C3.tif")

unit.bound <- st_read(paste0(getwd(), "/Derived_data/Shapefiles/unit_bound.shp"))

# scale raster
landscape.covs.s <- c(scale(landscape.covs.3$stem),
                      scale(landscape.covs.3$edge),
                      landscape.covs.3$mature)

#_______________________________________________________________________
# 3. Define simulation parameters ----
#_______________________________________________________________________
# 3a. Movement parameter distributions ----

# here we'll just approximate hare movement from preliminary data

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
# 3b. Habitat selection coefficients and identifiers ----
#_______________________________________________________________________

# mean coefficients
coef.stem <- 1.0          # selection for stem density
coef.stem.sl <- -0.3      # higher selection with shorter sl       
coef.edge <- -0.5         # avoidance of edge distance
coef.mature <- -1.5       # base avoidance of mature (start of step)
coef.mature.sl <- 0.5     # interaction with log(sl) (longer movements when starting in mature)

id.landscape <- "complex"
id.variability <- "high"
id.rep <- 3

#_______________________________________________________________________
# 3c. Home ranging parameters ----
#_______________________________________________________________________

# potential "home range" centroids
# allow this to be a random uniform draw from x and y within the unit boundary

# extract bounding box
unit.bbox <- st_bbox(unit.bound)

# we'll take a random draw from both distributions as each path is initialized

# home range size variance
# during initial runs it appears that 5,000 gives reasonable HRs (~ 2 ha)
e.var <- 5000
e.var.sd <- 100

# we'll draw this value from a normal distribution with mean = 5000 and sd = 100

# the Bx and By coefficients are then calculated from the "home range" centroid and the variance

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
                      e.var.sd = e.var.sd,    # sd of the draws for the bivariate normal variance
                      hrc = hrc)              # home range centroid as previously drawn
  
{
  
  # variance
  var.focal <- rnorm(n = 1, mean = e.var, sd = e.var.sd)
  
  # x2 + y2 coefficient
  b.x2y2 <- -1 / var.focal
  
  # solve for x and y coefficients
  b.x <- hrc[1] * 2 * -b.x2y2 
  b.y <- hrc[2] * 2 * -b.x2y2
  
  hr.params <- c(b.x, b.y, b.x2y2)
  
  return(hr.params)
  
}

#_______________________________________________________________________
# 5. Run simulations iteratively ----

# n of replicates
n.reps <- 100

#_______________________________________________________________________
# 5a. Create dfs to hold all sims ----
#_______________________________________________________________________

sims.df <- data.frame()

all.coef.draws <- data.frame()

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
  # draw from appropriate distributions
  coef.sd <- ifelse(id.variability == "low",
                    0.10,
                    0.75)
  
  coef.draws <- data.frame(coef.stem.1 = rnorm(1, coef.stem, coef.sd),
                           coef.edge.1 = rnorm(1, coef.edge, coef.sd),
                           coef.mature.1 = rnorm(1, coef.mature, coef.sd))
  
  # here the terms are important to get right so redistribution_kernel() works okay
  issf.model <- make_issf_model(coefs = c("stem_end" = coef.draws$coef.stem.1[1],  
                                          "stem_end:log(sl_)" = coef.stem.sl,
                                          "edge_end" = coef.draws$coef.edge.1[1],
                                          "mature_start" = coef.draws$coef.mature.1[1],
                                          "log(sl_):mature_start" = coef.mature.sl,
                                          x2_ = hr.params[1],
                                          y2_ = hr.params[2], 
                                          "I(x2_^2 + y2_^2)" = hr.params[3]),              
                                sl = sl.dist,
                                ta = ta.dist)
  
  # initialize redistribution kernel
  rk <- redistribution_kernel(x = issf.model,
                              start = start.step,
                              map = landscape.covs.s,
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
           variability = id.variability,
           rep = id.rep)
  
  # bind to df
  sims.df <- bind_rows(sims.df, sim.path.1)
  
  all.coef.draws <- rbind(all.coef.draws, coef.draws)
  
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
            alpha = 0.10) +
  
  # remove legend
  theme(legend.position = "none") +
  
  # coordinate system to make nice axis labels
  coord_sf(datum = st_crs(32611))

#_______________________________________________________________________
# 7. Write to .csvs ----
#_______________________________________________________________________

write.csv(sims.df, paste0(getwd(), "/Derived_data/Simulated data/sims_C3H.csv"))
write.csv(all.coef.draws, paste0(getwd(), "/Derived_data/Simulated data/coefs_C3H.csv"))
