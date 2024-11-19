# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 02a - Movement simulations (simple - weak HS)
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 15 Nov 2024
# Date completed: 
# Date last modified: 
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation
library(terra)           # rasters
library(amt)             # simulate tracks
library(lubridate)       # work with time

#_______________________________________________________________________
# 2. Read in rasters ----
#_______________________________________________________________________

landscape.covs <- rast("Rasters/simple.tif")

#_______________________________________________________________________
# 3. Define simulation parameters ----
#_______________________________________________________________________
# 3a. Movement parameter distributions ----

# here we'll just approximate hare movement from preliminary data

#_______________________________________________________________________

# step lengths (gamma)
sl.df <- data.frame(x = seq(0,
                            750,
                            length.out = 1000),
                    y = dgamma(x = seq(0,
                                       750,
                                       length.out = 1000),
                               shape = 1.2,
                               scale = 70))

# plot sl
ggplot() +
  
  theme_bw() +
  
  # gamma
  geom_line(data = sl.df,
            aes(x = x,
                y = y))

# make distribution
sl.dist <- make_gamma_distr(shape = 1.2,
                            scale = 70)

# turning angles (uniform)
# make distribution
ta.dist <- make_unif_distr(min = -pi,
                           max = pi)

#_______________________________________________________________________
# 3b. Make iSSF model and redistribution kernel ----
#_______________________________________________________________________

# make an iSSF model
issf.model <- make_issf_model(coefs = c("stem_end" = 0.5,               # must specify start and end here
                                        "edge_end" = -0.5,
                                        "mature_start" = -0.1,
                                        "mature_start:log(sl_)" = 0.5),
                              sl = sl.dist,
                              ta = ta.dist)

#_______________________________________________________________________
# 4. Initialize simulations ----
#_______________________________________________________________________

# make a starting step
start.step <- make_start(x = c(1119, 1119),
                         ta_ = 0,
                         time = ymd_hm("2024-09-01 18:00", tz = "America/Los_Angeles"),
                         dt = hours(2),
                         crs = crs("EPSG:32611"))

# make redistribution kernel
redist <- redistribution_kernel(x = issf.model,
                                start = start.step,
                                map = landscape.covs,
                                max.dist = get_max_dist(issf.model),
                                tolerance.outside = 0.0008)

#_______________________________________________________________________
# 5. Run a test simulation ----
#_______________________________________________________________________

test.path <- simulate_path(redist,
                           n.steps = 20,
                           start = start.step,
                           verbose = TRUE)

plot(as_track(test.path))

#_______________________________________________________________________
# 6. Include a distance to HR center covariate to incorporate ranging ----
#_______________________________________________________________________

# define the HR center from the start.step
hrc <- st_as_sf(start.step,
                coords = c("x_", "y_"),
                crs = 32611)

hrc <- vect(hrc)

hrc.dist <- distance(rasterize(hrc, landscape.covs$stem, mask = T))

plot(hrc.dist)

hrc.dist.s <- scale(hrc.dist)

# here we'll use this procedure to allow for home ranging before running all sims
landscape.covs.1 <- c(landscape.covs, hrc.dist.s)

names(landscape.covs.1)[4] <- "hrc"

issf.model.1 <- make_issf_model(coefs = c("stem_end" = 0.5,               
                                          "edge_end" = -0.5,
                                          "mature_start" = -0.1,
                                          "mature_start:log(sl_)" = 0.5,
                                          "hrc_end" = -2.0),
                                sl = sl.dist,
                                ta = ta.dist)

# make redistribution kernel
redist.1 <- redistribution_kernel(x = issf.model.1,
                                  start = start.step,
                                  map = landscape.covs.1,
                                  n.control = 4e5,   # approximate pixel number
                                  max.dist = get_max_dist(issf.model.1),
                                  tolerance.outside = 0.0008)

# run another sim
start.time <- Sys.time()
test.path.1 <- simulate_path(redist.1,
                             n.steps = 20,
                             start = start.step,
                             verbose = TRUE)

Sys.time() - start.time

# 31 secs for 20 steps @ 1,000 control steps

plot(as_ltraj(as_track(test.path.1)))

# secs per step
(sec.step <- 31 / 20)

(sec.step * 336 * 50) / (60 * 60)

# let's try the full month
# run another sim
start.time <- Sys.time()
test.path.2 <- simulate_path(redist.1,
                             n.steps = 336,
                             start = start.step,
                             verbose = TRUE)

Sys.time() - start.time

plot(as_ltraj(as_track(test.path.2)))

test.track.lines <- as_sf_lines(as_track(test.path.2))

plot(landscape.covs.1$edge)
plot(test.track.lines, add = T)

# home range is way too big