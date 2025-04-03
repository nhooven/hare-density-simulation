# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 04a - CTMM simulations (fold 1)
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 25 Mar 2025
# Date completed: 26 Mar 2025
# Date last modified: 03 Apr 2025
# R version: 4.4.3

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation
library(lubridate)            # work with dates
library(ctmm)                 # CTSP movement modeling
library(amt)                  # work with animal movement tracks
library(sf)                   # spatial operations

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

# movement parameters
Q1.cam <- read.csv(paste0(getwd(), "/Derived data/Parameters/camera_lo.csv"))
Q1.col <- read.csv(paste0(getwd(), "/Derived data/Parameters/collar_lo.csv"))
Q2.cam <- read.csv(paste0(getwd(), "/Derived data/Parameters/camera_hi.csv"))
Q2.col <- read.csv(paste0(getwd(), "/Derived data/Parameters/collar_hi.csv"))

# subset to the correct fold
Q1.cam <- Q1.cam %>% filter(fold == 1)
Q1.col <- Q1.col %>% filter(fold == 1)
Q2.cam <- Q2.cam %>% filter(fold == 1)
Q2.col <- Q2.col %>% filter(fold == 1)

# unit boundary
unit.bound <- st_read(paste0(getwd(), "/Derived data/Shapefiles/unit_bound.shp"))

# camera viewsheds
vs <- st_read(paste0(getwd(), "/Derived data/Shapefiles/cams_vs.shp"))

vs$cam.id <- 1:nrow(vs)

st_crs(vs) <- "epsg:32611"

#_______________________________________________________________________
# 3. Simulation parameters ----

# all camera contact sims will last 4 weeks and will be sampled at a 60-sec "fundamental step"

#_______________________________________________________________________

# simulation duration
sim.duration.wk <- 4

sim.duration.sec <- sim.duration.wk * 7 * 24 * 60 * 60 

# "sampling rate" for fundamental (i.e., straight line) steps
sim.samp.rate <- 60

# time steps to simulate on (from, to, by)
sim.timestep <- seq(1, sim.duration.sec, sim.samp.rate)

# collar sims ONLY
sim.duration.halfhr <- sim.duration.wk * 7 * 24 * 2

sim.timestep.halfhr <- seq(1, sim.duration.sec, 3600 / 2)

#_______________________________________________________________________
# 4. Loop ----

# number of iterations (let's stop calling them replicates)
n.iter <- 1

#_______________________________________________________________________
# 4a. Define functions ----
#_______________________________________________________________________

# simulations for camera contacts
cam_sim_contact <- function (focal.row) {
  
  # output list (will hold 3 dfs)
  focal.list <- list()
  
  # define OUF model
  ouf.mod <- ctmm(tau = c(8 %#% "hours",             # positional autocorr time
                          1 %#% "hours"),            # velocity autocorr time
                  isotropic = FALSE,                  # anisotropic
                  sigma = c(5000,      # variance along major axis
                            5000,      # variance along minor axis
                            0),         # angle of axis
                  mu = c(0 + rnorm(1, 0, 50),             # x coord
                         0 + rnorm(1, 0, 50)))            # y coord
  
  # run simulation
  ouf.sim <- simulate(object = ouf.mod, 
                      t = sim.timestep,
                      complete = T)
  
  # coerce telemetry object to sf
  telem.df <- data.frame(t = ouf.sim$t,
                         x = ouf.sim$x,
                         y = ouf.sim$y,
                         timestamp = as.POSIXct(ouf.sim$timestamp))
  
  telem.sf <- st_as_sf(telem.df,
                       coords = c("x", "y"),
                       crs = "epsg:32611")
  
  # tally camera contacts ("passes")
  # in this case, 60-second locations that fall within each viewshed
  # each "pass" through a viewshed will only count as one
  passes <- st_intersection(vs, 
                            telem.sf) %>%
    
    # drop the geometry
    st_drop_geometry() %>% 
    
    # group by camera
    group_by(cam.id) %>% 
    
    # add a time difference column
    #mutate(time.diff = as.numeric(t - lag(t))) %>%
    
    # keep intersections that are not part of the same "pass"
    #filter(time.diff > 60 |
    #       is.na(time.diff) == T) %>%
    
    # tally all intersections
    tally() %>%
    
    # and add identifiers
    mutate(iter = focal.row$iter,
           indiv = focal.row$indiv,
           Q = focal.row$Q,
           use = focal.row$use)
  
  # sample to a "GPS track"
  # this will be a regular 1-hr track, which can be subsampled
  telem.track <- telem.df %>%
    
    make_track(.x = x,
               .y = y,
               .t = timestamp)
  
  # calculate "true" speed (m/s)
  telem.track.speed <- telem.track %>%
    
    steps() 
  
  # 4 week
  speed.4wk <- sum(telem.track.speed$sl_) / sim.duration.sec
  
  # resample track for telemetering
  telem.track.1 <- telem.track %>%
    
    track_resample(rate = minutes(30),
                   tolerance = minutes(0)) %>%
    
    # remove "burst"
    dplyr::select(-burst_) %>%
    
    # and add identifiers
    mutate(iter = focal.row$iter,
           indiv = focal.row$indiv,
           Q = focal.row$Q,
           use = focal.row$use)
  
  # speeds df
  focal.speeds <- data.frame(iter = focal.row$iter,
                             indiv = focal.row$indiv,
                             Q = focal.row$Q,
                             use = focal.row$use,
                             speed.4wk = speed.4wk)
  
  # bind each into list
  focal.list[[1]] <- passes
  focal.list[[2]] <- telem.track.1
  focal.list[[3]] <- focal.speeds
  
  # return
  return(focal.list)
  
}

# simulations for collared individuals only
col_sim <- function (focal.row) {
  
  # output list (will hold 3 dfs)
  focal.list <- list()
  
  # define OUF model
  ouf.mod <- ctmm(tau = c(focal.row$tau1,             # positional autocorr time
                          focal.row$tau2),            # velocity autocorr time
                  isotropic = FALSE,                  # anisotropic
                  sigma = c(focal.row$sigma.maj,      # variance along major axis
                            focal.row$sigma.min,      # variance along minor axis
                            focal.row$angle),         # angle of axis
                  mu = c(focal.row$mean1,             # x coord
                         focal.row$mean2))            # y coord
  
  # run simulation
  ouf.sim <- simulate(object = ouf.mod, 
                      t = sim.timestep.halfhr,        # half-hourly only
                      complete = T)
  
  # coerce telemetry object to sf
  telem.df <- data.frame(t = ouf.sim$t,
                         x = ouf.sim$x,
                         y = ouf.sim$y,
                         timestamp = as.POSIXct(ouf.sim$timestamp))
  
  telem.sf <- st_as_sf(telem.df,
                       coords = c("x", "y"),
                       crs = "epsg:32611")
  
  # convert to a "GPS track"
  # this will be a regular 1-hr track
  telem.track <- telem.df %>%
    
    make_track(.x = x,
               .y = y,
               .t = timestamp) %>%
    
    # and add identifiers
    mutate(iter = focal.row$iter,
           indiv = focal.row$indiv,
           Q = focal.row$Q,
           use = focal.row$use)
  
  # return
  return(telem.track)
  
}

#_______________________________________________________________________
# 4b. Initialize dfs ----
#_______________________________________________________________________

Q1.cam.all.passes <- data.frame()
Q1.cam.all.relocs <- data.frame()
Q1.cam.all.speeds <- data.frame()

Q1.col.all.relocs <- data.frame()

Q2.cam.all.passes <- data.frame()
Q2.cam.all.relocs <- data.frame()
Q2.cam.all.speeds <- data.frame()

Q2.col.all.relocs <- data.frame()

#_______________________________________________________________________
# 4c. Run loop ----
#_______________________________________________________________________

# loop through all iterations (i)

# start time
start.time <- Sys.time()

for (i in 1:n.iter) {
  
  # subset iteration from each df
  focal.Q1.cam <- Q1.cam %>% filter(iter == i)
  focal.Q1.col <- Q1.col %>% filter(iter == i)
  focal.Q2.cam <- Q2.cam %>% filter(iter == i)
  focal.Q2.col <- Q2.col %>% filter(iter == i)
  
  # loop through all individuals, by dataset (w, x, y, z)
  # Q1.cam (w)
  for (w in 1:nrow(focal.Q1.cam)) {
    
    # subset row
    focal.row.w <- focal.Q1.cam %>% slice(w)
    
    # use function
    focal.sims.w <- cam_sim_contact(focal.row.w)
    
    # bind in
    Q1.cam.all.passes <- rbind(Q1.cam.all.passes, focal.sims.w[[1]])
    Q1.cam.all.relocs <- rbind(Q1.cam.all.relocs, focal.sims.w[[2]])
    Q1.cam.all.speeds <- rbind(Q1.cam.all.speeds, focal.sims.w[[3]])
    
  }
  
  # Q1.col (x)
  for (x in 1:nrow(focal.Q1.col)) {
    
    # subset row
    focal.row.x <- focal.Q1.col %>% slice(x)
    
    # use function
    focal.sims.x <- col_sim(focal.row.x)
    
    # bind in
    Q1.col.all.relocs <- rbind(Q1.col.all.relocs, focal.sims.x)
    
  }
  
  # Q2.cam (y)
  for (y in 1:nrow(focal.Q2.cam)) {
    
    # subset row
    focal.row.y <- focal.Q2.cam %>% slice(y)
    
    # use function
    focal.sims.y <- cam_sim_contact(focal.row.y)
    
    # bind in
    Q2.cam.all.passes <- rbind(Q2.cam.all.passes, focal.sims.y[[1]])
    Q2.cam.all.relocs <- rbind(Q2.cam.all.relocs, focal.sims.y[[2]])
    Q2.cam.all.speeds <- rbind(Q2.cam.all.speeds, focal.sims.y[[3]])
    
  }
  
  # Q2.col (z)
  for (z in 1:nrow(focal.Q2.col)) {
    
    # subset row
    focal.row.z <- focal.Q2.col %>% slice(z)
    
    # use function
    focal.sims.z <- col_sim(focal.row.z)
    
    # bind in
    Q2.col.all.relocs <- rbind(Q2.col.all.relocs, focal.sims.z)
    
  }
  
  # every 50 iterations:
    # save each to file for my sanity
    # print status message
  
  if (i %% 50 == 0) {
    
    # save each to file
    write.csv(Q1.cam.all.passes, paste0(getwd(), "/Derived data/Camera detections/detections_lo1.csv")) 
    write.csv(Q1.cam.all.relocs, paste0(getwd(), "/Derived data/Simulated tracks/tracks_lo1.csv"))
    write.csv(Q1.cam.all.speeds, paste0(getwd(), "/Derived data/Speeds/speeds_lo1.csv"))
    
    write.csv(Q1.col.all.relocs, paste0(getwd(), "/Derived data/Simulated tracks/tracks_col_lo1.csv"))
    
    write.csv(Q2.cam.all.passes, paste0(getwd(), "/Derived data/Camera detections/detections_hi1.csv"))
    write.csv(Q2.cam.all.relocs, paste0(getwd(), "/Derived data/Simulated tracks/tracks_hi1.csv"))
    write.csv(Q2.cam.all.speeds, paste0(getwd(), "/Derived data/Speeds/speeds_hi1.csv"))
    
    write.csv(Q2.col.all.relocs, paste0(getwd(), "/Derived data/Simulated tracks/tracks_col_hi1.csv"))
    
    # print status message
    elapsed.time <- round(as.numeric(difftime(Sys.time(), 
                                              start.time, 
                                              units = "mins")), 
                          digits = 1)
    
    print(paste0("Completed sim ", i, " of ", n.iter, " - ", elapsed.time, " mins"))
    
  }
  
}

#_______________________________________________________________________
# 5. Save all to file ----
#_______________________________________________________________________

write.csv(Q1.cam.all.passes, paste0(getwd(), "/Derived data/Camera detections/detections_lo1.csv")) 
write.csv(Q1.cam.all.relocs, paste0(getwd(), "/Derived data/Simulated tracks/tracks_lo1.csv"))
write.csv(Q1.cam.all.speeds, paste0(getwd(), "/Derived data/Speeds/speeds_lo1.csv"))

write.csv(Q1.col.all.relocs, paste0(getwd(), "/Derived data/Simulated tracks/tracks_col_lo1.csv"))

write.csv(Q2.cam.all.passes, paste0(getwd(), "/Derived data/Camera detections/detections_hi1.csv"))
write.csv(Q2.cam.all.relocs, paste0(getwd(), "/Derived data/Simulated tracks/tracks_hi1.csv"))
write.csv(Q2.cam.all.speeds, paste0(getwd(), "/Derived data/Speeds/speeds_hi1.csv"))

write.csv(Q2.col.all.relocs, paste0(getwd(), "/Derived data/Simulated tracks/tracks_col_hi1.csv"))
