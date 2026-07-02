# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 04 - Track simulations
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 25 Mar 2025
# Date completed: 26 Mar 2025
# Date last modified: 02 Jul 2026
# R version: 4.5.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation
library(ctmm)                 # CTSP movement modeling
library(amt)                  # work with animal movement tracks
library(sf)                   # spatial operations
library(mefa4)                # %notin%

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

# models/model parameters
model.TV1 <- readRDS("data_derived/Q1_models/tv1.rds")
model.TV2 <- readRDS("data_derived/Q1_models/tv2.rds")
model.TV3 <- readRDS("data_derived/Q1_models/tv3.rds")

Q2.params <- readRDS("data_derived/Q2_models/Q2_params.rds")

# sampled individuals
Q1.indivs <- readRDS("data_derived/sampled_indivs/Q1_indivs.rds")
Q2.indivs <- readRDS("data_derived/sampled_indivs/Q2_indivs.rds")

# camera viewsheds
vs <- st_read(paste0(getwd(), "/data_derived/Shapefiles/cams_vs.shp"))

vs$cam.id <- 1:nrow(vs)

st_crs(vs) <- "epsg:32611"

#_______________________________________________________________________
# 3. Simulation parameters ----

# all camera contact sims will last 4 weeks and will be sampled at a 60-sec "fundamental step"
# then, we'll save tracks subsampled to 0.5 h (the finest fix rate we'll test, ~ 0.5 tau_v)

#_______________________________________________________________________

# simulation duration (in seconds)
sim.duration.sec <- 28 * 24 * 60 * 60 

# "sampling rate" for fundamental (i.e., straight line) steps
sim.samp.rate <- 60

# time steps to simulate on (from, to, by)
sim.timestep <- seq(1, sim.duration.sec, sim.samp.rate)

# finest collar fix interval
sim.duration.halfhr <- 28 * 24 * 2

sim.timestep.halfhr <- seq(1, sim.duration.sec, 3600 / 2)

#_______________________________________________________________________
# 4. Simulation functions ----
#_______________________________________________________________________
# 4a. sim_tracks_Q1 ----

# by individual

# (1) Simulate a "fundamental" track
# (2) Calculate "true" speeds (total distance / time)
# (3) Tally camera contacts (both total points within camera viewsheds and separate detections)
# (4) Subset to 0.5 hr "fixes" for movement modeling

# this will yield a list in which:
# [[1]] = df of camera contacts (500 x 9 cams)
# [[2]] = df of true speeds (500)
# [[3]] = df of "GPS fixes" (500 x 1344 fixes)

#_______________________________________________________________________

sim_tracks_Q1 <- function (.df, .ctmm) {
  
  # model
  focal.model <- .ctmm
  
  # define list
  sim.track.list <- list()
  
  # define dfs
  all.contacts <- data.frame()
  all.speeds <- data.frame()
  all.fixes <- data.frame()
  
  # start time
  start.time <- Sys.time()
  
  # loop through all individuals
  for (i in 1:nrow(.df)) {
    
    # subset indiv
    focal.indiv <- .df |> slice(i)
    
    # modify movement model 
    # assign "home range centroid"
    focal.model$mu <- as.numeric(focal.indiv[ , c("hrc.x", "hrc.y")])
    
    # assign angle
    focal.model$sigma@par[3] <- focal.indiv$angle
    
    # remove UERE (assume that all relocations are measured without error)
    focal.model$UERE <- c(0, 0, 0)
    
    # simulate month-long track
    model.sim <- simulate(object = focal.model, 
                          t = sim.timestep,
                          complete = T)
    
    # coerce telemetry object to df, then sf
    telem.df <- data.frame(t = model.sim$t,
                           x = model.sim$x,
                           y = model.sim$y,
                           timestamp = as.POSIXct(model.sim$timestamp))
    
    telem.sf <- st_as_sf(telem.df,
                         coords = c("x", "y"),
                         crs = "epsg:32611")
    
    # tally camera contacts
    # suppress warnings 
    st_agr(vs) <- "constant"
    st_agr(telem.sf) <- "constant"
    
    # we'll tally both:
    
    # POINTS - total minute-by-minute locations within a camera viewshed
    cam.contacts.points <- st_intersection(vs, 
                                           telem.sf) |>
      
      # drop the geometry
      st_drop_geometry() |>
      
      # group by camera
      group_by(cam.id) |>
      
      # tally intersections
      summarize(points = n())
    
    # PASSES - separate passes of potentially multiple points
    cam.contacts.passes <- st_intersection(vs, 
                                           telem.sf) |>
      
      # drop the geometry
      st_drop_geometry() |>
      
      # group by camera
      group_by(cam.id) |>
      
      # add a time difference column
      mutate(time.diff = as.numeric(t - lag(t))) |>
      
      # keep intersections that are not part of the same "pass"
      filter(time.diff > 60 |
               is.na(time.diff) == T) |>
      
      # tally intersections
      summarize(passes = n())
    
    # join together
    cam.contacts <- cam.contacts.points |>
      
      left_join(cam.contacts.passes,
                by = "cam.id") 
    
    # add cameras with zero passes (thus also zero points)
    # BUT ONLY if there are cams with zero passes, should be the common case
    if (length(which(1:9 %notin% cam.contacts$cam.id)) > 0) {
      
      cam.contacts.1 <- cam.contacts |>
        
        bind_rows(data.frame(cam.id = which(1:9 %notin% cam.contacts$cam.id),
                             points = 0,
                             passes = 0)) |>
        
        arrange(cam.id) |>
        
        # and add identifiers
        mutate(question = focal.indiv$question,
               target = focal.indiv$target,
               indiv = focal.indiv$indiv)
      
    } else {
      
      cam.contacts.1 <- cam.contacts |>
        
        arrange(cam.id) |>
        
        # and add identifiers
        mutate(question = focal.indiv$question,
               target = focal.indiv$target,
               indiv = focal.indiv$indiv)
      
    }
    
    # sample to a "GPS track"
    # this will be a regular 0.5-hr track, which can be subsampled
    telem.track <- telem.df |>
      
      make_track(.x = x,
                 .y = y,
                 .t = timestamp)
    
    # calculate "true" speed (m/s)
    telem.track.speed <- telem.track |> steps() 
    
    true.speed <- sum(telem.track.speed$sl_) / sim.duration.sec
    
    # resample track for telemetering
    telem.track.1 <- telem.track |>
      
      track_resample(rate = minutes(30),
                     tolerance = minutes(0)) |>
      
      # remove "burst"
      dplyr::select(-burst_) |>
      
      # and add identifiers
      mutate(question = focal.indiv$question,
             target = focal.indiv$target,
             indiv = focal.indiv$indiv)
    
    # speeds df
    focal.speeds <- data.frame(question = focal.indiv$question,
                               target = focal.indiv$target,
                               indiv = focal.indiv$indiv,
                               true.speed = true.speed)
    
    # bind each into df
    all.contacts <- rbind(all.contacts, cam.contacts.1)
    all.speeds <- rbind(all.speeds, focal.speeds)
    all.fixes <- rbind(all.fixes, telem.track.1)
    
    # print status message after every 50 individuals
    if (i %% 50 == 0) {
      
      elapsed.time <- round(as.numeric(difftime(Sys.time(), 
                                                start.time, 
                                                units = "mins")), 
                            digits = 1)
      
      print(paste0("Completed indiv ", 
                   i, 
                   " of ", 
                   nrow(.df), 
                   " - ", 
                   elapsed.time, 
                   " mins"))
      
    }
    
  }
  
  # bind into list
  sim.track.list[[1]] <- all.contacts
  sim.track.list[[2]] <- all.speeds
  sim.track.list[[3]] <- all.fixes
  
  # return
  return(sim.track.list)
  
}

#_______________________________________________________________________
# 4b. sim_tracks_Q2 ----
#_______________________________________________________________________

sim_tracks_Q2 <- function (.df) {
  
  # define list
  sim.track.list <- list()
  
  # define dfs
  all.contacts <- data.frame()
  all.speeds <- data.frame()
  all.fixes <- data.frame()
  
  # start time
  start.time <- Sys.time()
  
  # loop through all individuals
  for (i in 1:nrow(.df)) {
    
    # subset indiv
    focal.indiv <- .df |> slice(i)
    
    # define movement model 
    focal.params <- Q2.params |> slice(i)
    
    focal.model <- ctmm(
      
      tau = c(focal.params$tau.p,
              focal.params$tau.v),
      isotropic = FALSE,
      range = T,
      error = FALSE,
      sigma = matrix(c(focal.params$sigma.major,
                       0,
                       0,
                       focal.params$sigma.minor),
                     nrow = 2),
      mu = as.numeric(focal.indiv[ , c("hrc.x", "hrc.y")])
      
    )
    
    # assign angle
    focal.model$sigma@par[3] <- focal.indiv$angle
    
    # remove UERE (assume that all relocations are measured without error)
    focal.model$UERE <- c(0, 0, 0)
    
    # simulate month-long track
    model.sim <- simulate(object = focal.model, 
                          t = sim.timestep,
                          complete = T)
    
    # coerce telemetry object to df, then sf
    telem.df <- data.frame(t = model.sim$t,
                           x = model.sim$x,
                           y = model.sim$y,
                           timestamp = as.POSIXct(model.sim$timestamp))
    
    telem.sf <- st_as_sf(telem.df,
                         coords = c("x", "y"),
                         crs = "epsg:32611")
    
    # tally camera contacts
    # suppress warnings 
    st_agr(vs) <- "constant"
    st_agr(telem.sf) <- "constant"
    
    # we'll tally both:
    
    # POINTS - total minute-by-minute locations within a camera viewshed
    cam.contacts.points <- st_intersection(vs, 
                                           telem.sf) |>
      
      # drop the geometry
      st_drop_geometry() |>
      
      # group by camera
      group_by(cam.id) |>
      
      # tally intersections
      summarize(points = n())
    
    # PASSES - separate passes of potentially multiple points
    cam.contacts.passes <- st_intersection(vs, 
                                           telem.sf) |>
      
      # drop the geometry
      st_drop_geometry() |>
      
      # group by camera
      group_by(cam.id) |>
      
      # add a time difference column
      mutate(time.diff = as.numeric(t - lag(t))) |>
      
      # keep intersections that are not part of the same "pass"
      filter(time.diff > 60 |
               is.na(time.diff) == T) |>
      
      # tally intersections
      summarize(passes = n())
    
    # join together
    cam.contacts <- cam.contacts.points |>
      
      left_join(cam.contacts.passes,
                by = "cam.id") 
    
    # add cameras with zero passes (thus also zero points)
    # BUT ONLY if there are cams with zero passes, should be the common case
    if (length(which(1:9 %notin% cam.contacts$cam.id)) > 0) {
      
      cam.contacts.1 <- cam.contacts |>
        
        bind_rows(data.frame(cam.id = which(1:9 %notin% cam.contacts$cam.id),
                             points = 0,
                             passes = 0)) |>
        
        arrange(cam.id) |>
        
        # and add identifiers
        mutate(question = focal.indiv$question,
               target = focal.indiv$target,
               indiv = focal.indiv$indiv)
      
    } else {
      
      cam.contacts.1 <- cam.contacts |>
        
        arrange(cam.id) |>
        
        # and add identifiers
        mutate(question = focal.indiv$question,
               target = focal.indiv$target,
               indiv = focal.indiv$indiv)
      
    }
    
    # sample to a "GPS track"
    # this will be a regular 0.5-hr track, which can be subsampled
    telem.track <- telem.df |>
      
      make_track(.x = x,
                 .y = y,
                 .t = timestamp)
    
    # calculate "true" speed (m/s)
    telem.track.speed <- telem.track |> steps() 
    
    true.speed <- sum(telem.track.speed$sl_) / sim.duration.sec
    
    # resample track for telemetering
    telem.track.1 <- telem.track |>
      
      track_resample(rate = minutes(30),
                     tolerance = minutes(0)) |>
      
      # remove "burst"
      dplyr::select(-burst_) |>
      
      # and add identifiers
      mutate(question = focal.indiv$question,
             target = focal.indiv$target,
             indiv = focal.indiv$indiv)
    
    # speeds df
    focal.speeds <- data.frame(question = focal.indiv$question,
                               target = focal.indiv$target,
                               indiv = focal.indiv$indiv,
                               true.speed = true.speed)
    
    # bind each into df
    all.contacts <- rbind(all.contacts, cam.contacts.1)
    all.speeds <- rbind(all.speeds, focal.speeds)
    all.fixes <- rbind(all.fixes, telem.track.1)
    
    # print status message after every 50 individuals
    if (i %% 50 == 0) {
      
      elapsed.time <- round(as.numeric(difftime(Sys.time(), 
                                                start.time, 
                                                units = "mins")), 
                            digits = 1)
      
      print(paste0("Completed indiv ", 
                   i, 
                   " of ", 
                   nrow(.df), 
                   " - ", 
                   elapsed.time, 
                   " mins"))
      
    }
    
  }
  
  # bind into list
  sim.track.list[[1]] <- all.contacts
  sim.track.list[[2]] <- all.speeds
  sim.track.list[[3]] <- all.fixes
  
  # return
  return(sim.track.list)
  
}
  
#_______________________________________________________________________
# 5. Use simulation function ----

# benchmark 09 Apr 2025: ~ 42 sec / 10 indivs (35 min / 500 indivs)

#_______________________________________________________________________
# 5a. Q1 ----
#_______________________________________________________________________

# subset dfs
Q1.1T.TV1    <- Q1.indivs |> filter(target == "T" & TV == 1)
Q1.1T.TV2    <- Q1.indivs |> filter(target == "T" & TV == 2)
Q1.1T.TV3    <- Q1.indivs |> filter(target == "T" & TV == 3)

Q1.1NT.TV1    <- Q1.indivs |> filter(target == "NT" & TV == 1)
Q1.1NT.TV2    <- Q1.indivs |> filter(target == "NT" & TV == 2)
Q1.1NT.TV3    <- Q1.indivs |> filter(target == "NT" & TV == 3)

# use function
sim.tracks.1T.TV1 <- sim_tracks_Q1(Q1.1T.TV1, model.TV1) # done
sim.tracks.1T.TV2 <- sim_tracks_Q1(Q1.1T.TV2, model.TV2) # done
sim.tracks.1T.TV3 <- sim_tracks_Q1(Q1.1T.TV3, model.TV3) # done

sim.tracks.1NT.TV1 <- sim_tracks_Q1(Q1.1NT.TV1, model.TV1) # done
sim.tracks.1NT.TV2 <- sim_tracks_Q1(Q1.1NT.TV2, model.TV2) # done
sim.tracks.1NT.TV3 <- sim_tracks_Q1(Q1.1NT.TV3, model.TV3) # done

#_______________________________________________________________________
# 5b. Q2 ----
#_______________________________________________________________________

Q2.T <- Q2.indivs |> filter(target == "T")
Q2.NT <- Q2.indivs |> filter(target == "NT")

sim.tracks.2T <- sim_tracks_Q2(Q2.T)    # done
sim.tracks.2NT <- sim_tracks_Q2(Q2.NT)  # done

#_______________________________________________________________________
# 6. Write to .csvs ----
#_______________________________________________________________________
# 6a. Q1 ----
#_______________________________________________________________________

# contacts
saveRDS(sim.tracks.1T.TV1[[1]], "data_derived/sampled_contacts/contacts_1T_TV1.rds")
saveRDS(sim.tracks.1T.TV2[[1]], "data_derived/sampled_contacts/contacts_1T_TV2.rds")
saveRDS(sim.tracks.1T.TV3[[1]], "data_derived/sampled_contacts/contacts_1T_TV3.rds")

saveRDS(sim.tracks.1NT.TV1[[1]], "data_derived/sampled_contacts/contacts_1NT_TV1.rds")
saveRDS(sim.tracks.1NT.TV2[[1]], "data_derived/sampled_contacts/contacts_1NT_TV2.rds")
saveRDS(sim.tracks.1NT.TV3[[1]], "data_derived/sampled_contacts/contacts_1NT_TV3.rds")

# speeds
saveRDS(sim.tracks.1T.TV1[[2]], "data_derived/sampled_speeds/speeds_1T_TV1.rds")
saveRDS(sim.tracks.1T.TV2[[2]], "data_derived/sampled_speeds/speeds_1T_TV2.rds")
saveRDS(sim.tracks.1T.TV3[[2]], "data_derived/sampled_speeds/speeds_1T_TV3.rds")

saveRDS(sim.tracks.1NT.TV1[[2]], "data_derived/sampled_speeds/speeds_1NT_TV1.rds")
saveRDS(sim.tracks.1NT.TV2[[2]], "data_derived/sampled_speeds/speeds_1NT_TV2.rds")
saveRDS(sim.tracks.1NT.TV3[[2]], "data_derived/sampled_speeds/speeds_1NT_TV3.rds")

# tracks
saveRDS(sim.tracks.1T.TV1[[3]], "data_derived/sampled_tracks/tracks_1T_TV1.rds")
saveRDS(sim.tracks.1T.TV2[[3]], "data_derived/sampled_tracks/tracks_1T_TV2.rds")
saveRDS(sim.tracks.1T.TV3[[3]], "data_derived/sampled_tracks/tracks_1T_TV3.rds")

saveRDS(sim.tracks.1NT.TV1[[3]], "data_derived/sampled_tracks/tracks_1NT_TV1.rds")
saveRDS(sim.tracks.1NT.TV2[[3]], "data_derived/sampled_tracks/tracks_1NT_TV2.rds")
saveRDS(sim.tracks.1NT.TV3[[3]], "data_derived/sampled_tracks/tracks_1NT_TV3.rds")

#_______________________________________________________________________
# 6b. Q2 ----
#_______________________________________________________________________

# contacts
saveRDS(sim.tracks.2T[[1]], "data_derived/sampled_contacts/contacts_2T.rds")
saveRDS(sim.tracks.2NT[[1]], "data_derived/sampled_contacts/contacts_2NT.rds")

# speeds
saveRDS(sim.tracks.2T[[2]], "data_derived/sampled_speeds/speeds_2T.rds")
saveRDS(sim.tracks.2NT[[2]], "data_derived/sampled_speeds/speeds_2NT.rds")

# tracks
saveRDS(sim.tracks.2T[[3]], "data_derived/sampled_tracks/tracks_2T.rds")
saveRDS(sim.tracks.2NT[[3]], "data_derived/sampled_tracks/tracks_2NT.rds")
