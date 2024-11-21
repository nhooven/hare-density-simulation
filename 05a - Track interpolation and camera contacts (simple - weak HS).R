# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 05a - Track interpolation and camera contacts (simple - weak HS)
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 21 Nov 2024
# Date completed: 
# Date last modified: 
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation
library(amt)             # work with tracks
library(sf)              # spatial data
library(ctmm)            # movement modeling

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

sims.df <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_simple_weak.csv"))

#_______________________________________________________________________
# 3. Create telemetry objects and fit CTSP models ----
#_______________________________________________________________________

for (i in unique(sims.df$indiv)) {
  
  # subset individual's data and create a track
  indiv.track <- sims.df %>% 
    
    filter(indiv == i) %>%
    
    # coerce t_ to posix
    mutate(t_ = as.POSIXct(t_)) %>%
    
    make_track(.x = x_, 
               .y = y_, 
               .t = t_, 
               crs = 32611)
  
  # convert to telemetry object
  indiv.telem <- as_telemetry(indiv.track,
                              timeformat = "auto",
                              timezone = "America/Los_Angeles",
                              projection = 32611,
                              keep = TRUE)
  
  # fit CTSP models
  # guesstimated model parameters from the variogram
  guess.param <- variogram.fit(variogram(indiv.telem), 
                               name = "guess.param", 
                               interactive = FALSE)
  
  # model selection
  fitted.mods <- ctmm.select(indiv.telem, 
                             CTMM = guess.param, 
                             verbose = TRUE)
  
  # baseline model (we'll use the Ornstein-Uhlenbeck Foraging process for simplicity)
  ctmm.model.1 <- ctmm(tau = guess.param$tau,
                       omega = guess.param$omega,
                       range = TRUE,
                       error = FALSE,
                       data = indiv.telem)
  
  # fit model
  ctmm.model.2 <- ctmm.fit(data = indiv.telem,
                           CTMM = ctmm.model.1)
  
  # interpolate track
  ctmm.interp <- simulate(ctmm.model.2,
                          data = indiv.telem,
                          res = 60,             # here we'll get a predicted location every ~2 minutes
                          complete = TRUE)
  
  # coerce to sf and transform to UTM
  ctmm.interp.sf <- ctmm.interp %>%
    
    as.sf() %>%
    
    st_transform(crs = "epsg:32611") %>%
    
    # mutate timestamp
    mutate(t = as.POSIXct(substr(timestamp, 1, 19),
                          tz = "America/Los_Angeles")) %>%
    
    # select only columns we need
    dplyr::select(t)
  

  
}



# cast to lines
interp.lines <- ctmm.interp.sf %>% 
  
  summarize(do_union = FALSE) %>%           # this is a critical step!
          
  st_cast(to = "LINESTRING")

  
  

ggplot() +
  
  theme_bw() +
  
  geom_point(data = indiv.track,
             aes(x = x_,
                 y = y_),
             size = 0.5) +
  
  geom_sf(data = interp.lines,
            alpha = 0.25) +
  
  coord_sf(datum = sf::st_crs(32611))

