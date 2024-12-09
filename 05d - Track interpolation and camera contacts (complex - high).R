# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 05d - Track interpolation and camera contacts (complex - high)
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 21 Nov 2024
# Date completed: 26 Nov 2024
# Date last modified: 09 Dec 2024
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

sims.df <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_complex_high.csv"))

# viewsheds
vs <- st_read(paste0(getwd(), "/Derived_data/Shapefiles/cams_9_vs.shp"))

# add camera ID
vs$cam.id <- 1:nrow(vs)

#_______________________________________________________________________
# 3. Create telemetry objects and fit CTSP models ----
#_______________________________________________________________________
# 3a. Initialize a list and a data.frame ----
#_______________________________________________________________________

all.ctmms <- vector("list", length(unique(sims.df$indiv)))

all.passes <- data.frame()

#_______________________________________________________________________
# 3b. Loop through all individuals ----
#_______________________________________________________________________

start.time <- Sys.time()

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
  
  # baseline model (we'll use the Ornstein-Uhlenbeck process)
  ctmm.model.1 <- ctmm(tau = guess.param$tau[1],
                       omega = FALSE,
                       range = TRUE,
                       error = FALSE,
                       isotropic = TRUE,
                       data = indiv.telem)
  
  # fit model
  ctmm.model.2 <- ctmm.fit(data = indiv.telem,
                           CTMM = ctmm.model.1)
  
  # save ctmm object into a list
  all.ctmms[[i]] <- ctmm.model.2
  
  # interpolate track
  ctmm.interp <- simulate(ctmm.model.2,
                          data = indiv.telem,
                          res = 30,             # here we'll get a predicted location every ~4 minutes
                          complete = TRUE)
  
  # coerce to sf and cast to lines
  interp.lines <- ctmm.interp %>%
    
    as.sf() %>%
    
    # transform to UTM
    st_transform(crs = "epsg:32611") %>%
    
    # mutate timestamp
    mutate(t = as.POSIXct(substr(timestamp, 1, 19),
                          tz = "America/Los_Angeles")) %>%
    
    # select only columns we need
    dplyr::select(t) %>%
    
    # duplicate each row so lines are connected
    slice(rep(1:n(), each = 2)) %>%
    
    # and now drop the first and last rows
    slice(-c(1, n()))
  
  # add pairings to group by
  interp.lines$id <- rep(1:(nrow(interp.lines) / 2), each = 2) 
  
  # group and create linestrings
  interp.lines.1 <- interp.lines %>%
    
    group_by(id) %>%
    
    summarize(do_union = FALSE) %>%
    
    st_cast("LINESTRING")
  
  # tallied intersections by camera
  passes <- st_intersection(vs, interp.lines.1) %>%
    
    # drop the geometry
    st_drop_geometry() %>% 
    
    # group by camera
    group_by(cam.id) %>% 
    
    # tally all intersections
    tally() %>%
    
    # add individual id
    mutate(indiv = i)
  
  # bind into df
  all.passes <- rbind(all.passes, passes)
  
  # status message
  elapsed.time <- round(as.numeric(difftime(Sys.time(), 
                                            start.time, 
                                            units = "mins")),
                        digits = 1)
  
  print(paste0("Completed passes ", 
               i, 
               " of ", 
               length(unique(sims.df$indiv)), 
               " - ", 
               elapsed.time, 
               " mins"))
  
}

#_______________________________________________________________________
# 4. Example plots ----
#_______________________________________________________________________

ggplot() +
  
  theme_bw() +
  
  geom_point(data = indiv.track,
             aes(x = x_,
                 y = y_),
             size = 0.5) +
  
  geom_sf(data = interp.lines.1,
          alpha = 0.15) +
  
  geom_sf(data = vs,
          color = "gold",
          fill = NA) +
  
  geom_sf_text(data = vs,
               aes(label = cam.id),
               color = "gold") +
  
  coord_sf(datum = sf::st_crs(32611))

#_______________________________________________________________________
# 5. Write data to files ----
#_______________________________________________________________________

# passes
write.csv(all.passes, paste0(getwd(), "/Derived_data/Passes/passes_complex_high.csv"))

# ctmms (do we need to do this if we're not running speed estimates?)
save(all.ctmms, file = paste0(getwd(), "/Derived_data/CTMMs/ctmms_complex_high.RData"))
