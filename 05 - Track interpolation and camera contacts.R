# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 05 - Track interpolation and camera contacts
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 21 Nov 2024
# Date completed: 26 Nov 2024
# Date last modified: 13 Dec 2024
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation
library(amt)             # work with tracks
library(sf)              # spatial data

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

sims.SL.df <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_simple_low.csv"))
sims.CL.df <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_complex_low.csv"))
sims.SH.df <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_simple_high.csv"))
sims.CH.df <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_complex_high.csv"))

# viewsheds
vs <- st_read(paste0(getwd(), "/Derived_data/Shapefiles/cams_9_vs.shp"))

# add camera ID
vs$cam.id <- 1:nrow(vs)

#_______________________________________________________________________
# 3. Create tracks, interpolate, and tabulate camera contacts ----
#_______________________________________________________________________
# 3a. Define function ----
#_______________________________________________________________________

track_contacts <- function (sims.df) {
  
  start.time <- Sys.time()
  
  # initialize df
  all.passes <- data.frame()
  
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
    
    # coerce to sf and cast to lines
    interp.lines <- indiv.track %>%
      
      as_sf() %>%
      
      # transform to UTM
      st_transform(crs = "epsg:32611") %>%
      
      # mutate timestamp
      mutate(t = as.POSIXct(substr(t_, 1, 19),
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
  
  # return
  return(all.passes)
  
}

#_______________________________________________________________________
# 3b. Use function ----
#_______________________________________________________________________

passes.SL <- track_contacts(sims.SL.df)
passes.CL <- track_contacts(sims.CL.df)
passes.SH <- track_contacts(sims.SH.df)
passes.CH <- track_contacts(sims.CH.df)

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
write.csv(passes.SL, paste0(getwd(), "/Derived_data/Passes/passes_simple_low.csv"))
write.csv(passes.CL, paste0(getwd(), "/Derived_data/Passes/passes_complex_low.csv"))
write.csv(passes.SH, paste0(getwd(), "/Derived_data/Passes/passes_simple_high.csv"))
write.csv(passes.CH, paste0(getwd(), "/Derived_data/Passes/passes_complex_high.csv"))
