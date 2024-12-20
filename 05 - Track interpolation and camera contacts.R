# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 05 - Track interpolation and camera contacts
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 21 Nov 2024
# Date completed: 26 Nov 2024
# Date last modified: 20 Dec 2024
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

sims.S1L.df <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_S1L.csv"))
sims.S2L.df <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_S2L.csv"))
sims.S3L.df <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_S3L.csv"))
sims.S1H.df <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_S1H.csv"))
sims.S2H.df <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_S2H.csv"))
sims.S3H.df <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_S3H.csv"))

sims.C1L.df <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_C1L.csv"))
sims.C2L.df <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_C2L.csv"))
sims.C3L.df <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_C3L.csv"))
sims.C1H.df <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_C1H.csv"))
sims.C2H.df <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_C2H.csv"))
sims.C3H.df <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_C3H.csv"))

# bind together
sims.df <- rbind(sims.S1L.df, sims.S2L.df, sims.S3L.df,
                 sims.S1H.df, sims.S2H.df, sims.S3H.df,
                 sims.C1L.df, sims.C2L.df, sims.C3L.df,
                 sims.C1H.df, sims.C2H.df, sims.C3H.df)

# viewsheds
vs.4 <- st_read(paste0(getwd(), "/Derived_data/Shapefiles/cams_4_vs.shp"))
vs.9 <- st_read(paste0(getwd(), "/Derived_data/Shapefiles/cams_9_vs.shp"))
vs.16 <- st_read(paste0(getwd(), "/Derived_data/Shapefiles/cams_16_vs.shp"))

# add camera ID
vs.4$cam.id <- 1:nrow(vs.4)
vs.9$cam.id <- 1:nrow(vs.9)
vs.16$cam.id <- 1:nrow(vs.16)

#_______________________________________________________________________
# 3. Create tracks, interpolate, and tabulate camera contacts ----
#_______________________________________________________________________
# 3a. Define function ----
#_______________________________________________________________________

track_contacts <- function (id.landscape,
                            id.variability,
                            id.rep,
                            n.cams) {
  
  # subset sims.df
  sims.df.1 <- sims.df %>%
    
    filter(landscape == id.landscape &
           variability == id.variability &
           rep == id.rep)
  
  # use correct vs
  if (n.cams == 4) {
    
    vs <- vs.4
    
  } else {
    
    if (n.cams == 9) {
      
      vs <- vs.9
      
    } else {
      
      vs <- vs.16
      
    }
    
  }
  
  # initialize df
  all.passes <- data.frame()
  
  for (i in unique(sims.df.1$indiv)) {
    
    # subset individual's data and create a track
    indiv.track <- sims.df.1 %>% 
      
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
    
  }
  
  # bind in ids
  all.passes <- all.passes %>%
    
    mutate(landscape = id.landscape,
           variability = id.variability,
           rep = id.rep,
           cams = n.cams)
  
  # return
  return(all.passes)
  
}

#_______________________________________________________________________
# 3b. Use function ----
#_______________________________________________________________________

# 4
passes.4 <- rbind(track_contacts("simple", "low", 1, 4),
                  track_contacts("simple", "low", 2, 4),
                  track_contacts("simple", "low", 3, 4),
                  track_contacts("simple", "high", 1, 4),
                  track_contacts("simple", "high", 2, 4),
                  track_contacts("simple", "high", 3, 4),
                  track_contacts("complex", "low", 1, 4),
                  track_contacts("complex", "low", 2, 4),
                  track_contacts("complex", "low", 3, 4),
                  track_contacts("complex", "high", 1, 4),
                  track_contacts("complex", "high", 2, 4),
                  track_contacts("complex", "high", 3, 4))

# 9
passes.9 <- rbind(track_contacts("simple", "low", 1, 9),
                  track_contacts("simple", "low", 2, 9),
                  track_contacts("simple", "low", 3, 9),
                  track_contacts("simple", "high", 1, 9),
                  track_contacts("simple", "high", 2, 9),
                  track_contacts("simple", "high", 3, 9),
                  track_contacts("complex", "low", 1, 9),
                  track_contacts("complex", "low", 2, 9),
                  track_contacts("complex", "low", 3, 9),
                  track_contacts("complex", "high", 1, 9),
                  track_contacts("complex", "high", 2, 9),
                  track_contacts("complex", "high", 3, 9))

# 16
passes.16 <- rbind(track_contacts("simple", "low", 1, 16),
                   track_contacts("simple", "low", 2, 16),
                   track_contacts("simple", "low", 3, 16),
                   track_contacts("simple", "high", 1, 16),
                   track_contacts("simple", "high", 2, 16),
                   track_contacts("simple", "high", 3, 16),
                   track_contacts("complex", "low", 1, 16),
                   track_contacts("complex", "low", 2, 16),
                   track_contacts("complex", "low", 3, 16),
                   track_contacts("complex", "high", 1, 16),
                   track_contacts("complex", "high", 2, 16),
                   track_contacts("complex", "high", 3, 16))

#_______________________________________________________________________
# 4. Example plot ----
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
write.csv(passes.4, paste0(getwd(), "/Derived_data/Passes/passes_4.csv"))
write.csv(passes.9, paste0(getwd(), "/Derived_data/Passes/passes_9.csv"))
write.csv(passes.16, paste0(getwd(), "/Derived_data/Passes/passes_16.csv"))
