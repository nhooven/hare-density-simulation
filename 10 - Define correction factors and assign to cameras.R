# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 10 - Define correction factors and assign to cameras
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 12 Dec 2024
# Date completed: 12 Dec 2024
# Date last modified: 13 Dec 2024
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation
library(sf)              # read in shapefiles
library(terra)           # rasters
library(cowplot)         # multiple plots

#_______________________________________________________________________
# 2. Read in spatial data ----
#_______________________________________________________________________

# viewsheds
vs <- st_read(paste0(getwd(), "/Derived_data/Shapefiles/cams_9_vs.shp"))

# add camera ID
vs$cam.id <- 1:nrow(vs)

# unit boundary
unit.bound <- st_read(paste0(getwd(), "/Derived_data/Shapefiles/unit_bound.shp"))

# SSF/iSSF predictions
ssf.rast <- exp(rast(paste0(getwd(), "/Rasters/SSF_pred.tif")))
issf.rast <- rast(paste0(getwd(), "/Rasters/kernel_all_sim.tif"))

#_______________________________________________________________________
# 3. Plot raw values ----
#_______________________________________________________________________

plot(ssf.rast)
plot(issf.rast)

#_______________________________________________________________________
# 4. Crop to unit boundary ----
#_______________________________________________________________________

ssf.rast.crop <- crop(ssf.rast, unit.bound)
issf.rast.crop <- crop(issf.rast, unit.bound)

#_______________________________________________________________________
# 5. Scale correctly ----
#_______________________________________________________________________
# 5a. SSF raster - relative selection strength

# RSS is interpreted as how many times more likely an animal is to be
# found at a location with the characteristics of the numerator compared
# to the denominator:

# RSS = w(x)2 / w(x)1

# we'll use the average w(x) across the camera grid as the baseline,
# and normalize everything to that - thus the average will be 1

# and then the correction factor will be 1 / RSS

#_______________________________________________________________________

ssf.rast.scale <- c(ssf.rast.crop$SL / mean(values(ssf.rast.crop$SL)),
                    ssf.rast.crop$CL / mean(values(ssf.rast.crop$CL)),
                    ssf.rast.crop$SH / mean(values(ssf.rast.crop$SH)),
                    ssf.rast.crop$CH / mean(values(ssf.rast.crop$CH)))



issf.rast.scale <- c(issf.rast.crop$SL / mean(values(issf.rast.crop$SL)),
                     issf.rast.crop$CL / mean(values(issf.rast.crop$CL)),
                     issf.rast.crop$SH / mean(values(issf.rast.crop$SH)),
                     issf.rast.crop$CH / mean(values(issf.rast.crop$CH)))

#_______________________________________________________________________
# 5b. iSSF utilization distribution ----
#_______________________________________________________________________

# plot again
plot(ssf.rast.scale)
plot(issf.rast.scale)

#_______________________________________________________________________
# 5. Extract mean values for each viewshed ----
#_______________________________________________________________________

# extract using the mean of each polygon
ssf.extract <- extract(ssf.rast.scale, vs, fun = "mean")[2:5]
issf.extract <- extract(issf.rast.scale, vs, fun = "mean")[2:5]

# rename
names(ssf.extract) <- c("ssf.SL", "ssf.CL", "ssf.SH", "ssf.CH")
names(issf.extract) <- c("issf.SL", "issf.CL", "issf.SH", "issf.CH")

# bind into vs
vs.df <- as.data.frame(vs$cam.id)

vs.df <- cbind(vs.df, ssf.extract, issf.extract)

#_______________________________________________________________________
# 6. Final faceted plot ----
#_______________________________________________________________________

ggplot() +
  
  theme_bw() +
  
  facet_wrap(~ lyr, nrow = 1) +
  
  tidyterra::geom_spatraster(data = ssf.rast.scale) +
  
  # add unit boundary
  geom_sf(data = unit.bound,
          color = "white",
          fill = NA) +
  
  # add camera viewsheds at a scale that makes them visible
  geom_sf(data = vs,
          fill = "white",
          color = "white",
          linewidth = 1.1) +
  
  theme(axis.text = element_blank(),
        legend.position = "none") + 
  
  scale_fill_viridis_c(option = "magma") -> ssf.final.plot

ggplot() +
  
  theme_bw() +
  
  facet_wrap(~ lyr, nrow = 1) +
  
  tidyterra::geom_spatraster(data = issf.rast.scale) +
  
  # add unit boundary
  geom_sf(data = unit.bound,
          color = "white",
          fill = NA) +
  
  # add camera viewsheds at a scale that makes them visible
  geom_sf(data = vs,
          fill = "white",
          color = "white",
          linewidth = 1.1) +
  
  theme(axis.text = element_blank(),
        legend.position = "none") + 
  
  scale_fill_viridis_c(option = "magma") -> issf.final.plot

cowplot::plot_grid(ssf.final.plot, issf.final.plot, nrow = 2)

#_______________________________________________________________________
# 7. Write to .csv ----
#_______________________________________________________________________

write.csv(vs.df, paste0(getwd(), "/Derived_data/Detections/cams.csv"))
