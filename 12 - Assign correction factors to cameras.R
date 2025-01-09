# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 12 - Define correction factors and assign to cameras
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 12 Dec 2024
# Date completed: 12 Dec 2024
# Date last modified: 08 Jan 2025
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation
library(sf)              # read in shapefiles
library(terra)           # rasters
library(cowplot)         # multiple plots

#_______________________________________________________________________
# 2. Read in viewsheds ----
#_______________________________________________________________________

# viewsheds
vs.4 <- st_read(paste0(getwd(), "/Derived_data/Shapefiles/cams_4_vs.shp"))
vs.9 <- st_read(paste0(getwd(), "/Derived_data/Shapefiles/cams_9_vs.shp"))
vs.16 <- st_read(paste0(getwd(), "/Derived_data/Shapefiles/cams_16_vs.shp"))

# add camera ID
vs.4$cam.id <- 1:nrow(vs.4)
vs.9$cam.id <- 1:nrow(vs.9)
vs.16$cam.id <- 1:nrow(vs.16)

#_______________________________________________________________________
# 3. Read in rasters ----
#_______________________________________________________________________
# 3a. Adjusted day range ----
#_______________________________________________________________________

dr.S1L <- rast(paste0(getwd(), "/Rasters/Adjusted day range/S1L.tif"))
dr.S2L <- rast(paste0(getwd(), "/Rasters/Adjusted day range/S2L.tif"))
dr.S3L <- rast(paste0(getwd(), "/Rasters/Adjusted day range/S3L.tif"))
dr.S1H <- rast(paste0(getwd(), "/Rasters/Adjusted day range/S1H.tif"))
dr.S2H <- rast(paste0(getwd(), "/Rasters/Adjusted day range/S2H.tif"))
dr.S3H <- rast(paste0(getwd(), "/Rasters/Adjusted day range/S3H.tif"))

dr.C1L <- rast(paste0(getwd(), "/Rasters/Adjusted day range/C1L.tif"))
dr.C2L <- rast(paste0(getwd(), "/Rasters/Adjusted day range/C2L.tif"))
dr.C3L <- rast(paste0(getwd(), "/Rasters/Adjusted day range/C3L.tif"))
dr.C1H <- rast(paste0(getwd(), "/Rasters/Adjusted day range/C1H.tif"))
dr.C2H <- rast(paste0(getwd(), "/Rasters/Adjusted day range/C2H.tif"))
dr.C3H <- rast(paste0(getwd(), "/Rasters/Adjusted day range/C3H.tif"))

#_______________________________________________________________________
# 3b. Naive SSF predictions ----
#_______________________________________________________________________

ssf.S1L <- rast(paste0(getwd(), "/Rasters/SSF predictions/S1L.tif"))
ssf.S2L <- rast(paste0(getwd(), "/Rasters/SSF predictions/S2L.tif"))
ssf.S3L <- rast(paste0(getwd(), "/Rasters/SSF predictions/S3L.tif"))
ssf.S1H <- rast(paste0(getwd(), "/Rasters/SSF predictions/S1H.tif"))
ssf.S2H <- rast(paste0(getwd(), "/Rasters/SSF predictions/S2H.tif"))
ssf.S3H <- rast(paste0(getwd(), "/Rasters/SSF predictions/S3H.tif"))

ssf.C1L <- rast(paste0(getwd(), "/Rasters/SSF predictions/C1L.tif"))
ssf.C2L <- rast(paste0(getwd(), "/Rasters/SSF predictions/C2L.tif"))
ssf.C3L <- rast(paste0(getwd(), "/Rasters/SSF predictions/C3L.tif"))
ssf.C1H <- rast(paste0(getwd(), "/Rasters/SSF predictions/C1H.tif"))
ssf.C2H <- rast(paste0(getwd(), "/Rasters/SSF predictions/C2H.tif"))
ssf.C3H <- rast(paste0(getwd(), "/Rasters/SSF predictions/C3H.tif"))

#_______________________________________________________________________
# 3c. Simulated iSSF predictions ----
#_______________________________________________________________________

issf.S1L <- rast(paste0(getwd(), "/Rasters/UD rasters/S1L.tif"))
issf.S2L <- rast(paste0(getwd(), "/Rasters/UD rasters/S2L.tif"))
issf.S3L <- rast(paste0(getwd(), "/Rasters/UD rasters/S3L.tif"))
issf.S1H <- rast(paste0(getwd(), "/Rasters/UD rasters/S1H.tif"))
issf.S2H <- rast(paste0(getwd(), "/Rasters/UD rasters/S2H.tif"))
issf.S3H <- rast(paste0(getwd(), "/Rasters/UD rasters/S3H.tif"))

issf.C1L <- rast(paste0(getwd(), "/Rasters/UD rasters/C1L.tif"))
issf.C2L <- rast(paste0(getwd(), "/Rasters/UD rasters/C2L.tif"))
issf.C3L <- rast(paste0(getwd(), "/Rasters/UD rasters/C3L.tif"))
issf.C1H <- rast(paste0(getwd(), "/Rasters/UD rasters/C1H.tif"))
issf.C2H <- rast(paste0(getwd(), "/Rasters/UD rasters/C2H.tif"))
issf.C3H <- rast(paste0(getwd(), "/Rasters/UD rasters/C3H.tif"))

#_______________________________________________________________________
# 4. Extract values for each viewshed ----
#_______________________________________________________________________
# 4a. Define function ----
#_______________________________________________________________________

extract_vs_vals <- function(n.cams,
                            id.landscape,
                            id.variability,
                            id.rep,
                            dr,
                            ssf,
                            issf) {
  
  # use correct viewshed
  if (n.cams == 4) {
    
    vs <- vs.4
    
  } else {
    
    if (n.cams == 9) {
      
      vs <- vs.9
      
    } else {
      
      vs <- vs.16
      
    }
    
  }
  
  # initialize a data.frame
  cam.values <- data.frame(cam.id = vs$cam.id,
                           n.cams = n.cams,
                           landscape = id.landscape,
                           variability = id.variability,
                           rep = id.rep,
                           dr.mean = extract(dr$mean, vs, fun = "mean")[ , 2],
                           dr.se = extract(dr$se, vs, fun = "mean")[ , 2],
                           ssf.mean = extract(ssf$mean, vs, fun = "mean")[ , 2],
                           ssf.se = extract(ssf$se, vs, fun = "mean")[ , 2],
                           issf.mean = extract(issf$mean, vs, fun = "mean")[ , 2],
                           issf.se = extract(issf$se, vs, fun = "mean")[ , 2])
  
  # return
  return(cam.values)
  
}

#_______________________________________________________________________
# 4b. Use function (bind into df) ----
#_______________________________________________________________________

all.cam.values <- rbind(extract_vs_vals(4, "simple", "low", 1, dr.S1L, ssf.S1L, issf.S1L),
                        extract_vs_vals(4, "simple", "low", 2, dr.S2L, ssf.S2L, issf.S2L),
                        extract_vs_vals(4, "simple", "low", 3, dr.S3L, ssf.S3L, issf.S3L),
                        extract_vs_vals(4, "simple", "high", 1, dr.S1H, ssf.S1H, issf.S1H),
                        extract_vs_vals(4, "simple", "high", 2, dr.S2H, ssf.S2H, issf.S2H),
                        extract_vs_vals(4, "simple", "high", 3, dr.S3H, ssf.S3H, issf.S3H),
                        extract_vs_vals(4, "complex", "low", 1, dr.C1L, ssf.C1L, issf.C1L),
                        extract_vs_vals(4, "complex", "low", 2, dr.C2L, ssf.C2L, issf.C2L),
                        extract_vs_vals(4, "complex", "low", 3, dr.C3L, ssf.C3L, issf.C3L),
                        extract_vs_vals(4, "complex", "high", 1, dr.C1H, ssf.C1H, issf.C1H),
                        extract_vs_vals(4, "complex", "high", 2, dr.C2H, ssf.C2H, issf.C2H),
                        extract_vs_vals(4, "complex", "high", 3, dr.C3H, ssf.C3H, issf.C3H),
                        extract_vs_vals(9, "simple", "low", 1, dr.S1L, ssf.S1L, issf.S1L),
                        extract_vs_vals(9, "simple", "low", 2, dr.S2L, ssf.S2L, issf.S2L),
                        extract_vs_vals(9, "simple", "low", 3, dr.S3L, ssf.S3L, issf.S3L),
                        extract_vs_vals(9, "simple", "high", 1, dr.S1H, ssf.S1H, issf.S1H),
                        extract_vs_vals(9, "simple", "high", 2, dr.S2H, ssf.S2H, issf.S2H),
                        extract_vs_vals(9, "simple", "high", 3, dr.S3H, ssf.S3H, issf.S3H),
                        extract_vs_vals(9, "complex", "low", 1, dr.C1L, ssf.C1L, issf.C1L),
                        extract_vs_vals(9, "complex", "low", 2, dr.C2L, ssf.C2L, issf.C2L),
                        extract_vs_vals(9, "complex", "low", 3, dr.C3L, ssf.C3L, issf.C3L),
                        extract_vs_vals(9, "complex", "high", 1, dr.C1H, ssf.C1H, issf.C1H),
                        extract_vs_vals(9, "complex", "high", 2, dr.C2H, ssf.C2H, issf.C2H),
                        extract_vs_vals(9, "complex", "high", 3, dr.C3H, ssf.C3H, issf.C3H),
                        extract_vs_vals(16, "simple", "low", 1, dr.S1L, ssf.S1L, issf.S1L),
                        extract_vs_vals(16, "simple", "low", 2, dr.S2L, ssf.S2L, issf.S2L),
                        extract_vs_vals(16, "simple", "low", 3, dr.S3L, ssf.S3L, issf.S3L),
                        extract_vs_vals(16, "simple", "high", 1, dr.S1H, ssf.S1H, issf.S1H),
                        extract_vs_vals(16, "simple", "high", 2, dr.S2H, ssf.S2H, issf.S2H),
                        extract_vs_vals(16, "simple", "high", 3, dr.S3H, ssf.S3H, issf.S3H),
                        extract_vs_vals(16, "complex", "low", 1, dr.C1L, ssf.C1L, issf.C1L),
                        extract_vs_vals(16, "complex", "low", 2, dr.C2L, ssf.C2L, issf.C2L),
                        extract_vs_vals(16, "complex", "low", 3, dr.C3L, ssf.C3L, issf.C3L),
                        extract_vs_vals(16, "complex", "high", 1, dr.C1H, ssf.C1H, issf.C1H),
                        extract_vs_vals(16, "complex", "high", 2, dr.C2H, ssf.C2H, issf.C2H),
                        extract_vs_vals(16, "complex", "high", 3, dr.C3H, ssf.C3H, issf.C3H))

#_______________________________________________________________________
# 5. Plot of RSS predictions ----
#_______________________________________________________________________

# mean
ggplot(data = all.cam.values,
       aes(x = ssf.mean,
           y = issf.mean,
           color = landscape,
           shape = variability)) +
  
  theme_bw() +
  
  geom_point() +
  
  geom_abline(slope = 1, 
              intercept = 0, 
              linetype = "dashed") +
  
  theme(panel.grid = element_blank()) +
  
  scale_color_manual(values = c("orange", "purple")) +
  
  scale_x_continuous(breaks = c(0.5, 0.75, 1.0, 1.25, 1.5)) +
  scale_y_continuous(breaks = c(0.5, 0.75, 1.0, 1.25, 1.5)) +
  
  coord_cartesian(xlim = c(0.5, 1.55),
                  ylim = c(0.5, 1.55)) +
  
  xlab("RSS from naive SSF prediction") +
  ylab("RSS from iSSF simulated prediction")

#_______________________________________________________________________
# 6. Write to .csv ----
#_______________________________________________________________________

write.csv(all.cam.values, paste0(getwd(), "/Derived_data/Passes/cam_data.csv"))
