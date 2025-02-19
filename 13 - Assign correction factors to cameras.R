# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 12 - Define correction factors and assign to cameras
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 12 Dec 2024
# Date completed: 12 Dec 2024
# Date last modified: 19 Feb 2025
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

dr.B1 <- rast(paste0(getwd(), "/Rasters/Adjusted day range/B1.tif"))
dr.B2 <- rast(paste0(getwd(), "/Rasters/Adjusted day range/B2.tif"))
dr.B3 <- rast(paste0(getwd(), "/Rasters/Adjusted day range/B3.tif"))
dr.A1 <- rast(paste0(getwd(), "/Rasters/Adjusted day range/A1.tif"))
dr.A2 <- rast(paste0(getwd(), "/Rasters/Adjusted day range/A2.tif"))
dr.A3 <- rast(paste0(getwd(), "/Rasters/Adjusted day range/A3.tif"))

#_______________________________________________________________________
# 3b. RSF predictions ----
#_______________________________________________________________________

rsf.B1 <- rast(paste0(getwd(), "/Rasters/RSF predictions/B1.tif"))
rsf.B2 <- rast(paste0(getwd(), "/Rasters/RSF predictions/B2.tif"))
rsf.B3 <- rast(paste0(getwd(), "/Rasters/RSF predictions/B3.tif"))
rsf.A1 <- rast(paste0(getwd(), "/Rasters/RSF predictions/A1.tif"))
rsf.A2 <- rast(paste0(getwd(), "/Rasters/RSF predictions/A2.tif"))
rsf.A3 <- rast(paste0(getwd(), "/Rasters/RSF predictions/A3.tif"))

#_______________________________________________________________________
# 3c. Simulated iSSF predictions ----
#_______________________________________________________________________

issf.B1 <- rast(paste0(getwd(), "/Rasters/UD predictions/B1.tif"))
issf.B2 <- rast(paste0(getwd(), "/Rasters/UD predictions/B2.tif"))
issf.B3 <- rast(paste0(getwd(), "/Rasters/UD predictions/B3.tif"))
issf.A1 <- rast(paste0(getwd(), "/Rasters/UD predictions/A1.tif"))
issf.A2 <- rast(paste0(getwd(), "/Rasters/UD predictions/A2.tif"))
issf.A3 <- rast(paste0(getwd(), "/Rasters/UD predictions/A3.tif"))

#_______________________________________________________________________
# 4. Extract values for each viewshed ----
#_______________________________________________________________________
# 4a. Define function ----
#_______________________________________________________________________

extract_vs_vals <- function(buffer = 0,
                            n.cams,
                            id.trt,
                            id.rep,
                            dr,
                            rsf,
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
  
  # determine if a buffer is used
  if (buffer == 0) {
    
    # initialize a data.frame (make sure the names are correct!)
    cam.values <- data.frame(cam.id = vs$cam.id,
                             n.cams = n.cams,
                             trt = id.trt,
                             rep = id.rep,
                             dr.mean = extract(dr$mean, vs, fun = "mean")[ , 2],
                             dr.se = extract(dr$se, vs, fun = "mean")[ , 2],
                             rsf.mean = extract(rsf$cf, vs, fun = "mean")[ , 2],       
                             rsf.se = extract(rsf$se, vs, fun = "mean")[ , 2],
                             issf.mean = extract(issf$mean, vs, fun = "mean")[ , 2],
                             issf.se = extract(issf$std, vs, fun = "mean")[ , 2])
    
    # else case - use a buffer to do the extraction
  } else {
    
    # buffer the viewsheds
    vs.buffer <- st_buffer(vs, dist = buffer)
    
    # initialize a data.frame (make sure the names are correct!)
    cam.values <- data.frame(cam.id = vs$cam.id,
                             n.cams = n.cams,
                             trt = id.trt,
                             rep = id.rep,
                             dr.mean = extract(dr$mean, vs.buffer, fun = "mean")[ , 2],
                             dr.se = extract(dr$se, vs.buffer, fun = "mean")[ , 2],
                             rsf.mean = extract(rsf$cf, vs.buffer, fun = "mean")[ , 2],       
                             rsf.se = extract(rsf$se, vs.buffer, fun = "mean")[ , 2],
                             issf.mean = extract(issf$mean, vs.buffer, fun = "mean")[ , 2],
                             issf.se = extract(issf$std, vs.buffer, fun = "mean")[ , 2])
    
  }
  
  # return
  return(cam.values)
  
}

#_______________________________________________________________________
# 4b. Use function (bind into df) ----
#_______________________________________________________________________

all.cam.values <- rbind(extract_vs_vals(buffer = 0, 4, "before", 1, dr.B1, rsf.B1, issf.B1),
                        extract_vs_vals(buffer = 0, 4, "before", 2, dr.B2, rsf.B2, issf.B2),
                        extract_vs_vals(buffer = 0, 4, "before", 3, dr.B3, rsf.B3, issf.B3),
                        extract_vs_vals(buffer = 0, 4, "after", 1, dr.A1, rsf.A1, issf.A1),
                        extract_vs_vals(buffer = 0, 4, "after", 2, dr.A2, rsf.A2, issf.A2),
                        extract_vs_vals(buffer = 0, 4, "after", 3, dr.A3, rsf.A3, issf.A3),
                        extract_vs_vals(buffer = 0, 9, "before", 1, dr.B1, rsf.B1, issf.B1),
                        extract_vs_vals(buffer = 0, 9, "before", 2, dr.B2, rsf.B2, issf.B2),
                        extract_vs_vals(buffer = 0, 9, "before", 3, dr.B3, rsf.B3, issf.B3),
                        extract_vs_vals(buffer = 0, 9, "after", 1, dr.A1, rsf.A1, issf.A1),
                        extract_vs_vals(buffer = 0, 9, "after", 2, dr.A2, rsf.A2, issf.A2),
                        extract_vs_vals(buffer = 0, 9, "after", 3, dr.A3, rsf.A3, issf.A3),
                        extract_vs_vals(buffer = 0, 16, "before", 1, dr.B1, rsf.B1, issf.B1),
                        extract_vs_vals(buffer = 0, 16, "before", 2, dr.B2, rsf.B2, issf.B2),
                        extract_vs_vals(buffer = 0, 16, "before", 3, dr.B3, rsf.B3, issf.B3),
                        extract_vs_vals(buffer = 0, 16, "after", 1, dr.A1, rsf.A1, issf.A1),
                        extract_vs_vals(buffer = 0, 16, "after", 2, dr.A2, rsf.A2, issf.A2),
                        extract_vs_vals(buffer = 0, 16, "after", 3, dr.A3, rsf.A3, issf.A3))

all.cam.values.buff <- rbind(extract_vs_vals(buffer = 10, 4, "before", 1, dr.B1, rsf.B1, issf.B1),
                             extract_vs_vals(buffer = 10, 4, "before", 2, dr.B2, rsf.B2, issf.B2),
                             extract_vs_vals(buffer = 10, 4, "before", 3, dr.B3, rsf.B3, issf.B3),
                             extract_vs_vals(buffer = 10, 4, "after", 1, dr.A1, rsf.A1, issf.A1),
                             extract_vs_vals(buffer = 10, 4, "after", 2, dr.A2, rsf.A2, issf.A2),
                             extract_vs_vals(buffer = 10, 4, "after", 3, dr.A3, rsf.A3, issf.A3),
                             extract_vs_vals(buffer = 10, 9, "before", 1, dr.B1, rsf.B1, issf.B1),
                             extract_vs_vals(buffer = 10, 9, "before", 2, dr.B2, rsf.B2, issf.B2),
                             extract_vs_vals(buffer = 10, 9, "before", 3, dr.B3, rsf.B3, issf.B3),
                             extract_vs_vals(buffer = 10, 9, "after", 1, dr.A1, rsf.A1, issf.A1),
                             extract_vs_vals(buffer = 10, 9, "after", 2, dr.A2, rsf.A2, issf.A2),
                             extract_vs_vals(buffer = 10, 9, "after", 3, dr.A3, rsf.A3, issf.A3),
                             extract_vs_vals(buffer = 10, 16, "before", 1, dr.B1, rsf.B1, issf.B1),
                             extract_vs_vals(buffer = 10, 16, "before", 2, dr.B2, rsf.B2, issf.B2),
                             extract_vs_vals(buffer = 10, 16, "before", 3, dr.B3, rsf.B3, issf.B3),
                             extract_vs_vals(buffer = 10, 16, "after", 1, dr.A1, rsf.A1, issf.A1),
                             extract_vs_vals(buffer = 10, 16, "after", 2, dr.A2, rsf.A2, issf.A2),
                             extract_vs_vals(buffer = 10, 16, "after", 3, dr.A3, rsf.A3, issf.A3))

#_______________________________________________________________________
# 5. Plot of CF predictions ----
#_______________________________________________________________________

# mean
ggplot(data = all.cam.values.buff,
       aes(x = rsf.mean,
           y = issf.mean,
           color = trt,
           shape = as.factor(rep))) +
  
  theme_bw() +
  
  geom_point() +
  
  geom_abline(slope = 1, 
              intercept = 0, 
              linetype = "dashed") +
  
  theme(panel.grid = element_blank()) +
  
  scale_color_manual(values = c("orange", "purple")) +
  
  scale_x_continuous(breaks = c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0)) +
  scale_y_continuous(breaks = c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0)) +
  
  coord_cartesian(xlim = c(0.5, 3.0),
                  ylim = c(0.5, 3.0)) +
  
  xlab("CF from RSF prediction") +
  ylab("CF from iSSF simulated prediction")

#_______________________________________________________________________
# 6. Write to .csv ----
#_______________________________________________________________________

write.csv(all.cam.values, paste0(getwd(), "/Derived_data/For REM/cam_data.csv"))
write.csv(all.cam.values.buff, paste0(getwd(), "/Derived_data/For REM/cam_data_buff.csv"))
