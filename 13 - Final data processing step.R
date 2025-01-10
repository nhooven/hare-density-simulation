# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 13 - Final data processing step
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 08 Jan 2025
# Date completed: 
# Date last modified: 09 Jan 2025
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation
library(terra)           # work with rasters

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

# camera passes
passes <- read.csv(paste0(getwd(), "/Derived_data/Passes/passes_gp_all.csv"))

# camera day range and corrections
cam.data <- read.csv(paste0(getwd(), "/Derived_data/Passes/cam_data.csv"))

#_______________________________________________________________________
# 3. Merge dataset ----
#_______________________________________________________________________

# drop unneeded columns and ensure corresponding variables are correctly named
passes <- passes %>% dplyr::select(-X) %>% rename(n.cams = cams)
cam.data <- cam.data %>% dplyr::select(-X)

# left join
passes.1 <- left_join(passes, cam.data)

#_______________________________________________________________________
# 4. Add mean day range by landscape/variability/rep ----
#_______________________________________________________________________
# 4a. Read in rasters ----
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
# 4b. Calculate means and SD ----
#_______________________________________________________________________

constant.dr <- data.frame(dr.c.mean = c(mean(values(dr.S1L$mean)),
                                         mean(values(dr.S2L$mean)),
                                         mean(values(dr.S3L$mean)),
                                         mean(values(dr.S1H$mean)),
                                         mean(values(dr.S2H$mean)),
                                         mean(values(dr.S3H$mean)),
                                         mean(values(dr.C1L$mean)),
                                         mean(values(dr.C2L$mean)),
                                         mean(values(dr.C3L$mean)),
                                         mean(values(dr.C1H$mean)),
                                         mean(values(dr.C2H$mean)),
                                         mean(values(dr.C3H$mean))),
                          dr.c.sd = c(sd(values(dr.S1L$mean)),
                                      sd(values(dr.S2L$mean)),
                                      sd(values(dr.S3L$mean)),
                                      sd(values(dr.S1H$mean)),
                                      sd(values(dr.S2H$mean)),
                                      sd(values(dr.S3H$mean)),
                                      sd(values(dr.C1L$mean)),
                                      sd(values(dr.C2L$mean)),
                                      sd(values(dr.C3L$mean)),
                                      sd(values(dr.C1H$mean)),
                                      sd(values(dr.C2H$mean)),
                                      sd(values(dr.C3H$mean))),
                          landscape = c("simple", "simple", "simple",
                                        "simple", "simple", "simple",
                                        "complex", "complex", "complex",
                                        "complex", "complex", "complex"),
                          variability = c("low", "low", "low",
                                          "high", "high", "high",
                                          "low", "low", "low",
                                          "high", "high", "high"),
                          rep = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3))

#_______________________________________________________________________
# 4c. Merge ----
#_______________________________________________________________________

passes.2 <- left_join(passes.1, constant.dr)

#_______________________________________________________________________
# 5. Add constant REM inputs (assume no variance) ----
#_______________________________________________________________________

passes.3 <- passes.2 %>%
  
  mutate(lens = (57.3 * pi) / 180,
         days = 28)

#_______________________________________________________________________
# 6. Write to .csv ----
#_______________________________________________________________________

write.csv(passes.3, paste0(getwd(), "/Derived_data/Passes/final_passes.csv"))
