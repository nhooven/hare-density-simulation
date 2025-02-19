# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 13 - Final data processing step
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 08 Jan 2025
# Date completed: 
# Date last modified: 19 Feb 2025
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
cam.data <- read.csv(paste0(getwd(), "/Derived_data/For REM/cam_data.csv"))
cam.data.buff <- read.csv(paste0(getwd(), "/Derived_data/For REM/cam_data_buff.csv"))

#_______________________________________________________________________
# 3. Merge datasets ----
#_______________________________________________________________________

# drop unneeded columns and ensure corresponding variables are correctly named
passes <- passes %>% dplyr::select(-X) %>% rename(n.cams = cams)
cam.data <- cam.data %>% dplyr::select(-X)
cam.data.buff <- cam.data.buff %>% dplyr::select(-X)

# left join
passes.1 <- left_join(passes, cam.data)
passes.1.buff <- left_join(passes, cam.data.buff)

#_______________________________________________________________________
# 4. Add mean day range ----
#_______________________________________________________________________
# 4a. Read in dataset ----
#_______________________________________________________________________

static.dr <- read.csv(paste0(getwd(), "/Derived_data/Model parameters/static_dr.csv"))

# rename columns
static.dr <- static.dr %>%
  
  rename(static.dr.mean = dr.mean,
         static.dr.se = dr.sd)

#_______________________________________________________________________
# 4b. Join ----
#_______________________________________________________________________

passes.2 <- left_join(passes.1, static.dr[ , -1])
passes.2.buff <- left_join(passes.1.buff, static.dr[ , -1])

#_______________________________________________________________________
# 5. Add constant REM inputs (assume no variance) ----
#_______________________________________________________________________

passes.3 <- passes.2 %>%
  
  mutate(lens = (57.3 * pi) / 180,
         days = 28)

passes.3.buff <- passes.2.buff %>%
  
  mutate(lens = (57.3 * pi) / 180,
         days = 28)

#_______________________________________________________________________
# 6. Write to .csv ----
#_______________________________________________________________________

write.csv(passes.3, paste0(getwd(), "/Derived_data/For REM/final_passes.csv"))
write.csv(passes.3.buff, paste0(getwd(), "/Derived_data/For REM/final_passes_buff.csv"))
