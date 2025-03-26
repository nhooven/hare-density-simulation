# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 05 - Prepare detection data
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 26 Mar 2025
# Date completed: 
# Date last modified: 26 Mar 2025
# R version: 4.4.3

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

detections.lo <- read.csv(paste0(getwd(), "/Derived data/Camera detections/detections_lo.csv"))

#_______________________________________________________________________
# 3. Aggregate detections by camera ----
#_______________________________________________________________________
# 3a. Q1 (low variability) ----
#_______________________________________________________________________

q1.agg <- detections.lo %>% 
  
  group_by(abund.group,
           rep,
           cam.id) %>%
  
  summarize(detections = n())

#_______________________________________________________________________
# 4. Add "zero" count rows where necessary ----
#_______________________________________________________________________
# 4a. Q1 (low variability) ----
#_______________________________________________________________________

# df of possibilities
q1.grid <- expand.grid(abund.group = unique(q1.agg$abund.group),
                       rep = 1:max(q1.agg$rep),
                       cam.id = 1:9)

# loop through aggregated df
q1.zero.rows <- data.frame()

for (i in 1:nrow(q1.grid)) {
  
  focal.row <- q1.grid %>% slice(i)
  
  focal.row.match <- plyr::match_df(q1.agg[ , -4], focal.row)
  
  # does this row exist?
  if (nrow(focal.row.match) == 0) {
    
    # then bind in a new row
    q1.zero.rows <- rbind(q1.zero.rows,
                          focal.row %>% mutate(detections = 0))
    
  }
  
}

# bind in and sort
q1.agg.1 <- q1.agg %>%
  
  bind_rows(q1.zero.rows) %>%
  
  arrange(abund.group, 
          rep,
          cam.id)
