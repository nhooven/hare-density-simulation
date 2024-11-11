# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 01 - Subset GPS data
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 11 Nov 2024
# Date completed: 11 Nov 2024
# Date last modified: 11 Nov 2024
# R version: 4.2.2

# NEXT - Subset "during" data further to treatment day

#_______________________________________________________________________________________________
# 1. Load required packages and directories ----
#_______________________________________________________________________________________________

library(tidyverse)
library(lubridate)

# processed GPS data directory
gps.dir <- "D:/Hare project/Data analysis/Spatial ecology/A - Process data/Derived_data/"

#_______________________________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________________________

# pre
data.pre <- read.csv(paste0(gps.dir, "all_data_pre_1.csv"))

# during
data.dur <- read.csv(paste0(gps.dir, "all_data_dur_1.csv"))

# post (for adding later if we retrieve fall 2024 GPS collars)
# data.post <- read.csv(paste0(gps.dir, "all_data_post_1.csv"))

#_______________________________________________________________________________________________
# 3. Keep relevant columns and ensure they're in the correct format ----
#_______________________________________________________________________________________________

# pre
data.pre.1 <- data.pre %>%
  
  dplyr::select(timestamp,
                longitude,
                latitude,
                easting,
                northing,
                Group,
                Site,
                Ear.tag,
                TrackID,
                Sex,
                days.cap) %>%
  
  # format date correctly
  mutate(timestamp = ymd_hms(timestamp, tz = "America/Los_Angeles"),
         Ear.tag = as.character(Ear.tag))

# dur
data.dur.1 <- data.dur %>%
  
  dplyr::select(timestamp,
                longitude,
                latitude,
                easting,
                northing,
                Group,
                Site,
                Ear.tag,
                TrackID,
                Sex,
                days.cap) %>%
  
  # format date and Ear.tag correctly
  mutate(timestamp = ymd_hms(timestamp, tz = "America/Los_Angeles"),
         Ear.tag = as.character(Ear.tag))

#_______________________________________________________________________________________________
# 4. Subset individual tracks (necessary for pre and post only) ----
#_______________________________________________________________________________________________

# pre 
data.pre.2 <- data.pre.1 %>%
  
  filter(TrackID %in% c("1876_1",
                        "1857_1",
                        "84_1",
                        "1877_1",
                        "1855_1",
                        "1885_1",
                        "1854_1",
                        "1880_1",
                        "1841_1",
                        "1842_1",
                        "1866_1",
                        "724_1",
                        "1880_2",
                        "1854_2",
                        "721_1",
                        "979_1",
                        "1877_2",
                        "1581_2",
                        "1580_2",
                        "396_2",
                        "446_1",
                        "932_1",
                        "828_1",
                        "803_1",
                        "1051_1",
                        "835_1",
                        "589_1",
                        "445_1",
                        "832_1",
                        "1067_1",
                        "539_1",
                        "986_1"))

#_______________________________________________________________________________________________
# 5. Bind together and filter out snow-on dates ----
#_______________________________________________________________________________________________

data.all <- bind_rows(data.pre.2, data.dur.1)

# define snow season cutoffs
start.2022 <- as.Date(ymd("2022-11-06", tz = "America/Los_Angeles"))
end.2023 <- as.Date(ymd("2023-04-30", tz = "America/Los_Angeles"))

start.2023 <- as.Date(ymd("2023-10-25", tz = "America/Los_Angeles"))
end.2024 <- as.Date(ymd("2024-04-24", tz = "America/Los_Angeles"))

#start.2024

data.all.1 <- data.all %>%
  
  filter(# before snow 2022
         (timestamp <= start.2022) |
         
         # after melt 2023 but before snow 2023
         (timestamp >= end.2023 & timestamp <= start.2023) | 
           
         # or after melt 2024
         (timestamp >= end.2024))

# quick visualization to see if this did what I hoped it would
ggplot(data.all.1,
       aes(x = month(timestamp),
           y = as.factor(TrackID))) +
  
  geom_point()

# looks good!

#_______________________________________________________________________________________________
# 6. Write to .csv ----
#_______________________________________________________________________________________________

write.csv(data.all.1, "Derived_data/all_data.csv")
