# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 02 - Process steps 
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 11 Nov 2024
# Date completed: 11 Nov 2024
# Date last modified: 11 Nov 2024
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages and define CRS ----
#_______________________________________________________________________________________________

library(tidyverse)
library(lubridate)
library(amt)
library(terra)

utm.epsg <- 32611

#_______________________________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________________________

data.all <- read.csv("Derived_data/all_data.csv")

# raster directory
rast.dir <- "D:/Hare project/Spatial data/Rasters/Derived_data/"

# read in rasters
# vegetation
cc <- rast(paste0(rast.dir, "cc_utm.tif"))
ch <- rast(paste0(rast.dir, "eth_ch_utm.tif"))

# topography
slope <- rast(paste0(rast.dir, "dtm_slope_deg.tif"))
tpi <- rast(paste0(rast.dir, "tpi_10.tif"))
twi <- rast(paste0(rast.dir, "twi_10.tif"))

#_______________________________________________________________________________________________
# 3. Ensure timestamp is in correct format ----
#_______________________________________________________________________________________________

data.all <- data.all %>%
  
  mutate(timestamp = ymd_hms(timestamp, tz = "America/Los_Angeles"))

#_______________________________________________________________________________________________
# 4. Make amt "tracks" by TrackID ----
#_______________________________________________________________________________________________

all.tracks <- tibble()

for (i in unique(data.all$TrackID)) {
  
  focal.id <- i
  
  # subset
  focal.data <- data.all %>% filter(TrackID == i)
  
  # make a track
  focal.track <- focal.data %>% 
    
    make_track(.x = easting,
               .y = northing,
               .t = timestamp,
               all_cols = TRUE,
               check_duplicates = TRUE,
               crs = utm.epsg) %>%
    
    # resample track
    track_resample(rate = hours(2),
                   tolerance = minutes(15)) %>%
    
    # steps by burst
    steps_by_burst(keep_cols = "start")
  
  # bind together
  all.tracks <- rbind(all.tracks, focal.track)
  
}

#_______________________________________________________________________________________________
# 5. Examine changes in movement post-collaring ----
#_______________________________________________________________________________________________

ggplot(data = all.tracks,
       aes(x = days.cap,
           y = log(sl_))) +
  
  theme_bw() +
  
  geom_point(alpha = 0.5) +
  
  geom_smooth(se = FALSE)

# we'll keep all of them for now and likely change this later

#_______________________________________________________________________________________________
# 6. Tentative SL and TA distributions ----
#_______________________________________________________________________________________________
# 6a. Fit  distributions ----
#_______________________________________________________________________________________________

# sl
# gamma
dist.gamma <- fit_distr(all.tracks$sl_, 
                        dist_name = "gamma")

# exponential
dist.exp <- fit_distr(all.tracks$sl_, 
                      dist_name = "exp")

# ta
# uniform
dist.unif <- fit_distr(all.tracks$ta_, 
                        dist_name = "unif")

# Von Mises
dist.vm <- fit_distr(all.tracks$ta_, 
                      dist_name = "vonmises")

#_______________________________________________________________________________________________
# 6b. Calculate implied distributions for plotting ----
#_______________________________________________________________________________________________

# sl
sl.dists.df <- data.frame(x = seq(min(all.tracks$sl_),
                                  max(all.tracks$sl_),
                                  length.out = 1000),
                          gamma = dgamma(x = seq(min(all.tracks$sl_),
                                                 max(all.tracks$sl_),
                                                 length.out = 1000),
                                         shape = dist.gamma$params$shape,
                                         scale = dist.gamma$params$scale),
                          exp = dexp(x = seq(min(all.tracks$sl_),
                                             max(all.tracks$sl_),
                                             length.out = 1000),
                                       rate = dist.exp$params$rate))

# ta
ta.dists.df <- data.frame(x = seq(dist.unif$params$min,
                                  dist.unif$params$max,
                                  length.out = 1000),
                          unif = dunif(x = seq(dist.unif$params$min,
                                               dist.unif$params$max,
                                               length.out = 1000),
                                       min = dist.unif$params$min,
                                       max = dist.unif$params$max),
                          vm = circular::dvonmises(x = seq(dist.unif$params$min,
                                                           dist.unif$params$max,
                                                           length.out = 1000),
                                                   mu = dist.vm$params$mu,
                                                   kappa = dist.vm$params$kappa))

#_______________________________________________________________________________________________
# 6c. Plot ----
#_______________________________________________________________________________________________

# sl
ggplot() +
  
  theme_bw() +
  
  # density
  geom_density(data = all.tracks,
                 aes(x = sl_),
                 fill = "lightgray",
                 color = "black") +
  
  # gamma
  geom_line(data = sl.dists.df,
            aes(x = x,
                y = gamma),
            linewidth = 1.05,
            color = "blue") +
  
  # exponential
  geom_line(data = sl.dists.df,
            aes(x = x,
                y = exp),
            linewidth = 1.05,
            color = "red")

# ta
ggplot() +
  
  theme_bw() +
  
  # density
  geom_density(data = all.tracks,
               aes(x = ta_),
               fill = "lightgray",
               color = "black") +
  
  # uniform
  geom_line(data = ta.dists.df,
            aes(x = x,
                y = unif),
            linewidth = 1.05,
            color = "blue") +
  
  # Von Mises
  geom_line(data = ta.dists.df,
            aes(x = x,
                y = vm),
            linewidth = 1.05,
            color = "red")

# from here, it looks like the gamma and the uniform do a decent job of approximating the dists

#_______________________________________________________________________________________________
# 7. Sample random steps ----
#_______________________________________________________________________________________________

all.tracks.1 <- all.tracks %>%
  
  random_steps(n_control = 25,
               sl_distr = dist.gamma,
               ta_distr = dist.unif)

#_______________________________________________________________________________________________
# 8. Extract covariates at the end of every step ----
#_______________________________________________________________________________________________

all.tracks.2 <- all.tracks.1 %>% 
  
  extract_covariates(cc) %>%
  extract_covariates(ch) %>%
  extract_covariates(slope) %>%
  extract_covariates(tpi) %>%
  extract_covariates(twi) 

# rename covariates
all.tracks.2 <- all.tracks.2 %>%
  
  rename(cc = Layer_1,
         ch = eth_ch_utm,
         slope = dtm_slope_deg,
         tpi = tpi_10,
         twi = twi_10)

#_______________________________________________________________________________________________
# 9. Write to csv ----
#_______________________________________________________________________________________________

save.image("Progress/11_11_2024.RData")
