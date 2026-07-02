# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 03 - Sample individuals
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 25 Mar 2025
# Date completed: 25 Mar 2025
# Date last modified: 02 Jul 2026
# R version: 4.5.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation
library(sf)                   # spatial operations
library(mefa4)
library(ctmm)                 # helper function for unit conversion

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________
# 2a. Shapefiles ----
#_______________________________________________________________________

unit.bound <- st_read(paste0(getwd(), "/data_derived/Shapefiles/unit_bound.shp"))
outer.bound <- st_read(paste0(getwd(), "/data_derived/Shapefiles/outer_bound.shp"))

st_crs(unit.bound) <- "epsg:32611"
st_crs(outer.bound) <- "epsg:32611"

# min and max x/y
unit.coord.min <- as.numeric(st_bbox(unit.bound)[1])
unit.coord.max <- as.numeric(st_bbox(unit.bound)[3])

#_______________________________________________________________________
# 3. Numbers for simulations ----

# for each question (Q1 and Q2), we'll simulate 1000 tracks (500 each for target and non-target)

# Notation:
# Question 1 (mean movement model), target individuals: 1T
# Question 1 (mean movement model), non-target individuals: 1NT
  # 3 groups of each for different hypothetical Tau v values

# Question 2 (all movement models), target individuals: 2T
# Question 2 (all movement models), non-target individuals: 2NT

#_______________________________________________________________________

n.sim.tracks <- 500

#_______________________________________________________________________
# 4. Construct dfs of all individuals ----

# here we'll:
# (1) Sample home range centroids
# (2) Sample movement models

# remember that for Q1, we'll use only the "mean" movement models
# and for Q2, it will be a random draw from the full suite

#_______________________________________________________________________
# 4a. 1T ----
#_______________________________________________________________________

df.1T.TV1 <- data.frame(question = 1,
                        target = "T",
                        TV = 1,
                        indiv = 1:n.sim.tracks,
                        hrc.x = runif(n.sim.tracks, unit.coord.min, unit.coord.max),
                        hrc.y = runif(n.sim.tracks, unit.coord.min, unit.coord.max),
                        angle = runif(n.sim.tracks, - pi / 2, pi / 2))

df.1T.TV2 <- data.frame(question = 1,
                        target = "T",
                        TV = 2,
                        indiv = 1:n.sim.tracks,
                        hrc.x = runif(n.sim.tracks, unit.coord.min, unit.coord.max),
                        hrc.y = runif(n.sim.tracks, unit.coord.min, unit.coord.max),
                        angle = runif(n.sim.tracks, - pi / 2, pi / 2))

df.1T.TV3 <- data.frame(question = 1,
                        target = "T",
                        TV = 3,
                        indiv = 1:n.sim.tracks,
                        hrc.x = runif(n.sim.tracks, unit.coord.min, unit.coord.max),
                        hrc.y = runif(n.sim.tracks, unit.coord.min, unit.coord.max),
                        angle = runif(n.sim.tracks, - pi / 2, pi / 2))

#_______________________________________________________________________
# 4b. 1NT ----

# here we'll need to do a spatial sample within the area
sample.1NT <- st_coordinates(st_sample(outer.bound, size = n.sim.tracks))

#_______________________________________________________________________

df.1NT.TV1 <- data.frame(question = 1,
                         target = "NT",
                         TV = 1,
                         indiv = 1:n.sim.tracks,
                         hrc.x = sample.1NT[ , 1],
                         hrc.y = sample.1NT[ , 2],
                         angle = runif(n.sim.tracks, - pi / 2, pi / 2))

df.1NT.TV2 <- data.frame(question = 1,
                         target = "NT",
                         TV = 2,
                         indiv = 1:n.sim.tracks,
                         hrc.x = sample.1NT[ , 1],
                         hrc.y = sample.1NT[ , 2],
                         angle = runif(n.sim.tracks, - pi / 2, pi / 2))

df.1NT.TV3 <- data.frame(question = 1,
                         target = "NT",
                         TV = 3,
                         indiv = 1:n.sim.tracks,
                         hrc.x = sample.1NT[ , 1],
                         hrc.y = sample.1NT[ , 2],
                         angle = runif(n.sim.tracks, - pi / 2, pi / 2))

#_______________________________________________________________________
# 4c. 2T ----
#_______________________________________________________________________

df.2T <- data.frame(question = 2,
                    target = "T",
                    indiv = 1:n.sim.tracks,
                    hrc.x = runif(n.sim.tracks, unit.coord.min, unit.coord.max),
                    hrc.y = runif(n.sim.tracks, unit.coord.min, unit.coord.max),
                    angle = runif(n.sim.tracks, - pi / 2, pi / 2),
                    model = sample(1:1000, size = n.sim.tracks, replace = TRUE))

#_______________________________________________________________________
# 4d. 2NT ----

# here we'll need to do a spatial sample within the area
sample.2NT <- st_coordinates(st_sample(outer.bound, size = n.sim.tracks))

#_______________________________________________________________________

df.2NT <- data.frame(question = 2,
                     target = "NT",
                     indiv = 1:n.sim.tracks,
                     hrc.x = sample.2NT[ , 1],
                     hrc.y = sample.2NT[ , 2],
                     angle = runif(n.sim.tracks, - pi / 2, pi / 2),
                     model = sample(1:1000, size = n.sim.tracks, replace = TRUE))

#_______________________________________________________________________
# 5. Plot home range centroids ----
#_______________________________________________________________________

# bind together
df.all.Q1 <- rbind(df.1T.TV1, df.1T.TV2, df.1T.TV3,
                   df.1NT.TV1, df.1NT.TV2, df.1NT.TV3)

df.all.Q2 <- rbind(df.2T, df.2NT)

# promote to spatial to make a nice map
df.all.Q1.sf <- st_as_sf(df.all.Q1, coords = c("hrc.x", "hrc.y"), crs = 32611)
df.all.Q2.sf <- st_as_sf(df.all.Q2, coords = c("hrc.x", "hrc.y"), crs = 32611)

# Q1
ggplot(data = df.all.Q1.sf) +
  
  theme_bw() +
  
  geom_sf(data = outer.bound,
          fill = NA) +
  
  geom_sf(data = unit.bound,
          fill = NA) +
  
  facet_wrap(~ TV) +
  
  geom_sf(aes(color = target),
          size = 0.65) +
  
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        legend.position = "none") +
  
  scale_color_manual(values = c("gray",
                                "purple"))

# Q1
ggplot(data = df.all.Q2.sf) +
  
  theme_bw() +
  
  geom_sf(data = outer.bound,
          fill = NA) +
  
  geom_sf(data = unit.bound,
          fill = NA) +
  
  geom_sf(aes(color = target),
          size = 0.65) +
  
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        legend.position = "none") +
  
  scale_color_manual(values = c("gray",
                                "purple"))

#_______________________________________________________________________
# 6. Write to file ----
#_______________________________________________________________________

saveRDS(df.all.Q1, "data_derived/sampled_indivs/Q1_indivs.rds")
saveRDS(df.all.Q2, "data_derived/sampled_indivs/Q2_indivs.rds")
