# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 12 - Figures for simulation design diagram
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 20 Jun 2025
# Date completed: 20 Jun 2025
# Date last modified: 20 Jun 2025
# R version: 4.4.3

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation
library(sf)                   # spatial operations
library(cowplot)

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

unit.bound <- st_read(paste0(getwd(), "/Derived data/Shapefiles/unit_bound.shp"))
outer.bound <- st_read(paste0(getwd(), "/Derived data/Shapefiles/outer_bound.shp"))
vs <- st_read(paste0(getwd(), "/Derived data/Shapefiles/cams_vs.shp"))

st_crs(unit.bound) <- "epsg:32611"
st_crs(outer.bound) <- "epsg:32611"
st_crs(vs) <- "epsg:32611"

# min and max x/y
unit.coord.min <- as.numeric(st_bbox(unit.bound)[1])
unit.coord.max <- as.numeric(st_bbox(unit.bound)[3])

# extract S vertex from each VS
vs.vert <- st_as_sf(as.data.frame(st_coordinates(vs)[seq(from = 1,
                                  to = nrow(st_coordinates(vs)),
                                  by = nrow(st_coordinates(vs)) / 9), ]),
                    coords = c("X", "Y"))

st_crs(vs.vert) <- "epsg:32611"

# plot together
plot(st_geometry(outer.bound), col = "blue")
plot(st_geometry(unit.bound), col = "green", add = T)
plot(st_geometry(vs.vert), col = "black", add = T)

#_______________________________________________________________________
# 3. Base activity center samples ----
#_______________________________________________________________________

sample.T <- st_as_sf(st_sample(unit.bound, size = 32))
sample.NT <- st_as_sf(st_sample(outer.bound, size = 64))

# add index variable
sample.T$index <- 1:nrow(sample.T)
sample.NT$index <- 1:nrow(sample.NT)

#_______________________________________________________________________
# 4. Resample ----
#_______________________________________________________________________

# sample targets
T.16 <- sample(sample.T$index, size = 16, replace = FALSE)
T.08 <- sample(T.16, size = 8, replace = FALSE)
T.04 <- sample(T.08, size = 4, replace = FALSE)
T.col <- T.04  # final four are "collared"

# sample non-targets (must be twice the targets)
NT.16 <- sample(sample.NT$index, size = 16 * 2, replace = FALSE)
NT.08 <- sample(NT.16, size = 8 * 2, replace = FALSE)
NT.04 <- sample(NT.08, size = 4 * 2, replace = FALSE)
NT.col <- sample(NT.04, size = 4, replace = FALSE)    # sample a final four for collaring

#_______________________________________________________________________
# 5. Plots ----
#_______________________________________________________________________
# 5a. All targets and non-targets ----
#_______________________________________________________________________

# targets
ggplot() +
  
  theme_bw() +
  
  geom_sf(data = outer.bound,
          fill = NA) +
  
  geom_sf(data = unit.bound,
          fill = NA) +
  
  geom_sf(data = vs.vert,
          shape = 3,
          size = 2) +
  
  geom_sf(data = sample.T,
          color = "purple3") +
  
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank()) -> target.plot

# non-targets
ggplot() +
  
  theme_bw() +
  
  geom_sf(data = outer.bound,
          fill = NA) +
  
  geom_sf(data = unit.bound,
          fill = NA) +
  
  geom_sf(data = vs.vert,
          shape = 3,
          size = 2) +
  
  geom_sf(data = sample.NT,
          color = "purple3") +
  
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank()) -> nontarget.plot

# plot together
plot_grid(target.plot, nontarget.plot, nrow = 2)

# 325 x 455

#_______________________________________________________________________
# 5b. Abundance range ----
#_______________________________________________________________________

# targets
ggplot() +
  
  theme_bw() +
  
  geom_sf(data = outer.bound,
          fill = NA) +
  
  geom_sf(data = unit.bound,
          fill = NA) +
  
  geom_sf(data = vs.vert,
          shape = 3,
          size = 2) +
  
  geom_sf(data = sample.T,
          color = "gray",
          size = 0.5) +
  
  geom_sf(data = sample.T[T.16, ],
          color = "purple3") +
  
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank()) -> target.16.plot

ggplot() +
  
  theme_bw() +
  
  geom_sf(data = outer.bound,
          fill = NA) +
  
  geom_sf(data = unit.bound,
          fill = NA) +
  
  geom_sf(data = vs.vert,
          shape = 3,
          size = 2) +
  
  geom_sf(data = sample.T,
          color = "gray",
          size = 0.5) +
  
  geom_sf(data = sample.T[T.04, ],
          color = "purple3") +
  
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank()) -> target.04.plot

# non-targets
ggplot() +
  
  theme_bw() +
  
  geom_sf(data = outer.bound,
          fill = NA) +
  
  geom_sf(data = unit.bound,
          fill = NA) +
  
  geom_sf(data = vs.vert,
          shape = 3,
          size = 2) +
  
  geom_sf(data = sample.NT,
          color = "gray",
          size = 0.5) +
  
  geom_sf(data = sample.NT[NT.16, ],
          color = "purple3") +
  
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank()) -> nontarget.16.plot

ggplot() +
  
  theme_bw() +
  
  geom_sf(data = outer.bound,
          fill = NA) +
  
  geom_sf(data = unit.bound,
          fill = NA) +
  
  geom_sf(data = vs.vert,
          shape = 3,
          size = 2) +
  
  geom_sf(data = sample.NT,
          color = "gray",
          size = 0.5) +
  
  geom_sf(data = sample.NT[NT.04, ],
          color = "purple3") +
  
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank()) -> nontarget.04.plot

# plot together
plot_grid(target.16.plot, target.04.plot,
          nontarget.16.plot, nontarget.04.plot,
          nrow = 2)

# 539 x 337

#_______________________________________________________________________
# 5c. Collar sample ----
#_______________________________________________________________________

# targets
ggplot() +
  
  theme_bw() +
  
  geom_sf(data = outer.bound,
          fill = NA) +
  
  geom_sf(data = unit.bound,
          fill = NA) +
  
  geom_sf(data = sample.T,
          color = "gray",
          size = 1) +
  
  geom_sf(data = sample.T[T.04, ],
          color = "green4",
          size = 2) +
  
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank()) -> target.collar.plot

# non-targets
ggplot() +
  
  theme_bw() +
  
  geom_sf(data = outer.bound,
          fill = NA) +
  
  geom_sf(data = unit.bound,
          fill = NA) +
  
  geom_sf(data = sample.NT,
          color = "gray",
          size = 1) +
  
  geom_sf(data = sample.NT[NT.col, ],
          color = "green4",
          size = 2) +
  
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank()) -> nontarget.collar.plot

# plot together
plot_grid(target.collar.plot,
          nontarget.collar.plot,
          nrow = 1)

# 487 x 315

#_______________________________________________________________________
# 6. Illustrative parameter distributions ----
#_______________________________________________________________________

# "typical" individual
typical.df <- data.frame(x = seq(0, 40, 0.01),
                         y = dlnorm(seq(0, 40, 0.01), 
                                    meanlog = 2.25, 
                                    sdlog = 0.75))

ggplot(data = typical.df,
       aes(x = x,
           y = y)) +
  
  theme_bw() +
  
  geom_line(color = "purple3") +
  
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) -> lo.param

# full population
full.df <- data.frame()

# loop
for (i in 1:15) {
  
  focal.df <- data.frame(x = seq(0, 40, 0.01),
                         y = dlnorm(seq(0, 40, 0.01), 
                                    meanlog = rnorm(1, 2.25, 0.50),
                                    sdlog = rnorm(1, 0.75, 0.10)),
                         i = i)
  
  full.df <- rbind(full.df, focal.df)
  
}

# weighted df
weighted.df <- data.frame()

# loop
for (i in 1:15) {
  
  focal.df <- data.frame(x = seq(0, 40, 0.01),
                         y = dlnorm(seq(0, 40, 0.01), 
                                    meanlog = rnorm(1, 3.25, 0.25),
                                    sdlog = rnorm(1, 0.25, 0.25)),
                         i = i)
  
  weighted.df <- rbind(weighted.df, focal.df)
  
}

# all together
ggplot() +
  
  theme_bw() +
  
  geom_line(data = weighted.df,
            aes(x = x,
                y = y,
                group = factor(i)),
            color = "purple",
            alpha = 0.5) +
  
  geom_line(data = full.df,
            aes(x = x,
                y = y,
                group = factor(i)),
            color = "purple",
            alpha = 0.5) +
  
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") -> hi.param

# weighted
ggplot() +
  
  theme_bw() +
  
  geom_line(data = full.df,
            aes(x = x,
                y = y,
                group = factor(i)),
            color = "gray",
            alpha = 0.5) +
  
  geom_line(data = weighted.df,
            aes(x = x,
                y = y,
                group = factor(i)),
            color = "purple",
            alpha = 0.5) +
  
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") -> weight.param

# plot all together
plot_grid(lo.param, lo.param,
          hi.param, hi.param,
          hi.param, weight.param,
          nrow = 3)

# 299 x 402