# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 01 - Landscape simulation
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 15 Nov 2024
# Date completed: 15 Nov 2024
# Date last modified: 15 Nov 2024
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation
library(terra)           # work with rasters
library(NLMR)            # simulate landscapes
library(landscapetools)  # work with simulated landscapes
library(sf)              # bounding polygons
library(tidyterra)       # plot in ggplot
library(cowplot)         # multiple plots

#_______________________________________________________________________
# 2. Define landscape size
#_______________________________________________________________________

# here we want a landscape large enough so simulated tracks in our focal area (10 ha)
# don't tend to move outside

# how many hectares?
hectares <- 500

# square meters
n.sqm <- hectares * 10000

# desired resolution (in m)
resol <- 1

# calculate required number of cells per dimension
columns <- as.integer(ceiling(sqrt(n.sqm / resol))) + 2
rows <- as.integer(ceiling(sqrt(n.sqm / resol))) + 2

#_______________________________________________________________________
# 3. Covariate 1 - "stem density" ----

# simulated hares will select for stem density
# this will vary across the entire landscape, not just in the "unit"
# we'll use the midpoint displacement neutral landscape model here

#_______________________________________________________________________
# 3a. Simple landscape
#_______________________________________________________________________

set.seed(84)

simple.cov1 <- nlm_mpd(ncol = columns, 
                       nrow = rows, 
                       resolution = resol, 
                       roughness = 0.2)

simple.cov1

plot(simple.cov1)

# this one looks decent

# add UTM CRS so spatial layers work together
crs(simple.cov1) <- crs("EPSG:32611")

#_______________________________________________________________________
# 3b. Complex landscape
#_______________________________________________________________________

set.seed(82)

complex.cov1 <- nlm_mpd(ncol = columns, 
                        nrow = rows, 
                        resolution = resol, 
                        roughness = 0.7)

complex.cov1

plot(complex.cov1)

# add UTM CRS so spatial layers work together
crs(complex.cov1) <- crs("EPSG:32611")

#_______________________________________________________________________
# 4. Covariate 2 - "edge distance" ----

# simulated hares will avoid unit edge distance

#_______________________________________________________________________
# 4a. Create a unit boundary polygon ----
#_______________________________________________________________________

# extract coordinates - find centroid of raster
rast.centroid <- c(simple.cov1@extent@xmax / 2,
                   simple.cov1@extent@ymax / 2)

# how many meters per side? (let's make our unit 10 ha)
m.side <- sqrt(10 * 10000)

# how many meters to add and subtract?
m.side.half <- m.side / 2

# create polygon feature
unit.bound <- st_polygon(list(cbind(c(rast.centroid[1] - m.side.half, 
                                      rast.centroid[1] + m.side.half,
                                      rast.centroid[1] + m.side.half,
                                      rast.centroid[1] - m.side.half,
                                      rast.centroid[1] - m.side.half), 
                                    c(rast.centroid[2] - m.side.half, 
                                      rast.centroid[2] - m.side.half,
                                      rast.centroid[2] + m.side.half,
                                      rast.centroid[2] + m.side.half,
                                      rast.centroid[2] - m.side.half))))

# plot
plot(complex.cov1)
plot(unit.bound, add = T)

# assign to sf object
unit.bound.sf <- st_as_sf(st_sfc(unit.bound))

st_crs(unit.bound.sf) <- crs("EPSG:32611")

# write to shapefile
st_write(unit.bound.sf,
         dsn = paste0(getwd(), "/Derived_data/Shapefiles/unit_bound.shp"),
         layer = "unit_bound.shp")

#_______________________________________________________________________
# 4b. Define an edge variable ----
#_______________________________________________________________________

# convert unit boundary to SpatialLines
unit.bound.lines <- as(as_Spatial(st_sfc(unit.bound)), "SpatialLines")

# rasterize the unit boundary
boundary.rast <- rasterize(x = unit.bound.lines, 
                           y = simple.cov1, 
                           mask = TRUE)

# create distance to edge of unit raster
cov2 <- distance(boundary.rast)

plot(cov2)

crs(cov2) <- crs("EPSG:32611")

#_______________________________________________________________________
# 5. Covariate 3 - "mature forest" ----

# simulated hares will weakly avoid mature forest in general
# they will select for it when steps that start in it are longer
# while avoid it when steps that start in it are shorter

# to do this, we'll classify a midpoint landscape

#_______________________________________________________________________
# 5a. Simple landscape
#_______________________________________________________________________

set.seed(65)

simple.cov3.start <- nlm_mpd(ncol = columns, 
                             nrow = rows, 
                             resolution = resol, 
                             roughness = 0.5,
                             rand_dev = 1.5)

plot(simple.cov3.start)

simple.cov3 <- util_classify(simple.cov3.start,
                             weighting = c(0.5, 0.5))
                             
                             #,
                             #level_names = c("other", "mature"))

plot(simple.cov3)
plot(unit.bound.sf, add = T)

crs(simple.cov3) <- crs("EPSG:32611")

# reclassify
simple.cov3 <- simple.cov3 - 1

#_______________________________________________________________________
# 5b. Complex landscape
#_______________________________________________________________________

set.seed(66)

complex.cov3.start <- nlm_mpd(ncol = columns, 
                              nrow = rows, 
                              resolution = resol, 
                              roughness = 0.9,
                              rand_dev = 1.5)

plot(complex.cov3.start)

complex.cov3 <- util_classify(complex.cov3.start,
                              weighting = c(0.5, 0.5))

                              #,
                              #level_names = c("other", "mature"))

plot(complex.cov3)
plot(unit.bound.sf, add = T)

crs(complex.cov3) <- crs("EPSG:32611")

# reclassify
complex.cov3 <- complex.cov3 - 1

#_______________________________________________________________________
# 5. Scale continuous landscapes ----
#_______________________________________________________________________

# cov1
simple.cov1.s <- scale(simple.cov1)
complex.cov1.s <- scale(complex.cov1)

# cov2
cov2.s <- scale(cov2)

#_______________________________________________________________________
# 6. Promote all to SpatRaster and bind ----
#_______________________________________________________________________

simple.cov1.spat <- rast(simple.cov1.s)
complex.cov1.spat <- rast(complex.cov1.s)
cov2.spat <- rast(cov2.s)
simple.cov3.spat <- rast(simple.cov3)
complex.cov3.spat <- rast(complex.cov3)

# bind together
simple <- c(simple.cov1.spat,
            cov2.spat,
            simple.cov3.spat)

names(simple) <- c("stem", "edge", "mature")

complex <- c(complex.cov1.spat,
             cov2.spat,
             complex.cov3.spat)

names(complex) <- c("stem", "edge", "mature")

#_______________________________________________________________________
# 7. Plot ----
#_______________________________________________________________________

# simple
ggplot() +
  
  theme_bw() +
  
  # add in raster
  geom_spatraster(data = simple) +

  # facet
  facet_wrap(~ lyr) +
  
  # viridis colors
  scale_fill_viridis_c() +
  
  # add unit boundary
  geom_sf(data = unit.bound.sf,
          fill = NA,
          color = "white") +
  
  # remove axis text
  theme(axis.text = element_blank()) -> simple.plot

# complex
ggplot() +
  
  theme_bw() +
  
  # add in raster
  geom_spatraster(data = complex) +
  
  # facet
  facet_wrap(~ lyr) +
  
  # viridis colors
  scale_fill_viridis_c() +
  
  # add unit boundary
  geom_sf(data = unit.bound.sf,
          fill = NA,
          color = "white") +
  
  # remove axis text
  theme(axis.text = element_blank()) -> complex.plot

# plot together
plot_grid(simple.plot, complex.plot, nrow = 2)

#_______________________________________________________________________
# 8. Save and write rasters ----
#_______________________________________________________________________

save.image("Progress/rasters_11_15_2024.RData")

writeRaster(simple, filename = "Rasters/simple.tif", overwrite = T)
writeRaster(complex, filename = "Rasters/complex.tif", overwrite = T)
