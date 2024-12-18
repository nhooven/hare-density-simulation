# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 01 - Landscape simulation
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 15 Nov 2024
# Date completed: 15 Nov 2024
# Date last modified: 18 Dec 2024
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
hectares <- 1000

# square meters
n.sqm <- hectares * 10000

# desired resolution (in m)
resol <- 10

# calculate required number of cells per dimension
columns <- as.integer(ceiling(sqrt(n.sqm / resol)))
rows <- as.integer(ceiling(sqrt(n.sqm / resol)))

#_______________________________________________________________________
# 3. Covariate 1 - "stem density" ----

# strong selection
# decreased speed with higher stem
# we'll use the planar gradient function here, with a random orientation

#_______________________________________________________________________
# 3a. Simple landscape ----
#_______________________________________________________________________

plot(nlm_planargradient(ncol = columns + 1,
                     nrow = rows + 1,
                     resolution = resol))


# S1
set.seed(74)
S1.cov1 <- nlm_planargradient(ncol = columns,
                              nrow = rows,
                              resolution = resol)

# S2
set.seed(539)
S2.cov1 <- nlm_planargradient(ncol = columns,
                              nrow = rows,
                              resolution = resol)

# S3
set.seed(1010)
S3.cov1 <- nlm_planargradient(ncol = columns,
                              nrow = rows,
                              resolution = resol)

# bind together
S.cov1 <- c(rast(S1.cov1), rast(S2.cov1), rast(S3.cov1))
names(S.cov1) <- c("S1", "S2", "S3")

# add UTM CRS so spatial layers work together
crs(S.cov1) <- crs("EPSG:32611")

# plot
plot(S.cov1)

#_______________________________________________________________________
# 3b. Complex landscape ----
#_______________________________________________________________________

# C1
set.seed(841)
C1.cov1 <- nlm_mpd(ncol = columns + 2,  # mpd will force these to be odd
                   nrow = rows + 2, 
                   resolution = resol, 
                   roughness = 0.7)

# C2
set.seed(45)
C2.cov1 <- nlm_mpd(ncol = columns + 2,  # mpd will force these to be odd
                   nrow = rows + 2, 
                   resolution = resol, 
                   roughness = 0.7)

# C3
set.seed(989)
C3.cov1 <- nlm_mpd(ncol = columns + 2,  # mpd will force these to be odd
                   nrow = rows + 2, 
                   resolution = resol, 
                   roughness = 0.7)

# bind together
C.cov1 <- c(rast(C1.cov1), rast(C2.cov1), rast(C3.cov1))
names(C.cov1) <- c("C1", "C2", "C3")

# add UTM CRS so spatial layers work together
crs(C.cov1) <- crs("EPSG:32611")

# plot
plot(C.cov1)

#_______________________________________________________________________
# 4. Covariate 2 - "edge distance" ----

# weak avoidance (weak selection for proximity)

#_______________________________________________________________________
# 4a. Create a unit boundary polygon ----
#_______________________________________________________________________

# extract coordinates - find centroid of raster
rast.centroid <- c(S1.cov1@extent@xmax / 2,
                   S1.cov1@extent@ymax / 2)

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

# assign to sf object
unit.bound.sf <- st_as_sf(st_sfc(unit.bound))

st_crs(unit.bound.sf) <- crs("EPSG:32611")

# write to shapefile
st_write(unit.bound.sf,
         dsn = paste0(getwd(), "/Derived_data/Shapefiles/unit_bound.shp"),
         layer = "unit_bound.shp",
         append = FALSE)

#_______________________________________________________________________
# 4b. Define an edge variable ----
#_______________________________________________________________________

# convert unit boundary to SpatialLines
unit.bound.lines <- as(as_Spatial(st_sfc(unit.bound)), "SpatialLines")

# rasterize the unit boundary
boundary.rast <- rasterize(x = unit.bound.lines, 
                           y = S1.cov1, 
                           mask = TRUE)

# create distance to edge of unit raster
cov2 <- distance(boundary.rast)

# promote to SpatRaster
cov2 <- rast(cov2)

# scale between 0 and 1 for easy comparison
cov2 <- cov2 / max(values(cov2))

# att UTM crs
crs(cov2) <- crs("EPSG:32611")

#_______________________________________________________________________
# 5. Covariate 3 - "mature forest" ----

# weak base avoidance
# longer steps that start in it

# to do this, we'll classify a midpoint landscape

# define function
create_discrete_ls <- function (seed = runif(1, 1, 1000),
                                rough = 0.5) {
  
  # set seed
  set.seed(seed)
  
  # generate landscape
  start.ls <- nlm_mpd(ncol = columns + 2, 
                      nrow = rows + 2, 
                      resolution = resol, 
                      roughness = rough,
                      rand_dev = 1.5)
  
  # classify into two discrete classes
  ls.class <- util_classify(start.ls,
                            weighting = c(0.5, 0.5))
  
  # reclassify
  ls.class.1 <- ls.class - 1
  
  # return
  return(ls.class.1)
  
}

#_______________________________________________________________________
# 5a. Simple landscape ----
#_______________________________________________________________________

# run function 
S1.cov3 <- create_discrete_ls(455, 0.5)
S2.cov3 <- create_discrete_ls(1889, 0.5)
S3.cov3 <- create_discrete_ls(101, 0.5)

# bind together
S.cov3 <- c(rast(S1.cov3), rast(S2.cov3), rast(S3.cov3))
names(S.cov3) <- c("S1", "S2", "S3")

crs(S.cov3) <- crs("EPSG:32611")

#_______________________________________________________________________
# 5b. Complex landscape
#_______________________________________________________________________

# run function 
C1.cov3 <- create_discrete_ls(1603, 0.9)
C2.cov3 <- create_discrete_ls(78, 0.9)
C3.cov3 <- create_discrete_ls(366, 0.9)

# bind together
C.cov3 <- c(rast(C1.cov3), rast(C2.cov3), rast(C3.cov3))
names(C.cov3) <- c("C1", "C2", "C3")

crs(C.cov3) <- crs("EPSG:32611")

#_______________________________________________________________________
# 5. Bind each landscape / replicate together ----
#_______________________________________________________________________
# 5a. Simple ----
#_______________________________________________________________________

# bind
S1.ls <- c(S.cov1$S1, cov2, S.cov3$S1)
S2.ls <- c(S.cov1$S2, cov2, S.cov3$S2)
S3.ls <- c(S.cov1$S3, cov2, S.cov3$S3)

# rename
names(S1.ls) <- c("stem", "edge", "mature")
names(S2.ls) <- c("stem", "edge", "mature")
names(S3.ls) <- c("stem", "edge", "mature")

#_______________________________________________________________________
# 5b. Complex ----
#_______________________________________________________________________

# bind
C1.ls <- c(C.cov1$C1, cov2, C.cov3$C1)
C2.ls <- c(C.cov1$C2, cov2, C.cov3$C2)
C3.ls <- c(C.cov1$C3, cov2, C.cov3$C3)

# rename
names(C1.ls) <- c("stem", "edge", "mature")
names(C2.ls) <- c("stem", "edge", "mature")
names(C3.ls) <- c("stem", "edge", "mature")

#_______________________________________________________________________
# 7. Plot ----
#_______________________________________________________________________
# 7a. Crop for visualization ----
#_______________________________________________________________________

# define extent (buffered unit)
unit.buff <- st_buffer(unit.bound.sf, dist = 100)

# crop
S1.ls.crop <- crop(S1.ls, unit.buff)
S2.ls.crop <- crop(S2.ls, unit.buff)
S3.ls.crop <- crop(S3.ls, unit.buff)

C1.ls.crop <- crop(C1.ls, unit.buff)
C2.ls.crop <- crop(C2.ls, unit.buff)
C3.ls.crop <- crop(C3.ls, unit.buff)

#_______________________________________________________________________
# 7b. Define function ----
#_______________________________________________________________________

plot_group_rast <- function (ls) {
  
  # stem
  ggplot() +
    
    theme_bw() +
    
    # add in raster
    geom_spatraster(data = ls$stem) +
    
    # viridis colors
    scale_fill_viridis_c() +
    
    # add unit boundary
    geom_sf(data = unit.bound.sf,
            fill = NA,
            color = "white") +
    
    # remove axis text
    theme(axis.text = element_blank(),
          legend.position = "none") -> stem
  
  # edge
  ggplot() +
    
    theme_bw() +
    
    # add in raster
    geom_spatraster(data = ls$edge) +
    
    # viridis colors
    scale_fill_viridis_c() +
    
    # add unit boundary
    geom_sf(data = unit.bound.sf,
            fill = NA,
            color = "white") +
    
    # remove axis text
    theme(axis.text = element_blank(),
          legend.position = "none") -> edge
  
  # mature
  ggplot() +
    
    theme_bw() +
    
    # add in raster
    geom_spatraster(data = ls$mature) +
    
    # viridis colors
    scale_fill_viridis_c() +
    
    # add unit boundary
    geom_sf(data = unit.bound.sf,
            fill = NA,
            color = "white") +
    
    # remove axis text
    theme(axis.text = element_blank(),
          legend.position = "none") -> mature
  
  # plot together
  focal.plot <- plot_grid(stem, edge, mature, nrow = 1)
  
  # return
  return(focal.plot)
  
}

#_______________________________________________________________________
# 7c. Assign plots ----
#_______________________________________________________________________

# simple
S1.plot <- plot_group_rast(S1.ls.crop)
S2.plot <- plot_group_rast(S2.ls.crop)
S3.plot <- plot_group_rast(S3.ls.crop)

# complex
C1.plot <- plot_group_rast(C1.ls.crop)
C2.plot <- plot_group_rast(C2.ls.crop)
C3.plot <- plot_group_rast(C3.ls.crop)

#_______________________________________________________________________
# 7d. Plot together ----
#_______________________________________________________________________

# simple
plot_grid(S1.plot, S2.plot, S3.plot, nrow = 3)

# complex
plot_grid(C1.plot, C2.plot, C3.plot, nrow = 3)

#_______________________________________________________________________
# 8. Save and write rasters ----
#_______________________________________________________________________

save.image("Progress/rasters_12_18_2024.RData")

writeRaster(S1.ls, filename = "Rasters/S1.tif", overwrite = T)
writeRaster(S2.ls, filename = "Rasters/S2.tif", overwrite = T)
writeRaster(S3.ls, filename = "Rasters/S3.tif", overwrite = T)

writeRaster(C1.ls, filename = "Rasters/C1.tif", overwrite = T)
writeRaster(C2.ls, filename = "Rasters/C2.tif", overwrite = T)
writeRaster(C3.ls, filename = "Rasters/C3.tif", overwrite = T)
