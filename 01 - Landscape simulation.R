# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 01 - Landscape simulation
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 15 Nov 2024
# Date completed: 15 Nov 2024
# Date last modified: 10 Feb 2025
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation
library(terra)           # work with rasters
library(NLMR)            # simulate landscapes
library(spatialEco)      # Gaussian blur
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
hectares <- 501

# square meters
n.sqm <- hectares * 10000

# desired resolution (in m)
resol <- 10

# calculate required number of cells per dimension
columns <- as.integer(ceiling(sqrt(n.sqm / resol)))
rows <- as.integer(ceiling(sqrt(n.sqm / resol)))

#_______________________________________________________________________
# 3. Covariate 1 - "forage" ----

# strong selection
# decreased speed with higher forage ("foraging behavior")
# we'll use the Gaussian random field function here

#_______________________________________________________________________
# 3a. Before landscape ----

# high autocorrelation in this continuous variable

# parameters
# autocorrelation range: maximum range of spatial autocorrelation
B.cov1.autocorr <- 50

# magnitude of variation
B.cov1.magvar <- 2

#_______________________________________________________________________

# B1
set.seed(89)
B1.cov1 <- nlm_gaussianfield(ncol = columns,
                             nrow = rows,
                             resolution = resol,
                             autocorr_range = B.cov1.autocorr,
                             mag_var = B.cov1.magvar)

# C2
set.seed(923)
B2.cov1 <- nlm_gaussianfield(ncol = columns,
                             nrow = rows,
                             resolution = resol,
                             autocorr_range = B.cov1.autocorr,
                             mag_var = B.cov1.magvar)

# C3
set.seed(12)
B3.cov1 <- nlm_gaussianfield(ncol = columns,
                             nrow = rows,
                             resolution = resol,
                             autocorr_range = B.cov1.autocorr,
                             mag_var = B.cov1.magvar)

# bind together
B.cov1 <- c(rast(B1.cov1), rast(B2.cov1), rast(B3.cov1))
names(B.cov1) <- c("B1", "B2", "B3")

# add UTM CRS so spatial layers work together
crs(B.cov1) <- crs("EPSG:32611")

# plot
plot(B.cov1)

#_______________________________________________________________________
# 4. Covariate 2 - "edge distance" ----

# weak avoidance (weak selection for proximity)

#_______________________________________________________________________
# 4a. Create a unit boundary polygon ----
#_______________________________________________________________________

# extract coordinates - find centroid of raster
rast.centroid <- c(B1.cov1@extent@xmax / 2,
                   B1.cov1@extent@ymax / 2)

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
                           y = B1.cov1, 
                           mask = TRUE)

# create distance to edge of unit raster
cov2 <- distance(boundary.rast)

# promote to SpatRaster
cov2 <- rast(cov2)

# scale between 0 and 1 for easy comparison
cov2 <- cov2 / max(values(cov2))

# att UTM crs
crs(cov2) <- crs("EPSG:32611")

# plot
plot(cov2)

 #_______________________________________________________________________
# 5. Covariate 3 - "cover type" ----

# this is what we'll manipulate
# strong base avoidance of cover type 1 ( we can think of it as open cover)
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
  
  # add crs
  crs(ls.class.1) <- crs("EPSG:32611")
  
  # crop to be the same size as the other variables
  ls.class.2 <- terra::crop(x = ls.class.1, 
                            y = B1.cov1)
  
  # return
  return(ls.class.2)
  
}

#_______________________________________________________________________
# 5a. Before landscape ----

# high "roughness" = 0.9

#_______________________________________________________________________

# run function 
B1.cov3 <- create_discrete_ls(1284, 0.9)
B2.cov3 <- create_discrete_ls(78, 0.9)
B3.cov3 <- create_discrete_ls(394, 0.9)

# bind together
B.cov3 <- c(rast(B1.cov3), rast(B2.cov3), rast(B3.cov3))
names(B.cov3) <- c("B1", "B2", "B3")

#_______________________________________________________________________
# 5b. After landscape ----
#_______________________________________________________________________

# first, we'll crop each one by the unit boundary to examine what we'll change
B.cov3.crop <- crop(B.cov3, unit.bound.sf)

plot(B.cov3.crop)

# now we'll simulate alteration by flipping zeroes ("closed") to ones ("open")
# each landscape has a different proportion of 0s and 1s, this will be good!
# write a function that does it
alter_cover <- function(landscape,
                        percent) {
  
  # focal landscape
  landscape.1 <- landscape
  
  # how many can we flip?
  n.0 <- length(values(landscape.1)[values(landscape.1) == 0])
  
  # randomly sample a percentage of them to flip
  to.flip <- sample(which(values(landscape.1) == 0), size = round(n.0 * percent))
  
  # flip them
  values(landscape.1)[to.flip] <- 1
  
  # return
  return(landscape.1)
  
}

# use function
A.cov3.crop <- c(alter_cover(B.cov3.crop$B1, 0.30),
                 alter_cover(B.cov3.crop$B2, 0.30),
                 alter_cover(B.cov3.crop$B3, 0.30))

plot(B.cov3.crop)
plot(A.cov3.crop)

# looks decent!
# now just to mosaic into the landscape
A1.cov3 <- merge(A.cov3.crop$B1, B.cov3$B1)
A2.cov3 <- merge(A.cov3.crop$B2, B.cov3$B2)
A3.cov3 <- merge(A.cov3.crop$B3, B.cov3$B3)

#_______________________________________________________________________
# 6. Bind each landscape / replicate together ----
#_______________________________________________________________________
# 6a. Before ----
#_______________________________________________________________________

# bind
B1.ls <- c(B.cov1$B1, cov2, B.cov3$B1)
B2.ls <- c(B.cov1$B2, cov2, B.cov3$B2)
B3.ls <- c(B.cov1$B3, cov2, B.cov3$B3)

# rename
names(B1.ls) <- c("forage", "edge", "open")
names(B2.ls) <- c("forage", "edge", "open")
names(B3.ls) <- c("forage", "edge", "open")

#_______________________________________________________________________
# 6b. After ----
#_______________________________________________________________________

# bind
A1.ls <- c(B.cov1$B1, cov2, A1.cov3)
A2.ls <- c(B.cov1$B2, cov2, A2.cov3)
A3.ls <- c(B.cov1$B3, cov2, A3.cov3)

# rename
names(A1.ls) <- c("forage", "edge", "open")
names(A2.ls) <- c("forage", "edge", "open")
names(A3.ls) <- c("forage", "edge", "open")

#_______________________________________________________________________
# 7. Plot ----
#_______________________________________________________________________
# 7a. Crop for visualization ----
#_______________________________________________________________________

# define extent (buffered unit)
unit.buff <- st_buffer(unit.bound.sf, dist = 100)

# crop
B1.ls.crop <- crop(B1.ls, unit.buff)
B2.ls.crop <- crop(B2.ls, unit.buff)
B3.ls.crop <- crop(B3.ls, unit.buff)

A1.ls.crop <- crop(A1.ls, unit.buff)
A2.ls.crop <- crop(A2.ls, unit.buff)
A3.ls.crop <- crop(A3.ls, unit.buff)

#_______________________________________________________________________
# 7b. Define function ----
#_______________________________________________________________________

plot_group_rast <- function (ls) {
  
  # forage
  ggplot() +
    
    theme_bw() +
    
    # add in raster
    geom_spatraster(data = ls$forage) +
    
    # viridis colors
    scale_fill_viridis_c() +
    
    # add unit boundary
    geom_sf(data = unit.bound.sf,
            fill = NA,
            color = "white") +
    
    # remove axis text
    theme(axis.text = element_blank(),
          legend.position = "none") -> forage
  
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
    geom_spatraster(data = ls$open) +
    
    # viridis colors
    scale_fill_viridis_c() +
    
    # add unit boundary
    geom_sf(data = unit.bound.sf,
            fill = NA,
            color = "white") +
    
    # remove axis text
    theme(axis.text = element_blank(),
          legend.position = "none") -> open
  
  # plot together
  focal.plot <- plot_grid(forage, edge, open, nrow = 1)
  
  # return
  return(focal.plot)
  
}

#_______________________________________________________________________
# 7c. Assign plots ----
#_______________________________________________________________________

# simple
B1.plot <- plot_group_rast(B1.ls.crop)
B2.plot <- plot_group_rast(B2.ls.crop)
B3.plot <- plot_group_rast(B3.ls.crop)

# complex
A1.plot <- plot_group_rast(A1.ls.crop)
A2.plot <- plot_group_rast(A2.ls.crop)
A3.plot <- plot_group_rast(A3.ls.crop)

#_______________________________________________________________________
# 7d. Plot together ----
#_______________________________________________________________________

# simple
plot_grid(B1.plot, B2.plot, B3.plot, nrow = 3)

# complex
plot_grid(A1.plot, A2.plot, A3.plot, nrow = 3)

#_______________________________________________________________________
# 8. Save and write rasters ----
#_______________________________________________________________________

save.image("Progress/rasters_02_10_2025.RData")

writeRaster(B1.ls, filename = "Rasters/B1.tif", overwrite = T)
writeRaster(B2.ls, filename = "Rasters/B2.tif", overwrite = T)
writeRaster(B3.ls, filename = "Rasters/B3.tif", overwrite = T)

writeRaster(A1.ls, filename = "Rasters/A1.tif", overwrite = T)
writeRaster(A2.ls, filename = "Rasters/A2.tif", overwrite = T)
writeRaster(A3.ls, filename = "Rasters/A3.tif", overwrite = T)

