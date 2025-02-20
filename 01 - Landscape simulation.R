# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 01 - Landscape simulation
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 15 Nov 2024
# Date completed: 15 Nov 2024
# Date last modified: 20 Feb 2025
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
# 3. Covariate 1 - "forage" - fora ----

# Start: Faster movements when starting in in low forage
# End: More likely to switch direction to get to high forage
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
# 4. Create a unit boundary polygon ----
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
# 5. Covariate 2 - "elevation" - elev ----

# intermediate (i.e., quadratic) selection
# we'll use a planar gradient

#_______________________________________________________________________

# B1
set.seed(558)
B1.cov2 <- nlm_planargradient(ncol = columns,
                              nrow = rows,
                              resolution = resol)

# C2
set.seed(1072)
B2.cov2 <- nlm_planargradient(ncol = columns,
                              nrow = rows,
                              resolution = resol)

# C3
set.seed(223)
B3.cov2 <- nlm_planargradient(ncol = columns,
                              nrow = rows,
                              resolution = resol)

# bind together
B.cov2 <- c(rast(B1.cov2), rast(B2.cov2), rast(B3.cov2))
names(B.cov2) <- c("B1", "B2", "B3")

# add UTM CRS so spatial layers work together
crs(B.cov2) <- crs("EPSG:32611")

# plot
plot(B.cov2)

#_______________________________________________________________________
# 5. Covariate 3 - "openness" - open ----

# this is what we'll manipulate
# Start: Faster movements when starting in higher openness
# End: Negative selection

# midpoint landscape, classify, the calculate percent open within a buffer

#_______________________________________________________________________
# 5a. Define function to create discrete landscape ----
#_______________________________________________________________________

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
  
  crs(ls.class.1) <- crs("EPSG:32611")
  
  # crop to be the same size as the other variables
  ls.class.2 <- terra::crop(x = ls.class.1, 
                            y = B1.cov1)
  
  
  # return
  return(ls.class.2)

}

#_______________________________________________________________________
# 5a. Before landscape - discrete ----
#_______________________________________________________________________

B1.cov3.dis <- create_discrete_ls(1987, 0.9)
B2.cov3.dis <- create_discrete_ls(78, 0.9)
B3.cov3.dis <- create_discrete_ls(394, 0.9)

# bind together
B.cov3.dis <- c(rast(B1.cov3.dis), rast(B2.cov3.dis), rast(B3.cov3.dis))
names(B.cov3.dis) <- c("B1", "B2", "B3")

#_______________________________________________________________________
# 5b. Define function to calculate percent open in a buffer ----
#_______________________________________________________________________

# define function
perc_open <- function (landscape,
                       buffer) {
  
  # define focal matrix
  focal.mat <- focalMat(x = landscape,
                        d = buffer,
                        type = "circle")
  
  # replace all weights with 1
  focal.mat[focal.mat > 0] <- 1
  
  # sum these to determine how many cells we're dealing with
  sum.focal.mat <- sum(focal.mat)
  
  # use focal to calculate percent 1
  # define function 
  # https://gis.stackexchange.com/a/410811
  pclass <- function (x) {
    
    return(length(which(x == 1)) / length(x))
    
    }
  
  # run
  focal.rast <- terra::focal(x = landscape,
                             w = focal.mat,
                             fun = pclass)
  
  # return
  return(focal.rast)
  
}

#_______________________________________________________________________
# 5c. Before - calculate percent open ----
#_______________________________________________________________________

B.cov3 <- perc_open(B.cov3.dis, 50)

B.cov3.crop <- crop(B.cov3, unit.bound.sf)

plot(B.cov3.crop)

#_______________________________________________________________________
# 5b. After landscape ----
#_______________________________________________________________________

# now we'll simulate alteration by:
# (1) splitting each landscape into four quadrants
# (2) choosing two diagonal quadrants at random
# (3) flipping 30% of non-open pixels to open

# create quadrants - BL, BR, TR, TL
# initialize coordinates
quad.1.df <- data.frame(x = c(st_bbox(unit.bound.sf)[1],
                              st_bbox(unit.bound.sf)[1] + sqrt(25000),
                              st_bbox(unit.bound.sf)[1] + sqrt(25000),
                              st_bbox(unit.bound.sf)[1]),
                        y = c(st_bbox(unit.bound.sf)[4] - sqrt(25000),
                              st_bbox(unit.bound.sf)[4] - sqrt(25000),
                              st_bbox(unit.bound.sf)[4],
                              st_bbox(unit.bound.sf)[4]))

quad.2.df <- data.frame(x = c(st_bbox(unit.bound.sf)[1] + sqrt(25000),
                              st_bbox(unit.bound.sf)[3],
                              st_bbox(unit.bound.sf)[3],
                              st_bbox(unit.bound.sf)[1] + sqrt(25000)),
                        y = c(st_bbox(unit.bound.sf)[4] - sqrt(25000),
                              st_bbox(unit.bound.sf)[4] - sqrt(25000),
                              st_bbox(unit.bound.sf)[4],
                              st_bbox(unit.bound.sf)[4]))

quad.3.df <- data.frame(x = c(st_bbox(unit.bound.sf)[1],
                              st_bbox(unit.bound.sf)[1] + sqrt(25000),
                              st_bbox(unit.bound.sf)[1] + sqrt(25000),
                              st_bbox(unit.bound.sf)[1]),
                        y = c(st_bbox(unit.bound.sf)[2],
                              st_bbox(unit.bound.sf)[2],
                              st_bbox(unit.bound.sf)[2] + sqrt(25000),
                              st_bbox(unit.bound.sf)[2] + sqrt(25000)))

quad.4.df <- data.frame(x = c(st_bbox(unit.bound.sf)[1] + sqrt(25000),
                              st_bbox(unit.bound.sf)[3],
                              st_bbox(unit.bound.sf)[3],
                              st_bbox(unit.bound.sf)[1] + sqrt(25000)),
                        y = c(st_bbox(unit.bound.sf)[2],
                              st_bbox(unit.bound.sf)[2],
                              st_bbox(unit.bound.sf)[2] + sqrt(25000),
                              st_bbox(unit.bound.sf)[2] + sqrt(25000)))

# convert to polygons
quad.sf.list <- list()

quad.sf.list[[1]] <- quad.1.df %>% 
  
  st_as_sf(coords = c("x", "y"), crs = 32611) %>% 
  
  summarise(geometry = st_combine(geometry)) %>%
  
  st_cast("POLYGON")

quad.sf.list[[2]] <- quad.2.df %>% 
  
  st_as_sf(coords = c("x", "y"), crs = 32611) %>% 
  
  summarise(geometry = st_combine(geometry)) %>%
  
  st_cast("POLYGON")

quad.sf.list[[3]] <- quad.3.df %>% 
  
  st_as_sf(coords = c("x", "y"), crs = 32611) %>% 
  
  summarise(geometry = st_combine(geometry)) %>%
  
  st_cast("POLYGON")

quad.sf.list[[4]] <- quad.4.df %>% 
  
  st_as_sf(coords = c("x", "y"), crs = 32611) %>% 
  
  summarise(geometry = st_combine(geometry)) %>%
  
  st_cast("POLYGON")

# plot
plot(st_geometry(unit.bound.sf))
plot(st_geometry(quad.sf.list[[1]]), add = T)
plot(st_geometry(quad.sf.list[[2]]), add = T)
plot(st_geometry(quad.sf.list[[3]]), add = T)
plot(st_geometry(quad.sf.list[[4]]), add = T)

# determine which quadrants to treat
to.treat <- data.frame(rep = 1:3,
                       first.quad = sample(1:4, size = 3)) 

to.treat$second.quad <- case_when(to.treat$first.quad == 1 ~ 4,
                                  to.treat$first.quad == 2 ~ 3,
                                  to.treat$first.quad == 3 ~ 2,
                                  to.treat$first.quad == 4 ~ 1)
  
# write a function that does it
alter_cover <- function(landscape,
                        rep,
                        percent) {
  
  # focal landscape
  landscape.1 <- landscape
  
  # crop correct quadrants
  first.quad <- to.treat$first.quad[to.treat$rep == rep]
  second.quad <- to.treat$second.quad[to.treat$rep == rep]
  
  first.quad.crop <- crop(landscape.1, quad.sf.list[[first.quad]])
  second.quad.crop <- crop(landscape.1, quad.sf.list[[second.quad]])
  
  # how many can we flip?
  n.0.first <- length(values(first.quad.crop)[values(first.quad.crop) == 0])
  n.0.second <- length(values(second.quad.crop)[values(second.quad.crop) == 0])
  
  # randomly sample a percentage of them to flip
  to.flip.first <- sample(which(values(first.quad.crop) == 0), size = round(n.0.first * percent))
  to.flip.second <- sample(which(values(second.quad.crop) == 0), size = round(n.0.second * percent))
  
  # flip them
  values(first.quad.crop)[to.flip.first] <- 1
  values(second.quad.crop)[to.flip.second] <- 1
  
  # merge back in
  merge.1 <- merge(first.quad.crop, landscape.1)
  merge.2 <- merge(second.quad.crop, merge.1)
  
  # ensure that crs is correct
  crs(merge.2) <- crs("EPSG:32611")
  
  # return
  return(merge.2)
  
}

# use function
A.cov3.dis.crop <- c(alter_cover(B.cov3.dis.crop$B1, 1, 0.50),
                     alter_cover(B.cov3.dis.crop$B2, 2, 0.50),
                     alter_cover(B.cov3.dis.crop$B3, 3, 0.50))

plot(B.cov3.dis.crop)
plot(A.cov3.dis.crop)

# now to calculate percent open
# merge into landscape
A.cov3.dis <- merge(A.cov3.dis.crop, B.cov3.dis)

A.cov3 <- perc_open(A.cov3.dis, 50)

A.cov3.crop <- crop(A.cov3, unit.bound.sf)

# compare
plot(B.cov3.crop)
plot(A.cov3.crop)

#_______________________________________________________________________
# 6. Bind each landscape / replicate together ----
#_______________________________________________________________________
# 6a. Before ----
#_______________________________________________________________________

# bind
B1.ls <- c(B.cov1$B1, B.cov2$B1, B.cov3$B1)
B2.ls <- c(B.cov1$B2, B.cov2$B2, B.cov3$B2)
B3.ls <- c(B.cov1$B3, B.cov2$B3, B.cov3$B3)

# rename
names(B1.ls) <- c("fora", "elev", "open")
names(B2.ls) <- c("fora", "elev", "open")
names(B3.ls) <- c("fora", "elev", "open")

#_______________________________________________________________________
# 6b. After ----
#_______________________________________________________________________

# bind
A1.ls <- c(B.cov1$B1, B.cov2$B1, A1.cov3)
A2.ls <- c(B.cov1$B2, B.cov2$B2, A2.cov3)
A3.ls <- c(B.cov1$B3, B.cov2$B3, A3.cov3)

# rename
names(A1.ls) <- c("fora", "elev", "open")
names(A2.ls) <- c("fora", "elev", "open")
names(A3.ls) <- c("fora", "elev", "open")

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
    geom_spatraster(data = ls$fora) +
    
    # viridis colors
    scale_fill_viridis_c() +
    
    # add unit boundary
    geom_sf(data = unit.bound.sf,
            fill = NA,
            color = "white") +
    
    # remove axis text
    theme(axis.text = element_blank(),
          legend.position = "none") -> fora
  
  # edge
  ggplot() +
    
    theme_bw() +
    
    # add in raster
    geom_spatraster(data = ls$elev) +
    
    # viridis colors
    scale_fill_viridis_c() +
    
    # add unit boundary
    geom_sf(data = unit.bound.sf,
            fill = NA,
            color = "white") +
    
    # remove axis text
    theme(axis.text = element_blank(),
          legend.position = "none") -> elev
  
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
  focal.plot <- plot_grid(fora, elev, open, nrow = 1)
  
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

# BEFORE
plot_grid(B1.plot, B2.plot, B3.plot, nrow = 3)

# AFTER
plot_grid(A1.plot, A2.plot, A3.plot, nrow = 3)

#_______________________________________________________________________
# 8. Save and write rasters ----
#_______________________________________________________________________

save.image("Progress/rasters_02_19_2025.RData")

writeRaster(B1.ls, filename = "Rasters/B1.tif", overwrite = T)
writeRaster(B2.ls, filename = "Rasters/B2.tif", overwrite = T)
writeRaster(B3.ls, filename = "Rasters/B3.tif", overwrite = T)

writeRaster(A1.ls, filename = "Rasters/A1.tif", overwrite = T)
writeRaster(A2.ls, filename = "Rasters/A2.tif", overwrite = T)
writeRaster(A3.ls, filename = "Rasters/A3.tif", overwrite = T)

