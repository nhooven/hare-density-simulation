# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 04 - Set up virtual camera arrays
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 22 Nov 2024
# Date completed: 22 Nov 2024
# Date last modified: 20 Dec 2024
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation
library(sf)              # spatial data
library(secr)            # convenient camera array function

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

unit.bound <- st_read(paste0(getwd(), "/Derived_data/Shapefiles/unit_bound.shp"))

#_______________________________________________________________________
# 3. Define viewshed parameters ----
#_______________________________________________________________________

# (assume all cameras face north) create triangle
# max distance - 3.5 m
max.dist <- 3.5

# lens angle (degrees)
lens.angle.deg <- 57.3

#_______________________________________________________________________
# 4. Write function to generate viewsheds ----
#_______________________________________________________________________

make_viewshed <- function(cams) {
  
  # loop through all points
  all.polys <- data.frame()
  
  n.cams <- nrow(cams)
  
  for (i in 1:n.cams) {
    
    # subset
    focal.point <- cams %>% slice(i)
    
    # extract coordinates 
    focal.coords <- st_coordinates(focal.point)
    
    # calculate lens angle in rad
    lens.angle <- lens.angle.deg * (pi / 180)
    
    # y axis adjustment
    y.adjust <- max.dist
    
    # x axis adjustment
    x.adjust <- y.adjust * tan(lens.angle / 2)
    
    # create triangle vertices
    vert.w <- c(focal.coords[1] - x.adjust,
                focal.coords[2] + y.adjust)
    
    vert.e <- c(focal.coords[1] + x.adjust,
                focal.coords[2] + y.adjust)
    
    # bind vertices together
    tri.verts <- as.data.frame(matrix(c(focal.coords, 
                                        vert.w, 
                                        vert.e),
                                      nrow = 3,
                                      ncol = 2,
                                      byrow = TRUE))
    
    # convert to sf polygon
    tri.polygon <- tri.verts %>%
      
      st_as_sf(coords = c("V1", 
                          "V2"), 
               crs = "epsg:32611") %>%
      
      summarise(geometry = st_combine(geometry)) %>%
      
      st_cast("POLYGON")
    
    # bind into df
    all.polys <- rbind(all.polys, tri.polygon)
    
  }
  
  # return
  return(all.polys)
  
}

#_______________________________________________________________________
# 5. Sample camera locations ----
#_______________________________________________________________________
# 5a. Define sampling area ----
#_______________________________________________________________________

# define 10-ha sampling area
# must be 100,000 sq m
# vertical and horizontal distance from centroid
dist.adjust <- sqrt(100000) / 2

centroid <- st_coordinates(st_centroid(unit.bound))

#_______________________________________________________________________
# 5a. 4 cameras ----
#_______________________________________________________________________

# equal spacing between each camera and the sampling area edge
cams.4.grid <- make.grid(nx = 2,
                         ny = 2,
                         spacing = sqrt(100000) / 3, 
                         originxy = c(centroid[1] - dist.adjust + (sqrt(100000) / 3), 
                                      centroid[2] - dist.adjust + (sqrt(100000) / 3)))

cams.4 <- st_as_sf(cams.4.grid,
                   coords = c("x", "y"),
                   crs = "epsg:32611")

plot(st_geometry(unit.bound))
plot(st_geometry(cams.4), add = T)

#_______________________________________________________________________
# 5b. 9 cameras ----
#_______________________________________________________________________

# equal spacing between each camera and the sampling area edge
cams.9.grid <- make.grid(nx = 3,
                         ny = 3,
                         spacing = sqrt(100000) / 4, 
                         originxy = c(centroid[1] - dist.adjust + (sqrt(100000) / 4), 
                                      centroid[2] - dist.adjust + (sqrt(100000) / 4)))

cams.9 <- st_as_sf(cams.9.grid,
                   coords = c("x", "y"),
                   crs = "epsg:32611")

plot(st_geometry(unit.bound))
plot(st_geometry(cams.9), add = T)

#_______________________________________________________________________
# 5c. 16 cameras ----
#_______________________________________________________________________

# equal spacing between each camera and the sampling area edge
cams.16.grid <- make.grid(nx = 4,
                          ny = 4,
                          spacing = sqrt(100000) / 5, 
                          originxy = c(centroid[1] - dist.adjust + (sqrt(100000) / 5), 
                                       centroid[2] - dist.adjust + (sqrt(100000) / 5)))

cams.16 <- st_as_sf(cams.16.grid,
                    coords = c("x", "y"),
                    crs = "epsg:32611")

plot(st_geometry(unit.bound))
plot(st_geometry(cams.16), add = T)

#_______________________________________________________________________
# 6. Make viewshed polygons ----
#_______________________________________________________________________

cams.4.viewsheds <- make_viewshed(cams.4)
cams.9.viewsheds <- make_viewshed(cams.9)
cams.16.viewsheds <- make_viewshed(cams.16)

#_______________________________________________________________________
# 7. Write to shapefiles ----
#_______________________________________________________________________

st_write(cams.4.viewsheds, 
         dsn = paste0(getwd(), "/Derived_data/Shapefiles/cams_4_vs.shp"))

st_write(cams.9.viewsheds, 
         dsn = paste0(getwd(), "/Derived_data/Shapefiles/cams_9_vs.shp"))

st_write(cams.16.viewsheds, 
         dsn = paste0(getwd(), "/Derived_data/Shapefiles/cams_16_vs.shp"))
