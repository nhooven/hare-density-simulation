# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 02 - Initialize simulation stage
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 25 Mar 2025
# Date completed: 25 Mar 2025
# Date last modified: 02 Apr 2025
# R version: 4.4.3

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation
library(sf)                   # spatial operations
library(secr)                 # convenient camera array function

#_______________________________________________________________________
# 2. Define unit boundary ----

# this will be a 10-ha area centered in the middle of the landscape (i.e., c(0, 0))

#_______________________________________________________________________

# centroid
unit.centroid <- c(0, 0)

# how many meters per side? (let's make our unit 10 ha)
m.side <- sqrt(10 * 10000)

# how many meters to add and subtract?
m.side.half <- m.side / 2

# create polygon feature
unit.bound <- st_polygon(list(cbind(c(unit.centroid[1] - m.side.half, 
                                      unit.centroid[1] + m.side.half,
                                      unit.centroid[1] + m.side.half,
                                      unit.centroid[1] - m.side.half,
                                      unit.centroid[1] - m.side.half), 
                                    c(unit.centroid[2] - m.side.half, 
                                      unit.centroid[2] - m.side.half,
                                      unit.centroid[2] + m.side.half,
                                      unit.centroid[2] + m.side.half,
                                      unit.centroid[2] - m.side.half))))

# assign to sf object
unit.bound.sf <- st_as_sf(st_sfc(unit.bound))

# plot 
plot(unit.bound.sf) 

# write to shapefile
st_write(unit.bound.sf,
         dsn = paste0(getwd(), "/Derived data/Shapefiles/unit_bound.shp"),
         layer = "unit_bound.shp",
         append = FALSE)

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

# let's make this a true "pie slice" rather than an approximated triangle

# https://stackoverflow.com/questions/59328707/how-do-i-partition-a-circle-into-equal-multipolygon-slices-with-sf-and-r

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
    
    # theta (sequence of angles)
    theta <- seq(-lens.angle / 2, lens.angle / 2, length.out = 50)
    
    # x coords of the arc
    xarc <- focal.coords[1] + max.dist * sin(theta)
    
    # y coords of the arc
    yarc <- focal.coords[2] + max.dist * cos(theta)
    
    # final coords
    xc <- c(focal.coords[1], xarc, focal.coords[1])
    yc <- c(focal.coords[2], yarc, focal.coords[2])
    
    pie.df <- data.frame(x = xc,
                         y = yc)
      
    # convert to sf
    pie.sf <- pie.df %>% 
      
      st_as_sf(coords = c("x",
                          "y")) %>%
      
      summarize(geometry = st_combine(geometry)) %>%
      
      st_cast("POLYGON")
    
    # bind into df
    all.polys <- rbind(all.polys, pie.sf)
    
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

#_______________________________________________________________________
# 5b. 9 cameras ----
#_______________________________________________________________________

# equal spacing between each camera and the sampling area edge
cams.9.grid <- make.grid(nx = 3,
                         ny = 3,
                         spacing = sqrt(100000) / 4, 
                         originxy = c(unit.centroid[1] - dist.adjust + (sqrt(100000) / 4), 
                                      unit.centroid[2] - dist.adjust + (sqrt(100000) / 4)))

cams.9 <- st_as_sf(cams.9.grid,
                   coords = c("x", "y"))

plot(st_geometry(unit.bound.sf))
plot(st_geometry(cams.9), add = T)

#_______________________________________________________________________
# 6. Make viewshed polygons ----
#_______________________________________________________________________

cams.viewsheds <- make_viewshed(cams.9)

#_______________________________________________________________________
# 7. Plot them ----
#_______________________________________________________________________

plot(st_geometry(cams.viewsheds))

# area for sanity check
st_area(cams.viewsheds)

(pi * 3.5^2) * (57.3 / 360)

#_______________________________________________________________________
# 8. Write to shapefile ----
#_______________________________________________________________________

st_write(cams.viewsheds, 
         dsn = paste0(getwd(), "/Derived data/Shapefiles/cams_vs.shp"),
         append = FALSE)

