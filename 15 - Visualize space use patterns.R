# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 15 - Visualize space use patterns
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 02 Aug 2025
# Date completed: 03 Aug 2025
# Date last modified: 02 Aug 2025
# R version: 4.4.3

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation
library(lubridate)            # work with dates
library(ctmm)                 # CTSP movement modeling
library(terra)
library(sf)
library(tidyterra)
library(cowplot)

#_______________________________________________________________________
# 2. Read in data ----

# error model
load("E:/Hare project/Data analysis/GPS processing/Derived data/error_model.RData")

#_______________________________________________________________________
# 2a. Location files ---- 
#_______________________________________________________________________

# directories
dir.pre <- "E:/Hare project/Data analysis/GPS processing/Derived data/Cleaned data 2/PRE/"
dir.dur <- "E:/Hare project/Data analysis/GPS processing/Derived data/Cleaned data 2/DUR/"
dir.post <- "E:/Hare project/Data analysis/GPS processing/Derived data/Cleaned data 2/POST/"

# location file list
location.names <- data.frame("file" = c("001_1863_1.csv", "002_528_1.csv",
                                        "006_1737_1.csv", "008_553_1.csv",
                                        "016_1730_1.csv", "020_456_1.csv",
                                        "004_1885_2.csv", "009_801_1.csv",
                                        "013_1020_1.csv", "014_432_2.csv",
                                        "017_986_2.csv", "021_1253_1.csv",
                                        "023_432_3.csv", "026_1766_1.csv",
                                        "034_1631_1.csv", "002_1876_1.csv",
                                        "010_1841_1.csv", "011_1842_1.csv",
                                        "019_1600_1.csv", "021_1897_1.csv",
                                        "039_446_1.csv", "043_803_1.csv",
                                        "044_1051_1.csv", "055_510_1.csv"),
                            "Group" = c("DUR", "DUR", "DUR", "DUR",
                                        "DUR", "DUR", "POST", "POST",
                                        "POST", "POST", "POST", "POST",
                                        "POST", "POST", "POST", "PRE",
                                        "PRE", "PRE", "PRE", "PRE",
                                        "PRE", "PRE", "PRE", "PRE"))
# loop
# blank list
all.telem <- list()

for (i in 1:nrow(location.names)) {
  
  # which location file
  focal.location <- location.names[i, ]
  
  # read.csv
  if (focal.location$Group == "PRE") {
    
    focal.data <- read.csv(paste0(dir.pre,
                                  focal.location$file))
    
  } else {
    
    if (focal.location$Group == "DUR") {
      
      focal.data <- read.csv(paste0(dir.dur,
                                    focal.location$file))
    
    } else {
    
      if (focal.location$Group == "POST") {
        
        focal.data <- read.csv(paste0(dir.post,
                                      focal.location$file)) 
      
      }
      
    }
    
  }
  
  # convert to Movebank format
  hare.movebank <- data.frame("timestamp" = focal.data$timestamp,
                              "location.lat" = focal.data$location.lat,
                              "location.long" = focal.data$location.long,
                              "height above mean sea level" = focal.data$height.above.mean.sea.level,
                              "GPS satellite count" = focal.data$GPS.satellite.count,
                              "GPS HDOP" = focal.data$GPS.HDOP)
  
  # convert to telemetry object
  hare.telem <- as.telemetry(object = hare.movebank,
                             timeformat = "auto",
                             timezone = "America/Los_Angeles",
                             keep = TRUE)
  
  # add a "class" variable for error model
  hare.telem$class <- as.factor(ifelse(hare.telem$GPS.satellite.count > 3,
                                       "3D",
                                       "2D"))
  
  # add in error model
  uere(hare.telem) <- best.uere.HDOP.class
  
  # bind into list
  all.telem[[i]] <- hare.telem
  
}
       
#_______________________________________________________________________
# 2b. Model files ---- 
#_______________________________________________________________________

all.files <- list.files("Derived data/Hares - Fitted models")

all.models <- list()

for (i in 1:length(all.files)) {
  
  load(paste0("Derived data/Hares - Fitted models/", all.files[i]))
  
  # bind into list
  all.models[[i]] <- top.model
  
  rm(top.model)
  
}

names(all.models) <- all.files

#_______________________________________________________________________
# 2c. PCT unit boundaries ---- 
#_______________________________________________________________________

pct.units <- st_read("E:/Hare project/Spatial data/Units/final_units_ground_poly/final_units_ground_poly.shp")

pct.units.utm <- st_transform(pct.units, crs = "epsg:32611")

#_______________________________________________________________________
# 3. Write AKDE visualization function ----
#_______________________________________________________________________

akde_viz <- function (index,
                      buffer = 350) {
  
  # telem file
  focal.telem <- all.telem[[index]]
  
  # model file
  focal.model <- all.models[[index]]
  
  # fit AKDE
  focal.akde <- akde(data = focal.telem,
                     CTMM = focal.model,
                     debias = TRUE,
                     weights = TRUE,
                     res = 25)
  
  # export as raster and convert to rast
  focal.rast <- rast(raster(focal.akde,
                            DF = "PDF"))
  
  # project to UTM
  focal.rast.utm <- project(focal.rast,
                            "epsg:32611")
  
  # export 95% contour (point estimate), transform
  focal.95 <- st_transform(as.sf(focal.akde,
                           level.UD = 0.95)[2, ],
                           "epsg:32611")
  
  # extract extent coordinates and centroid
  focal.extent <- ext(focal.rast.utm)
  focal.centroid <- as.numeric(unlist(st_centroid(focal.95))[2:3])
  
  # crop and mask raster
  focal.rast.mask <- mask(crop(focal.rast.utm,
                               focal.95),
                          focal.95)
  
  # ggplot
  ggplot() +
    
    theme_bw() +
  
    # PCT units
    geom_sf(data = pct.units.utm,
            fill = "black",
            color = "black",
            alpha = 0.45) +
    
    geom_spatraster(data = focal.rast.mask,
                    alpha = 0.85) +
    
    geom_sf(data = focal.95,
            fill = NA,
            color = "black",
            linewidth = 0.75) +
    
    coord_sf(datum = "epsg:32611",
             xlim = c(focal.centroid[1] - buffer,
                      focal.centroid[1] + buffer),
             ylim = c(focal.centroid[2] - buffer,
                      focal.centroid[2] + buffer)) +
    
    scale_fill_viridis_c(option = "plasma",
                         na.value = NA) +
    
    ggtitle(paste0(location.names$Group[index],
                   "_",
                   substr(location.names$file[index], 1, 3))) +
    
    theme(legend.position = "none",
          panel.grid = element_line(color = "gray"),
          axis.text = element_blank(),
          plot.title = element_text(size = 8),
          axis.ticks = element_blank()) -> focal.plot
    
  # return
  return(focal.plot)
    
}

#_______________________________________________________________________
# 4. Apply function to create 6 x 4 grid ----
#_______________________________________________________________________

plot_grid(plotlist = list(akde_viz(1, 200), akde_viz(2, 250), akde_viz(3, 400), akde_viz(4), akde_viz(5, 200), akde_viz(6, 300),
                          akde_viz(7, 250), akde_viz(8, 150), akde_viz(9, 300), akde_viz(10), akde_viz(11, 200), akde_viz(12, 250),
                          akde_viz(13, 550), akde_viz(14), akde_viz(15), akde_viz(16, 375), akde_viz(17), akde_viz(18, 300),
                          akde_viz(19, 200), akde_viz(20, 250), akde_viz(21), akde_viz(22), akde_viz(23, 375), akde_viz(24, 300)),
          ncol = 6,
          nrow = 4)

# 827 x 627