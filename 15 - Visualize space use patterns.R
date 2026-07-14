# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 15 - Visualize space use patterns
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 02 Aug 2025
# Date completed: 03 Aug 2025
# Date last modified: 14 Aug 2026
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
# 2. Read in models ----
#_______________________________________________________________________

# directory
dir.model <- "D:/hare_project/data_analysis/General/hare-gps-processing-new/data_cleaned/"

# read in
# AKDEs
all.akde <- readRDS(paste0(dir.model, "all_akde.rds"))

# model selection results
model.select <- readRDS(paste0(dir.model, "all_model_select.rds"))

# shift rownames to a column
model.select$model <- rownames(model.select) |>
  
  # remove digits
  gsub('[[:digit:]]+', '', x = _)

model.select.top <- model.select |> filter(mod == 1) |>
  
  mutate(model = case_when(
    
    model %in% c("IID", "IID anisotropic") ~ "IID",
    model %in% c("OU", "OU anisotropic") ~ "OU",
    model %in% c("OUF", "OUF anisotropic") ~ "OUF",
    model %in% c("OUf", "OUf anisotropic") ~ "OUf"
    
  )
  
  )

model.select.tau.v <- model.select.top |> filter(model %in% c("OUF", "OUf"))

# indices
indices.tau.v <- model.select.tau.v$i

# subset top model lists
models.tau.v <- top.models[indices.tau.v]

#_______________________________________________________________________
# 2c. PCT unit boundaries ---- 
#_______________________________________________________________________

pct.units <- st_read("D:/hare_project/data_spatial/Units/final_units_ground_poly/final_units_ground_poly.shp")

pct.units.utm <- st_transform(pct.units, crs = "epsg:32611")

#_______________________________________________________________________
# 3. Write AKDE visualization function ----
#_______________________________________________________________________

akde_viz <- function (index,
                      buffer = 350) {
  
  # AKDE
  focal.akde <- all.akde[[index]]
  
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
    
    ggtitle(model.select.tau.v$track_season_post[model.select.tau.v$i == index]) +
    
    theme(legend.position = "none",
          panel.grid = element_line(color = "gray"),
          axis.text = element_blank(),
          plot.title = element_text(size = 8),
          axis.ticks = element_blank()) -> focal.plot
    
  # return
  return(focal.plot)
    
}

#_______________________________________________________________________
# 4. Apply function to create 7 x 7 grid ----
#_______________________________________________________________________

plot_grid(plotlist = list(akde_viz(8, 250), akde_viz(12, 350), akde_viz(25, 350), akde_viz(30, 250), akde_viz(34, 400), akde_viz(36, 250), akde_viz(38, 350),
                          akde_viz(41, 250), akde_viz(43, 200), akde_viz(46, 250), akde_viz(47, 350), akde_viz(48, 350), akde_viz(52, 450), akde_viz(54, 350),
                          akde_viz(57, 350), akde_viz(58, 350), akde_viz(60, 350), akde_viz(63, 250), akde_viz(68, 250), akde_viz(72, 400), akde_viz(75, 200),
                          akde_viz(79, 250), akde_viz(81, 350), akde_viz(83, 250), akde_viz(87, 250), akde_viz(94, 250), akde_viz(104, 250), akde_viz(106, 150),
                          akde_viz(110, 350), akde_viz(113, 250), akde_viz(118, 250), akde_viz(120, 350), akde_viz(123, 250), akde_viz(124, 250), akde_viz(127, 150),
                          akde_viz(128, 250), akde_viz(132, 150), akde_viz(141, 250), akde_viz(143, 250), akde_viz(145, 250), akde_viz(156, 350), akde_viz(160, 250),
                          akde_viz(168, 250), akde_viz(170, 350), akde_viz(171, 450), akde_viz(172, 475), akde_viz(181, 150), akde_viz(184, 350), akde_viz(200, 150)),
          ncol = 7,
          nrow = 7)

# 900 x 900