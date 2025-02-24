# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 06 - Sample individuals for detection and movement modeling
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 26 Nov 2024
# Date completed: 09 Dec 2024
# Date last modified: 24 Feb 2025
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation
library(terra)
library(sf)

#_______________________________________________________________________
# 2. Read in naive RSF rasters and home range centers ----
#_______________________________________________________________________

# rasters 
B.rsf <- rast(paste0(getwd(), "/Rasters/B_rsf.tif"))
A.rsf <- rast(paste0(getwd(), "/Rasters/A_rsf.tif"))

# shapefiles
hrc.1 <- st_read(paste0(getwd(), "/Derived_data/Shapefiles/B1_hrc.shp"))
hrc.2 <- st_read(paste0(getwd(), "/Derived_data/Shapefiles/B2_hrc.shp"))
hrc.3 <- st_read(paste0(getwd(), "/Derived_data/Shapefiles/B3_hrc.shp"))

#_______________________________________________________________________
# 3. Sample individuals for before and after ----
#_______________________________________________________________________

# here, we need the same individuals to be present before and after treatment
# we have 100 individuals on each landscape, so let's just subsample iteratively

abundances <- c(100, 75, 50, 25, 15, 10, 5, 2)

#_______________________________________________________________________
# 3a. Extract "before" RSF score by HRC ----
#_______________________________________________________________________

hrc.1$rsf <- extract(B.rsf$B1, hrc.1)[ , 2]
hrc.2$rsf <- extract(B.rsf$B2, hrc.2)[ , 2]
hrc.3$rsf <- extract(B.rsf$B3, hrc.3)[ , 2]

#_______________________________________________________________________
# 3a. Function to subsample ----
#_______________________________________________________________________

subsample_indivs <- function(hrc) {
  
  samp.75 <- sample(1:100, size = 75, replace = FALSE, prob = hrc$rsf)
  samp.50 <- sample(samp.75, size = 50, replace = FALSE, prob = hrc[samp.75, ]$rsf)
  samp.25 <- sample(samp.50, size = 25, replace = FALSE, prob = hrc[samp.50, ]$rsf)
  samp.15 <- sample(samp.25, size = 15, replace = FALSE, prob = hrc[samp.25, ]$rsf)
  samp.10 <- sample(samp.15, size = 10, replace = FALSE, prob = hrc[samp.15, ]$rsf)
  samp.05 <- sample(samp.10, size = 5, replace = FALSE, prob = hrc[samp.10, ]$rsf)
  samp.02 <- sample(samp.05, size = 2, replace = FALSE, prob = hrc[samp.05, ]$rsf)
  
  # list
  sampled.indivs <- list(samp.75,
                         samp.50,
                         samp.25,
                         samp.15,
                         samp.10,
                         samp.05,
                         samp.02)
  
  # return
  return(sampled.indivs)
  
}

#_______________________________________________________________________
# 2b. Use function ----
#_______________________________________________________________________

sampled.indivs.1 <- subsample_indivs(hrc.1)
sampled.indivs.2 <- subsample_indivs(hrc.2)
sampled.indivs.3 <- subsample_indivs(hrc.3)

#_______________________________________________________________________
# 3. Sample "collared" individuals for RSF/iSSF modeling  ----
#_______________________________________________________________________

# these will be the 10 individuals sampled above, for simplicity
# they occur in all of the higher densities, and have individuals in the two lower densities
collared.1 <- sampled.indivs.1[[5]]
collared.2 <- sampled.indivs.2[[5]]
collared.3 <- sampled.indivs.3[[5]]

#_______________________________________________________________________
# 5. Write to files ----
#_______________________________________________________________________

# detection
save(sampled.indivs.1, file = paste0(getwd(), "/Derived_data/Lookup/detection_lookup_1.RData"))
save(sampled.indivs.2, file = paste0(getwd(), "/Derived_data/Lookup/detection_lookup_2.RData"))
save(sampled.indivs.3, file = paste0(getwd(), "/Derived_data/Lookup/detection_lookup_3.RData"))

# collared
save(collared.1, file = paste0(getwd(), "/Derived_data/Lookup/collared_lookup_1.RData"))
save(collared.2, file = paste0(getwd(), "/Derived_data/Lookup/collared_lookup_2.RData"))
save(collared.3, file = paste0(getwd(), "/Derived_data/Lookup/collared_lookup_3.RData"))
