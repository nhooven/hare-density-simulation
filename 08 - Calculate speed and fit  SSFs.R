# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 08 - Calculate speed and fit SSFs
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 26 Nov 2024
# Date completed: 
# Date last modified: 
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation
library(ctmm)
library(amt)
library(glmmTMB)

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

# lookup table
#lookup.movement <- read.csv()

# simulated tracks
sims.simple.weak <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_simple_weak.csv"))
#sims.simple.weak <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_simple_weak.csv"))
#sims.simple.weak <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_simple_weak.csv"))
#sims.simple.weak <- read.csv(paste0(getwd(), "/Derived_data/Simulated data/sims_simple_weak.csv"))

# CTMMs (if we assume no fix acquisition failure)
load(paste0(getwd(), "/Derived_data/CTMMs/ctmms_simple_weak.RData"))
ctmm.simple.weak <- all.ctmms
rm(all.ctmms)

# and so forth

#_______________________________________________________________________
# 3. Filter data (sampled individuals from script 06) ----
#_______________________________________________________________________
# 3a. Simple / weak ----
#_______________________________________________________________________

# pull sampled individuals from lookup
sampled.simple.weak <- lookup.movement %>%
  
  filter(landscape == "simple",
         hs = "weak") %>%
  
  dplyr::select(indiv) %>%
  
  # pull out as a vector (neat!)
  pull(indiv)

# subset
sims.simple.weak.sub <- sims.simple.weak %>%
  
  filter(indiv %in% sampled.simple.weak)

#_______________________________________________________________________
# 4. Subset CTMMs and calculate speeds ----
#_______________________________________________________________________

ctmm.sampled.simple.weak <- ctmm.simple.weak[sampled.simple.weak]

# loop through all
for (i in 1:20) {
  
  # subset
  focal.ctmm <- ctmm.sampled.simple.weak[[i]]
  
  ctmm::speed(focal.ctmm,
              data = focal.ctmm$data)
  
  
}
