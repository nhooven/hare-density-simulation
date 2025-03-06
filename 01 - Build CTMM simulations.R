# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 01 - Build CTMM simulations
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 06 Mar 2025
# Date completed: 
# Date last modified: 
# R version: 4.4.3

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)
library(ctmm)
library(sf)

#_______________________________________________________________________
# 2. Orientation ----

# this is essentially Table 1 from Calabrese et al. 2016

#_______________________________________________________________________

# MODEL    Pos    Vel    RR    Tau
# ________________________________
# BM       Y      N      N     Inf
# OU       Y      N      Y     Tr
# IOU      Y      Y      N     {Inf, Tv}
# OUF      Y      Y      Y     {Tr, Tv}

#_______________________________________________________________________
# 3. Brownian motion model ----

# BM assumes an infinitely diffusing process with no velocity autocorrelation
# or restricted space use

# it would be nice to see how this model generates data initially

# is this what the REM is assuming? Not really because the animal must be 
# "resident" within the sampling grid of the cameras

#_______________________________________________________________________

# initialize a CTMM object for simulation
ctmm.bm <- ctmm(
  
  tau = Inf,            # tau is infinite, no range residency
  sigma = 50,            # Brownian motion variance
  mu = c(0, 0)          # starting location
  
  )

# time
t.bm <- 1:1000

# simulate
sim.bm <- simulate(object = ctmm.bm,
                   t = t.bm)

# plot
ggplot(sim.bm,
       aes(x = x,
           y = y)) +
  
  theme_bw() +
  
  geom_path() +
  
  geom_point()
