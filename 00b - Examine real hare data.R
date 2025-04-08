# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 00b - Examine real hare data
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 04 Apr 2025 
# Date completed:  
# Date last modified: 
# R version: 4.4.3

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation
library(lubridate)            # work with dates
library(ctmm)                 # CTSP movement modeling
library(amt)                  # work with animal movement tracks
library(sf)                   # spatial operations
library(cowplot)              # multiple plots

#_______________________________________________________________________
# 2. Directories ----
#_______________________________________________________________________

dir.pre <- "E:/Hare project/Data analysis/GPS processing/Derived data/Cleaned data 2/PRE/"

# file lists
list.pre <- list.files(dir.pre)

#_______________________________________________________________________
# 3. Read in data ----
#_______________________________________________________________________

# error model
load("E:/Hare project/Data analysis/GPS processing/Derived data/error_model.RData")

# tracks to import
# PRE - 002
# PRE - 010
# PRE - 011
# PRE - 019

# read in .csvs 
focal.csv.1 <- list.pre[2]
focal.csv.2 <- list.pre[10]
focal.csv.3 <- list.pre[11]
focal.csv.4 <- list.pre[19]

hare.data.1 <- read.csv(paste0(dir.pre, focal.csv.1))
hare.data.2 <- read.csv(paste0(dir.pre, focal.csv.2)) 
hare.data.3 <- read.csv(paste0(dir.pre, focal.csv.3)) 
hare.data.4 <- read.csv(paste0(dir.pre, focal.csv.4)) 

#_______________________________________________________________________
# 4. Convert to telemetry objects ----
#_______________________________________________________________________

convert_telem <- function(hare.data) {
  
  # convert to Movebank format
hare.movebank <- data.frame("timestamp" = hare.data$timestamp,
                            "location.lat" = hare.data$location.lat,
                            "location.long" = hare.data$location.long,
                            "height above mean sea level" = hare.data$height.above.mean.sea.level,
                            "GPS satellite count" = hare.data$GPS.satellite.count,
                            "GPS HDOP" = hare.data$GPS.HDOP)

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
uere(hare.telem) <- uere.HDOP.class

# return
return(hare.telem)
  
}

hare.telem.1 <- convert_telem(hare.data.1)
hare.telem.2 <- convert_telem(hare.data.2)
hare.telem.3 <- convert_telem(hare.data.3)
hare.telem.4 <- convert_telem(hare.data.4)

#_______________________________________________________________________
# 5. CTMM model selection ----
#_______________________________________________________________________

hare.telem <- hare.telem.4
  
# account for irregular sampling
samp.sched <- c(2, 12) %#% "hour"

# plot variogram
plot(variogram(hare.telem,
               dt  = samp.sched))

# guesstimated model parameters from the variogram
guess.param <- variogram.fit(variogram(hare.telem,
                                       dt  = samp.sched), 
                             name = "guess.param", 
                             interactive = FALSE)

# turn error on
guess.param$error <- TRUE

# model selection
fitted.mods <- ctmm.select(hare.telem, 
                           CTMM = guess.param, 
                           verbose = TRUE)

summary(fitted.mods)

top.model <- fitted.mods[[2]]

# parameters
top.model$tau
top.model$sigma

summary(top.model)

ctmm::speed(top.model)
ctmm::speed(object = hare.telem, CTMM = top.model, fast = T)
  
# parameters to use
# tau1: 15550, tau2: 15550, sigma major: 49744, sigma minor: 298

# speed to recover: 0.998 km/day

# what does a simulated track look like?

# we'll look at track differences between 1 and 0.5 hours
ctmm.ouf.real <- ctmm(tau = c(15550, 
                              15550), 
                      isotropic = FALSE, 
                      sigma = c(49744,
                                298,
                                0), 
                      mu = c(0, 0))

# simulate
sim.ouf.real <- simulate(object = ctmm.ouf.real, 
                         t = seq(1, 4 * 7 * 24 * 60 * 60, 60), 
                         complete = TRUE)


# plot
ggplot(data = sim.ouf.real) + 
  
  theme_bw() +
  
  geom_path(aes(x = x,
                y = y,
                color = timestamp)) +
  
  scale_color_viridis_c() +
  
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank())


