# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 01d - Fit models to hare data
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 07 Apr 2025
# Date completed: 08 Apr 2025
# Date last modified: 08 Apr 2025
# R version: 4.4.3

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation
library(lubridate)            # work with dates
library(ctmm)                 # CTSP movement modeling
library(sf)                   # spatial operations
library(cowplot)              # multiple plots

#_______________________________________________________________________
# 2. Directories ----
#_______________________________________________________________________

dir.pre <- "E:/Hare project/Data analysis/GPS processing/Derived data/Cleaned data 2/PRE/"
dir.dur <- "E:/Hare project/Data analysis/GPS processing/Derived data/Cleaned data 2/DUR/"
dir.post <- "E:/Hare project/Data analysis/GPS processing/Derived data/Cleaned data 2/POST/"

# file lists
list.pre <- list.files(dir.pre)
list.dur <- list.files(dir.dur)
list.post <- list.files(dir.post)

#_______________________________________________________________________
# 3. Read in data ----
#_______________________________________________________________________

# error model
load("E:/Hare project/Data analysis/GPS processing/Derived data/error_model.RData")

# read in .csv
(focal.csv <- list.pre[19])

hare.data <- read.csv(paste0(dir.pre,
                             focal.csv))        

#_______________________________________________________________________
# 4. Convert to telemetry object ----
#_______________________________________________________________________

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
uere(hare.telem) <- best.uere.HDOP.class

#_______________________________________________________________________
# 5. CTMM model selection ----
#_______________________________________________________________________

# account for irregular sampling
samp.sched <- c(2, 12) %#% "hour"

# plot periodogram
plot(periodogram(hare.telem))

# plot variogram
plot(variogram(hare.telem,
               dt  = samp.sched))

# guesstimated model parameters from the variogram
guess.param <- variogram.fit(variogram(hare.telem,
                                       dt  = samp.sched), 
                             name = "guess.param", 
                             interactive = FALSE)

# model selection
fitted.mods <- ctmm.select(hare.telem, 
                           CTMM = guess.param, 
                           verbose = TRUE)

summary(fitted.mods)

top.model <- fitted.mods[[1]]

# plot diagnostics
plot(variogram(hare.telem,
               dt  = samp.sched),
     CTMM = top.model)

# parameters
top.model$tau
top.model$sigma

summary(top.model)

# save to file
save(top.model, file = "Derived data/Hares - Fitted models/PRE_019.RData")
