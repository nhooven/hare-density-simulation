# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 03 - Sample movement parameters
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 25 Mar 2025
# Date completed: 25 Mar 2025
# Date last modified: 03 Apr 2025
# R version: 4.4.3

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation
library(sf)                   # spatial operations
library(mefa4)
library(ctmm)                 # helper function for unit conversion

#_______________________________________________________________________
# 2. Define numbers for simulations ----
#_______________________________________________________________________

# number of iterations
n.iter <- 500

# total individuals to simulate per iteration
# camera contacts for Q1
n.indiv.cam.Q1 <- 50

# additional out of sample collars for Q1
n.indiv.collar.Q1 <- 10

# high variability camera contacts for Q2
n.indiv.cam.Q2 <- 10

# high variability additional out of sample collars for Q2
n.indiv.collar.Q2 <- 10

#_______________________________________________________________________
# 3. Sample individuals ----
#_______________________________________________________________________
# 3a. Q1 camera ----
#_______________________________________________________________________

indivs.cam.Q1 <- expand.grid(iter = 1:n.iter,
                             indiv = 1:n.indiv.cam.Q1,
                             Q = "Q1",
                             use = "camera")

# add in folds
indivs.cam.Q1$fold <- case_when(indivs.cam.Q1$iter %in% 1:100 ~ 1,
                                indivs.cam.Q1$iter %in% 101:200 ~ 2,
                                indivs.cam.Q1$iter %in% 201:300 ~ 3,
                                indivs.cam.Q1$iter %in% 301:400 ~ 4,
                                indivs.cam.Q1$iter %in% 401:500 ~ 5)

#_______________________________________________________________________
# 3b. Q1 collar ----
#_______________________________________________________________________

indivs.collar.Q1 <- expand.grid(iter = 1:n.iter,
                                indiv = 1:n.indiv.collar.Q1,
                                Q = "Q1",
                                use = "collar")

# add in folds
indivs.collar.Q1$fold <- case_when(indivs.collar.Q1$iter %in% 1:100 ~ 1,
                                   indivs.collar.Q1$iter %in% 101:200 ~ 2,
                                   indivs.collar.Q1$iter %in% 201:300 ~ 3,
                                   indivs.collar.Q1$iter %in% 301:400 ~ 4,
                                   indivs.collar.Q1$iter %in% 401:500 ~ 5)

#_______________________________________________________________________
# 3c. Q2 camera ----
#_______________________________________________________________________

indivs.cam.Q2 <- expand.grid(iter = 1:n.iter,
                             indiv = 1:n.indiv.cam.Q2,
                             Q = "Q2",
                             use = "camera")

# add in folds
indivs.cam.Q2$fold <- case_when(indivs.cam.Q2$iter %in% 1:100 ~ 1,
                                indivs.cam.Q2$iter %in% 101:200 ~ 2,
                                indivs.cam.Q2$iter %in% 201:300 ~ 3,
                                indivs.cam.Q2$iter %in% 301:400 ~ 4,
                                indivs.cam.Q2$iter %in% 401:500 ~ 5)

#_______________________________________________________________________
# 3d. Q2 collar ----
#_______________________________________________________________________

indivs.collar.Q2 <- expand.grid(iter = 1:n.iter,
                                indiv = 1:n.indiv.collar.Q2,
                                Q = "Q2",
                                use = "collar")

# add in folds
indivs.collar.Q2$fold <- case_when(indivs.collar.Q2$iter %in% 1:100 ~ 1,
                                   indivs.collar.Q2$iter %in% 101:200 ~ 2,
                                   indivs.collar.Q2$iter %in% 201:300 ~ 3,
                                   indivs.collar.Q2$iter %in% 301:400 ~ 4,
                                   indivs.collar.Q2$iter %in% 401:500 ~ 5)

#_______________________________________________________________________
# 4. Define distributions of movement parameters ----

# log-normal parameter function
log_norm_params <- function (
    
  vals    # must be a length 2 vector with the desired mean and SD
  
) {
  
  # median = geometric mean
  lnorm.med <- log(vals[1])
  
  lnorm.mean <- log((vals[1]^2) / sqrt((vals[1]^2) + (vals[2]^2)))
  
  lnorm.sd <- sqrt(log(1 + ((vals[2]^2) / (vals[1]^2))))
  
  # bind together
  lnorm.params <- c(lnorm.med, lnorm.mean, lnorm.sd)
  
  # return
  return(lnorm.params)
  
}

#_______________________________________________________________________
# 4a. Home range centroids ----

# here, our definition of population density will be n home range centroids
# (activity centers, etc.) per a 10-ha area

# is this the most realistic? It's an assumption the SECR uses, and we can 
# model it pretty intuitively with the mean parameters in ctmm

# each individual will get a random uniform draw within the unit boundary, for each x and y

#_______________________________________________________________________

unit.bound <- st_read(paste0(getwd(), "/Derived data/Shapefiles/unit_bound.shp"))

# min and max x/y
coord.min <- as.numeric(st_bbox(unit.bound)[1])
coord.max <- as.numeric(st_bbox(unit.bound)[3])

# Q1
indivs.cam.Q1$mean1  <- runif(n = nrow(indivs.cam.Q1), coord.min, coord.max)
indivs.cam.Q1$mean2  <- runif(n = nrow(indivs.cam.Q1), coord.min, coord.max)

indivs.collar.Q1$mean1 <- runif(n = nrow(indivs.collar.Q1), coord.min, coord.max)
indivs.collar.Q1$mean2 <- runif(n = nrow(indivs.collar.Q1), coord.min, coord.max)

# Q2
indivs.cam.Q2$mean1  <- runif(n = nrow(indivs.cam.Q2), coord.min, coord.max)
indivs.cam.Q2$mean2  <- runif(n = nrow(indivs.cam.Q2), coord.min, coord.max)

indivs.collar.Q2$mean1 <- runif(n = nrow(indivs.collar.Q2), coord.min, coord.max)
indivs.collar.Q2$mean2 <- runif(n = nrow(indivs.collar.Q2), coord.min, coord.max)

#_______________________________________________________________________
# 4b. Tau 1 ----

# position autocorrelation / home range crossing time parameter

# geometric mean = 8 hr
# SDs = 2, 8 hr

lnorm.tau1.lo <- log_norm_params(c(8 %#% "hours", 2 %#% "hours"))
lnorm.tau1.hi <- log_norm_params(c(8 %#% "hours", 8 %#% "hours"))

#_______________________________________________________________________

# Q1
indivs.cam.Q1$tau1 <- rlnorm(nrow(indivs.cam.Q1), 
                             meanlog = lnorm.tau1.lo[1],
                             sdlog = lnorm.tau1.lo[3])

indivs.collar.Q1$tau1 <- rlnorm(nrow(indivs.collar.Q1), 
                                meanlog = lnorm.tau1.lo[1],
                                sdlog = lnorm.tau1.lo[3])

# Q2
indivs.cam.Q2$tau1 <- rlnorm(nrow(indivs.cam.Q2), 
                             meanlog = lnorm.tau1.lo[1],
                             sdlog = lnorm.tau1.lo[3])

indivs.collar.Q2$tau1 <- rlnorm(nrow(indivs.collar.Q2), 
                                meanlog = lnorm.tau1.lo[1],
                                sdlog = lnorm.tau1.lo[3])

#_______________________________________________________________________
# 4c. Tau 2 ----

# geometric mean = 1 hr
# SDs = 0.25, 1 hr

lnorm.tau2.lo <- log_norm_params(c(1 %#% "hours", 0.25 %#% "hours"))
lnorm.tau2.hi <- log_norm_params(c(1 %#% "hours", 1 %#% "hours"))

#_______________________________________________________________________

# Q1
indivs.cam.Q1$tau2 <- rlnorm(nrow(indivs.cam.Q1), 
                             meanlog = lnorm.tau2.lo[1],
                             sdlog = lnorm.tau2.lo[3])

indivs.collar.Q1$tau2 <- rlnorm(nrow(indivs.collar.Q1), 
                                meanlog = lnorm.tau2.lo[1],
                                sdlog = lnorm.tau2.lo[3])

# Q2
indivs.cam.Q2$tau2 <- rlnorm(nrow(indivs.cam.Q2), 
                             meanlog = lnorm.tau2.hi[1],
                             sdlog = lnorm.tau2.hi[3])

indivs.collar.Q2$tau2 <- rlnorm(nrow(indivs.collar.Q2), 
                                meanlog = lnorm.tau2.hi[1],
                                sdlog = lnorm.tau2.hi[3])

#_______________________________________________________________________
# 4d. Sigma ----

# sigma major
# geometric mean = 5,000 m
# SDs = 1,750, 5,000 m

lnorm.sigma.maj.lo <- log_norm_params(c(5000, 1750))
lnorm.sigma.maj.hi <- log_norm_params(c(5000, 5000))

# "aspect ratio"
# synthetic parameter that relates the major to minor axis sigmas
# these will also be log-normal draws, and then we'll calculate sigma minor
lnorm.aspect.lo <- log_norm_params(c(7, 1.75))
lnorm.aspect.hi <- log_norm_params(c(7, 7))

# angle
# this determine which axis the long dimension of the track aligns with
# it varies between - pi / 2 and pi / 2

#_______________________________________________________________________

# sigma major
# Q1
indivs.cam.Q1$sigma.maj <- rlnorm(nrow(indivs.cam.Q1), 
                                  meanlog = lnorm.sigma.maj.lo[1],
                                  sdlog = lnorm.sigma.maj.lo[3])

indivs.collar.Q1$sigma.maj <- rlnorm(nrow(indivs.collar.Q1), 
                                     meanlog = lnorm.sigma.maj.lo[1],
                                     sdlog = lnorm.sigma.maj.lo[3])

# Q2
indivs.cam.Q2$sigma.maj <- rlnorm(nrow(indivs.cam.Q2), 
                                  meanlog = lnorm.sigma.maj.hi[1],
                                  sdlog = lnorm.sigma.maj.hi[3])

indivs.collar.Q2$sigma.maj <- rlnorm(nrow(indivs.collar.Q2), 
                                     meanlog = lnorm.sigma.maj.hi[1],
                                     sdlog = lnorm.sigma.maj.hi[3])

# aspect ratio
# Q1
indivs.cam.Q1$aspect <- rlnorm(nrow(indivs.cam.Q1), 
                               meanlog = lnorm.aspect.lo[1],
                               sdlog = lnorm.aspect.lo[3])

indivs.collar.Q1$aspect <- rlnorm(nrow(indivs.collar.Q1), 
                                  meanlog = lnorm.aspect.lo[1],
                                  sdlog = lnorm.aspect.lo[3])

# Q2
indivs.cam.Q2$aspect <- rlnorm(nrow(indivs.cam.Q2), 
                               meanlog = lnorm.aspect.hi[1],
                               sdlog = lnorm.aspect.hi[3])

indivs.collar.Q2$aspect <- rlnorm(nrow(indivs.collar.Q2), 
                                  meanlog = lnorm.aspect.hi[1],
                                  sdlog = lnorm.aspect.hi[3])

# calculate sigma minor
indivs.cam.Q1$sigma.min         <- indivs.cam.Q1$sigma.maj / indivs.cam.Q1$aspect
indivs.collar.Q1$sigma.min      <- indivs.collar.Q1$sigma.maj / indivs.collar.Q1$aspect
indivs.cam.Q2$sigma.min         <- indivs.cam.Q2$sigma.maj / indivs.cam.Q2$aspect
indivs.collar.Q2$sigma.min      <- indivs.collar.Q2$sigma.maj / indivs.collar.Q2$aspect

# angle
indivs.cam.Q1$angle            <- runif(n = nrow(indivs.cam.Q1), -pi / 2, pi /2)
indivs.collar.Q1$angle         <- runif(n = nrow(indivs.collar.Q1), -pi / 2, pi /2)
indivs.cam.Q2$angle            <- runif(n = nrow(indivs.cam.Q2), -pi / 2, pi /2)
indivs.collar.Q2$angle         <- runif(n = nrow(indivs.collar.Q2), -pi / 2, pi /2)

#_______________________________________________________________________
# 5. Final check ----

# here we should make absolute certain that Tau2 is never > Tau1
# a while loop might be necessary here

#_______________________________________________________________________

# how many in each?
sum(indivs.cam.Q1$tau2 > indivs.cam.Q1$tau1)
sum(indivs.collar.Q1$tau2 > indivs.collar.Q1$tau1)
sum(indivs.cam.Q2$tau2 > indivs.cam.Q2$tau1)
sum(indivs.collar.Q2$tau2 > indivs.collar.Q2$tau1)

# replace with a new draw while there are still cases
fix_tau2 <- function(df) {
  
  while (sum(df$tau2 > df$tau1) > 0) {
    
    df$tau2[which(df$tau2 > df$tau1)] <- rlnorm(length(df$tau2[which(df$tau2 > df$tau1)]), 
                                                meanlog = lnorm.tau2.lo[1],
                                                sdlog = lnorm.tau2.lo[3])
    
  } 
    
  return(df)
    
}

# use function
indivs.cam.Q1 <- fix_tau2(indivs.cam.Q1)
indivs.collar.Q1 <- fix_tau2(indivs.collar.Q1)
indivs.cam.Q2 <- fix_tau2(indivs.cam.Q2)
indivs.collar.Q2 <- fix_tau2(indivs.collar.Q2)

# did I fix it?
sum(indivs.cam.Q1$tau2 > indivs.cam.Q1$tau1)
sum(indivs.collar.Q1$tau2 > indivs.collar.Q1$tau1)
sum(indivs.cam.Q2$tau2 > indivs.cam.Q2$tau1)
sum(indivs.collar.Q2$tau2 > indivs.collar.Q2$tau1)

#_______________________________________________________________________
# 5b. Examine distributions ----
#_______________________________________________________________________

hist(indivs.cam.Q1$tau1)
hist(indivs.cam.Q1$tau2)
hist(indivs.cam.Q1$sigma.maj)
hist(indivs.cam.Q1$sigma.min)

hist(indivs.cam.Q2$tau1)
hist(indivs.cam.Q2$tau2)
hist(indivs.cam.Q2$sigma.maj)
hist(indivs.cam.Q2$sigma.min)

#_______________________________________________________________________
# 7. Simulate a typical track ----
#_______________________________________________________________________

# simulation duration
sim.duration.wk <- 4

sim.duration.sec <- sim.duration.wk * 7 * 24 * 60 * 60 

# "sampling rate" for fundamental (i.e., straight line) steps
sim.samp.rate <- 60

# time steps to simulate on (from, to, by)
sim.timestep <- seq(1, sim.duration.sec, sim.samp.rate)

typical.track.ctmm <- ctmm(tau = c(2 %#% "hours",             # positional autocorr time
                                   1 %#% "hours"),            # velocity autocorr time
                           isotropic = FALSE,                 # anisotropic
                           sigma = c(5000,                   # variance along major axis
                                     714,                    # variance along minor axis
                                     0),                      # angle of axis
                           mu = c(0,                          # x coord
                                  0))                         # y coord

# run simulation
typical.track.sim <- simulate(object = typical.track.ctmm, 
                              t = sim.timestep,
                              complete = T)

# plot
ggplot(data = typical.track.sim) + 
  
  theme_bw() +
  
  geom_path(aes(x = x,
                y = y,
                color = timestamp)) +
  
  scale_color_viridis_c() +
  
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank())

#_______________________________________________________________________
# 6. Write to .csvs ----
#_______________________________________________________________________

write.csv(indivs.cam.Q1, paste0(getwd(), "/Derived data/Parameters/camera_lo.csv"))
write.csv(indivs.collar.Q1, paste0(getwd(), "/Derived data/Parameters/collar_lo.csv"))
write.csv(indivs.cam.Q2, paste0(getwd(), "/Derived data/Parameters/camera_hi.csv"))
write.csv(indivs.collar.Q2, paste0(getwd(), "/Derived data/Parameters/collar_hi.csv"))
