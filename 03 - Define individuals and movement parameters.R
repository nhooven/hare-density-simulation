# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 03 - Define individuals and movement parameters
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 25 Mar 2025
# Date completed: 25 Mar 2025
# Date last modified: 25 Mar 2025
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

# number of total replicates for each scenario
n.reps <- 6

#_______________________________________________________________________
# 2a. Question 1 ----

# For Q1, we have 5 different abundances totaling 92 individuals
# we can just take draws from these 92 for Qs 2 and 3

# importantly, these will all be ran for the full 8 week period even though
# for Q1 we only need half that

#_______________________________________________________________________

# abundances
abundances <- c(2, 5, 10, 25, 50)

# total indivs
n.indivs <- sum(abundances)

# expand grid by replicate
all.indivs <- expand.grid(indiv = 1:n.indivs,
                          rep = 1:n.reps)

# add in identifier to signify that these will all be used for camera contacts
all.indivs$use <- "camera"

# define function to sample without replacement and assign abundance IDs
subsample_indivs <- function(
    
  n.reps = n.reps,
  all.indivs = all.indivs
  
  ) {
  
  # define possible indivs
  possible.indivs <- 1:n.indivs
  
  # loop through iterations
  all.reps <- list()
  
  for (i in 1:n.reps) {
    
    # permute numbers
    permuted.indivs <- sample(possible.indivs)
    
    samp.1 <- permuted.indivs[1:50]
    samp.2 <- permuted.indivs[51:75]
    samp.3 <- permuted.indivs[76:85]
    samp.4 <- permuted.indivs[86:90]
    samp.5 <- permuted.indivs[91:92]
    
    # list
    sampled.indivs <- list(samp.1,
                           samp.2,
                           samp.3,
                           samp.4,
                           samp.5) 
    
    # bind in
    all.reps[[i]] <- sampled.indivs
    
  }
  
  # define sub-function to assign abundance IDs to each individual
  all.rows <- data.frame()
  
  for (j in 1:nrow(all.indivs)) {
    
    focal.row <- all.indivs[j, ]
    
    # subset list
    focal.list <- all.reps[[focal.row$rep]]
    
    # which abundance group is it in?
    focal.row$abund.group <- case_when(focal.row$indiv %in% focal.list[[1]] ~ abundances[5],
                                       focal.row$indiv %in% focal.list[[2]] ~ abundances[4],
                                       focal.row$indiv %in% focal.list[[3]] ~ abundances[3],
                                       focal.row$indiv %in% focal.list[[4]] ~ abundances[2],
                                       focal.row$indiv %in% focal.list[[5]] ~ abundances[1])
    
    # bind in
    all.rows <- rbind(all.rows, focal.row)
    
  }
  
  # return
  return(all.rows)
  
}

# use function
sampled.indivs <- subsample_indivs(n.reps, all.indivs)

#_______________________________________________________________________
# 2b. Question 2 ----

# For Q2, we'll draw 10 individuals from each replicate, at random
q2.size <- 10

#_______________________________________________________________________

sampled.indivs.1 <- data.frame()

for (i in 1:n.reps) {
  
  # focal rep
  focal.rep <- sampled.indivs %>% filter(rep == i)
  
  # sample 10 at random
  q2.samp <- sample(1:nrow(focal.rep), size = q2.size)
  
  # add in indicator
  focal.rep$q2 <- ifelse(focal.rep$indiv %in% q2.samp,
                         1,
                         0)
  
  # bind in
  sampled.indivs.1 <- rbind(sampled.indivs.1, focal.rep)
  
}

# sanity check
sum(sampled.indivs.1$q2) == q2.size * n.reps

#_______________________________________________________________________
# 2c. "Collared" individuals ----

# we'll simulate 10 tracks per replicate to add here, and draw 5 from each
# replicate as a supplement (perhaps mirroring how we do this in real life)
collared.size <- 5

add.collared.size <- 10

#_______________________________________________________________________

# sample 5 collared individuals from each replicate
sampled.indivs.2 <- data.frame()

for (i in 1:n.reps) {
  
  # focal rep
  focal.rep <- sampled.indivs %>% filter(rep == i)
  
  # sample 10 at random
  collared.samp <- sample(1:nrow(focal.rep), size = collared.size)
  
  # add in indicator
  focal.rep$collared <- ifelse(focal.rep$indiv %in% collared.samp,
                               1,
                               0)
  
  # bind in
  sampled.indivs.2 <- rbind(sampled.indivs.2, focal.rep)
  
}

# sanity check
sum(sampled.indivs.2$collared) == collared.size * n.reps

# create data frame for additional individuals (Qs 1 and 2 only)
add.collared.indivs <- expand.grid(indiv = 1:add.collared.size,
                                   rep = 1:n.reps,
                                   use = "collar",
                                   q = c("q1", "q2"))

#_______________________________________________________________________
# 2d. Question 3 ----

# these are the individuals with the high variability in their
# OUF parameters - we need 60 camera and 60 collar critters

#_______________________________________________________________________

# create data frame for additional individuals (Qs 1 and 2 only)
q3.indivs <- expand.grid(indiv = 1:10,
                         rep = 1:n.reps,
                         use = "camera",
                         q = "q3")

# sample 5 collared individuals from each replicate
q3.indivs.1 <- data.frame()

for (i in 1:n.reps) {
  
  # focal rep
  focal.rep <- q3.indivs %>% filter(rep == i)
  
  # sample 10 at random
  collared.samp <- sample(1:nrow(focal.rep), size = collared.size)
  
  # add in indicator
  focal.rep$collared <- ifelse(focal.rep$indiv %in% collared.samp,
                               1,
                               0)
  
  # bind in
  q3.indivs.1 <- rbind(q3.indivs.1, focal.rep)
  
}

# create data frame for additional individuals (Q3)
q3.add.collared.indivs <- expand.grid(indiv = 1:add.collared.size,
                                      rep = 1:n.reps,
                                      use = "collar",
                                      q = "q3")

#_______________________________________________________________________
# 3. Define distributions of movement parameters ----

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
# 3a. Home range centroids ----

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

# Q1 and Q2
sampled.indivs.2$mean1 <- runif(n = nrow(sampled.indivs.2), coord.min, coord.max)
sampled.indivs.2$mean2 <- runif(n = nrow(sampled.indivs.2), coord.min, coord.max)

add.collared.indivs$mean1 <- runif(n = nrow(add.collared.indivs), coord.min, coord.max)
add.collared.indivs$mean2 <- runif(n = nrow(add.collared.indivs), coord.min, coord.max)

# Q3
q3.indivs.1$mean1 <- runif(n = nrow(q3.indivs.1), coord.min, coord.max)
q3.indivs.1$mean2 <- runif(n = nrow(q3.indivs.1), coord.min, coord.max)

q3.add.collared.indivs$mean1 <- runif(n = nrow(q3.add.collared.indivs), coord.min, coord.max)
q3.add.collared.indivs$mean2 <- runif(n = nrow(q3.add.collared.indivs), coord.min, coord.max)

#_______________________________________________________________________
# 3b. Tau 1 ----

# geometric mean = 6 hr
# SDs = 1.5, 6 hr

lnorm.tau1.lo <- log_norm_params(c(6 %#% "hours", 1.5 %#% "hours"))
lnorm.tau1.hi <- log_norm_params(c(6 %#% "hours", 6 %#% "hours"))

#_______________________________________________________________________

# Q1 and Q2
sampled.indivs.2$tau1 <- rlnorm(nrow(sampled.indivs.2), 
                                meanlog = lnorm.tau1.lo[1],
                                sdlog = lnorm.tau1.lo[3])

add.collared.indivs$tau1 <- rlnorm(nrow(add.collared.indivs), 
                                   meanlog = lnorm.tau1.lo[1],
                                   sdlog = lnorm.tau1.lo[3])

# Q3
q3.indivs.1$tau1 <- rlnorm(nrow(q3.indivs.1), 
                           meanlog = lnorm.tau1.hi[1],
                           sdlog = lnorm.tau1.hi[3])

q3.add.collared.indivs$tau1 <- rlnorm(nrow(q3.add.collared.indivs), 
                                      meanlog = lnorm.tau1.hi[1],
                                      sdlog = lnorm.tau1.hi[3])

#_______________________________________________________________________
# 3c. Tau 2 ----

# geometric mean = 2 hr
# SDs = 0.5, 2 hr

lnorm.tau2.lo <- log_norm_params(c(2 %#% "hours", 0.5 %#% "hours"))
lnorm.tau2.hi <- log_norm_params(c(2 %#% "hours", 2 %#% "hours"))

#_______________________________________________________________________

# Q1 and Q2
sampled.indivs.2$tau2 <- rlnorm(nrow(sampled.indivs.2), 
                                meanlog = lnorm.tau2.lo[1],
                                sdlog = lnorm.tau2.lo[3])

add.collared.indivs$tau2 <- rlnorm(nrow(add.collared.indivs), 
                                   meanlog = lnorm.tau2.lo[1],
                                   sdlog = lnorm.tau2.lo[3])

# Q3
q3.indivs.1$tau2 <- rlnorm(nrow(q3.indivs.1), 
                           meanlog = lnorm.tau2.hi[1],
                           sdlog = lnorm.tau2.hi[3])

q3.add.collared.indivs$tau2 <- rlnorm(nrow(q3.add.collared.indivs), 
                                      meanlog = lnorm.tau2.hi[1],
                                      sdlog = lnorm.tau2.hi[3])

#_______________________________________________________________________
# 3d. Sigma ----

# sigma major
# geometric mean = 10,000 m
# SDs = 2,500, 10,000 m

lnorm.sigma.maj.lo <- log_norm_params(c(10000, 2500))
lnorm.sigma.maj.hi <- log_norm_params(c(10000, 10000))

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
# Q1 and Q2
sampled.indivs.2$sigma.maj <- rlnorm(nrow(sampled.indivs.2), 
                                     meanlog = lnorm.sigma.maj.lo[1],
                                     sdlog = lnorm.sigma.maj.lo[3])

add.collared.indivs$sigma.maj <- rlnorm(nrow(add.collared.indivs), 
                                        meanlog = lnorm.sigma.maj.lo[1],
                                        sdlog = lnorm.sigma.maj.lo[3])

# Q3
q3.indivs.1$sigma.maj <- rlnorm(nrow(q3.indivs.1), 
                                meanlog = lnorm.sigma.maj.hi[1],
                                sdlog = lnorm.sigma.maj.hi[3])

q3.add.collared.indivs$sigma.maj <- rlnorm(nrow(q3.add.collared.indivs), 
                                           meanlog = lnorm.sigma.maj.hi[1],
                                           sdlog = lnorm.sigma.maj.hi[3])

# aspect ratio
# Q1 and Q2
sampled.indivs.2$aspect <- rlnorm(nrow(sampled.indivs.2), 
                                  meanlog = lnorm.aspect.lo[1],
                                  sdlog = lnorm.aspect.lo[3])

add.collared.indivs$aspect <- rlnorm(nrow(add.collared.indivs), 
                                     meanlog = lnorm.aspect.lo[1],
                                     sdlog = lnorm.aspect.lo[3])

# Q3
q3.indivs.1$aspect <- rlnorm(nrow(q3.indivs.1), 
                             meanlog = lnorm.aspect.hi[1],
                             sdlog = lnorm.aspect.hi[3])

q3.add.collared.indivs$aspect <- rlnorm(nrow(q3.add.collared.indivs), 
                                        meanlog = lnorm.aspect.hi[1],
                                        sdlog = lnorm.aspect.hi[3])

# calculate sigma minor
sampled.indivs.2$sigma.min        <- sampled.indivs.2$sigma.maj / sampled.indivs.2$aspect
add.collared.indivs$sigma.min     <- add.collared.indivs$sigma.maj / add.collared.indivs$aspect
q3.indivs.1$sigma.min             <- q3.indivs.1$sigma.maj / q3.indivs.1$aspect
q3.add.collared.indivs$sigma.min  <- q3.add.collared.indivs$sigma.maj / q3.add.collared.indivs$aspect

# angle
sampled.indivs.2$angle            <- runif(n = nrow(sampled.indivs.2), -pi / 2, pi /2)
add.collared.indivs$angle         <- runif(n = nrow(add.collared.indivs), -pi / 2, pi /2)
q3.indivs.1$angle                 <- runif(n = nrow(q3.indivs.1), -pi / 2, pi /2)
q3.add.collared.indivs$angle      <- runif(n = nrow(q3.add.collared.indivs), -pi / 2, pi /2)

#_______________________________________________________________________
# 4. Final check ----

# let's make sure we know how many tracks we're simulating

#_______________________________________________________________________

sum(
  
  nrow(sampled.indivs.2),          # 552
  nrow(add.collared.indivs),       # 120
  nrow(q3.indivs.1),               # 60
  nrow(q3.add.collared.indivs)     # 60
  
  )

# total expected time (simulations only)
(23 * 552) / 3600           # 3.5 hr
(23 * 120) / 3600           # ~ 45 mins
(23 * 60) / 3600            # ~ 20 mins 
(23 * 60) / 3600            # ~ 20 mins 

#_______________________________________________________________________
# 5. Write to .csvs ----
#_______________________________________________________________________

write.csv(sampled.indivs.2, paste0(getwd(), "/Derived data/Individuals and parameters/camera_lo.csv"))
write.csv(add.collared.indivs, paste0(getwd(), "/Derived data/Individuals and parameters/collar_lo.csv"))
write.csv(q3.indivs.1, paste0(getwd(), "/Derived data/Individuals and parameters/camera_hi.csv"))
write.csv(q3.add.collared.indivs, paste0(getwd(), "/Derived data/Individuals and parameters/collar_hi.csv"))
