# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 01 - Build CTMM simulations
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 06 Mar 2025
# Date completed: 
# Date last modified: 07 Mar 2025 
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
# or restricted space use.

# It would be nice to see how this model generates data initially.

# Is this what the REM is assuming? Not really because the animal must be 
# "resident" within the sampling grid of the cameras.

#_______________________________________________________________________

# initialize a CTMM object for simulation
ctmm.bm <- ctmm(
  
  tau = Inf,            # tau is infinite, no range residency
  sigma = 50,            # Brownian motion variance
  mu = c(0, 0)          # starting location
  
  )

# time
t.bm <- 1:10000

# simulate
sim.bm <- simulate(object = ctmm.bm,
                   t = t.bm)

# plot
ggplot(sim.bm,
       aes(x = x,
           y = y)) +
  
  theme_bw() +
  
  geom_path()

#_______________________________________________________________________
# 4. Ornstein-Uhlenbeck Foraging model ----

# OUF accounts for both range residency and velocity autocorrelation-
# two key components of real animal movement

# In order to fit this to real tracking data, they need to be pretty finely-sampled.
# However, this is the only model that can allow speed estimation and thus
# is the only CTSP that I know of that is appropriate for estimating
# day range in the REM

# Our approach will be to assume animal movement can be well-approximated with
# one of these models and test various departures from an "ideal" scenario,
# e.g., high variance across individuals, short tracking durations, fix success, etc.

#_______________________________________________________________________

# initial simulations just to see what this looks like
# initialize a CTMM object for simulation
ctmm.ouf <- ctmm(
  
  tau = c(100, 10),     # two tau parameters (residency and velo autocorr)
  sigma = 50,           # asymptotic variance
  mu = c(0, 0)          # starting location
  
)

# time
t.ouf <- 1:10000

# simulate
sim.ouf <- simulate(object = ctmm.ouf,
                    t = t.ouf,
                    complete = TRUE)

# plot
ggplot(sim.ouf,
       aes(x = x,
           y = y)) +
  
  theme_bw() +
  
  geom_path(aes(color = timestamp)) +
  
  scale_color_viridis_c()

# let's compare different parameter values

#_______________________________________________________________________
# 4a. Define plotting function ----
#_______________________________________________________________________

plot_ouf <- function(sim,         # simulated telemetry object
                     param,       # tau.r, tau.v, sigma
                     ramp) {      # low, mid, or hig
  
  plot.out <- ggplot(data = sim) + 
    
    theme_bw() +
    
    geom_path(aes(x = x,
                  y = y,
                  color = timestamp)) +
    
    scale_color_viridis_c() +
    
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.title = element_blank()) +
    
    labs(title = paste0(param, " - ", ramp))
  
  return(plot.out)
    
}        

#_______________________________________________________________________
# NOTE 1 ----
#_______________________________________________________________________

# According to https://ctmm-initiative.github.io/ctmm/articles/variogram.html
# all units for parameters are in m and seconds

# You can use the utility function "%#%" to convert these, for example
6 %#% "day"

13 %#% "km^2"

#_______________________________________________________________________
# 4b. tau_r ----

# This parameter controls the timescale for range residency (first tau parameter)
# note that this will always be relative to the REAL time

# From Calabrese et al. 2016: position correlations die off over some period of time tau_r

#_______________________________________________________________________

# different tau_r values
tau.r.low <- c(1 %#% "minutes", 10)
tau.r.mid <- c(5 %#% "minutes", 10)
tau.r.hig <- c(10 %#% "minutes", 10)

# fit models
ctmm.ouf.tau.r.low <- ctmm(tau = tau.r.low, sigma = 50, mu = c(0, 0))
ctmm.ouf.tau.r.mid <- ctmm(tau = tau.r.mid, sigma = 50, mu = c(0, 0))
ctmm.ouf.tau.r.hig <- ctmm(tau = tau.r.hig, sigma = 50, mu = c(0, 0))

# simulate 10000 points
sim.ouf.tau.r.low <- simulate(object = ctmm.ouf.tau.r.low, t = t.ouf, complete = TRUE)
sim.ouf.tau.r.mid <- simulate(object = ctmm.ouf.tau.r.mid, t = t.ouf, complete = TRUE)
sim.ouf.tau.r.hig <- simulate(object = ctmm.ouf.tau.r.hig, t = t.ouf, complete = TRUE)

# plot together
plot_grid(plot_ouf(sim.ouf.tau.r.low,
                   "tau.r",
                   "low"),
          plot_ouf(sim.ouf.tau.r.mid,
                   "tau.r",
                   "mid"),
          plot_ouf(sim.ouf.tau.r.hig,
                   "tau.r",
                   "hig"))

# examine variograms
plot(variogram(sim.ouf.tau.r.low))
plot(variogram(sim.ouf.tau.r.mid))
plot(variogram(sim.ouf.tau.r.hig))

# these show how larger tau_r values lead to less concentrated movements - cool!

#_______________________________________________________________________
# 4c. tau_v ----

# This parameter controls the timescale for velocity autocorrelation (second tau parameter).
# This makes most sense if it is much smaller than tau_r

#_______________________________________________________________________

# different tau_v values
tau.v.low <- c(5 %#% "minutes", 1 %#% "seconds")
tau.v.mid <- c(5 %#% "minutes", 5 %#% "seconds")
tau.v.hig <- c(5 %#% "minutes", 10 %#% "seconds")

# fit models
ctmm.ouf.tau.v.low <- ctmm(tau = tau.v.low, sigma = 50, mu = c(0, 0))
ctmm.ouf.tau.v.mid <- ctmm(tau = tau.v.mid, sigma = 50, mu = c(0, 0))
ctmm.ouf.tau.v.hig <- ctmm(tau = tau.v.hig, sigma = 50, mu = c(0, 0))

# simulate 10000 points
sim.ouf.tau.v.low <- simulate(object = ctmm.ouf.tau.v.low, t = t.ouf, complete = TRUE)
sim.ouf.tau.v.mid <- simulate(object = ctmm.ouf.tau.v.mid, t = t.ouf, complete = TRUE)
sim.ouf.tau.v.hig <- simulate(object = ctmm.ouf.tau.v.hig, t = t.ouf, complete = TRUE)

# plot together
plot_grid(plot_ouf(sim.ouf.tau.v.low,
                   "tau.v",
                   "low"),
          plot_ouf(sim.ouf.tau.v.mid,
                   "tau.v",
                   "mid"),
          plot_ouf(sim.ouf.tau.v.hig,
                   "tau.v",
                   "hig"))

# examine variograms
plot(variogram(sim.ouf.tau.v.low))
plot(variogram(sim.ouf.tau.v.mid))
plot(variogram(sim.ouf.tau.v.hig))

# larger tau_v values allow for more directional persistence in parts of the track
# and less concentrated, tortuous movements overall

#_______________________________________________________________________
# 4d. Sigma ----

# This is the asymptotic variance (presumably controlling the total area used)

#_______________________________________________________________________

# different sigma values
sigma.low <- 10 %#% "ha"
sigma.mid <- 50 %#% "ha"
sigma.hig <- 100 %#% "ha"

# fit models
ctmm.ouf.sigma.low <- ctmm(tau = c(5 %#% "minutes", 5 %#% "seconds"), sigma = sigma.low, mu = c(0, 0))
ctmm.ouf.sigma.mid <- ctmm(tau = c(5 %#% "minutes", 5 %#% "seconds"), sigma = sigma.mid, mu = c(0, 0))
ctmm.ouf.sigma.hig <- ctmm(tau = c(5 %#% "minutes", 5 %#% "seconds"), sigma = sigma.hig, mu = c(0, 0))

# simulate 10000 points
sim.ouf.sigma.low <- simulate(object = ctmm.ouf.sigma.low, t = t.ouf, complete = TRUE)
sim.ouf.sigma.mid <- simulate(object = ctmm.ouf.sigma.mid, t = t.ouf, complete = TRUE)
sim.ouf.sigma.hig <- simulate(object = ctmm.ouf.sigma.hig, t = t.ouf, complete = TRUE)

# plot together
plot_grid(plot_ouf(sim.ouf.sigma.low,
                   "sigma",
                   "low"),
          plot_ouf(sim.ouf.sigma.mid,
                   "sigma",
                   "mid"),
          plot_ouf(sim.ouf.sigma.hig,
                   "sigma",
                   "hig"))

# examine variograms
plot(variogram(sim.ouf.sigma.low))
plot(variogram(sim.ouf.sigma.mid))
plot(variogram(sim.ouf.sigma.hig))

# larger sigma values change the overall scale of space use and should be in unit area
# I cranked these up to ha

# this would be a good parameter to vary to allow for different "home range sizes"

#_______________________________________________________________________
# 4e. Isotropy ----

# This controls whether "home ranges" are circular or elliptical

#_______________________________________________________________________

# fit models
ctmm.ouf.iso.f <- ctmm(tau = c(5 %#% "minutes", 5 %#% "seconds"), isotropic = FALSE, sigma = sigma.low, mu = c(0, 0))
ctmm.ouf.iso.t <- ctmm(tau = c(5 %#% "minutes", 5 %#% "seconds"), isotropic = TRUE, sigma = sigma.low, mu = c(0, 0))

# simulate 10000 points
sim.ouf.iso.f <- simulate(object = ctmm.ouf.iso.f, t = t.ouf, complete = TRUE)
sim.ouf.iso.t <- simulate(object = ctmm.ouf.iso.t, t = t.ouf, complete = TRUE)

# plot together
plot_grid(plot_ouf(sim.ouf.sigma.low,
                   "iso",
                   "f"),
          plot_ouf(sim.ouf.sigma.mid,
                   "iso",
                   "t"))

# when isotropic == TRUE, the resultant home range is more "circular" but
# is also much larger
# maybe this isn't the best parameter to vary - we'll look at some real data

#_______________________________________________________________________
# NOTE 2 ----
#_______________________________________________________________________

# These parameters are important for controlling movement behavior but also
# the quality of data needed to recover them with a ctmm fit
# ALWAYS keep that in mind when running the simulations and analysis for this study

#_______________________________________________________________________
# 5. Examine some real hare data ----

# 420's data will probably be the best to look at here

#_______________________________________________________________________
# 5a. Clean data and prepare for modeling ----
#_______________________________________________________________________

# read in .csv for 420 
hare.data <- read.csv("E:/Hare project/GPS data/alldata/DUR/005_420_2.csv", 
                      sep = "",           # required to separate columns
                      fill = TRUE)        # required to separate columns

# examine
summary(hare.data)

# clean
hare.data.1 <- hare.data %>%
  
  # keep only lat-long entries
  mutate(location.lon = as.numeric(location.lon),
         location.lat = as.numeric(location.lat)) %>%
  
  filter(location.lon > -121 & location.lon < -117 &
         location.lat > 47 & location.lat < 50) %>%
  
  # create timestamp and coerce to POSIXct
  mutate(timestamp = dmy_hms(paste0(Date, " ", Time))) %>%
  
  # keep columns we want
  dplyr::select(location.lon, 
                location.lat,
                timestamp,
                height.msl,
                satellites,
                hdop) %>%
  
  # coerce required columns
  mutate(height.msl = as.numeric(height.msl)) %>%
  
  # drop any NAs in the key columns
  drop_na(location.lon,
          location.lat,
          timestamp)

# examine erroneous locations
ggplot(data = hare.data.1,
       aes(x = location.lon,
           y = location.lat,
           color = timestamp)) +
  
  theme_bw() +
  
  geom_point()

# looks like we've got some, let's also look at hdop
ggplot(data = hare.data.1,
       aes(x = hdop,
           y = height.msl,
           color = timestamp)) +
  
  theme_bw() +
  
  geom_point()

# drop erroneous locations
hare.data.2 <- hare.data.1 %>%
  
  filter(location.lon < 119.7 &
         location.lat < 48.8 & location.lat > 48.75 &
         height.msl > 1400 & height.msl < 1600 &
         hdop < 10)

# examine again
ggplot(data = hare.data.2,
       aes(x = location.lon,
           y = location.lat,
           color = timestamp)) +
  
  theme_bw() +
  
  geom_point()

ggplot(data = hare.data.2,
       aes(x = hdop,
           y = height.msl,
           color = timestamp)) +
  
  theme_bw() +
  
  geom_point()

# these look good for running a model now!

#_______________________________________________________________________
# 5b. Spatialize and coerce to telemetry object ----
#_______________________________________________________________________

# promote to sf object
hare.sf <- st_as_sf(hare.data.2,
                    coords = c("location.lon",
                               "location.lat"),
                    crs = "EPSG:4326")

# transform to projected coords (UTM 11N)
hare.sf.utm <- st_transform(hare.sf, crs = "EPSG:32611")

# add in coords as columns
hare.sf.utm$x <- st_coordinates(hare.sf.utm)[ , 1]
hare.sf.utm$y <- st_coordinates(hare.sf.utm)[ , 2]

# plot for sanity check
plot(hare.sf.utm)

# make a track
hare.track <- make_track(hare.sf.utm,
                         .x = x,
                         .y = y,
                         .t = timestamp,
                         crs = "EPSG:32611",
                         all_cols = TRUE)

# convert to telemetry object
hare.telem <- as_telemetry(hare.track)

# plot 
plot(hare.telem, error = FALSE)

#_______________________________________________________________________
# 5c. Fit CTMMs ----
#_______________________________________________________________________