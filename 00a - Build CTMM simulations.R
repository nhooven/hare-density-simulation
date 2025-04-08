# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 00a - Build CTMM simulations
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 06 Mar 2025
# Date completed: 25 Mar 2025 
# Date last modified: 07 Apr 2025 
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

# OUF accounts for range residency, positional autocorrelation, and velocity autocorrelation-
# three key components of real animal movement

# In order to fit this to real tracking data, they need to be pretty finely-sampled.
# However, this is the only range-resident model that can allow speed estimation and thus
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

# This parameter controls the timescale for positional autocorrelation (first tau parameter)
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

# This is the asymptotic variance, defined by both the major and minor axes
# (assuming anisotropy, which we likely will)

# below we'll assume isotropy

# I'm not convinced that this is in area, I think these are linear (m)
# these tend to extend the "home range" away from zero by ~ sigma / 10

#_______________________________________________________________________

# different sigma values
sigma.low <- 100
sigma.mid <- 500
sigma.hig <- 1000

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
# 5. Tune simulation parameters (and time them) ----

# now the rubber meets the road - we want to examine how to make reasonable,
# realistic tracks that can be used both to tally camera contacts AND
# subsample to fit and simulate more models for

#_______________________________________________________________________
# 5a. Track duration ----
#_______________________________________________________________________

# how many seconds are in 4 weeks?
(sec.4wk <- 4 * 7 * 24 * 60 * 60)

# each simulation must run for 2,419,200 steps (if we want a location every second)
# this seems extravagant; let's assume that our animal's fundamental step length 
# is straight over a period of 60 seconds
# let's time this

seq.4wk <- seq(1, sec.4wk, 60)

length(seq.4wk)

ctmm.ouf.duration <- ctmm(tau = c(6 %#% "hours", 
                                  2 %#% "hours"), 
                          isotropic = FALSE, 
                          sigma = c(700,
                                    100,
                                    -0.05), 
                          mu = c(0, 0))

# simulate 10000 points
start.time <- Sys.time()

sim.ouf.duration <- simulate(object = ctmm.ouf.duration, 
                             t = seq.4wk, 
                             complete = TRUE)

# difftime
Sys.time() - start.time

# this took ~ 4.98 seconds. If we wanted 1000 simulated tracks this would take:
(4.98 * 1000) / 3600

# ~ 1.38 h

# plot the resultant track
plot_ouf(sim.ouf.duration,
         "4 weeks",
         "60-sec fundamental step")

#_______________________________________________________________________
# 5b. Tau 1 - positional autocorrelation ----

# We want this to be longer than any of our fix rates, so 6 hr seems reasonable

#_______________________________________________________________________

# we'll look at track differences between 4 and 6 hours
ctmm.ouf.tau1.1 <- ctmm(tau = c(6 %#% "hours", 
                                2 %#% "hours"), 
                        isotropic = FALSE, 
                        sigma = c(7000,
                                  1000,
                                  -0.05), 
                        mu = c(0, 0))

ctmm.ouf.tau1.2 <- ctmm(tau = c(12 %#% "hours", 
                                2 %#% "hours"), 
                        isotropic = FALSE, 
                        sigma = c(7000,
                                  1000,
                                  -0.05), 
                        mu = c(0, 0))

# simulate
sim.ouf.tau1.1 <- simulate(object = ctmm.ouf.tau1.1, 
                           t = seq.4wk, 
                           complete = TRUE)

sim.ouf.tau1.2 <- simulate(object = ctmm.ouf.tau1.2, 
                           t = seq.4wk, 
                           complete = TRUE)

# plot together
plot_grid(plot_ouf(sim.ouf.tau1.1,
                   "tau.p",
                   "6 hr"),
          plot_ouf(sim.ouf.tau1.2,
                   "tau.p",
                   "12 hr"))

# maybe we could do 3 h for this, ~ 3x Tau2 at 1 h

# occurrence distributions (subsampled to hourly)
occur1 <- occurrence(sim.ouf.tau1.1[seq(1, 40320, 60), ], CTMM = ctmm.ouf.tau1.1)
occur2 <- occurrence(sim.ouf.tau1.2[seq(1, 40320, 60), ], CTMM = ctmm.ouf.tau1.2)

plot(occur1)
plot(occur2)

# convert to sf
occur1.sf <- as(SpatialPolygonsDataFrame.UD(occur1, level.UD = 0.75), "sf")
occur2.sf <- as(SpatialPolygonsDataFrame.UD(occur2, level.UD = 0.75), "sf")

plot(st_geometry(occur1.sf))
plot(st_geometry(occur2.sf))

st_area(occur1.sf)
st_area(occur2.sf)

# is there too much core area in the simulations? is space use too concentrated???

#_______________________________________________________________________
# 5c. Tau 2 - velocity autocorrelation timescale ----

# this will be longer than our finest fix rate, equal to our medium fix rate, and
# shorter than our coarsest fix rate

# 1 hours seems reasonable and similar to real data
# then 3 h would be 3tau which should not allow for OUF fitting

#_______________________________________________________________________

# we'll look at track differences between 1 and 0.5 hours
ctmm.ouf.tau2.1 <- ctmm(tau = c(12 %#% "hours", 
                                1 %#% "hours"), 
                        isotropic = FALSE, 
                        sigma = c(1000,
                                  100,
                                  -0.05), 
                        mu = c(0, 0))

ctmm.ouf.tau2.2 <- ctmm(tau = c(12 %#% "hours", 
                                0.5 %#% "hours"), 
                        isotropic = FALSE, 
                        sigma = c(1000,
                                  100,
                                  -0.05), 
                        mu = c(0, 0))

# simulate
sim.ouf.tau2.1 <- simulate(object = ctmm.ouf.tau2.1, 
                           t = seq.4wk, 
                           complete = TRUE)

sim.ouf.tau2.2 <- simulate(object = ctmm.ouf.tau2.2, 
                           t = seq.4wk, 
                           complete = TRUE)

# plot together
plot_grid(plot_ouf(sim.ouf.tau2.1,
                   "tau.v",
                   "2 hr"),
          plot_ouf(sim.ouf.tau2.2,
                   "tau.v",
                   "0.5 hr"))

# it's easy to see here why you need fine fix rates to resolve tau v
# we'll work with them near to the grain of our fix rate,
# recognizing that this may not be the case with our real data

#_______________________________________________________________________
# 5d. Sigmas - asymptotic variance ----

# since we'll simulate anisotropic home ranges, we'll need to tune these
# to give the elongate, elliptical space use areas we see
# we can probably allow these to vary slightly for Qs 1 and 2

# angles vary from pi / 2 to - pi /2 because they represent
# deviations from 0 (i.e., aligned N-S, or in the case of ctmm's projection, E-W)

#_______________________________________________________________________

# let's do three sizes for the major axis:
# 20,000 m, 10,000 m, 5,000 m

# mean aspect ratio for real hare data ~ 7.4,
# so minor axes will be major / 7
sigmas.major <- c(20000, 10000, 5000)
sigmas.minor <- sigmas.major / 7

ctmm.ouf.sigma.1 <- ctmm(tau = c(6 %#% "hours", 
                                 2 %#% "hours"), 
                         isotropic = FALSE, 
                         sigma = c(sigmas.major[1],
                                  sigmas.minor[1],
                                  runif(1, 
                                        -pi / 2, 
                                        pi / 2)), 
                         mu = c(0, 0))

ctmm.ouf.sigma.2 <- ctmm(tau = c(6 %#% "hours", 
                                 2 %#% "hours"), 
                         isotropic = FALSE, 
                         sigma = c(sigmas.major[2],
                                   sigmas.minor[2],
                                   runif(1, 
                                         -pi / 2, 
                                         pi / 2)), 
                         mu = c(0, 0))

ctmm.ouf.sigma.3 <- ctmm(tau = c(6 %#% "hours", 
                                 2 %#% "hours"), 
                         isotropic = FALSE, 
                         sigma = c(sigmas.major[3],
                                   sigmas.minor[3],
                                   runif(1, 
                                         -pi / 2, 
                                         pi / 2)), 
                         mu = c(0, 0))

# simulate
sim.ouf.sigma.1 <- simulate(object = ctmm.ouf.sigma.1, 
                            t = seq.4wk, 
                            complete = TRUE)

sim.ouf.sigma.2 <- simulate(object = ctmm.ouf.sigma.2, 
                            t = seq.4wk, 
                            complete = TRUE)

sim.ouf.sigma.3 <- simulate(object = ctmm.ouf.sigma.3, 
                            t = seq.4wk, 
                            complete = TRUE)

# plot together
plot_grid(plot_ouf(sim.ouf.sigma.1,
                   "Major = 2e4",
                   "Minor = 2e4 / 7"),
          plot_ouf(sim.ouf.sigma.2,
                   "Major = 1e4",
                   "Minor = 1e4 / 7"),
          plot_ouf(sim.ouf.sigma.3,
                   "Major = 5e3",
                   "Minor = 5e3 / 7"))

# these change the spatial scale of movement, so for our purposes, this is related
# to the total space use area, or range size

# we'll want to tune the variability here for both the "low" and "high" scenarios

#_______________________________________________________________________
# 6. Movement parameter variability ----

# for Q 1, we'll keep these low and centered on the mean
# for Q2, allowing them to vary widely is central to our question

# these should be drawn from a strictly positive distribution, so the log-normal makes sense

# recall that the two parameters of the log-normal distribution are the 
# mean and SD of the LOGARITHM, not the parameter itself

# ln(X) ~ N(mu, sigma)

# From Wikipedia: https://en.wikipedia.org/wiki/Log-normal_distribution

# "In order to produce a distribution with desired mean mu_x and variance sigma^2_x:

# mu = ln(mu_x^2 / sqrt(mu_x^2 + sigma^2_x))

# sigma^2 = ln(1 + (sigma^2_x / mu_x^2))

# so we need a helper function:
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

# HOWEVER - the mean and SD might be less informative than the median and the "scatter"

#_______________________________________________________________________
# 6a. Tau 1 ---- 

# in our real data, Tau 1 varies wildly (2,200 to 80,000, corresponding to
# 0.61 to 22 hr, with a mean of ~ 11 hr).

# our mean will be 6 hours
# SDs range will be a quarter of the mean to the full mean

#_______________________________________________________________________

# define parameters
lnorm.tau1.lo <- log_norm_params(c(6 %#% "hours", 1.5 %#% "hours"))
lnorm.tau1.hi <- log_norm_params(c(6 %#% "hours", 6 %#% "hours"))

# both means should be 6 hours
tau.1.lo <- dlnorm(x = seq(0, 
                           24 %#% "hours", 
                           length.out = 500),
                   meanlog = lnorm.tau1.lo[1],       # geometric mean here
                   sdlog = lnorm.tau1.lo[3]) 

tau.1.hi <- dlnorm(x = seq(0, 
                           24 %#% "hours", 
                           length.out = 500),
                   meanlog = lnorm.tau1.hi[1],
                   sdlog = lnorm.tau1.hi[3]) 

# df
tau.1.df <- data.frame(x = seq(0, 
                               24 %#% "hours", 
                               length.out = 500),
                       y = c(tau.1.lo,
                             tau.1.hi),
                       variability = rep(c("low",
                                           "high"),
                                         each = 500))

# plot
ggplot(data = tau.1.df,
       aes(x = x,
           y = y,
           color = variability)) +
  
  theme_bw() +
  
  geom_line(linewidth = 1.1) +
  
  scale_x_continuous(breaks = c(0, 6, 12, 18, 24) %#% "hours",
                     labels = c(0, 6, 12, 18, 24)) +
  
  xlab("Positional autocorrelation timescale (hr)")

# look at means
mean(rlnorm(n = 500, meanlog = lnorm.tau1.lo[1], sdlog = lnorm.tau1.lo[3])) / 3600
mean(rlnorm(n = 500, meanlog = lnorm.tau1.hi[1], sdlog = lnorm.tau1.hi[3])) / 3600

# look at medians
median(rlnorm(n = 500, meanlog = lnorm.tau1.lo[1], sdlog = lnorm.tau1.lo[3])) / 3600
median(rlnorm(n = 500, meanlog = lnorm.tau1.hi[1], sdlog = lnorm.tau1.hi[3])) / 3600

# and SDs
sd(rlnorm(n = 500, meanlog = lnorm.tau1.lo[1], sdlog = lnorm.tau1.lo[3])) / 3600
sd(rlnorm(n = 500, meanlog = lnorm.tau1.hi[1], sdlog = lnorm.tau1.hi[3])) / 3600

#_______________________________________________________________________
# 6b. Tau 2 ---- 

# we have very little information about tau 2 (velocity)
# from our limited real data (20 OUf and 4 OUF)

# I suspect that these OUf parameters may not be representative
# of the population, but who knows

# the mean across both model types is 3,350
# OUF = 1,083 (under a half hour)
# OUf = 3,808 (just over an hour)

# our mean will be 2 hours
# again, SDs vary from half the mean, to the mean

#_______________________________________________________________________

# define parameters
lnorm.tau2.lo <- log_norm_params(c(2 %#% "hours", 0.5 %#% "hours"))
lnorm.tau2.hi <- log_norm_params(c(2 %#% "hours", 2 %#% "hours"))

# both means should be 6 hours
tau.2.lo <- dlnorm(x = seq(0, 
                           12 %#% "hours", 
                           length.out = 500),
                   meanlog = lnorm.tau2.lo[1],
                   sdlog = lnorm.tau2.lo[3]) 

tau.2.hi <- dlnorm(x = seq(0, 
                           12 %#% "hours", 
                           length.out = 500),
                   meanlog = lnorm.tau2.hi[1],
                   sdlog = lnorm.tau2.hi[3]) 

# df
tau.2.df <- data.frame(x = seq(0, 
                               12 %#% "hours", 
                               length.out = 500),
                       y = c(tau.2.lo,
                             tau.2.hi),
                       variability = rep(c("low",
                                           "high"),
                                         each = 500))

# plot
ggplot(data = tau.2.df,
       aes(x = x,
           y = y,
           color = variability)) +
  
  theme_bw() +
  
  geom_line(linewidth = 1.1) +
  
  scale_x_continuous(breaks = c(0, 3, 6, 9, 12) %#% "hours",
                     labels = c(0, 3, 6, 9, 12)) +
  
  xlab("Velocity autocorrelation timescale (hr)")

# look at means
mean(rlnorm(n = 500, meanlog = lnorm.tau2.lo[1], sdlog = lnorm.tau2.lo[3])) / 3600
mean(rlnorm(n = 500, meanlog = lnorm.tau2.hi[1], sdlog = lnorm.tau2.hi[3])) / 3600

# look at medians
median(rlnorm(n = 500, meanlog = lnorm.tau2.lo[1], sdlog = lnorm.tau2.lo[3])) / 3600
median(rlnorm(n = 500, meanlog = lnorm.tau2.hi[1], sdlog = lnorm.tau2.hi[3])) / 3600

# and SDs
sd(rlnorm(n = 500, meanlog = lnorm.tau2.lo[1], sdlog = lnorm.tau2.lo[3])) / 3600
sd(rlnorm(n = 500, meanlog = lnorm.tau2.hi[1], sdlog = lnorm.tau2.hi[3])) / 3600

# look like the geometric mean parameterization is a bit more logical

#_______________________________________________________________________
# 6c. Sigma ---- 

# there is a LOT of variability in our real dataset here
# the sigmas control the size of the area used, as well as the aspect ratio

# as mentioned earlier, the mean aspect ratio (major / minor) for our data
# is ~ 7.4 

# it seems reasonable to vary this as well

# the angle should just be a random draw

# we'll define the mean major to be 10,000 m
# SDs will again be half and full mean
# aspect ratio will also be lognormal, mean = 7 and half and full mean SDs

#_______________________________________________________________________

# define parameters
# sigma major
lnorm.sigma.maj.lo <- log_norm_params(c(10000, 2500))
lnorm.sigma.maj.hi <- log_norm_params(c(10000, 10000))

# aspect ratio
lnorm.aspect.lo <- log_norm_params(c(7, 1.75))
lnorm.aspect.hi <- log_norm_params(c(7, 7))

# major
sigma.maj.lo <- dlnorm(x = seq(0, 
                               100000, 
                               length.out = 500),
                       meanlog = lnorm.sigma.maj.lo[1],
                       sdlog = lnorm.sigma.maj.lo[3]) 

sigma.maj.hi <- dlnorm(x = seq(0, 
                               100000, 
                               length.out = 500),
                       meanlog = lnorm.sigma.maj.hi[1],
                       sdlog = lnorm.sigma.maj.hi[3]) 

# df
sigma.maj.df <- data.frame(x = seq(0, 
                                   100000, 
                                   length.out = 500),
                       y = c(sigma.maj.lo,
                             sigma.maj.hi),
                       variability = rep(c("low",
                                           "high"),
                                         each = 500))

# plot
ggplot(data = sigma.maj.df,
       aes(x = x,
           y = y,
           color = variability)) +
  
  theme_bw() +
  
  geom_line(linewidth = 1.1) +
  
  scale_x_continuous(breaks = c(0, 25000, 50000, 75000, 100000)) +
  
  xlab("Asymptotic variance (major axis; m)")

#_______________________________________________________________________
# 7. OU-omega ----

# this model is oscillatory, and perhaps can allow us to model
# multiple home range centroids

#_______________________________________________________________________

# initial simulations just to see what this looks like
# initialize a CTMM object for simulation
ctmm.omega <- ctmm(
  
  tau = c(6 %#% "hours", 6 %#% "hours"), 
  sigma = c(14000, 1000, 0), 
  #omega = (2 * pi) / 36 %#% "hours",
  isotropic = FALSE,
  range = TRUE,
  #sigma = c(1400,
  #          100,
  #          0),
  mu = c(0, 0)         # mean location
  
)

# time
t.omega <- seq(1, 28 * 24 %#% "hours", length.out = 672)

# simulate
sim.omega <- simulate(object = ctmm.omega,
                      t = t.omega,
                      complete = TRUE)

# plot
ggplot(sim.omega,
       aes(x = x,
           y = y)) +
  
  theme_bw() +

  #geom_path(aes(color = timestamp)) +
  
  geom_point(size = 0.05) +
  
  scale_color_viridis_c() +
  
  coord_cartesian(xlim = c(-300, 300),
                  ylim = c(-300, 300))

