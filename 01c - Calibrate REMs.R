# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 01c - Calibrate REMs
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
# 2. Purpose ----
#_______________________________________________________________________

# Upon initial simulation runs, our expected counts given the true densities
# needed for the REM are vastly undersetimated, by between ~60-90% so

# I tried everything- quadruple-checking the unit conversion, trying different
# movement model parameters, restricting where home range centers could be,
# all to no avail. Clearly there is a disconnect between the model inputs and
# how the model is function (i.e., what it's measuring)

# Here I will attempt (to the best of my ability) to figure out what's going on
# by using 10 individuals in a 10-ha unit

# density should be somewhat insensitive to sampling area and camera placement
# as long as they're random, right? That's the assumption that SCR works off of

#_______________________________________________________________________
# 3. Read in data ----
#_______________________________________________________________________

# unit boundary
unit.bound <- st_read(paste0(getwd(), "/Derived data/Shapefiles/unit_bound.shp"))

st_crs(unit.bound) <- "epsg:32611"

# camera viewsheds
vs <- st_read(paste0(getwd(), "/Derived data/Shapefiles/cams_vs.shp"))

vs$cam.id <- 1:nrow(vs)

st_crs(vs) <- "epsg:32611"

#_______________________________________________________________________
# 4. Simulation parameters ----

# all camera contact sims will last 2 weeks and will be sampled at a 60-sec "fundamental step"

#_______________________________________________________________________

# simulation duration (days)
sim.duration.day <- 14

sim.duration.sec <- sim.duration.day * 24 * 60 * 60 

# "sampling rate" for fundamental (i.e., straight line) steps
sim.samp.rate <- 60

# time steps to simulate on (from, to, by)
sim.timestep <- seq(1, sim.duration.sec, sim.samp.rate)

# number of iterations (let's stop calling them replicates)
n.iter <- 5

# number of indivs
n.indiv <- 10

# define all individuals and their HRCs
all.indivs <- data.frame(iter = rep(1:n.iter, each = n.indiv),
                         indiv = rep(1:n.indiv, times = n.iter))

all.indivs$mean1 <- runif(nrow(all.indivs), as.numeric(st_bbox(unit.bound)[1]), as.numeric(st_bbox(unit.bound)[3]))
all.indivs$mean2 <- runif(nrow(all.indivs), as.numeric(st_bbox(unit.bound)[1]), as.numeric(st_bbox(unit.bound)[3]))

# ensure that approximate home range area stays the same
# 1:1 (isotropic)
iso.area <- pi * 5000^2

# area of an oval = pi * (major / 2) * (minor / 2)
# 7:1
(maj.7 <- sqrt((iso.area / pi) * 28))
(min.7 <- sqrt((iso.area / pi) * 28) / 7)

(maj.14 <- sqrt((iso.area / pi) * 28 * 2)) 
(min.14 <- sqrt((iso.area / pi) * 28 * 2) / 14) 

# HR crossing time (tau1) should be correlated to sigma major
# from initial simulations, it appears that 10000 sig maj and 12-h tau1 seem to work well together
tau1.base <- 12

tau1.2 <- tau1.base * (maj.7 - 10000) / (10000)     # 64% larger
tau1.3 <- tau1.base * (maj.14 - 10000) / (10000)     # 270% larger

plot(c(10000, maj.7, maj.14), c(tau1.base, tau1.2, tau1.3))

#_______________________________________________________________________
# 5. Function ----
#_______________________________________________________________________

sim_contacts <- function (
    
  tau = c(12 %#% "hours", 1 %#% "hours"),
  sigma = c(10000, 10000)
  
) {
  
  # initialize dfs
all.passes <- data.frame()
all.speeds <- data.frame()
all.relocs <- data.frame()

# loop through iterations
for (i in 1:n.iter) {
  
  for (j in 1:n.indiv) {
    
    # focal row
    focal.row <- all.indivs %>%
      
      filter(iter == i &
             indiv == j)
    
    # define model
    if (sigma[1] == sigma[2]) {
      
      ouf.mod <- ctmm(tau = c(tau[1],             
                            tau[2]),    
                    isotropic = TRUE,                 
                    sigma = c(sigma[1],
                              sigma[2],
                              runif(1, - pi / 2, pi / 2)),                  
                    mu = c(focal.row$mean1, 
                           focal.row$mean2)) 
      
    } else {
      
      ouf.mod <- ctmm(tau = c(tau[1],             
                              tau[2]),    
                      isotropic = TRUE,                 
                      sigma = c(sigma[1],
                                sigma[2],
                                runif(1, - pi / 2, pi / 2)),                  
                      mu = c(focal.row$mean1, 
                             focal.row$mean2)) 
      
    }
    
    # run simulation
    ouf.sim <- simulate(object = ouf.mod, 
                        t = sim.timestep,
                        complete = T)
    
    # coerce telemetry object to sf
    telem.df <- data.frame(t = ouf.sim$t,
                           x = ouf.sim$x,
                           y = ouf.sim$y,
                           timestamp = as.POSIXct(ouf.sim$timestamp))
    
    telem.sf <- st_as_sf(telem.df,
                         coords = c("x", "y"),
                         crs = "epsg:32611")
    
    # tally pass locations
    passes <- st_intersection(vs, 
                              telem.sf) %>%
      
      # drop the geometry
      st_drop_geometry() %>% 
      
      # group by camera
      group_by(cam.id) %>% 
      
      # tally all intersections
      tally() %>%
      
      # and add identifiers
      mutate(iter = focal.row$iter,
             indiv = focal.row$indiv)
    
    # calculate "true" speed (m/s)
    telem.track <- telem.df %>%
      
      make_track(.x = x,
                 .y = y,
                 .t = timestamp) %>%
      
      steps() 
    
    # total speed (m/s)
    total.speed <- sum(telem.track$sl_) / sim.duration.sec
    
    # speeds
    focal.speeds <- data.frame(iter = focal.row$iter,
                               indiv = focal.row$indiv,
                               total.speed = total.speed)
    
    # relocs
    focal.relocs <- telem.df %>%
      
      make_track(.x = x,
                 .y = y,
                 .t = timestamp) %>%
      
      track_resample(rates = hours(1),
                     tolerance = minutes(0)) %>%
      
      mutate(iter = focal.row$iter,
             indiv = focal.row$indiv)
    
    # bind in
    all.passes <- rbind(all.passes, passes)
    all.speeds <- rbind(all.speeds, focal.speeds)
    all.relocs <- rbind(all.relocs, focal.relocs)
    
  }
  
}

# bind into list and return
all.list <- list(all.passes,
                 all.speeds,
                 all.relocs)

return(all.list)
  
}

#_______________________________________________________________________
# 5b. Use function ----
#_______________________________________________________________________

all.passes.1 <- sim_contacts(tau = c(12 %#% "hours", 1 %#% "hours"), sigma = c(10000, 10000))
all.passes.2 <- sim_contacts(tau = c(tau1.2 %#% "hours", 1 %#% "hours"), sigma = c(maj.7, min.7))
all.passes.3 <- sim_contacts(tau = c(tau1.3 %#% "hours", 1 %#% "hours"), sigma = c(maj.14, min.14))
all.passes.4 <- sim_contacts(tau = c(12 %#% "hours", 2 %#% "hours"), sigma = c(10000, 10000))
all.passes.5 <- sim_contacts(tau = c(tau1.2 %#% "hours", 2 %#% "hours"), sigma = c(maj.7, min.7))
all.passes.6 <- sim_contacts(tau = c(tau1.3 %#% "hours", 2 %#% "hours"), sigma = c(maj.14, min.14))
all.passes.7 <- sim_contacts(tau = c(12 %#% "hours", 4 %#% "hours"), sigma = c(10000, 10000))
all.passes.8 <- sim_contacts(tau = c(tau1.2 %#% "hours", 4 %#% "hours"), sigma = c(maj.7, min.7))
all.passes.9 <- sim_contacts(tau = c(tau1.3 %#% "hours", 4 %#% "hours"), sigma = c(maj.14, min.14))

#_______________________________________________________________________
# 6. Aggregate passes by camera ----

# prepare data
all.passes.1[[1]]$scenario <- 1
all.passes.2[[1]]$scenario <- 2
all.passes.3[[1]]$scenario <- 3
all.passes.4[[1]]$scenario <- 4
all.passes.5[[1]]$scenario <- 5
all.passes.6[[1]]$scenario <- 6
all.passes.7[[1]]$scenario <- 7
all.passes.8[[1]]$scenario <- 8
all.passes.9[[1]]$scenario <- 9

# bind together
all.passes.bind <- rbind(all.passes.1[[1]],
                         all.passes.2[[1]],
                         all.passes.3[[1]],
                         all.passes.4[[1]],
                         all.passes.5[[1]],
                         all.passes.6[[1]],
                         all.passes.7[[1]],
                         all.passes.8[[1]],
                         all.passes.9[[1]])

#_______________________________________________________________________

all.passes.bind.1 <- all.passes.bind %>% 
  
  group_by(scenario,
           iter,
           cam.id) %>%
  
  summarize(passes = sum(n))

#_______________________________________________________________________
# 7. Add "zero" count rows where necessary ----
#_______________________________________________________________________

# write function
add_zero <- function (df) {
  
    possible.grid <- expand.grid(scenario = 1:9,
                                 iter = 1:max(df$iter),
                                 cam.id = 1:9)
    
    # loop through aggregated df
    zero.rows <- data.frame()
    
    for (i in 1:nrow(possible.grid)) {
      
      focal.row <- possible.grid %>% slice(i)
      
      focal.row.match <- plyr::match_df(df, focal.row)
      
      # does this row exist?
      if (nrow(focal.row.match) == 0) {
        
        # then bind in a new row
        zero.rows <- rbind(zero.rows,
                           focal.row %>% mutate(passes = 0))
        
      }
      
    }
    
    # bind in and sort
    df.1 <- df %>%
      
      bind_rows(zero.rows) %>%
      
      arrange(iter,
              cam.id)
  
  # return
  return(df.1)
  
}

# use function
all.passes.bind.2 <- add_zero(all.passes.bind.1)

#_______________________________________________________________________
# 8. Mean speeds ----
#_______________________________________________________________________

# prepare data
all.passes.1[[2]]$scenario <- 1
all.passes.2[[2]]$scenario <- 2
all.passes.3[[2]]$scenario <- 3
all.passes.4[[2]]$scenario <- 4
all.passes.5[[2]]$scenario <- 5
all.passes.6[[2]]$scenario <- 6
all.passes.7[[2]]$scenario <- 7
all.passes.8[[2]]$scenario <- 8
all.passes.9[[2]]$scenario <- 9

# bind together
all.speeds.bind <- rbind(all.passes.1[[2]],
                         all.passes.2[[2]],
                         all.passes.3[[2]],
                         all.passes.4[[2]],
                         all.passes.5[[2]],
                         all.passes.6[[2]],
                         all.passes.7[[2]],
                         all.passes.8[[2]],
                         all.passes.9[[2]])

# mean and SD speeds
all.speeds.bind.1 <- all.speeds.bind %>% 
  
  group_by(scenario,
           iter) %>%
  
  summarize(mean.speed = mean(total.speed),
            sd.speed = sd(total.speed))

#_______________________________________________________________________
# 9. Examine detections to see if they seem reasonable ----

# here we'll see what the expected counts would be, given a mean speed

# equation 3 from Rowcliffe et al. (2008):

# y = ((2 + theta) / pi) * r * t * v * D 

# all in km and days

# lens in radians
lens.rad <- (57.3 * (pi / 180))

#_______________________________________________________________________

all.passes.bind.3 <- all.passes.bind.2 %>%
  
  group_by(scenario,
           iter) %>%
  
  summarize(mean.passes = mean(passes)) %>%
  
  left_join(all.speeds.bind.1) %>%
  
  # REM point estimate
  mutate(mean.REM = (mean.passes / 28) * (pi / ((mean.speed * (86400 / 1000)) * (3.5 / 1000) * (2 + lens.rad))) / 100) %>%
  
  # percent bias and space use parameter IDs
  mutate(perc.bias = (mean.REM - 1.0) / 1.0,
         hr.crossing.time = factor(case_when(scenario %in% c(1, 4, 7) ~ "12 h",
                                             scenario %in% c(2, 5, 8) ~ "20 h",
                                             scenario %in% c(3, 6, 9) ~ "33 h"),
                                   levels = c("12 h", "20 h", "33 h")),
         velo.autocorr.t = factor(case_when(scenario %in% c(1, 2, 3) ~ "1 h",
                                            scenario %in% c(4, 5, 6) ~ "2 h",
                                            scenario %in% c(7, 8, 9) ~ "4 h"),
                                   levels = c("1 h", "2 h", "4 h")),
         aspect.ratio = factor(case_when(scenario %in% c(1, 4, 7) ~ "1:1",
                                         scenario %in% c(2, 5, 8) ~ "7:1",
                                         scenario %in% c(3, 6, 9) ~ "14:1"),
                               levels = c("1:1", "7:1", "14:1")))

# plot
ggplot(data = all.passes.bind.3,
       aes(x = velo.autocorr.t,
           y = perc.bias)) +
  
  facet_wrap(~ aspect.ratio) +
  
  theme_bw() +
  
  geom_hline(yintercept = 0) +
  
  geom_point(size = 2,
             shape = 21,
             color = "black",
             aes(fill = perc.bias)) +
  
  scale_fill_gradient2(low = "darkblue",
                       high = "yellow") +
  
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  
  xlab("Velocity autocorrelation timescale") +
  ylab("Percent bias in density") +
  
  scale_y_continuous(breaks = seq(-1, 3.5, 0.5))

#_______________________________________________________________________
# 9. Plot space use areas ----
#_______________________________________________________________________
# 9a. Relocations ----
#_______________________________________________________________________

# prepare data
all.passes.1[[3]]$scenario <- 1
all.passes.2[[3]]$scenario <- 2
all.passes.3[[3]]$scenario <- 3
all.passes.4[[3]]$scenario <- 4
all.passes.5[[3]]$scenario <- 5
all.passes.6[[3]]$scenario <- 6
all.passes.7[[3]]$scenario <- 7
all.passes.8[[3]]$scenario <- 8
all.passes.9[[3]]$scenario <- 9

# bind together
all.relocs.bind <- rbind(all.passes.1[[3]],
                         all.passes.2[[3]],
                         all.passes.3[[3]],
                         all.passes.4[[3]],
                         all.passes.5[[3]],
                         all.passes.6[[3]],
                         all.passes.7[[3]],
                         all.passes.8[[3]],
                         all.passes.9[[3]])

# add IDs
all.relocs.bind.1 <- all.relocs.bind %>%
  
mutate(scenario = factor(scenario,
                         levels = 1:9,
                         labels = c("Tau2 = 1h, AR = 1:1",
                                    "Tau2 = 1h, AR = 7:1",
                                    "Tau2 = 1h, AR = 14:1",
                                    "Tau2 = 2h, AR = 1:1",
                                    "Tau2 = 2h, AR = 7:1",
                                    "Tau2 = 2h, AR = 14:1",
                                    "Tau2 = 4h, AR = 1:1",
                                    "Tau2 = 4h, AR = 7:1",
                                    "Tau2 = 4h, AR = 14:1")))
  
  

#_______________________________________________________________________
# 9b. Convert to sf ----
#_______________________________________________________________________

all.indivs.sf <- all.indivs %>%
  
  st_as_sf(coords = c("mean1",
                      "mean2"),
           crs = "EPSG:32611")

all.relocs.sf <- all.relocs.bind.1 %>%
  
  st_as_sf(coords = c("x_",
                      "y_"),
           crs = "EPSG:32611")

#_______________________________________________________________________
# 9c. Facetted plot ----
#_______________________________________________________________________

ggplot() +
  
  theme_bw() +
  
  geom_sf(data = unit.bound) +
  
  geom_sf(data = vs) +
  
  geom_sf(data = all.relocs.sf,
          alpha = 0.15,
          size = 0.005,
          aes(color = as.factor(indiv))) +
  
  geom_sf(data = all.indivs.sf,
          size = 0.75) +
  
  facet_grid(iter ~ scenario) +
  
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 8, 
                                  color = "white"),
        strip.background.x = element_rect(fill = "black",
                                          color = "white"),
        strip.background.y = element_blank(),
        panel.spacing = unit(0, "in")) +
  
  scale_color_viridis_d()
