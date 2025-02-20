# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 03 - Initial movement simulations
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 15 Nov 2024
# Date completed: 21 Nov 2024
# Date last modified: 20 Feb 2025
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation
library(sf)              # read in shapefiles
library(terra)           # rasters
library(amt)             # simulate tracks
library(lubridate)       # work with time

#_______________________________________________________________________
# 2. Read in rasters and unit boundary ----
#_______________________________________________________________________

# covariate rasters
B1 <- rast(paste0(getwd(), "/Rasters/B1.tif"))
B2 <- rast(paste0(getwd(), "/Rasters/B2.tif"))
B3 <- rast(paste0(getwd(), "/Rasters/B3.tif"))

A1 <- rast(paste0(getwd(), "/Rasters/A1.tif"))
A2 <- rast(paste0(getwd(), "/Rasters/A2.tif"))
A3 <- rast(paste0(getwd(), "/Rasters/A3.tif"))

# unit bound
unit.bound <- st_read(paste0(getwd(), "/Derived_data/Shapefiles/unit_bound.shp"))

# scale rasters
landscape.covs.B <- list(scale(B1),
                         scale(B2),
                         scale(B3))

landscape.covs.A <- list(scale(A1),
                         scale(A2),
                         scale(A3))

# naive RSF surfaces
B.rsf <- rast(paste0(getwd(), "/Rasters/B_rsf.tif"))
A.rsf <- rast(paste0(getwd(), "/Rasters/A_rsf.tif"))

#_______________________________________________________________________
# 3. Define simulation parameters ----
#_______________________________________________________________________
# 3a. Movement parameter distributions ----
#_______________________________________________________________________

# step lengths (gamma)
sl.dist <- make_gamma_distr(shape = 1,
                            scale = 50)

# turning angles (uniform)
ta.dist <- make_vonmises_distr(kappa = 1.5)

#_______________________________________________________________________
# 3b. Habitat selection coefficients ----
#_______________________________________________________________________

# standardized!!!

coef.fora.sl <- -0.05       # β1 - (start) fora and speed interaction = negative
coef.fora <- 1.5            # β2 = (end) fora selection = positive
coef.fora.ta <- -0.5        # β3 = (end) fora and concentration interaction = negative       
coef.elev <- 0.75           # β4 = (end) linear elev selection = positive
coef.elev2 <- -0.75         # β5 = (end) squared elev selection = negative
coef.open.sl <- 0.25        # β6 = (start) open and speed interaction = positive
coef.open <- -1.0           # β7 = (end) open selection = negative

#_______________________________________________________________________
# 3c. Home ranging parameters ----
#_______________________________________________________________________

# bivariate normal variance
# we want this to reflect a decently large home range (~ 2 ha)
library(MASS)
library(adehabitatHR)

# simulate bivariate normal data (assume no covariance)
bivariate_data <- as.data.frame(mvrnorm(n = 337,
                                        mu = c(0, 0),
                                        Sigma = matrix(c(2000, 0, 
                                                         0, 2000), ncol = 2)))

bv.sp <- SpatialPoints(coords = data.frame(bivariate_data$V1,
                                           bivariate_data$V2),
                       proj4string = CRS("EPSG:32611"))

bv.mcp <- mcp(bv.sp, percent = 90, unin = "m", unout = "ha")

bv.mcp

plot(bv.sp)
plot(bv.mcp, add = T)

# 2000 seems like a good compromise

# the Bx and By coefficients are then calculated from the "home range" centroid and the variance
e.var <- 2000

#_______________________________________________________________________
# 3d. Redistribution kernel parameters ----
#_______________________________________________________________________

# control steps
rk.control <- 100

# tolerance outside the landscape (hopefully we won't have to deal with this much)
rk.tolerance <- 0.01    # 1%

#_______________________________________________________________________
# 4. Home range centers ----

# n of replicates for each landscape replicate
n.reps <- 100

# here we'll draw 100 locations from our MFHK BEFORE rasters

#_______________________________________________________________________
# 4a. Draw in proportion to the "RSF" ----
#_______________________________________________________________________

# list so we can easily subset
hrc <- list()

hrc[[1]] <- st_as_sf(spatSample(B.rsf$B1, 
                                size = n.reps, 
                                method = "weights",
                                as.points = TRUE,
                                values = FALSE))

hrc[[2]] <- st_as_sf(spatSample(B.rsf$B2, 
                                size = n.reps, 
                                method = "weights",
                                as.points = TRUE,
                                values = FALSE))

hrc[[3]] <- st_as_sf(spatSample(B.rsf$B3, 
                                size = n.reps, 
                                method = "weights",
                                as.points = TRUE,
                                values = FALSE))

#_______________________________________________________________________
# 4b. Plot ----
#_______________________________________________________________________

plot(B.rsf$B1)
plot(st_geometry(hrc[[1]]), add = T)

plot(B.rsf$B2)
plot(st_geometry(hrc[[2]]), add = T)

plot(B.rsf$B3)
plot(st_geometry(hrc[[3]]), add = T)

# 500 x 500

#_______________________________________________________________________
# 4c. Save point shapefiles ----
#_______________________________________________________________________

st_write(hrc[[1]],
         dsn = paste0(getwd(), "/Derived_data/Shapefiles/B1_hrc.shp"),
         layer = "B1_hrc.shp",
         append = FALSE)

st_write(hrc[[2]],
         dsn = paste0(getwd(), "/Derived_data/Shapefiles/B2_hrc.shp"),
         layer = "B2_hrc.shp",
         append = FALSE)

st_write(hrc[[3]],
         dsn = paste0(getwd(), "/Derived_data/Shapefiles/B3_hrc.shp"),
         layer = "B3_hrc.shp",
         append = FALSE)

#_______________________________________________________________________
# 5. Define fuction that calculates home ranging parameters ----
#_______________________________________________________________________

hr_params <- function(e.var = e.var,          # expected variance of the bivariate normal
                      focal.hrc = focal.hrc)  # home range centroid as previously drawn
  
  {
  
  # x2 + y2 coefficient
  b.x2y2 <- -1 / e.var
  
  # solve for x and y coefficients
  b.x <- focal.hrc[1] * 2 * -b.x2y2 
  b.y <- focal.hrc[2] * 2 * -b.x2y2
  
  hr.params <- c(b.x, b.y, b.x2y2)
  
  return(hr.params)
  
}

#_______________________________________________________________________
# 6. Run simulations iteratively ----
#_______________________________________________________________________
# 6a. Define function ----
#_______________________________________________________________________

init_sim <- function(id.rep) {
  
  # extract landscapes
  focal.landscape.B <- landscape.covs.B[[id.rep]]
  focal.landscape.A <- landscape.covs.A[[id.rep]]
  
  # choose correct AFTER rsf surface
  if (id.rep == 1) {
    
    focal.A.rsf <- A.rsf$A1
    
  } else {
    
    if (id.rep == 2) {
      
      focal.A.rsf <- A.rsf$A2
      
    } else {
      
      focal.A.rsf <- A.rsf$A3
      
    }
    
  }
   
  # loop through each replicate
  sims.df <- data.frame()
  
  start.time <- Sys.time()
  
  for (i in 1:n.reps) {
    
    # define hrc - pull coordinates from the first HRC
    focal.hrc <- st_coordinates(hrc[[id.rep]][i, ])
    
    # define start step (add a little noise to the HRC)
    start.step <- make_start(x = c(focal.hrc[1] + rnorm(n = 1, mean = 0, sd = 50),
                                   focal.hrc[2] + rnorm(n = 1, mean = 0, sd = 50)),
                             ta_ = 0,
                             time = ymd_hm("2024-09-01 18:00", 
                                           tz = "America/Los_Angeles"),
                             dt = hours(2),
                             crs = crs("EPSG:32611"))
    
    # calculate home ranging parameters
    hr.params <- hr_params(e.var, focal.hrc)
    
    # make iSSF model
    # here the terms are important to get right so redistribution_kernel() works okay
    issf.model <- make_issf_model(coefs = c("fora_start:log(sl_)" = coef.fora.sl, 
                                            "fora_end" = coef.fora,
                                            "elev_end" = coef.elev,
                                            "I(elev_end^2)" = coef.elev2,
                                            "log(sl_):open_start" = coef.open.sl,
                                            "open_end" = coef.open,
                                            x2_ = hr.params[1],
                                            y2_ = hr.params[2], 
                                            "I(x2_^2 + y2_^2)" = hr.params[3]),              
                                  sl = sl.dist,
                                  ta = ta.dist)
    
    # BEFORE simulations
    # initialize redistribution kernel
    rk.B <- redistribution_kernel(x = issf.model,
                                  start = start.step,
                                  map = focal.landscape.B,     # use correct landscape here
                                  n.control = rk.control,   
                                  max.dist = get_max_dist(issf.model),
                                  tolerance.outside = rk.tolerance)
    
    # run simulation
    sim.path.B <- simulate_path(rk.B,
                                n.steps = 336,
                                start = start.step,
                                verbose = TRUE)
    
    # extract endpoint for starting point of AFTER simulation
    end.step <- make_start(x = c(sim.path.B$x_[337],
                                 sim.path.B$y_[337]),
                           ta_ = 0,
                           time = sim.path.B$t_[337],
                           dt = hours(2),
                           crs = crs("EPSG:32611"))
    
    # shift HRC from 25 nearby locations (i.e., within a 100-m buffer)
    nearby.samp <- st_as_sf(st_sample(st_buffer(hrc[[id.rep]][i, ], 
                                                dist = 100),
                                      size = 25))
    
    nearby.rsf <- extract(focal.A.rsf, 
                          nearby.samp)
    
    # change names
    names(nearby.rsf) <- c("ID", "RSF")
    
    # remove any rows with NAs
    nearby.rsf <- nearby.rsf[is.na(nearby.rsf$RSF) == FALSE, ]
    
    tentative.new.hrc <- sample(nearby.rsf$ID, 1, prob = nearby.rsf$RSF)
    
    # only use the new one if the RSF score is higher
    if (nearby.rsf$RSF[nearby.rsf$ID == tentative.new.hrc] > extract(focal.A.rsf, hrc[[id.rep]][i, ])[ , 2]) {
      
      new.hrc <- st_coordinates(nearby.samp[tentative.new.hrc, ])
      
    } else {
      
      new.hrc <- focal.hrc
      
    }
    
    # AFTER simulations
    # calculate home ranging parameters
    hr.params.A <- hr_params(e.var, new.hrc)
    
    # make iSSF model
    # here the terms are important to get right so redistribution_kernel() works okay
    issf.model.A <- make_issf_model(coefs = c("fora_start:log(sl_)" = coef.fora.sl, 
                                              "fora_end" = coef.fora,
                                              "elev_end" = coef.elev,
                                              "I(elev_end^2)" = coef.elev2,
                                              "log(sl_):open_start" = coef.open.sl,
                                              "open_end" = coef.open,
                                              x2_ = hr.params.A[1],
                                              y2_ = hr.params.A[2], 
                                              "I(x2_^2 + y2_^2)" = hr.params.A[3]),              
                                  sl = sl.dist,
                                  ta = ta.dist)
  
    # initialize redistribution kernel
    rk.A <- redistribution_kernel(x = issf.model.A,            # updated issf parameters
                                  start = end.step,            # use end step
                                  map = focal.landscape.A,     # use correct landscape here
                                  n.control = rk.control,   
                                  max.dist = get_max_dist(issf.model.A),
                                  tolerance.outside = rk.tolerance)
    
    # run simulation
    sim.path.A <- simulate_path(rk.A,
                                n.steps = 336,
                                start = end.step,
                                verbose = TRUE)
    
    # bind together
    sim.path.B$trt <- "before"
    sim.path.A$trt <- "after"
    
    sim.path <- rbind(sim.path.B,
                      sim.path.A)
    
    sim.path.1 <- sim.path %>% 
      
      mutate(indiv = i,
             rep = id.rep)
    
    # bind to df
    sims.df <- rbind(sims.df, sim.path.1)
    
    # status message (every 10 iterations)
    if (i %% 10 == 0) {
    
    elapsed.time <- round(as.numeric(difftime(Sys.time(), 
                                              start.time, 
                                              units = "mins")), 
                          digits = 1)
    
    print(paste0("Completed path ", i, " of ", n.reps, " - ", elapsed.time, " mins"))
    
    }
  
  }
  
  # return
  return(sims.df)

}

#_______________________________________________________________________
# 5c. Use function ----
#_______________________________________________________________________

init.sim.1 <- init_sim(1)
init.sim.2 <- init_sim(2)
init.sim.3 <- init_sim(3)

#_______________________________________________________________________
# 6. Plot tracks ----
#_______________________________________________________________________

# bind together
init.sim.all <- rbind(init.sim.1, init.sim.2, init.sim.3)

# factor level
init.sim.all$trt <- factor(init.sim.all$trt,
                           levels = c("before", "after"))

# plot paths
ggplot() +
  
  theme_bw() +
  
  # unit boundary as a black square
  geom_sf(data = unit.bound,
          fill = NA) +
  
  facet_grid(rep ~ trt) +
  
  # paths
  geom_point(data = init.sim.all,
            aes(x = x_,
                y = y_,
                color = trt),
            alpha = 0.05,
            size = 0.5) +
  
  scale_color_manual(values = c("purple", "orange")) +
  
  theme(legend.position = "none") +
  
  # coordinate system to make nice axis labels
  coord_sf(datum = st_crs(32611))

#_______________________________________________________________________
# 7. Write to file ----
#_______________________________________________________________________

# all simulated locations
write.csv(init.sim.all, paste0(getwd(), "/Derived_data/Simulated data/init_sims.csv"))
