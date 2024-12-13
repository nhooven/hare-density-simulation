# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 11 - Fit Bayesian REMs
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 12 Dec 2024
# Date completed: 
# Date last modified: 13 Dec 2024
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation
library(rstan)           # MCMC

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

detec <- read.csv(paste0(getwd(), "/Derived_data/Detections/detec_all.csv"))

cams <- read.csv(paste0(getwd(), "/Derived_data/Detections/cams.csv"))

day.range <- read.csv(paste0(getwd(), "/Derived_data/Model parameters/day_range.csv"))

#_______________________________________________________________________
# 3. Format data correctly ----
#_______________________________________________________________________
# 3a. Define function ----
#_______________________________________________________________________

add_cam_data <- function (ls, var) {
  
  # subset data
  x.1 <- detec %>% 
    
    filter(landscape == ls, variability == var) %>%
    
    dplyr::select(-X)
  
  # identifier
  ls.var.id <- case_when(ls == "simple" & var == "low" ~ "SL",
                         ls == "complex" & var == "low" ~ "CL",
                         ls == "simple" & var == "high" ~ "SH",
                         ls == "complex" & var == "high" ~ "CH")
  
  # add correct camera data
  # ssf (naive)
  ssf.data <- cams %>% dplyr::select(paste0("ssf.", ls.var.id)) %>% rename("ssf" = paste0("ssf.", ls.var.id))
  
  ssf.data <- rbind(ssf.data, ssf.data, ssf.data)
  
  x.1$ssf.data <- ssf.data$ssf
  
  # issf
  issf.data <- cams %>% dplyr::select(paste0("issf.", ls.var.id)) %>% rename("issf" = paste0("issf.", ls.var.id))
  
  issf.data <- rbind(issf.data, issf.data, issf.data)
  
  x.1$issf.data <- issf.data$issf
  
  # add correct day range
  day.range.focal <- day.range %>%
    
    filter(landscape == ls, 
           variability == var) %>%
    
    dplyr::select(day.range)
  
  # add unique "session" identifier (convert to integer later for Stan)
  x.1 <- x.1 %>% 
    
    mutate(session.unique = paste0(ls.var.id, ".", x.1$n.indiv),
           day.range = day.range.focal$day.range[1]) %>%
    
    # drop unused columns
    dplyr::select(-c(n.indiv, landscape, variability))
  
  # return
  return(x.1)
  
}

#_______________________________________________________________________
# 3b. Use function ----
#_______________________________________________________________________

detec.1 <- rbind(add_cam_data("simple", "low"),
                 add_cam_data("complex", "low"),
                 add_cam_data("simple", "high"),
                 add_cam_data("complex", "high"))

#_______________________________________________________________________
# 4. Format for Stan ----
#_______________________________________________________________________

data.stan <- list(
  
  # constants / data lengths
  N_session_station = nrow(detec.1),         # n cams
  N_session = length(unique(detec.1$session.unique)),

  # data - detection
  session = as.integer(factor(detec.1$session.unique,
                              levels = unique(detec.1$session.unique))),       # session integer
  session_station = 1:nrow(detec.1),              # camera identifiers
  
  # data - for REM
  lens = (57.3 * pi) / 180,
  days = 28,
  day_range = detec.1$day.range,  
  
  # detection probability for EDD
  det_prob = 0.30,
  max_dist = 3.5,
  
  # scaled predicted values for use as correction factors (matrix that matches detections)
  ssf_pred = matrix(detec.1$ssf.data,
                    nrow = 9,
                    ncol = 12),
  issf_pred = matrix(detec.1$issf.data,
                     nrow = 9,
                     ncol = 12),
  
  # response variables
  # count of indepedent detections
  count = detec.1$total.detections
  
)

#______________________________________________________________________________
# 4. Tune prior on density ----

# this should be a gamma distribution that's somewhat flat

# keep in mind that our densities range from 0.2 to 1.0 in this first example!

#______________________________________________________________________________

test.dist <- data.frame(x = seq(0, 5, length.out = 1000),
                        y = dgamma(x = seq(0, 5, length.out = 1000),
                                   shape = 1.25,
                                   rate = 1.5))

ggplot(test.dist,
       aes(x = x,
           y = y)) +
  
  theme_bw() +
  
  geom_line(color = "darkblue",
            linewidth = 1.05) +
  
  theme(panel.grid = element_blank())

# let's try this pretty flat version
# ~ gamma(1.25, 1.50)

#_____________________________________________________________________________
# 5. Model fitting ----
#_____________________________________________________________________________
# 5a. Model without correction factors ----
#_____________________________________________________________________________

m1 <- rstan::stan(
  
  file = "model_1.stan",
  data = data.stan,
  chains = 4,
  warmup = 1000,
  iter = 2000
  
)

print(m1, pars = c("D_site",
                   "CV_site"))

rstan::traceplot(m1, pars = c("D_site"))

print(m1, pars = c("D_station"))

#_____________________________________________________________________________
# 5b. Model with post-hoc correction factors ----
#_____________________________________________________________________________

m2 <- rstan::stan(
  
  file = "model_2.stan",
  data = data.stan,
  chains = 4,
  warmup = 1000,
  iter = 2000
  
)

print(m2, pars = c("D_site",
                   "D_site_ssf",
                   "D_site_issf"))

print(m2, pars = c("CV_site",
                   "CV_site_ssf",
                   "CV_site_issf"))

rstan::traceplot(m2, pars = c("D_site_ssf"))

#_____________________________________________________________________________
# 5c. Model with integrated correction factors ----
#_____________________________________________________________________________

m3 <- rstan::stan(
  
  file = "model_3.stan",
  data = data.stan,
  chains = 1,
  warmup = 5000,
  iter = 10000
  
)

print(m3, pars = c("D_site",
                   "D_site_ssf",
                   "D_site_issf"))

print(m3, pars = c("CV_site",
                   "CV_site_ssf",
                   "CV_site_issf"))

rstan::traceplot(m3, pars = c("D_site_issf"))

# multiplying the log(lambda counts) term by the CF seems to lead to a lot of
# divergent transitions

# for now, we'll shelve the Bayesian REMs (they are data-hungry!)
