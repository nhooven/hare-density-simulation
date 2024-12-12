# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 11 - Fit Bayesian REMs
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 12 Dec 2024
# Date completed: 
# Date last modified: 12 Dec 2024
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
  day_range = detec.1$day.range * 1000,  
  
  # detection probability for EDD
  det_prob = 0.30,
  max_dist = 3.5,
  
  # scaled predicted values for use as correction factors
  ssf_pred = detec.1$ssf.data,
  issf_pred = detec.1$issf.data,
  
  # response variables
  # count of indepedent detections
  count = detec.1$total.detections
  
)


#_______________________________________________________________________________________________
# 4. Model fitting ----
#_______________________________________________________________________________________________

# 12 Dec 2024
# The D prior is very important since there isn't much data here!

m1 <- rstan::stan(
  
  file = "model.stan",
  data = data.stan,
  chains = 4,
  warmup = 1000,
  iter = 2000
  
)

print(m1, pars = c("D_site",
                   "CV_site"))

rstan::traceplot(m1, pars = c("D_site"))

print(m1, pars = c("D_station"))
