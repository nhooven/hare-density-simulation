# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 07 - Speed distributions
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 02 Apr 2025
# Date completed: 
# Date last modified: 02 Apr 2025
# R version: 4.4.3

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

speeds.lo1 <- read.csv(paste0(getwd(), "/Derived data/CTMM speeds/speeds_lo1.csv"))

#_______________________________________________________________________
# 3. Sampling distributions ----

# here we'll assume a normal sampling distribution (in the future we're do this non-parametrically)

#_______________________________________________________________________

speeds.lo1.1 <- speeds.lo1 %>%
  
  # back-calculate the SD of the sampling dist. (i.e., the SE)
  # we'll use the mean absolute difference from the mean to the 95% CIs
  mutate(parametric.se = (((speed.mean.hi - speed.mean.est) +
                          (speed.mean.est - speed.mean.lo)) / 2) / 1.96)

#_______________________________________________________________________
# 4. Create aggregated speed distributions ----

# here we want to create a distribution of many (say 1500) observations that
# could have come from the sampling distributions of each scenario
# each individual will get 100 samples for now

# this will yield a matrix with 1500 rows and 12 cols

#_______________________________________________________________________

# create blank matrix
speed.samples <- matrix(data = NA,
                        nrow = 1500,
                        ncol = 12) 

# loop through each scenario
for (i in 1:12) {
  
  speeds.scenario <- speeds.lo1.1 %>%
    
    filter(scenario == i)
  
  # loop through each individual
  speeds.scenario.dist <- matrix(data = NA,
                                 nrow = 1500,
                                 ncol = 1)
  
  for (j in 1:length(unique(speeds.scenario$ID))) {
    
    speeds.indiv <- speeds.scenario %>%
      
      filter(ID == unique(speeds.scenario$ID)[j])
    
    # bind into matrix
    speeds.scenario.dist[1:100 + 100 * (j - 1), ] <- rnorm(n = 100,
                                                           mean = speeds.indiv$speed.mean.est,
                                                           sd = speeds.indiv$parametric.se)
    
  }
    
    # bind into full matrix
    speed.samples[ , i] <- speeds.scenario.dist
  
}

#_______________________________________________________________________
# 5. Examine distributions ----
#_______________________________________________________________________

# pivot
speed.samples.df <- as.data.frame(speed.samples) %>%
  
  pivot_longer(cols = V1:V12) %>%
  
  mutate(name = factor(name,
                       levels = c("V1", "V2", "V3", "V4",
                                  "V5", "V6", "V7", "V8",
                                  "V9", "V10", "V11", "V12"),
                       labels = c("1 hr, 100%, 4 wk",
                                  "2 hr, 100%, 4 wk",
                                  "4 hr, 100%, 4 wk",
                                  "1 hr, 80%, 4 wk",
                                  "2 hr, 80%, 4 wk",
                                  "4 hr, 80%, 4 wk",
                                  "1 hr, 100%, 2 wk",
                                  "2 hr, 100%, 2 wk",
                                  "4 hr, 100%, 2 wk",
                                  "1 hr, 80%, 2 wk",
                                  "2 hr, 80%, 2 wk",
                                  "4 hr, 80%, 2 wk")))

# plot
ggplot(data = speed.samples.df,
       aes(x = value)) +
  
  theme_bw() +
  
  facet_wrap(~ name,
             nrow = 4) +
  
  geom_histogram(color = "black",
                 fill = "white") +
  
  theme(panel.grid = element_blank())

#_______________________________________________________________________
# 6. Convert from m/s to km/day ----
#_______________________________________________________________________

# conversion factor is 1000 / 86400 ~ 0.0115574
speed.samples.day <- speed.samples / (1000 / 86400)

#_______________________________________________________________________
# 7. Write to file ----
#_______________________________________________________________________

write.csv(speed.samples.day, paste0(getwd(), "/Derived data/Speed distributions/speeds_lo1.csv"))
