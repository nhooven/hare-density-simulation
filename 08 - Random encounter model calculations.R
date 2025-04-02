# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 08 - Random encounter model calculations
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

# detections
detec.lo1 <- read.csv(paste0(getwd(), "/Derived data/Aggregated detections/detections_lo1.csv"))

# speeds
speeds.lo1 <- read.csv(paste0(getwd(), "/Derived data/Speed distributions/speeds_lo1.csv"))

#_______________________________________________________________________
# 3. Combinations to loop through ----
#_______________________________________________________________________

all.combos <- expand.grid(iter = 1,
                          group = unique(detec.lo1$group),
                          scenario = 1:12)

#_______________________________________________________________________
# 4. Monte Carlo REM calculation ----

# here, for each combination, we'll take a random draw from the cameras for detection
# and a random draw from the speed distribution, calculating the mean and quantiles for CIs
# on the sampled distribution of REM densities

#_______________________________________________________________________

all.rem <- data.frame()

for (i in 1:nrow(all.combos)) {
  
  focal.combo <- all.combos[i, ]
  
  # subset detec
  detec.i <- detec.lo1 %>% 
    
    filter(iter == focal.combo$iter,
           group == focal.combo$group)
  
  # subset speeds
  speeds.i <- speeds.lo1[ , focal.combo$scenario + 1]  # make sure the first row from the csv isn't counted
  
  # draw MC samples (5000 is fine for now)
  mc.rem.draws <- matrix(data = NA,
                         nrow = 5000)
  
  for (j in 1:5000) {
    
    # take random draws
    rand.count <- sample(detec.i$detections, size = 1)
    rand.speed <- sample(speeds.i, size = 1)
    
    # calculate and bind in REM
    mc.rem.draws[j, ] <- (rand.count / 28) * # constant for now
                         (pi / (rand.speed * (3.5 / 1000) * (2.0 + ((57.3 * pi) / 180)))) /
                         100   # per hectare
    
  }
  
  # calculate mean, SD, CV, and quantiles
  mean.D <- mean(mc.rem.draws, na.rm = T)
  sd.D <- sd(mc.rem.draws, na.rm = T)
  cv.D <- mean.D / sd.D
  lo.D <- quantile(mc.rem.draws, probs = 0.0275, na.rm = T)
  hi.D <- quantile(mc.rem.draws, probs = 0.975, na.rm = T)
  
  # bind in
  all.rem <- rbind(all.rem,
                   data.frame(iter = focal.combo$iter,
                              group = focal.combo$group,
                              scenario = focal.combo$scenario,
                              mean.D = mean.D,
                              sd.D = sd.D,
                              cv.D = cv.D,
                              lo.D = lo.D,
                              hi.D = hi.D))
  
}

#_______________________________________________________________________
# 5. Examine bias ----
#_______________________________________________________________________

# calculate bias
all.rem.1 <- all.rem %>%
  
  mutate(true.D = group / 10) %>%
  
  mutate(perc.bias = ((mean.D - true.D) / true.D) * 100) %>%
  
  # filter out group 2
  filter(group != 2) %>%
  
  # factor labels
  mutate(scenario = factor(scenario,
                       levels = c("1", "2", "3", "4",
                                  "5", "6", "7", "8",
                                  "9", "10", "11", "12"),
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
ggplot(all.rem.1,
       aes(y = as.factor(true.D),
           x = perc.bias)) +
  
  theme_bw() +
  
  facet_wrap(~ scenario,
             nrow = 4) +
  
  geom_point(size = 2,
             aes(color = perc.bias)) +
  
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank()) +
  
  scale_color_viridis_c() +
  
  xlab("Percent bias")

# everything is underestimated - guessing that the detection threshold is messing me up here