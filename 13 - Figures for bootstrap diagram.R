# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 13 - Figures for bootstrap diagram
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 29 Jun 2025
# Date completed: 
# Date last modified: 
# R version: 4.4.3

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation
library(cowplot)

#_______________________________________________________________________
# 2. Simulate datasets ----
#_______________________________________________________________________

sim.data <- rbind(
  
  data.frame("indiv" = 1,
             "speed" = rlnorm(n = 100,
                              meanlog = 1, 
                              sdlog = 0.25)),
  
  data.frame("indiv" = 2,
             "speed" = rlnorm(n = 100,
                              meanlog = 1, 
                              sdlog = 0.25)),
  
  data.frame("indiv" = "all",
             "speed" = rlnorm(n = 800,
                              meanlog = 2.25, 
                              sdlog = 0.25))
  
)

#_______________________________________________________________________
# 3. Create plots ----
#_______________________________________________________________________
# 3a. Indiv 1 ----
#_______________________________________________________________________

ggplot(data = sim.data %>% filter(indiv == "1")) +
  
  theme_minimal() +
  
  geom_histogram(aes(x = speed),
                 bins = 10,
                 fill = "purple3",
                 color = "white",
                 linewidth = 1.25,
                 alpha = 0.65) +
  
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()) -> plot.1

#_______________________________________________________________________
# 3b. Indiv 2 ----
#_______________________________________________________________________

ggplot(data = sim.data %>% filter(indiv == "2")) +
  
  theme_minimal() +
  
  geom_histogram(aes(x = speed),
                 bins = 10,
                 fill = "purple3",
                 color = "white",
                 linewidth = 1.25,
                 alpha = 0.75) +
  
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()) -> plot.2

#_______________________________________________________________________
# 3c. All ----
#_______________________________________________________________________

ggplot(data = sim.data) +
  
  theme_minimal() +
  
  geom_histogram(aes(x = speed),
                 bins = 10,
                 fill = "purple3",
                 color = "white",
                 linewidth = 1.25) +
  
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()) -> plot.all

#_______________________________________________________________________
# 4. Plot together ----
#_______________________________________________________________________

# extra plot to throw in
ggplot(data = sim.data) + theme_minimal() -> extra.plot

plot_grid(plot.1, extra.plot, plot.2,
          nrow = 3,
          rel_heights = c(1, 0.5, 1)) -> together.1

plot_grid(together.1, extra.plot, plot.all,
          nrow = 1,
          rel_widths = c(1, 0.75, 2))

# 690 x 323