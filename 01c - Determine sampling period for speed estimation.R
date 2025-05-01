# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 01c - Determine sampling period for speed estimation
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 01 May 2025
# Date completed: 01 May 2025
# Date last modified: 01 May 2025
# R version: 4.4.3

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation
library(lubridate)            # work with dates

#_______________________________________________________________________
# 2. Directories ----
#_______________________________________________________________________

dir.pre <- "E:/Hare project/Data analysis/GPS processing/Derived data/Cleaned data 2/PRE/"
dir.dur <- "E:/Hare project/Data analysis/GPS processing/Derived data/Cleaned data 2/DUR/"
dir.post <- "E:/Hare project/Data analysis/GPS processing/Derived data/Cleaned data 2/POST/"

# file lists
list.pre <- list.files(dir.pre)
list.dur <- list.files(dir.dur)
list.post <- list.files(dir.post)

# subset file lists
list.pre.1 <- list.pre[c(2, 10, 11, 19, 21, 39, 43, 44, 54)]
list.dur.1 <- list.dur[c(1, 2, 6, 8, 16, 20)]
list.post.1 <- list.post[c(4, 9, 13, 14, 17, 21, 23, 26, 34)]

#_______________________________________________________________________
# 3. Read in data ----

# Here we'll read everything in and pack into one df for summary and visualization

#_______________________________________________________________________

# blank df
all.relocs <- data.frame()

for (i in 1:3) {
  
  if (i == 1) {
    
    focal.dir <- dir.pre
    focal.list <- list.pre.1
    
  }
  
  if (i == 2) {
    
    focal.dir <- dir.dur
    focal.list <- list.dur.1
    
  }
  
  if (i == 3) {
    
    focal.dir <- dir.post
    focal.list <- list.post.1
    
  }
  
  # loop through each one
  for (j in 1:length(focal.list)) {
    
    focal.data <- read.csv(paste0(focal.dir,
                                  focal.list[j]))
    
    # bind in
    all.relocs <- rbind(all.relocs, focal.data)
    
  }
  
}

#_______________________________________________________________________
# 4. Visualize ----

# here we'll aggregate by hour of the day
all.relocs.1 <- all.relocs %>%
  
  mutate(hour = as.integer(substr(timestamp, 12, 13)),
         indiv.deploy = paste0(all.relocs.1$indiv, "_", all.relocs.1$deploy))

#_______________________________________________________________________

# pooled histogram
ggplot(data = all.relocs.1,
       aes(x = hour)) +
  
  theme_bw() +
  
  geom_histogram(bins = 24,
                 color = "black",
                 aes(fill = after_stat(count))) +
  
  geom_vline(xintercept = 8,
             linetype = "dashed") +
  
  geom_vline(xintercept = 15,
             linetype = "dashed") +
  
  scale_fill_viridis_c(option = "plasma") +
  
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  
  xlab("Hour of day") +
  
  scale_x_continuous(breaks = seq(0, 23, 2))

# facetted histogram
ggplot(data = all.relocs.1,
       aes(x = hour)) +
  
  facet_wrap(~ indiv.deploy,
             nrow = 4) +
  
  theme_bw() +
  
  geom_histogram(bins = 24,
                 aes(fill = after_stat(count),
                     color = after_stat(count))) +
  
  geom_vline(xintercept = 8,
             linetype = "dashed") +
  
  geom_vline(xintercept = 15,
             linetype = "dashed") +
  
  scale_color_viridis_c(option = "plasma") +
  scale_fill_viridis_c(option = "plasma") +
  
  theme(panel.grid = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(0, "lines"),
        axis.text = element_text(size = 6,
                                 color = "black")) +
  
  xlab("Hour of day") +
  
  scale_x_continuous(breaks = seq(0, 23, 4))

#_______________________________________________________________________
# 5. Summarize ----
#_______________________________________________________________________

# percent of total fixes
all.relocs.summary <- all.relocs.1 %>%
  
  group_by(hour) %>%
  
  summarize(perc.total = n() / nrow(all.relocs.1))

all.relocs.summary

1 - sum(all.relocs.summary$perc.total[all.relocs.summary$hour %in% c(8:15)])

