# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 05 - Prepare detection data
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 26 Mar 2025
# Date completed: 
# Date last modified: 01 Apr 2025
# R version: 4.4.3

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

# low
detections.lo1 <- read.csv(paste0(getwd(), "/Derived data/Camera detections/detections_lo1.csv"))

# hi
detections.hi1 <- read.csv(paste0(getwd(), "/Derived data/Camera detections/detections_hi1.csv"))

#_______________________________________________________________________
# 3. Split into each abundance group ----

# individuals in each group
group.25 <- 1:25
group.10 <- 1:10
group.05 <- 1:5
group.02 <- 1:2

# write function
split_abund <- function (df) {
  
  # pack subsetted dfs into a df
  each.abund.list <- rbind(df %>% mutate(group = 50),
                           df %>% filter(indiv %in% group.25) %>% mutate(group = 25),
                           df %>% filter(indiv %in% group.10) %>% mutate(group = 10),
                           df %>% filter(indiv %in% group.05) %>% mutate(group = 05),
                           df %>% filter(indiv %in% group.02) %>% mutate(group = 02))
  
  return(each.abund.list)
  
}

#_______________________________________________________________________
# 3a. Use function ----
#_______________________________________________________________________

lo1.group <- split_abund(detections.lo1)

#_______________________________________________________________________
# 4. Aggregate detections by camera ----
#_______________________________________________________________________
# 4a. Q1 (low variability) ----
#_______________________________________________________________________

# fold 1
q1.agg1 <- lo1.group %>% 
  
  group_by(iter,
           group,
           cam.id) %>%
  
  summarize(detections = n())

#_______________________________________________________________________
# 4b. Q2 (high variability) ----
#_______________________________________________________________________

# fold 1
q2.agg1 <- detections.hi1 %>% 
  
  group_by(iter,
           cam.id) %>%
  
  summarize(detections = n())

#_______________________________________________________________________
# 5. Add "zero" count rows where necessary ----

# write function
add_zero <- function (
    
  df,
  type
  
  ) {
  
  # define df of "possible" cameras
  if(type == "lo") {
    
    possible.grid <- expand.grid(group = c(2, 5, 10, 25, 50),
                                 iter = 1:max(df$iter),
                                 cam.id = 1:9)
    
    # loop through aggregated df
    zero.rows <- data.frame()
    
    for (i in 1:nrow(possible.grid)) {
      
      focal.row <- possible.grid %>% slice(i)
      
      focal.row.match <- plyr::match_df(df[ , -4], focal.row)
      
      # does this row exist?
      if (nrow(focal.row.match) == 0) {
        
        # then bind in a new row
        zero.rows <- rbind(zero.rows,
                           focal.row %>% mutate(detections = 0))
        
      }
      
    }
    
    # bind in and sort
    df.1 <- df %>%
      
      bind_rows(zero.rows) %>%
      
      arrange(group, 
              iter,
              cam.id)
    
  } else {
    
    possible.grid <- expand.grid(iter = 1:max(df$iter),
                                 cam.id = 1:9)
    
    # loop through aggregated df
    zero.rows <- data.frame()
    
    for (i in 1:nrow(possible.grid)) {
      
      focal.row <- possible.grid %>% slice(i)
      
      focal.row.match <- plyr::match_df(df[ , -3], focal.row)
      
      # does this row exist?
      if (nrow(focal.row.match) == 0) {
        
        # then bind in a new row
        zero.rows <- rbind(zero.rows,
                           focal.row %>% mutate(detections = 0))
        
      }
      
    }
    
    # bind in and sort
    df.1 <- df %>%
      
      bind_rows(zero.rows) %>%
      
      arrange(iter,
              cam.id)
    
  }
  
 # return
 return(df.1)
  
}

#_______________________________________________________________________
# 5a. Q1 (low variability) ----
#_______________________________________________________________________

q1.agg1.1 <- add_zero(q1.agg1, "lo")

#_______________________________________________________________________
# 5b. Q2 (high variability) ----
#_______________________________________________________________________

q2.agg1.1 <- add_zero(q2.agg1, "hi")

#_______________________________________________________________________
# 6. Examine detections to see if they seem reasonable ----

# here we'll see what the expected counts would be, given a mean speed

# equation 3 from Rowcliffe et al. (2008):

# y = ((2 + theta) / pi) * r * t * v * D 

# all in km and days

# lens in radians
lens.rad <- (57.3 * (pi / 180))

#_______________________________________________________________________

mean(Q1.cam.all.speeds$speed.4wk)
mean(Q2.cam.all.speeds$speed.4wk)

(e.count.02 <- (((2 + lens.rad) / pi) * (3.5 / 1000) * 28 * (0.00879 / (1000 / 86400)) * (0.2 * 100)))
(e.count.05 <- (((2 + lens.rad) / pi) * (3.5 / 1000) * 28 * (0.00879 / (1000 / 86400)) * (0.5 * 100)))
(e.count.10 <- (((2 + lens.rad) / pi) * (3.5 / 1000) * 28 * (0.00879 / (1000 / 86400)) * (1.0 * 100)))
(e.count.25 <- (((2 + lens.rad) / pi) * (3.5 / 1000) * 28 * (0.00879 / (1000 / 86400)) * (2.5 * 100)))
(e.count.50 <- (((2 + lens.rad) / pi) * (3.5 / 1000) * 28 * (0.00879 / (1000 / 86400)) * (5.0 * 100)))

q1.agg1.1 %>% group_by(group) %>% summarize(mean.detec = mean(detections)) %>% 
  
  mutate(true.count = c(e.count.02,
                        e.count.05,
                        e.count.10,
                        e.count.25,
                        e.count.50)) %>%
  
  mutate(perc.diff = (mean.detec - true.count) / true.count)

q2.agg1.1 %>% summarize(mean.detec = mean(detections))

# seems like I'm underestimating everything by 3-4 fold
# is this because I'm removing relocations in the same pass? 
# I think movement here is too linear for too long

# what if there's just a mismatch in total area?

#_______________________________________________________________________
# 7. Write to .csvs ----
#_______________________________________________________________________

write.csv(q1.agg1.1, paste0(getwd(), "/Derived data/Aggregated detections/detections_lo1.csv"))


write.csv(q2.agg1.1, paste0(getwd(), "/Derived data/Aggregated detections/detections_hi1.csv"))
