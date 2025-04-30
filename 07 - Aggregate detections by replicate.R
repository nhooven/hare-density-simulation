# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 07 - Aggregate detections by replicate
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 28 Apr 2025
# Date completed: 28 Apr 2025
# Date last modified: 29 Apr 2025
# R version: 4.4.3

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation

#_______________________________________________________________________
# 2. Purpose ----
#_______________________________________________________________________

# Now that we have our sampled individuals, we'll need to aggregate camera
# contacts for each replicate (1000) by abundance (4) combination to allow for
# easy sampling for the REM Monte Carlo/bootstrapping procedure in the next script

# since this is camera-specific, each question will get 1000 x 4 x 9 (36000) observations

#_______________________________________________________________________
# 3. Read in data ----
#_______________________________________________________________________

# individuals for subsetting
indivs.1 <- read.csv(paste0(getwd(), "/Derived data/Sampled reps/contacts_1.csv"))
indivs.2 <- read.csv(paste0(getwd(), "/Derived data/Sampled reps/contacts_2.csv"))
indivs.3 <- read.csv(paste0(getwd(), "/Derived data/Sampled reps/contacts_3.csv"))

# camera contacts (tallied points + passes for all 500 indivs)
contacts.1T <- read.csv(paste0(getwd(), "/Derived data/Sampled - Camera contacts/contacts_1T.csv"))
contacts.1NT <- read.csv(paste0(getwd(), "/Derived data/Sampled - Camera contacts/contacts_1NT.csv"))
contacts.2T <- read.csv(paste0(getwd(), "/Derived data/Sampled - Camera contacts/contacts_2T.csv"))
contacts.2NT <- read.csv(paste0(getwd(), "/Derived data/Sampled - Camera contacts/contacts_2NT.csv"))

#_______________________________________________________________________
# 4. Define function ----
#_______________________________________________________________________

aggregate_contacts <- function (q = 1) {
  
  # use correct dfs
  if (q == 1) {
    
    indivs <- indivs.1
    contacts.T <- contacts.1T
    contacts.NT <- contacts.1NT
    
  } 
  
  if (q == 2) {
    
    indivs <- indivs.2
    contacts.T <- contacts.2T
    contacts.NT <- contacts.2NT
    
  }
  
  if (q == 3) {
    
    indivs <- indivs.3
    contacts.T <- contacts.2T
    contacts.NT <- contacts.2NT
    
  }
  
  # loop through replicates
  all.contacts.agg.df <- data.frame()
  
  for (i in 1:1000) {
    
    # loop through abundances
    contacts.agg.df.i <- data.frame()
    
    for (j in c(32, 16, 8, 4)) {
      
     # subset indivs
     focal.indivs.T <- indivs %>% filter(rep == i, purpose == "T", abund == j)
     focal.indivs.NT <- indivs %>% filter(rep == i, purpose == "NT", abund == j)
     
     # subset contacts.T
     focal.contacts.T <- contacts.T %>% filter(indiv %in% focal.indivs.T$indiv)
     
     # subset contacts.NT
     focal.contacts.NT <- contacts.NT %>% filter(indiv %in% focal.indivs.NT$indiv) 
     
     # bind together
     focal.contacts <- rbind(focal.contacts.T, focal.contacts.NT)
     
     # aggregate by camera
     focal.contacts.agg <- focal.contacts %>%
       
       group_by(cam.id) %>%
       
       summarize(points = sum(points),
                 passes = sum(passes)) %>%
       
       # add in id columns
       mutate(rep = i,
              abund = j)
     
     # bind in
     contacts.agg.df.i <- rbind(contacts.agg.df.i, focal.contacts.agg) 
     
    }
    
    # bind in 
    all.contacts.agg.df <- rbind(all.contacts.agg.df, contacts.agg.df.i)
    
  }
 
  # return
  return(all.contacts.agg.df)
  
}
  
#_______________________________________________________________________
# 5. Use function ----
#_______________________________________________________________________

agg.contacts.1 <- aggregate_contacts(q = 1)
agg.contacts.2 <- aggregate_contacts(q = 2)
agg.contacts.3 <- aggregate_contacts(q = 3)

#_______________________________________________________________________
# 6. Plots for sanity check ----
#_______________________________________________________________________

# question 1 (points)
ggplot(data = agg.contacts.1,
       aes(x = as.factor(abund),
           y = points)) +
  
  theme_bw() +
  
  geom_jitter(alpha = 0.05)

# question 1 (passes)
ggplot(data = agg.contacts.1,
       aes(x = as.factor(abund),
           y = passes)) +
  
  theme_bw() +
  
  geom_jitter(alpha = 0.05)

# question 2 (points)
ggplot(data = agg.contacts.2,
       aes(x = as.factor(abund),
           y = points)) +
  
  theme_bw() +
  
  geom_jitter(alpha = 0.05)

# question 2 (passes)
ggplot(data = agg.contacts.2,
       aes(x = as.factor(abund),
           y = passes)) +
  
  theme_bw() +
  
  geom_jitter(alpha = 0.05)

#_______________________________________________________________________
# 5. Write to files ----
#_______________________________________________________________________

write.csv(agg.contacts.1, file = paste0(getwd(), "/Derived data/Aggregated contacts/agg_contacts_1.csv"))
write.csv(agg.contacts.2, file = paste0(getwd(), "/Derived data/Aggregated contacts/agg_contacts_2.csv"))
write.csv(agg.contacts.3, file = paste0(getwd(), "/Derived data/Aggregated contacts/agg_contacts_3.csv"))
