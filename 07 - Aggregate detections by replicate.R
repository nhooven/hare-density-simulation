# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 07 - Aggregate detections by replicate
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 28 Apr 2025
# Date completed: 28 Apr 2025
# Date last modified: 06 Jul 2026
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
indivs.1 <- readRDS(paste0(getwd(), "/data_derived/sampled_reps/contacts_1.rds"))
indivs.2 <- readRDS(paste0(getwd(), "/data_derived/sampled_reps/contacts_2.rds"))
indivs.3 <- readRDS(paste0(getwd(), "/data_derived/sampled_reps/contacts_3.rds"))

# camera contacts (tallied points + passes for all 500 indivs)
contacts.1T.TV1 <- readRDS(paste0(getwd(), "/data_derived/sampled_contacts/contacts_1T_TV1.rds"))
contacts.1T.TV2 <- readRDS(paste0(getwd(), "/data_derived/sampled_contacts/contacts_1T_TV2.rds"))
contacts.1T.TV3 <- readRDS(paste0(getwd(), "/data_derived/sampled_contacts/contacts_1T_TV3.rds"))

contacts.1NT.TV1 <- readRDS(paste0(getwd(), "/data_derived/sampled_contacts/contacts_1NT_TV1.rds"))
contacts.1NT.TV2 <- readRDS(paste0(getwd(), "/data_derived/sampled_contacts/contacts_1NT_TV2.rds"))
contacts.1NT.TV3 <- readRDS(paste0(getwd(), "/data_derived/sampled_contacts/contacts_1NT_TV3.rds"))

contacts.2T <- readRDS(paste0(getwd(), "/data_derived/sampled_contacts/contacts_2T.rds"))
contacts.2NT <- readRDS(paste0(getwd(), "/data_derived/sampled_contacts/contacts_2NT.rds"))

#_______________________________________________________________________
# 4. Define function ----
#_______________________________________________________________________

aggregate_contacts <- function (.indivs, .contacts.T, .contacts.NT) {
  
  # loop through replicates
  all.contacts.agg.df <- data.frame()
  
  for (i in 1:1000) {
    
    # loop through abundances
    contacts.agg.df.i <- data.frame()
    
    for (j in c(32, 16, 8, 4)) {
      
     # subset indivs
     focal.indivs.T <- .indivs %>% filter(rep == i, purpose == "T", abund == j)
     focal.indivs.NT <- .indivs %>% filter(rep == i, purpose == "NT", abund == j)
     
     # subset contacts.T
     focal.contacts.T <- .contacts.T %>% filter(indiv %in% focal.indivs.T$indiv)
     
     # subset contacts.NT
     focal.contacts.NT <- .contacts.NT %>% filter(indiv %in% focal.indivs.NT$indiv) 
     
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

agg.contacts.1.TV1 <- aggregate_contacts(indivs.1, contacts.1T.TV1, contacts.1NT.TV1)
agg.contacts.1.TV2 <- aggregate_contacts(indivs.1, contacts.1T.TV2, contacts.1NT.TV2)
agg.contacts.1.TV3 <- aggregate_contacts(indivs.1, contacts.1T.TV3, contacts.1NT.TV3)

agg.contacts.2 <- aggregate_contacts(indivs.2, contacts.2T, contacts.2NT)

agg.contacts.3 <- aggregate_contacts(indivs.3, contacts.2T, contacts.2NT)

#_______________________________________________________________________
# 6. Plots for sanity check ----
#_______________________________________________________________________

agg.contacts.3 |>

ggplot(aes(x = as.factor(abund),
           y = points)) +
  
  theme_bw() +
  
  geom_jitter(alpha = 0.05)

#_______________________________________________________________________
# 5. Write to files ----
#_______________________________________________________________________

saveRDS(agg.contacts.1.TV1, paste0(getwd(), "/data_derived/aggregated_contacts/agg_contacts_1_TV1.rds"))
saveRDS(agg.contacts.1.TV2, paste0(getwd(), "/data_derived/aggregated_contacts/agg_contacts_1_TV2.rds"))
saveRDS(agg.contacts.1.TV3, paste0(getwd(), "/data_derived/aggregated_contacts/agg_contacts_1_TV3.rds"))

saveRDS(agg.contacts.2, paste0(getwd(), "/data_derived/aggregated_contacts/agg_contacts_2.rds"))
saveRDS(agg.contacts.3, paste0(getwd(), "/data_derived/aggregated_contacts/agg_contacts_3.rds"))
