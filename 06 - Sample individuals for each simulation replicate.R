# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 06 - Sample individuals for each simulation replicate
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 28 Apr 2025
# Date completed: 28 Apr 2025
# Date last modified: 30 Apr 2025
# R version: 4.4.3

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation
library(mefa4)                # %notin%

#_______________________________________________________________________
# 2. Purpose ----
#_______________________________________________________________________

# This script's main purpose is to define which individuals will be included for each
# REM replicate. We have already simulated tracks and camera contacts as well as 
# fit CTMMs to estimate speed for all 2,000 individuals (1T, 1NT, 2T, 2NT), so this step
# will create lookup tables for each question.

# Each replicate will require:
  # Target individuals (i.e., those will activity centers within the 10-ha "unit")
    # Either 4, 8, 16, or 32
  # Non-target individuals (i.e., those with activity centers outside the unit in the 20-ha "matrix")
    # Double the number of targets since this area is double the unit
    # Reminder that this is to reduce the probability of edge effects
  # Extra collared individuals (i.e., those not included as target or non-target individuals)
    # These are drawn from the remaining pool, and we'll draw 4 for each replicate

  # Collared individuals
  # From the three pools above, 4 individuals will be drawn each, totaling 12 "collared" animals
  # We will then use the same individuals for each of the 12 collar scenarios
  # Note that each abundance level will be a subsample of the 32, where the 4 individuals will be collared

#_______________________________________________________________________
# 3. Write function ----

# QUESTIONS 1 AND 2
# since the structure is the same here, we can just write a generic function and 
# apply it twice

#_______________________________________________________________________

sample_indiv_reps <- function (n.rep = 1000,
                               q = 1) {
  
  # loop
  all.contact.df <- data.frame()
  all.col.df <- data.frame()
  
  for (i in 1:n.rep) {
    
    # available individuals
    indiv.ids <- 1:500
    
    # sample targets
    T.32 <- sample(indiv.ids, size = 32, replace = FALSE)
    T.16 <- sample(T.32, size = 16, replace = FALSE)
    T.08 <- sample(T.16, size = 8, replace = FALSE)
    T.04 <- sample(T.08, size = 4, replace = FALSE)
    T.col <- T.04  # final four are "collared"
    
    # sample non-targets (must be twice the targets)
    NT.32 <- sample(indiv.ids, size = 32 * 2, replace = FALSE)
    NT.16 <- sample(NT.32, size = 16 * 2, replace = FALSE)
    NT.08 <- sample(NT.16, size = 8 * 2, replace = FALSE)
    NT.04 <- sample(NT.08, size = 4 * 2, replace = FALSE)
    NT.col <- sample(NT.04, size = 4, replace = FALSE)    # sample a final four for collaring
    
    # sample extra collared individuals
    # we'll just draw these 4 from the targets 
    extra.col <- sample(indiv.ids[indiv.ids %notin% T.32], size = 4, replace = FALSE)
    
    # create dfs
    # targets
    T.df <- data.frame(q = q,
                       rep = i,
                       purpose = "T",
                       indiv = c(T.32,
                                 T.16,
                                 T.08,
                                 T.04),
                       abund = c(rep(32, times = 32),
                                 rep(16, times = 16),
                                 rep(8, times = 8),
                                 rep(4, times = 4)))
    
    # non-targets
    NT.df <- data.frame(q = q,
                        rep = i,
                        purpose = "NT",
                        indiv = c(NT.32,
                                  NT.16,
                                  NT.08,
                                  NT.04),
                        abund = c(rep(32, times = 64),
                                  rep(16, times = 32),
                                  rep(8, times = 16),
                                  rep(4, times = 8)))
    
    # collared individuals
    col.df <- data.frame(q = q,
                         rep = i,
                         purpose = "col",
                         indiv = c(T.col,
                                   NT.col,
                                   extra.col),
                         source = c(rep("T", times = 4),
                                    rep("NT", times = 4),
                                    rep("T", times = 4)))
    
    # bind contact individuals together
    contact.df <- rbind(T.df, NT.df)
    
    # bind both into master df
    all.contact.df <- rbind(all.contact.df, contact.df)
    all.col.df <- rbind(all.col.df, col.df)
    
  }
  
  # bind into list and return
  both.df <- list(all.contact.df, all.col.df)
  
  return(both.df)

}

#_______________________________________________________________________
# 4. Use function ----
#_______________________________________________________________________

sampled.indiv.reps.1 <- sample_indiv_reps(1000, 1)
sampled.indiv.reps.2 <- sample_indiv_reps(1000, 2)

#_______________________________________________________________________
# 5. Question 3 ----

# Here we'll sample camera individuals as usual, but
# sample collared individuals weighted by their speed
# this emulates collaring individuals with a more restless movement syndrome
# than the population on average

# read in "true speeds" - must be the highly variable subset
speeds.2T <- read.csv(paste0(getwd(), "/Derived data/Sampled - Speeds/speeds_2T.csv"))
speeds.2NT <- read.csv(paste0(getwd(), "/Derived data/Sampled - Speeds/speeds_2NT.csv"))

hist(speeds.2T$true.speed)
hist(speeds.2NT$true.speed)

#_______________________________________________________________________
# 5a. Write function ----
#_______________________________________________________________________

sample_indiv_reps_q3 <- function (n.rep = 1000) {
  
  # loop
  all.contact.df <- data.frame()
  all.col.df <- data.frame()
  
  for (i in 1:n.rep) {
    
    # available individuals
    indiv.ids.T <- speeds.2T$indiv
    
    # sample targets
    T.32 <- sample(indiv.ids.T, size = 32, replace = FALSE)
    T.16 <- sample(T.32, size = 16, replace = FALSE)
    T.08 <- sample(T.16, size = 8, replace = FALSE)
    T.04 <- sample(T.08, size = 4, replace = FALSE)
    
    # sample non-targets (must be twice the targets)
    NT.32 <- sample(indiv.ids.T, size = 32 * 2, replace = FALSE)
    NT.16 <- sample(NT.32, size = 16 * 2, replace = FALSE)
    NT.08 <- sample(NT.16, size = 8 * 2, replace = FALSE)
    NT.04 <- sample(NT.08, size = 4 * 2, replace = FALSE)
    
    # sample  collared individuals
    # we'll take six from the targets and 6 from the non-targets, weighted by speed
    extra.col.T <- sample(indiv.ids.T, size = 6, replace = FALSE, prob = speeds.2T$true.speed)
    extra.col.NT <- sample(indiv.ids.T, size = 6, replace = FALSE, prob = speeds.2NT$true.speed)
    
    # create dfs
    # targets
    T.df <- data.frame(q = 3,
                       rep = i,
                       purpose = "T",
                       indiv = c(T.32,
                                 T.16,
                                 T.08,
                                 T.04),
                       abund = c(rep(32, times = 32),
                                 rep(16, times = 16),
                                 rep(8, times = 8),
                                 rep(4, times = 4)))
    
    # non-targets
    NT.df <- data.frame(q = 3,
                        rep = i,
                        purpose = "NT",
                        indiv = c(NT.32,
                                  NT.16,
                                  NT.08,
                                  NT.04),
                        abund = c(rep(32, times = 64),
                                  rep(16, times = 32),
                                  rep(8, times = 16),
                                  rep(4, times = 8)))
    
    # collared individuals
    col.df <- data.frame(q = 3,
                         rep = i,
                         purpose = "col",
                         indiv = c(extra.col.T,
                                   extra.col.NT),
                         source = c(rep("T", times = 6),
                                    rep("NT", times = 6)))
    
    # bind contact individuals together
    contact.df <- rbind(T.df, NT.df)
    
    # bind both into master df
    all.contact.df <- rbind(all.contact.df, contact.df)
    all.col.df <- rbind(all.col.df, col.df)
    
  }
  
  # bind into list and return
  both.df <- list(all.contact.df, all.col.df)
  
  return(both.df)
  
}

#_______________________________________________________________________
# 5b. Use function ----
#_______________________________________________________________________

sampled.indiv.reps.3 <- sample_indiv_reps_q3(1000)

#_______________________________________________________________________
# 6. Write to files ----
#_______________________________________________________________________

write.csv(sampled.indiv.reps.1[[1]], file = paste0(getwd(), "/Derived data/Sampled reps/contacts_1.csv"))
write.csv(sampled.indiv.reps.1[[2]], file = paste0(getwd(), "/Derived data/Sampled reps/collared_1.csv"))

write.csv(sampled.indiv.reps.2[[1]], file = paste0(getwd(), "/Derived data/Sampled reps/contacts_2.csv"))
write.csv(sampled.indiv.reps.2[[2]], file = paste0(getwd(), "/Derived data/Sampled reps/collared_2.csv"))

write.csv(sampled.indiv.reps.3[[1]], file = paste0(getwd(), "/Derived data/Sampled reps/contacts_3.csv"))
write.csv(sampled.indiv.reps.3[[2]], file = paste0(getwd(), "/Derived data/Sampled reps/collared_3.csv"))
