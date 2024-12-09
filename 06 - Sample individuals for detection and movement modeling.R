# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 06 - Sample individuals for detection and movement modeling
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 26 Nov 2024
# Date completed: 09 Dec 2024
# Date last modified: 09 Dec 2024
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation

#_______________________________________________________________________
# 2. Read in passes data ----
#_______________________________________________________________________

passes.simple.low <- read.csv(paste0(getwd(), "/Derived_data/Passes/passes_simple_low.csv"))
passes.complex.low <- read.csv(paste0(getwd(), "/Derived_data/Passes/passes_complex_low.csv"))
passes.simple.high <- read.csv(paste0(getwd(), "/Derived_data/Passes/passes_simple_high.csv"))
passes.complex.high <- read.csv(paste0(getwd(), "/Derived_data/Passes/passes_complex_high.csv"))

#_______________________________________________________________________
# 3. Sample individuals from each scenario without replacement ----
#_______________________________________________________________________

# how are abundances set up?
abundances <- c(2, 5, 10)
                
                #, 20, 50)

# how will individuals be allocated to each abundance?
abundances.allocation <- c(rep(abundances[1], abundances[1]),
                           rep(abundances[2], abundances[2]),
                           rep(abundances[3], abundances[3]))
                           
                           #,
                           #rep(abundances[4], abundances[4]),
                           #rep(abundances[5], abundances[5]))

# how many individuals are available? (should be something like 100, whatever we simulated)
possible.indivs <- max(passes.simple.low$indiv)

# how many individuals do we need to sample?
total.indivs <- 2 + 5 + 10
#total.indivs <- 2 + 5 + 10 + 20 + 50     # should be however many we need

#_______________________________________________________________________
# 3a. Sample total.indivs from all simulated tracks  ----
#_______________________________________________________________________

set.seed(867)

sampled.simple.low <- sample(1:possible.indivs, size = total.indivs)
sampled.complex.low <- sample(1:possible.indivs, size = total.indivs)
sampled.simple.high <- sample(1:possible.indivs, size = total.indivs)
sampled.complex.high <- sample(1:possible.indivs, size = total.indivs)

#_______________________________________________________________________
# 3b. Bind into a lookup table  ----
#_______________________________________________________________________

lookup.simple.low <- data.frame("landscape" = "simple",
                                "variability" = "low",
                                "indiv" = sampled.simple.low,
                                "n" = abundances.allocation)

lookup.complex.low <- data.frame("landscape" = "complex",
                                 "variability" = "low",
                                 "indiv" = sampled.complex.low,
                                 "n" = abundances.allocation)

lookup.simple.high <- data.frame("landscape" = "simple",
                                 "variability" = "high",
                                 "indiv" = sampled.simple.high,
                                 "n" = abundances.allocation)

lookup.complex.high <- data.frame("landscape" = "complex",
                                  "variability" = "high",
                                  "indiv" = sampled.complex.high,
                                  "n" = abundances.allocation)

# and bind together
lookup.all <- rbind(lookup.simple.low,
                    lookup.complex.low,
                    lookup.simple.high,
                    lookup.complex.high)

#_______________________________________________________________________
# 4. Sample "collared" individuals  ----

# let's take a random sample of n individuals
collared.n <- 10

#_______________________________________________________________________
# 4a. Sample ----
#_______________________________________________________________________

set.seed(223)

# draw from each sample
collared.simple.low <- sample(1:possible.indivs, size = collared.n)
collared.complex.low <- sample(1:possible.indivs, size = collared.n)
collared.simple.high <- sample(1:possible.indivs, size = collared.n)
collared.complex.high <- sample(1:possible.indivs, size = collared.n)

#_______________________________________________________________________
# 4b. Bind into a lookup table ----
#_______________________________________________________________________

lookup.collared.simple.low <- data.frame("landscape" = "simple",
                                         "variability" = "low",
                                         "indiv" = collared.simple.low)

lookup.collared.complex.low <- data.frame("landscape" = "complex",
                                         "variability" = "low",
                                         "indiv" = collared.complex.low)

lookup.collared.simple.high <- data.frame("landscape" = "simple",
                                          "variability" = "high",
                                          "indiv" = collared.simple.high)

lookup.collared.complex.high <- data.frame("landscape" = "complex",
                                           "variability" = "high",
                                           "indiv" = collared.complex.high)

# and bind together
lookup.collared.all <- rbind(lookup.collared.simple.low,
                             lookup.collared.complex.low,
                             lookup.collared.simple.high,
                             lookup.collared.complex.high)

#_______________________________________________________________________
# 5. Write to .csvs ----
#_______________________________________________________________________

write.csv(lookup.all, paste0(getwd(), "/Derived_data/Lookup/lookup_passes.csv"))
write.csv(lookup.collared.all, paste0(getwd(), "/Derived_data/Lookup/lookup_collared.csv"))
