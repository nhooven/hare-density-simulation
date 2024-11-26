# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 06 - Sample individuals for detection and movement modeling
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 26 Nov 2024
# Date completed: 
# Date last modified: 26 Nov 2024
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation

#_______________________________________________________________________
# 2. Read in passes data ----
#_______________________________________________________________________

passes.simple.weak <- read.csv(paste0(getwd(), "/Derived_data/Passes/passes_simple_weak.csv"))
#passes.simple.strong <- read.csv(paste0(getwd(), "/Derived_data/Passes/passes_simple_strong.csv"))
#passes.complex.weak <- read.csv(paste0(getwd(), "/Derived_data/Passes/passes_complex_weak.csv"))
#passes.complex.strong <- read.csv(paste0(getwd(), "/Derived_data/Passes/passes_complex_strong.csv"))

#_______________________________________________________________________
# 3. Sample individuals from each scenario without replacement ----
#_______________________________________________________________________

# how are abundances set up?
abundances <- c(2, 5, 10, 20, 50)

# how will individuals be allocated to each abundance?
abundances.allocation <- c(rep(abundances[1], abundances[1]),
                           rep(abundances[2], abundances[2]),
                           rep(abundances[3], abundances[3]),
                           rep(abundances[4], abundances[4]),
                           rep(abundances[5], abundances[5]))

# how many individuals are available? (should be something like 100, whatever we simulated)
possible.indivs <- max(passes.simple.weak$indiv)

# how many individuals do we need to sample?
total.indivs <- max(passes.simple.weak$indiv)
#total.indivs <- 2 + 5 + 10 + 20 + 50     # should be however many we need

#_______________________________________________________________________
# 3a. Sample total.indivs from all simulated tracks  ----
#_______________________________________________________________________

set.seed(867)

sampled.simple.weak <- sample(1:possible.indivs, size = total.indivs)
#sampled.simple.strong <- sample(1:possible.indivs, size = total.indivs)
#sampled.complex.weak <- sample(1:possible.indivs, size = total.indivs)
#sampled.complex.strong <- sample(1:possible.indivs, size = total.indivs)

#_______________________________________________________________________
# 3b. Bind into a lookup table  ----
#_______________________________________________________________________

lookup.simple.weak <- data.frame("landscape" = "simple",
                                 "hs" = "weak",
                                 "indiv" = 1:87,  # change to sampled.simple.weak
                                 "n" = abundances.allocation)

# and so forth

# and bind together
# rbind

#_______________________________________________________________________
# 4. Sample "collared" individuals  ----

# let's take a random sample of n individuals
collared.n <- 20

#_______________________________________________________________________
# 4a. Sample ----
#_______________________________________________________________________

set.seed(223)

# draw from each sample
collared.simple.weak <- sample(sampled.simple.weak, size = collared.n)
#collared.simple.strong <- sample(sampled.simple.strong, size = collared.n)
#collared.complex.weak <- sample(sampled.complex.weak, size = collared.n)
#collared.complex.strong <- sample(sampled.complex.strong, size = collared.n)

#_______________________________________________________________________
# 4b. Bind into a lookup table ----
#_______________________________________________________________________

lookup.collared.simple.weak <- data.frame("landscape" = "simple",
                                          "hs" = "weak",
                                          "indiv" = collared.simple.weak)

# and so forth

# and bind together
# rbind

#_______________________________________________________________________
# 5. Write to .csvs ----
#_______________________________________________________________________

# lookup folder