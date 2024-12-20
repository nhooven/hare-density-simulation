# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 06 - Sample individuals for detection and movement modeling
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 26 Nov 2024
# Date completed: 09 Dec 2024
# Date last modified: 20 Dec 2024
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
abundances <- c(2, 5, 10, 15, 25, 40)

# how will individuals be allocated to each abundance?
abundances.allocation <- c(rep(abundances[1], abundances[1]),
                           rep(abundances[2], abundances[2]),
                           rep(abundances[3], abundances[3]),
                           rep(abundances[4], abundances[4]),
                           rep(abundances[5], abundances[5]),
                           rep(abundances[6], abundances[6]))

# how many individuals are available?
possible.indivs <- 100

# how many individuals do we need to sample?
total.indivs <- sum(abundances)

#_______________________________________________________________________
# 3. Sample individuals for detections and bind into lookup table  ----
#_______________________________________________________________________

all.combos <- expand.grid(landscape = c("simple", "complex"),
                          variability = c("low", "high"),
                          rep = 1:3)

# loop through all combos
detection.lookup <- data.frame()

for (i in 1:nrow(all.combos)) {
  
  focal.combo <- all.combos[i, ]
  
  # create df
  focal.df <- data.frame("landscape" = focal.combo$landscape,
                         "variability" = focal.combo$variability,
                         "rep" = focal.combo$rep,
                         "indiv" = sample(1:possible.indivs, size = total.indivs),
                         "n" = abundances.allocation)
  
  # bind into full df
  detection.lookup <- rbind(detection.lookup, focal.df)
  
}

#_______________________________________________________________________
# 4. Sample "collared" individuals for SSF/iSSF modeling  ----
#_______________________________________________________________________

# let's take a random sample of n individuals
collared.n <- 10

# loop through all combos
collared.lookup <- data.frame()

for (i in 1:nrow(all.combos)) {
  
  focal.combo <- all.combos[i, ]
  
  # create df
  focal.df <- data.frame("landscape" = focal.combo$landscape,
                         "variability" = focal.combo$variability,
                         "rep" = focal.combo$rep,
                         "indiv" = sample(1:possible.indivs, size = collared.n))
  
  # bind into full df
  collared.lookup <- rbind(collared.lookup, focal.df)
  
}

#_______________________________________________________________________
# 5. Write to .csvs ----
#_______________________________________________________________________

write.csv(detection.lookup, paste0(getwd(), "/Derived_data/Lookup/detection_lookup.csv"))
write.csv(collared.lookup, paste0(getwd(), "/Derived_data/Lookup/collared_lookup.csv"))
