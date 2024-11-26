# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 07 - Prepare detection data for REMs
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 26 Nov 2024
# Date completed: 
# Date last modified: 
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation

#_______________________________________________________________________
# 2. Read in passes data ----
#_______________________________________________________________________

passes.simple.weak <- read.csv(paste0(getwd(), "/Derived_data/Passes/passes_simple_weak.csv"))

#_______________________________________________________________________
# 3. Sample individuals for each density WITHOUT replacement ----

possible.indivs <- max(passes.simple.weak$indiv)

total.indivs <- max(passes.simple.weak$indiv)
#total.indivs <- 2 + 5 + 10 + 20 + 50

#_______________________________________________________________________
# 3a. Simple / weak ----
#_______________________________________________________________________

# sample total.indivs from all simulated tracks 
sampled.indivs <- sample(1:possible.indivs,
                         size = total.indivs)

passes.simple.weak.02 <- passes.simple.weak %>% filter(indiv %in% sampled.indivs[1:2])
#passes.simple.weak.05 <- passes.simple.weak %>% filter(indiv %in% sampled.indivs[3:8])

# and so forth

#_______________________________________________________________________
# 4. Aggregate together and add in cameras without passes ----
#______________________________________________________________________
# 4a. Define function ----
#_______________________________________________________________________

group_passes <- function (x) {
  
  # group by cam and sum all passes
  x.1 <- x %>% group_by(cam.id) %>%
    
    summarize(total.passes = sum(n))
  
  # determine if any cams have no passes
  no.passes <- which(!c(1:9) %in% x.1$cam.id)
  
  # add a zero entry for each one
  for (i in no.passes) {
    
    new.row <- data.frame(cam.id = i,
                          total.passes = 0)
    
    x.1 <- rbind(x.1, new.row)
    
  }
  
  # sort by cam.id
  x.2 <- x.1 %>% arrange(cam.id)
  
  # return
  return(x.2)
  
}

#_______________________________________________________________________
# 4b. Use function ----
#_______________________________________________________________________

passes.simple.weak.gp <- group_passes(passes.simple.weak)

#_______________________________________________________________________
# 5. Draw passes based on detection probability ----

# define possible detection probabilities
# to keep this simple we'll assume a constant p for each camera

det.p <- 0.30

#_______________________________________________________________________

# work by row and bind in
passes.simple.weak.gp$total.passes.p <- as.vector(by(passes.simple.weak.gp, 
                                                     passes.simple.weak.gp$cam.id, 
                                                     function (x) rbinom(1, 
                                                                         size = x$total.passes, 
                                                                         prob = 0.30)))
