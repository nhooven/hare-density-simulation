# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 07 - Prepare detection data for REMs
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
# 2. Read in data and subset lookup table ----
#_______________________________________________________________________

# passes
passes.simple.low <- read.csv(paste0(getwd(), "/Derived_data/Passes/passes_simple_low.csv"))
passes.complex.low <- read.csv(paste0(getwd(), "/Derived_data/Passes/passes_complex_low.csv"))
passes.simple.high <- read.csv(paste0(getwd(), "/Derived_data/Passes/passes_simple_high.csv"))
passes.complex.high <- read.csv(paste0(getwd(), "/Derived_data/Passes/passes_complex_high.csv"))

# lookup table
lookup.passes <- read.csv(paste0(getwd(), "/Derived_data/Lookup/lookup_passes.csv"))

# subset
lookup.simple.low <- lookup.passes %>% filter(landscape == "simple" & variability == "low")
lookup.complex.low <- lookup.passes %>% filter(landscape == "complex" & variability == "low")
lookup.simple.high <- lookup.passes %>% filter(landscape == "simple" & variability == "high")
lookup.complex.high <- lookup.passes %>% filter(landscape == "complex" & variability == "high")

#_______________________________________________________________________
# 3. Extract individual data at each density ----
#_______________________________________________________________________
# 3a. Define function ----
#_______________________________________________________________________

extract_passes <- function(landscape,     # "simple" or "complex"
                           variability,   # "low" or "high"
                           densities) {   # vector of densities, assume indivs/ha
  
  # assign correct passes and lookup tables
  if (landscape == "simple" & variability == "low") {
    
    focal.passes <- passes.simple.low
    focal.lookup <- lookup.simple.low
    
  } else {
    
    if (landscape == "complex" & variability == "low") {
      
      focal.passes <- passes.complex.low
      focal.lookup <- lookup.complex.low
      
    } else {
      
      if (landscape == "simple" & variability == "high") {
        
        focal.passes <- passes.simple.high
        focal.lookup <- lookup.simple.high
        
      } else {
        
        focal.passes <- passes.complex.high
        focal.lookup <- lookup.complex.high
        
      }
      
    }
    
  }
  
  # loop through densities
  # convert densities to n.indiv 
  indivs <- densities * 10
  
  # df to hold all
  all.passes <- data.frame()
  
  for (i in indivs) {
    
    # extract as a vector
    focal.indivs <- focal.lookup %>% filter(n == i) %>% dplyr::select(indiv) %>% pull()
    
    # extract passes and add identifiers
    focal.passes.1 <- focal.passes %>% filter(indiv %in% focal.indivs) %>%
      
                                       mutate(landscape = landscape,
                                              variability = variability,
                                              n.indiv = i) %>%
      
                                       # drop X
                                       dplyr::select(-X)
    
    # bind in
    all.passes <- rbind(all.passes, focal.passes.1)
    
  }
  
  # return
  return(all.passes)
  
}     

#_______________________________________________________________________
# 3b. Run through ----
#_______________________________________________________________________

passes.extracted.simple.low <- extract_passes("simple", "low", c(0.2, 0.5, 1.0))
passes.extracted.complex.low <- extract_passes("complex", "low", c(0.2, 0.5, 1.0))
passes.extracted.simple.high <- extract_passes("simple", "high", c(0.2, 0.5, 1.0))
passes.extracted.complex.high <- extract_passes("complex", "high", c(0.2, 0.5, 1.0))

#_______________________________________________________________________
# 4. Aggregate together and add in cameras without passes ----
#______________________________________________________________________
# 4a. Define function ----
#_______________________________________________________________________

group_passes <- function (x) {
  
  # loop through all n.indiv
  # blank df
  all.passes.grouped <- data.frame()
  
  for (i in unique(x$n.indiv)) {
    
    # subset
    x.1 <- x %>% filter(n.indiv == i)
    
    # group by cam and sum all passes
    x.2 <- x.1 %>% group_by(cam.id) %>%
      
      summarize(total.passes = sum(n))
    
    # determine if any cams have no passes
    no.passes <- which(!c(1:9) %in% x.2$cam.id)
    
    # add a zero entry for each one
    for (j in no.passes) {
      
      new.row <- data.frame(cam.id = j,
                            total.passes = 0)
      
      x.2 <- rbind(x.2, new.row)
      
    }
    
    # sort by cam.id and add in identifiers
    x.3 <- x.2 %>% 
      
      arrange(cam.id) %>%
      
      mutate(n.indiv = x.1$n.indiv[1],
             landscape = x.1$landscape[1],
             variability = x.1$variability[1])
    
    # bind
    all.passes.grouped <- rbind(all.passes.grouped, x.3)
    
  }
  
  # return
  return(all.passes.grouped)
  
}

#_______________________________________________________________________
# 4b. Use function ----
#_______________________________________________________________________

passes.simple.low.gp <- group_passes(passes.extracted.simple.low)
passes.complex.low.gp <- group_passes(passes.extracted.complex.low)
passes.simple.high.gp <- group_passes(passes.extracted.simple.high)
passes.complex.high.gp <- group_passes(passes.extracted.complex.high)

#_______________________________________________________________________
# 5. Draw passes based on detection probability ----

# define possible detection probabilities
# to keep this simple we'll assume a constant p for each camera

det.p <- 0.30

#_______________________________________________________________________
# 5a. Define function ----
#_______________________________________________________________________

sample_detections <- function (x) {
  
  # loop through n.indiv
  all.detec <- data.frame()
  
  for (i in unique(x$n.indiv)) {
    
    x.1 <- x %>% filter(n.indiv == i)
    
    # work by row and bind in
    total.detections <- as.vector(by(x.1, 
                                     x.1$cam.id, 
                                     function (y) rbinom(1, 
                                                         size = y$total.passes, 
                                                         prob = det.p)))
    
    x.1$total.detections <- total.detections
    
    # bind
    all.detec <- rbind(all.detec, x.1)
    
  }
  
  # return
  return(all.detec)
  
}

#_______________________________________________________________________
# 5b. Use function ----
#_______________________________________________________________________

detec.simple.low <- sample_detections(passes.simple.low.gp)
detec.complex.low <- sample_detections(passes.complex.low.gp)
detec.simple.high <- sample_detections(passes.simple.high.gp)
detec.complex.high <- sample_detections(passes.complex.high.gp)

#_______________________________________________________________________
# 5c. Bind together ----
#_______________________________________________________________________

detec.all <- rbind(detec.simple.low,
                   detec.complex.low,
                   detec.simple.high,
                   detec.complex.high)

#_______________________________________________________________________
# 6. Plot (sanity check) ----
#_______________________________________________________________________

ggplot(detec.all) +
  
  theme_bw() +
  
  facet_grid(landscape ~ variability) +
  
  geom_point(aes(x = as.factor(n.indiv),
                 y = total.detections)) +
  
  theme(panel.grid = element_blank())

#_______________________________________________________________________
# 7. Write to .csv ----
#_______________________________________________________________________
  
write.csv(detec.all, paste0(getwd(), "/Derived_data/Detections/detec_all.csv"))
