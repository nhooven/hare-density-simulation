# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 07 - Prepare detection data for REMs
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
# 2. Read in data and subset lookup table ----
#_______________________________________________________________________

# passes
passes.4 <- read.csv(paste0(getwd(), "/Derived_data/Passes/passes_4.csv"))
passes.9 <- read.csv(paste0(getwd(), "/Derived_data/Passes/passes_9.csv"))
passes.16 <- read.csv(paste0(getwd(), "/Derived_data/Passes/passes_16.csv"))

# lookup table
lookup.passes <- read.csv(paste0(getwd(), "/Derived_data/Lookup/detection_lookup.csv"))

#_______________________________________________________________________
# 3. Extract individual data at each density ----
#_______________________________________________________________________
# 3a. Define function ----
#_______________________________________________________________________

# expand grid
all.combos <- expand.grid(landscape = c("simple", "complex"),
                          variability = c("low", "high"),
                          rep = 1:3,
                          n = c(2, 5, 10, 15, 25, 40))

extract_passes <- function(df) {
  
  # loop through all combos
  all.focal.passes <- data.frame()
  
  for (i in 1:nrow(all.combos)) {
    
    focal.combo <- all.combos[i, ]
    
    # subset lookup table
    focal.lookup <- lookup.passes %>%
      
      filter(landscape == focal.combo$landscape &
             variability == focal.combo$variability &
             rep == focal.combo$rep &
             n == focal.combo$n)
    
    # subset correct individuals from passes df
    focal.passes <- df %>%
      
      filter(landscape == focal.combo$landscape &
             variability == focal.combo$variability &
             rep == focal.combo$rep &
             indiv %in% focal.lookup$indiv) %>%
      
      # ensure that I'm not confused by what "n" means
      mutate(n.indiv = focal.combo$n) %>%
      rename(passes = n)
    
    # bind into df
    all.focal.passes <- rbind(all.focal.passes, focal.passes)
  
  }
  
  # return
  return(all.focal.passes)
  
}

#_______________________________________________________________________
# 3b. Run through ----
#_______________________________________________________________________

passes.extracted.4 <- extract_passes(passes.4)
passes.extracted.9 <- extract_passes(passes.9)
passes.extracted.16 <- extract_passes(passes.16)

#_______________________________________________________________________
# 4. Aggregate together and add in cameras without passes ----
#______________________________________________________________________
# 4a. Define function ----
#_______________________________________________________________________

group_passes <- function (df,
                          n.cams) {
  
  # loop through all combos
  all.passes.grouped <- data.frame()
  
  for (i in 1:nrow(all.combos)) {
    
    focal.combo <- all.combos[i, ]
    
    # subset extracted passes
    focal.passes <- df %>%
      
      filter(landscape == focal.combo$landscape &
             variability == focal.combo$variability &
             rep == focal.combo$rep &
             n.indiv == focal.combo$n)
    
    # group by cam and sum all passes
    focal.passes.1 <- focal.passes %>% 
      
      group_by(cam.id) %>%
      
      summarize(total.passes = sum(passes))
    
    # determine if any cams have no passes
    no.passes <- which(!c(1:n.cams) %in% focal.passes.1$cam.id)
    
    # add a zero entry for each one
    for (j in no.passes) {
      
      new.row <- data.frame(cam.id = j,
                            total.passes = 0)
      
      focal.passes.1 <- rbind(focal.passes.1, new.row)
      
    }
    
    # sort by cam.id and add in identifiers
    focal.passes.2 <- focal.passes.1 %>% 
      
      arrange(cam.id) %>%
      
      mutate(n.indiv = focal.passes$n.indiv[1],
             landscape = focal.passes$landscape[1],
             variability = focal.passes$variability[1],
             rep = focal.passes$rep[1],
             cams = focal.passes$cams[1])
    
    # bind
    all.passes.grouped <- rbind(all.passes.grouped, focal.passes.2)
    
  }
  
  # return
  return(all.passes.grouped)
  
}

#_______________________________________________________________________
# 4b. Use function ----
#_______________________________________________________________________

passes.gp.4 <- group_passes(passes.extracted.4, 4)
passes.gp.9 <- group_passes(passes.extracted.9, 9)
passes.gp.16 <- group_passes(passes.extracted.16, 16)

#_______________________________________________________________________
# 5. Bind together ----
#_______________________________________________________________________

passes.gp.all <- rbind(passes.gp.4, passes.gp.9, passes.gp.16)

#_______________________________________________________________________
# 6. Plot (sanity check) ----
#_______________________________________________________________________

ggplot(passes.gp.all) +
  
  theme_bw() +
  
  facet_grid(landscape ~ variability) +
  
  geom_point(aes(x = as.factor(n.indiv),
                 y = total.passes,
                 shape = as.factor(rep),
                 color = as.factor(rep)),
             alpha = 0.5) +
  
  theme(panel.grid = element_blank(),
        legend.position = "none")

#_______________________________________________________________________
# 7. Write to .csv ----
#_______________________________________________________________________
  
write.csv(passes.gp.all, paste0(getwd(), "/Derived_data/Passes/passes_gp_all.csv"))
