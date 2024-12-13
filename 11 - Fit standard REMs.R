# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 11 - Fit standard REMs
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 13 Dec 2024
# Date completed: 
# Date last modified: 13 Dec 2024
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation
library(boot)            # nonparametric bootstrap

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

detec <- read.csv(paste0(getwd(), "/Derived_data/Detections/detec_all.csv"))

cams <- read.csv(paste0(getwd(), "/Derived_data/Detections/cams.csv"))

day.range <- read.csv(paste0(getwd(), "/Derived_data/Model parameters/day_range.csv"))

# add zeroes for SL.2
detec <- rbind(detec,
               data.frame(X = NA,
                          cam.id = 1:9, 
                          total.passes = 0,
                          n.indiv = 2,
                          landscape = "simple",
                          variability = "low"))

#_______________________________________________________________________
# 3. Format data correctly ----
#_______________________________________________________________________
# 3a. Define function ----
#_______________________________________________________________________

add_cam_data <- function (ls, var) {
  
  # subset data
  x.1 <- detec %>% 
    
    filter(landscape == ls, variability == var) %>%
    
    dplyr::select(-X)
  
  # identifier
  ls.var.id <- case_when(ls == "simple" & var == "low" ~ "SL",
                         ls == "complex" & var == "low" ~ "CL",
                         ls == "simple" & var == "high" ~ "SH",
                         ls == "complex" & var == "high" ~ "CH")
  
  # add correct camera data
  # ssf (naive)
  ssf.data <- cams %>% dplyr::select(paste0("ssf.", ls.var.id)) %>% rename("ssf" = paste0("ssf.", ls.var.id))
  
  ssf.data <- rbind(ssf.data, ssf.data, ssf.data)
  
  x.1$ssf.data <- ssf.data$ssf
  
  # issf
  issf.data <- cams %>% dplyr::select(paste0("issf.", ls.var.id)) %>% rename("issf" = paste0("issf.", ls.var.id))
  
  issf.data <- rbind(issf.data, issf.data, issf.data)
  
  x.1$issf.data <- issf.data$issf
  
  # add correct day range
  day.range.focal <- day.range %>%
    
    filter(landscape == ls, 
           variability == var) %>%
    
    dplyr::select(day.range)
  
  # add unique "session" identifier (convert to integer later for Stan)
  x.1 <- x.1 %>% 
    
    mutate(session.unique = paste0(ls.var.id, ".", x.1$n.indiv),
           day.range = day.range.focal$day.range[1])
  
  # return
  return(x.1)
  
}

#_______________________________________________________________________
# 3b. Use function ----
#_______________________________________________________________________

detec.1 <- rbind(add_cam_data("simple", "low"),
                 add_cam_data("complex", "low"),
                 add_cam_data("simple", "high"),
                 add_cam_data("complex", "high"))

#_______________________________________________________________________
# 4. Prepare dataset for modeling ----
#_______________________________________________________________________

detec.2 <- detec.1 %>%
  
  mutate(lens = (57.3 * pi) / 180,
         days = 28,
         r.star = 3.5) %>%
  
  # rename
  rename(total.detections = total.passes)

#_______________________________________________________________________
# 5. Fit REMs with bootstrapped confidence intervals ----
#_______________________________________________________________________
# 5a. Define function ----
#_______________________________________________________________________

REM_boot <- function (df, scenario) {
  
  # calculate D by camera (in individuals/ha)
  if (scenario == "naive") {
    
    df$D <- (df$total.detections / df$days) *
             (pi / ((df$day.range * 1000) * (df$r.star) * (2.0 + df$lens))) *
             10000
    
  } else {
    
    if (scenario == "ssf") {
      
      df$D <- (df$total.detections / df$days) *
              (pi / ((df$day.range * 1000) * (df$r.star) * (2.0 + df$lens))) *
              10000 *
              (1 / df$ssf.data)
      
    } else {
      
      df$D <- (df$total.detections / df$days) *
              (pi / ((df$day.range * 1000) * (df$r.star) * (2.0 + df$lens))) *
              10000 *
              (1 / df$issf.data)
      
    }
    
  }
  
  # calculate means and CVs
  df.grouped <- df %>%
    
    group_by(landscape, variability, n.indiv) %>%
    
    summarize(D.mean = mean(D),
              D.CV = sd(D) / mean(D)) %>%
    
    # add scenario id
    mutate(scenario = scenario) %>%
    
    ungroup()
  
  # bootstrap
  df.summary <- data.frame()
  
  for (i in 1:nrow(df.grouped)) {
    
    # subset
    x <- df.grouped %>% slice(i)
    
    # define data for resampling
    y <- detec.2 %>% 
      
      filter(landscape == x$landscape[1] &
             variability == x$variability[1] &
             n.indiv == x$n.indiv[1])
    
    # determine function to use
    if (scenario == "naive") {
      
      rem_func <- function (data, indices) {
        
        df <- data[indices, ]
        
        D <- (df$total.detections / df$days) *
             (pi / ((df$day.range * 1000) * (df$r.star) * (2.0 + df$lens))) *
             10000
        
        return(mean(D))
        
      }
      
    } else {
      
      if (scenario == "ssf") {
        
        rem_func <- function (data, indices) {
          
          df <- data[indices, ]
          
          D <- (df$total.detections / df$days) *
               (pi / ((df$day.range * 1000) * (df$r.star) * (2.0 + df$lens))) *
               10000 *
               (1 / df$ssf.data)
          
          return(mean(D))
          
        }
        
      } else {
        
        rem_func <- function (data, indices) {
          
          df <- data[indices, ]
          
          D <- (df$total.detections / df$days) *
               (pi / ((df$day.range * 1000) * (df$r.star) * (2.0 + df$lens))) *
               10000 *
               (1 / df$issf.data)
          
          return(mean(D))
          
        }
        
      }
      
    }
    
  # use function within bootstrap
  boot.obj <- boot.ci(boot.out = boot(data = y, statistic = rem_func, R = 5000), type = "perc")
  
  # add into df
  x$low.95 <- boot.obj$percent[4]
  x$upp.95 <- boot.obj$percent[5]
  
  # bind in to df
  df.summary <- rbind(df.summary, x)
  
  }
  
  # return
  return(df.summary)
  
}

#_______________________________________________________________________
# 5b. Use function ----
#_______________________________________________________________________

rem.naive <- REM_boot(detec.2, "naive")
rem.ssf <- REM_boot(detec.2, "ssf")
rem.issf <- REM_boot(detec.2, "issf")

# bind together
rem.all <- rbind(rem.naive, rem.ssf, rem.issf)

#_______________________________________________________________________
# 6. Plot results ----
#_______________________________________________________________________

# reorder factors
rem.all <- rem.all %>%
  
  mutate(scenario = factor(scenario,
                           levels = c("naive", "ssf", "issf")),
         landscape = factor(landscape,
                            levels = c("simple", "complex")),
         variability = factor(variability,
                              levels = c("low", "high")))

ggplot(data = rem.all) +
  
  theme_bw() +
  
  facet_grid(scenario ~ landscape*variability) +
  
  coord_cartesian(xlim = c(0, 2),
                  ylim = c(0, 2)) +
  
  geom_abline(slope = 1, 
              intercept = 0,
              linetype = "dashed") +
  
  geom_errorbarh(aes(y = n.indiv / 10,
                     xmin = low.95,
                     xmax = upp.95),
                 height = 0) +
  
  geom_point(aes(x = D.mean,
                 y = n.indiv / 10,
                 fill = as.factor(n.indiv)),
             shape = 21,
             size = 1.5) +
  
  scale_fill_viridis_d(option = "magma") +
  
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  
  xlab("Estimated density (individuals/ha)") +
  ylab("True density (individuals/ha)")



# NEXT STEP
# Add bias/coverage metrics to really compare these

