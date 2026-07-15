# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 14 - Simulated speed performance
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 16 Jul 2025
# Date completed: 16 Jul 2025
# Date last modified: 16 Jul 2025
# R version: 4.4.3

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

# "true" speeds
true.speeds.1T.TV1 <- readRDS(paste0(getwd(), "/data_derived/sampled_speeds/speeds_1T_TV1.rds"))
true.speeds.1T.TV2 <- readRDS(paste0(getwd(), "/data_derived/sampled_speeds/speeds_1T_TV2.rds"))
true.speeds.1T.TV3 <- readRDS(paste0(getwd(), "/data_derived/sampled_speeds/speeds_1T_TV3.rds"))

true.speeds.1NT.TV1 <- readRDS(paste0(getwd(), "/data_derived/sampled_speeds/speeds_1NT_TV1.rds"))
true.speeds.1NT.TV2 <- readRDS(paste0(getwd(), "/data_derived/sampled_speeds/speeds_1NT_TV2.rds"))
true.speeds.1NT.TV3 <- readRDS(paste0(getwd(), "/data_derived/sampled_speeds/speeds_1NT_TV3.rds"))

true.speeds.2T <- readRDS(paste0(getwd(), "/data_derived/sampled_speeds/speeds_2T.rds"))
true.speeds.2NT <- readRDS(paste0(getwd(), "/data_derived/sampled_speeds/speeds_2NT.rds"))

# simulated speeds
speeds.1T.TV1 <- readRDS(paste0(getwd(), "/data_derived/sampled_ctmm_speeds/1T_TV1.rds"))
speeds.1T.TV2 <- readRDS(paste0(getwd(), "/data_derived/sampled_ctmm_speeds/1T_TV2.rds"))
speeds.1T.TV3 <- readRDS(paste0(getwd(), "/data_derived/sampled_ctmm_speeds/1T_TV3.rds"))

speeds.1NT.TV1 <- readRDS(paste0(getwd(), "/data_derived/sampled_ctmm_speeds/1NT_TV1.rds"))
speeds.1NT.TV2 <- readRDS(paste0(getwd(), "/data_derived/sampled_ctmm_speeds/1NT_TV2.rds"))
speeds.1NT.TV3 <- readRDS(paste0(getwd(), "/data_derived/sampled_ctmm_speeds/1NT_TV3.rds"))

speeds.2T <- readRDS(paste0(getwd(), "/data_derived/sampled_ctmm_speeds/2T.rds"))
speeds.2NT <- readRDS(paste0(getwd(), "/data_derived/sampled_ctmm_speeds/2NT.rds"))

#_______________________________________________________________________
# 3. Clean data ----

# we'll need to join, add identifiers, and bind

#_______________________________________________________________________

# function
clean_speeds <- function (.trueT, .trueNT, .estT, .estNT) {
  
  # start with target individuals
  all.speeds <- .trueT |>
    
    # drop unneeded columns
    dplyr::select(-c(question, target)) |>
    
    # join in estimates
    left_join(
      
      .estT |>
        
        # keep only mean estimates for simplicity
        dplyr::select(indiv, scenario, speed.mean.est, sld.speed.mean.est, rms.speed.est)
      
    ) |>
    
    # drop NAs (we can't use them anyway)
    drop_na() |>
    
    # bind in NT individuals
    bind_rows(
      
      .trueNT |>
        
        # drop unneeded columns
        dplyr::select(-c(question, target)) |>
        
        # join in estimates
        left_join(
          
          .estNT |>
            
            # keep only mean estimates for simplicity
            dplyr::select(indiv, scenario, speed.mean.est, sld.speed.mean.est, rms.speed.est)
          
        ) |>
        
        # drop NAs (we can't use them anyway)
        drop_na()
      
    ) |>
    
    mutate(
      
      # fix rate
      fix.rate = case_when(scenario %in% 1:4 ~ "0.5",
                           scenario %in% 5:8 ~ "1",
                           scenario %in% 9:12 ~ "4"),
      
      # fix success
      fix.success = factor(case_when(scenario %in% c(1, 2, 5, 6, 9, 10) ~ "100",
                                     scenario %in% c(3, 4, 7, 8, 11, 12) ~ "60"),
                           levels = c("100", "60")),
      
      # tracking duration
      duration = factor(case_when(scenario %in% c(1, 3, 5, 7, 9, 11) ~ "4",
                                  scenario %in% c(2, 4, 6, 8, 10, 12) ~ "2"),
                        levels = c("4", "2"))
      
    )
  
  return(all.speeds)
  
}

# use function
speeds.TV1 <- clean_speeds(true.speeds.1T.TV1,
                           true.speeds.1NT.TV1,
                           speeds.1T.TV1,
                           speeds.1NT.TV1) |>
  
  mutate(TV = "0.5 mins")

speeds.TV2 <- clean_speeds(true.speeds.1T.TV2,
                           true.speeds.1NT.TV2,
                           speeds.1T.TV2,
                           speeds.1NT.TV2) |>
  
  mutate(TV = "30 mins")

speeds.TV3 <- clean_speeds(true.speeds.1T.TV3,
                           true.speeds.1NT.TV3,
                           speeds.1T.TV3,
                           speeds.1NT.TV3) |>
  
  mutate(TV = "120 mins")

# bind together
speeds.1 <- rbind(speeds.TV1, speeds.TV2, speeds.TV3)

speeds.2 <- clean_speeds(true.speeds.2T,
                         true.speeds.2NT,
                         speeds.2T,
                         speeds.2NT)

#_______________________________________________________________________
# 4. Visualization ----
#_______________________________________________________________________
# 4. Bias plots ----

# these will be 1:1 line plots

# x panels: fix rate
# y panels: fix success within duration
# colors/shapes: estimated vs SLD 

# we should center speeds based on each tau v
speeds.center.TV <- speeds.1 |> 
  
  group_by(TV) |>
  
  summarize(speed.center = mean(true.speed))

#_______________________________________________________________________
# Q1 ----
#_______________________________________________________________________

# pivot
speeds.1.pivot <- speeds.1 |> 
  
  # join in centers
  left_join(speeds.center.TV) |>
  
  # center
  mutate(
    
    true.speed = true.speed - speed.center,
    speed.mean.est = speed.mean.est - speed.center,
    rms.speed.est = rms.speed.est - speed.center,
    sld.speed.mean.est = sld.speed.mean.est - speed.center 
    
  ) |>
  
  pivot_longer(cols = c(speed.mean.est,
                        rms.speed.est,
                        sld.speed.mean.est)) |>
  
  mutate(name = factor(name,
                       levels = c("speed.mean.est",
                                  "rms.speed.est",
                                  "sld.speed.mean.est"),
                       labels = c("CVM", "RMS", "SLD"))) |>
  
  # TV factor
  mutate(TV = factor(TV, levels = c("0.5 mins", "30 mins", "120 mins")))

# Q1
ggplot(data = speeds.1.pivot) +
  
  theme_bw() +
  
  # facet
  ggh4x::facet_nested(TV + duration + fix.success ~ fix.rate,
                      labeller = labeller(fix.rate = as_labeller(c("0.5" = "0.5 h",
                                                                   "1" = "1 h",
                                                                   "4" = "4 h")),
                                          duration = as_labeller(c("4" = "4 weeks",
                                                                   "2" = "2 weeks")),
                                          fix.success = as_labeller(c("100" = "100%",
                                                                      "60" = "60%")))) +
  
  # 1:1 line
  geom_abline(slope = 1,
              intercept = 0,
              linetype = "dashed") +
  
  # points
  geom_point(aes(x = true.speed,
                 y = value,
                 color = name,
                 shape = name),
             size = 0.45) +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_text(hjust = 0.01,
                                    face = "bold",
                                    size = 8),
        strip.text.y = element_text(hjust = 0.01,
                                    size = 8),
        strip.background = element_rect(fill = "white"),
        plot.margin = unit(c(0, 0, 0, 0), "pt"),
        axis.text = element_blank()) +
  
  # legend guide
  guides(color = guide_legend(override.aes = list(size = 2.5))) +
  
  scale_color_manual(values = c("palegreen4", "palegreen3", "palegreen2")) +
  
  xlab("True travel speed") +
  ylab("Estimated travel speed")

# 625 x 530

#_______________________________________________________________________
# Q1 ----
#_______________________________________________________________________

# pivot
speeds.2.pivot <- speeds.2 |> 
  
  pivot_longer(cols = c(speed.mean.est,
                        rms.speed.est,
                        sld.speed.mean.est)) |>
  
  mutate(name = factor(name,
                       levels = c("speed.mean.est",
                                  "rms.speed.est",
                                  "sld.speed.mean.est"),
                       labels = c("CVM", "RMS", "SLD")))

# Q1
ggplot(data = speeds.2.pivot) +
  
  theme_bw() +
  
  # facet
  ggh4x::facet_nested(duration + fix.success ~ fix.rate,
                      labeller = labeller(fix.rate = as_labeller(c("0.5" = "0.5 h",
                                                                   "1" = "1 h",
                                                                   "4" = "4 h")),
                                          duration = as_labeller(c("4" = "4 weeks",
                                                                   "2" = "2 weeks")),
                                          fix.success = as_labeller(c("100" = "100%",
                                                                      "60" = "60%")))) +
  
  # 1:1 line
  geom_abline(slope = 1,
              intercept = 0,
              linetype = "dashed") +
  
  # points
  geom_point(aes(x = true.speed,
                 y = value,
                 color = name,
                 shape = name),
             size = 0.45) +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_text(hjust = 0.01,
                                    face = "bold",
                                    size = 8),
        strip.text.y = element_text(hjust = 0.01,
                                    size = 8),
        strip.background = element_rect(fill = "white"),
        plot.margin = unit(c(0, 0, 0, 0), "pt"),
        axis.text = element_blank()) +
  
  # legend guide
  guides(color = guide_legend(override.aes = list(size = 2.5))) +
  
  scale_color_manual(values = c("salmon4", "salmon3", "salmon1")) +
  
  xlab("True travel speed") +
  ylab("Estimated travel speed")

# 471 x 444