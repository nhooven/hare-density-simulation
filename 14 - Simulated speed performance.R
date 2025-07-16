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
speeds.1T <- read.csv(paste0(getwd(), "/Derived data/Sampled - Speeds/speeds_1T.csv"))
speeds.1NT <- read.csv(paste0(getwd(), "/Derived data/Sampled - Speeds/speeds_1NT.csv"))
speeds.2T <- read.csv(paste0(getwd(), "/Derived data/Sampled - Speeds/speeds_2T.csv"))
speeds.2NT <- read.csv(paste0(getwd(), "/Derived data/Sampled - Speeds/speeds_2NT.csv"))

# ctmm speed csvs
ctmm.speeds.1T <- read.csv(paste0(getwd(), "/Derived data/Sampled - CTMM speeds/ctmm_speeds_1T.csv"))
ctmm.speeds.1NT <- read.csv(paste0(getwd(), "/Derived data/Sampled - CTMM speeds/ctmm_speeds_1NT.csv"))
ctmm.speeds.2T <- read.csv(paste0(getwd(), "/Derived data/Sampled - CTMM speeds/ctmm_speeds_2T.csv"))
ctmm.speeds.2T.2 <- read.csv(paste0(getwd(), "/Derived data/Sampled - CTMM speeds/ctmm_speeds_2T_2.csv"))
ctmm.speeds.2NT <- read.csv(paste0(getwd(), "/Derived data/Sampled - CTMM speeds/ctmm_speeds_2NT.csv"))

#_______________________________________________________________________
# 3. Clean data ----
#_______________________________________________________________________
# 3a. Add identifiers ----
#_______________________________________________________________________

# Q1
ctmm.speeds.1T <- ctmm.speeds.1T %>%
  
  mutate(dist = "typical",
         target = "T")

ctmm.speeds.1NT <- ctmm.speeds.1NT %>%
  
  mutate(dist = "typical",
         target = "NT")

# Q2
ctmm.speeds.2T <- ctmm.speeds.2T %>%
  
  mutate(dist = "full",
         target = "T")

ctmm.speeds.2T.2 <- ctmm.speeds.2T.2 %>%
  
  mutate(dist = "full",
         target = "T")

ctmm.speeds.2NT <- ctmm.speeds.2NT %>%
  
  mutate(dist = "full",
         target = "NT")

# bind together
ctmm.speeds.2T <- rbind(ctmm.speeds.2T, ctmm.speeds.2T.2)

#_______________________________________________________________________
# 3b. Join in true speeds ----
#_______________________________________________________________________

# Q1 - T
ctmm.speeds.1T <- speeds.1T %>%
  
  dplyr::select(indiv,
                true.speed) %>%
  
  right_join(ctmm.speeds.1T)

# Q1 - NT
ctmm.speeds.1NT <- speeds.1NT %>%
  
  dplyr::select(indiv,
                true.speed) %>%
  
  right_join(ctmm.speeds.1NT)

# Q2 - T
ctmm.speeds.2T <- speeds.2T %>%
  
  dplyr::select(indiv,
                true.speed) %>%
  
  right_join(ctmm.speeds.2T)

# Q2 - NT
ctmm.speeds.2NT <- speeds.2NT %>%
  
  dplyr::select(indiv,
                true.speed) %>%
  
  right_join(ctmm.speeds.2NT)

#_______________________________________________________________________
# 3c. Bind together ----
#_______________________________________________________________________

ctmm.speeds.1 <- rbind(ctmm.speeds.1T, ctmm.speeds.1NT)
ctmm.speeds.2 <- rbind(ctmm.speeds.2T, ctmm.speeds.2NT)

#_______________________________________________________________________
# 3. Add scenario identifiers ----
#_______________________________________________________________________

ctmm.speeds.1 <- ctmm.speeds.1 %>%
  
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

ctmm.speeds.2 <- ctmm.speeds.2 %>%
  
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

#_______________________________________________________________________
# 4. Visualization ----
#_______________________________________________________________________
# 4. Bias plots ----

# these will be 1:1 line plots

# x panels: fix rate
# y panels: fix success within duration
# colors/shapes: estimated vs SLD 

#_______________________________________________________________________

# Q1
ggplot(data = ctmm.speeds.1) +
  
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
                 y = speed.mean.est),
             color = "palegreen4",
             alpha = 0.25,
             size = 0.75) +
  
  # points
  geom_point(aes(x = true.speed,
                 y = sld.speed.mean.est),
             color = "palegreen3",
             alpha = 0.25,
             size = 0.55,
             shape = 2) +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_text(hjust = 0.01,
                                    face = "bold",
                                    size = 8),
        strip.text.y = element_text(hjust = 0.01,
                                    size = 8),
        strip.background = element_rect(fill = "white"),
        axis.text = element_text(size = 7,
                                color = "black"),
        plot.margin = unit(c(0, 0, 0, 0), "pt")) +
  
  xlab("True travel speed (m/s)") +
  ylab("Estimated travel speed (m/s)") -> agreement.Q1.plot

agreement.Q1.plot

# 600 x 600

# Q2
ggplot(data = ctmm.speeds.2) +
  
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
                 y = speed.mean.est),
             color = "salmon4",
             alpha = 0.25,
             size = 0.75) +
  
  # points
  geom_point(aes(x = true.speed,
                 y = sld.speed.mean.est),
             color = "salmon3",
             alpha = 0.25,
             size = 0.55,
             shape = 2) +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_text(hjust = 0.01,
                                    face = "bold",
                                    size = 8),
        strip.text.y = element_text(hjust = 0.01,
                                    size = 8),
        strip.background = element_rect(fill = "white"),
        axis.text = element_text(size = 7,
                                 color = "black"),
        plot.margin = unit(c(0, 0, 0, 0), "pt")) +
  
  xlab("True travel speed (m/s)") +
  ylab("Estimated travel speed (m/s)") -> agreement.Q2.plot

agreement.Q2.plot

# 600 x 600

#_______________________________________________________________________
# 4b. Count of each top model ----

# sensible labels
ctmm.speeds.1$top.model <- factor(ctmm.speeds.1$top.model,
                                  levels = c("IID anisotropic",
                                             "OU anisotropic" ,
                                             "OUF anisotropic",
                                             "OUf anisotropic"),
                                  labels = c("IID", "OU", "OUF", "OUf"))

# sensible labels
ctmm.speeds.2$top.model <- factor(ctmm.speeds.2$top.model,
                                  levels = c("IID",
                                             "IID anisotropic",
                                             "OU",
                                             "OU anisotropic" ,
                                             "OUF",
                                             "OUF anisotropic",
                                             "OUf",
                                             "OUf anisotropic",
                                             "OUΩ anisotropic"),
                                  labels = c("IID",
                                             "IID",
                                             "OU", 
                                             "OU",
                                             "OUF",
                                             "OUF", 
                                             "OUf",
                                             "OUf",
                                             "OUΩ"))

#_______________________________________________________________________

# Q1
ggplot(data = ctmm.speeds.1) +
  
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
  
  # points
  geom_bar(aes(x = top.model),
             fill = "palegreen4",
           color = "black",
             stat = "count",
             size = 0.75) +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_text(hjust = 0.01,
                                    face = "bold",
                                    size = 8),
        strip.text.y = element_text(hjust = 0.01,
                                    size = 8),
        strip.background = element_rect(fill = "white"),
        axis.text.y = element_text(size = 7,
                                   color = "black"),
        axis.text.x = element_text(size = 7,
                                   color = "black",
                                   angle = 35,
                                   hjust = 1),
        plot.margin = unit(c(0, 0, 0, 0), "pt"),
        axis.title = element_blank())

# 600 x 600

# Q2
ggplot(data = ctmm.speeds.2) +
  
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
  
  # points
  geom_bar(aes(x = top.model),
           fill = "salmon4",
           color = "black",
           stat = "count",
           size = 0.75) +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_text(hjust = 0.01,
                                    face = "bold",
                                    size = 8),
        strip.text.y = element_text(hjust = 0.01,
                                    size = 8),
        strip.background = element_rect(fill = "white"),
        axis.text.y = element_text(size = 7,
                                   color = "black"),
        axis.text.x = element_text(size = 7,
                                   color = "black",
                                   angle = 35,
                                   hjust = 1),
        plot.margin = unit(c(0, 0, 0, 0), "pt"),
        axis.title = element_blank())

# 600 x 600

#_______________________________________________________________________
# 5. Summaries ----
#_______________________________________________________________________
# 5a. Bias from true speeds and counts (proportions) of each model type ----
#_______________________________________________________________________

# Q1 
Q1.summary.bias <- ctmm.speeds.1 %>%
  
  # calculate bias
  mutate(percent.bias.ctmm = ((speed.mean.est - true.speed) / true.speed * 100),
         percent.bias.sld = ((sld.speed.mean.est - true.speed) / true.speed * 100)) %>%
  
  group_by(scenario) %>% 
  
  summarize(mean.percent.bias.ctmm = mean(percent.bias.ctmm, na.rm = T),
            sd.percent.bias.ctmm = sd(percent.bias.ctmm, na.rm = T),
            mean.percent.bias.sld = mean(percent.bias.sld, na.rm = T),
            sd.percent.bias.sld = sd(percent.bias.sld, na.rm = T))

Q1.summary.bias

# Q12
Q2.summary.bias <- ctmm.speeds.2 %>%
  
  # calculate bias
  mutate(percent.bias.ctmm = ((speed.mean.est - true.speed) / true.speed * 100),
         percent.bias.sld = ((sld.speed.mean.est - true.speed) / true.speed * 100)) %>%
  
  group_by(scenario) %>% 
  
  summarize(mean.percent.bias.ctmm = mean(percent.bias.ctmm, na.rm = T),
            sd.percent.bias.ctmm = sd(percent.bias.ctmm, na.rm = T),
            mean.percent.bias.sld = mean(percent.bias.sld, na.rm = T),
            sd.percent.bias.sld = sd(percent.bias.sld, na.rm = T))

Q2.summary.bias

# write tables
write.table(Q1.summary.bias, "clipboard", sep = "\t")
write.table(Q2.summary.bias, "clipboard", sep = "\t")

#_______________________________________________________________________
# 5b. Countsof each model type ----
#_______________________________________________________________________

# Q1 
Q1.summary.count <- ctmm.speeds.1 %>%
  
  group_by(scenario) %>% 
  
  count(top.model) %>%
  
  pivot_wider(names_from = top.model,
              values_from = n)

Q1.summary.count

# Q2
Q2.summary.count <- ctmm.speeds.2 %>%
  
  group_by(scenario) %>% 
  
  count(top.model) %>%
  
  pivot_wider(names_from = top.model,
              values_from = n)

Q2.summary.count

# write tables
write.table(Q1.summary.count, "clipboard", sep = "\t")
write.table(Q2.summary.count, "clipboard", sep = "\t")

# these are a lot nicer than the plots (plots are probably unneccessary)