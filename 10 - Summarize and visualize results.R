# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 10 - Summarize and visualize results
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 29 Apr 2025
# Date completed:
# Date last modified: 29 Apr 2025
# R version: 4.4.3

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation
library(glmmTMB)
library(performance)

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

rem.1 <- read.csv("Derived data/REM results/rem_1.csv")
rem.2 <- read.csv("Derived data/REM results/rem_2.csv")

#_______________________________________________________________________
# 3. Calculate bias metrics ----

# here we'll calculate:
# percent bias: ((estimated - true) / true) * 100
# absolute percent bias: abs(((estimated - true) / true) * 100)
# CI coverage: is the true value within the 95% CI?

#_______________________________________________________________________

rem.1 <- rem.1 %>%
  
  mutate(perc.bias = ((mean.REM.D - true.D) / true.D * 100),
         abs.perc.bias = abs((mean.REM.D - true.D) / true.D * 100),
         ci.cov = ifelse(true.D >= l95.REM.D & true.D <= u95.REM.D,
                         1,
                         0))

rem.2 <- rem.2 %>%
  
  mutate(perc.bias = ((mean.REM.D - true.D) / true.D * 100),
         abs.perc.bias = abs((mean.REM.D - true.D) / true.D * 100),
         ci.cov = ifelse(true.D >= l95.REM.D & true.D <= u95.REM.D,
                         1,
                         0))

#_______________________________________________________________________
# 4. Scenario variables ----

# We'll separate out each level for comparison (modeling?) later
# recall that:
# 1: "0.5 hr, 100% , 4 wk",
# 2: "0.5 hr, 100% , 2 wk",
# 3: "0.5 hr, 60% , 4 wk",
# 4: "0.5 hr, 60% , 2 wk",
# 5: "1 hr, 100% , 4 wk",
# 6: "1 hr, 100% , 2 wk",
# 7: "1 hr, 60% , 4 wk",
# 8: "1 hr, 60% , 2 wk",
# 9: "4 hr, 100% , 4 wk",
# 10: "4 hr, 100% , 2 wk",
# 11: "4 hr, 60% , 4 wk",
# 12: "4 hr, 60% , 2 wk"

#_______________________________________________________________________

rem.1 <- rem.1 %>%
  
  mutate(
    
    # rename factor levels
    scenario.fac = factor(scenario,
                          labels = c("0.5 hr, 100% , 4 wk",
                                     "0.5 hr, 100% , 2 wk",
                                     "0.5 hr, 60% , 4 wk",
                                     "0.5 hr, 60% , 2 wk",
                                     "1 hr, 100% , 4 wk",
                                     "1 hr, 100% , 2 wk",
                                     "1 hr, 60% , 4 wk",
                                     "1 hr, 60% , 2 wk",
                                     "4 hr, 100% , 4 wk",
                                     "4 hr, 100% , 2 wk",
                                     "4 hr, 60% , 4 wk",
                                     "4 hr, 60% , 2 wk")),
    
    # fix rate
    fix.rate = case_when(scenario %in% 1:4 ~ "0.5",
                         scenario %in% 5:8 ~ "1",
                         scenario %in% 9:12 ~ "4"),
    
    # fix success
    fix.success = case_when(scenario %in% c(1, 2, 5, 6, 9, 10) ~ "100",
                            scenario %in% c(3, 4, 7, 8, 11, 12) ~ "60"),
    
    # tracking duration
    duration = case_when(scenario %in% c(1, 3, 5, 7, 9, 11) ~ "4",
                         scenario %in% c(2, 4, 6, 8, 10, 12) ~ "2"))
    
rem.2 <- rem.2 %>%
  
  mutate(
    
    # rename factor levels
    scenario.fac = factor(scenario,
                          labels = c("0.5 hr, 100%, 4 wk",
                                     "0.5 hr, 100%, 2 wk",
                                     "0.5 hr, 60%, 4 wk",
                                     "0.5 hr, 60%, 2 wk",
                                     "1 hr, 100%, 4 wk",
                                     "1 hr, 100%, 2 wk",
                                     "1 hr, 60%, 4 wk",
                                     "1 hr, 60%, 2 wk",
                                     "4 hr, 100%, 4 wk",
                                     "4 hr, 100%, 2 wk",
                                     "4 hr, 60%, 4 wk",
                                     "4 hr, 60%, 2 wk")),
    
    # fix rate
    fix.rate = case_when(scenario %in% 1:4 ~ "0.5",
                         scenario %in% 5:8 ~ "1",
                         scenario %in% 9:12 ~ "4"),
    
    # fix success
    fix.success = case_when(scenario %in% c(1, 2, 5, 6, 9, 10) ~ "100",
                            scenario %in% c(3, 4, 7, 8, 11, 12) ~ "60"),
    
    # tracking duration
    duration = case_when(scenario %in% c(1, 3, 5, 7, 9, 11) ~ "4",
                         scenario %in% c(2, 4, 6, 8, 10, 12) ~ "2"))

#_______________________________________________________________________
# 5. Explanatory models ----
#_______________________________________________________________________
# 5a. Absolute percent bias ----
#_______________________________________________________________________

# Q1
Q1.abs.bias <- glmmTMB(abs.perc.bias ~
                         fix.rate +
                         fix.success +
                         duration +
                         (1 | rep),
                       data = rem.1,
                       family = gaussian)

summary(Q1.abs.bias)
performance(Q1.abs.bias)
plot(simulate_residuals(Q1.abs.bias))

# Q2
Q2.abs.bias <- glmmTMB(abs.perc.bias ~
                         fix.rate +
                         fix.success +
                         duration +
                         (1 | rep),
                       data = rem.2,
                       family = gaussian)

summary(Q2.abs.bias)
performance(Q2.abs.bias)
plot(simulate_residuals(Q2.abs.bias))

# this is a good start! Perhaps the best thing to do is just summarize each of these
# quantities by scenario x abundance

#_______________________________________________________________________
# 5. Summary statistics ----
#_______________________________________________________________________

rem.1.summary <- rem.1 %>% 
  
  group_by(scenario, abund) %>%
  
  summarize(mean.abs.bias = mean(abs.perc.bias, na.rm = T),
            sd.abs.bias = sd(abs.perc.bias, na.rm = T),
            mean.cv = mean(cv.REM.D, na.rm = T),
            sd.cv = sd(cv.REM.D, na.rm = T))

rem.2.summary <- rem.2 %>% 
  
  group_by(scenario, abund) %>%
  
  summarize(mean.abs.bias = mean(abs.perc.bias, na.rm = T),
            sd.abs.bias = sd(abs.perc.bias, na.rm = T),
            mean.cv = mean(cv.REM.D, na.rm = T),
            sd.cv = sd(cv.REM.D, na.rm = T))

#_______________________________________________________________________
# 5. Estimated vs. true ----

# Here we'll look at each scenario in turn

#_______________________________________________________________________
# 5a. Q1 ----
#_______________________________________________________________________

ggplot(data = rem.1) +
  
  theme_bw() +
  
  facet_wrap(~ scenario.fac) +
  
  # 1:1 line
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dashed") +
  
  geom_point(aes(x = true.D,
                 y = mean.REM.D),
             alpha = 0.05,
             size = 0.5) +
  
  # coordinates
  coord_cartesian(xlim = c(0, 6),
                  ylim = c(0, 6)) +
  
  theme(panel.grid = element_blank()) +
  
  xlab("True density (individuals/ha") +
  ylab("Estimated density (individuals/ha)")

#_______________________________________________________________________
# 5b. Q2 ----
#_______________________________________________________________________

ggplot(data = rem.2) +
  
  theme_bw() +
  
  facet_wrap(~ scenario.fac) +
  
  # 1:1 line
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dashed") +
  
  geom_point(aes(x = true.D,
                 y = mean.REM.D),
             alpha = 0.05,
             size = 0.5) +
  
  # coordinates
  coord_cartesian(xlim = c(0, 7),
                  ylim = c(0, 7)) +
  
  theme(panel.grid = element_blank()) +
  
  xlab("True density (individuals/ha") +
  ylab("Estimated density (individuals/ha)")

#_______________________________________________________________________
# 6. Percent bias ----
#_______________________________________________________________________
# 6a. Q1 ----
#_______________________________________________________________________

ggplot(data = rem.1,
       aes(x = perc.bias)) +
  
  theme_bw() +
  
  facet_grid(duration ~ fix.success) +
  
  geom_vline(xintercept = 0) +
  
  geom_density(aes(color = fix.rate,
                   fill = fix.rate),
               alpha = 0.25) +
  
  theme(panel.grid = element_blank()) +
  
  scale_x_continuous(breaks = c(-50, 0, 50, 100, 150)) +
  
  xlab("Percent bias")

#_______________________________________________________________________
# 6b. Q2 ----
#_______________________________________________________________________

ggplot(data = rem.2,
       aes(x = perc.bias)) +
  
  theme_bw() +
  
  facet_grid(duration ~ fix.success) +
  
  geom_vline(xintercept = 0) +
  
  geom_density(aes(color = fix.rate,
                   fill = fix.rate),
               alpha = 0.25) +
  
  theme(panel.grid = element_blank()) +
  
  scale_x_continuous(breaks = c(-50, 0, 50, 100, 150)) +
  
  xlab("Percent bias")

#_______________________________________________________________________
# 7. Coefficients of variation ----
#_______________________________________________________________________
# 7a. Q1 ----
#_______________________________________________________________________

ggplot(data = rem.1,
       aes(x = cv.REM.D)) +
  
  theme_bw() +
  
  facet_grid(duration ~ fix.success) +
  
  geom_histogram(aes(group = fix.rate,
    color = fix.rate,
                   fill = fix.rate),
               alpha = 0.25) +
  
  theme(panel.grid = element_blank()) +
  
  xlab("Coefficient of variation")

#_______________________________________________________________________
# 7b. Q2 ----
#_______________________________________________________________________

ggplot(data = rem.2,
       aes(x = cv.REM.D)) +
  
  theme_bw() +
  
  facet_grid(duration ~ fix.success) +
  
  geom_density(aes(color = fix.rate,
                   fill = fix.rate),
               alpha = 0.25) +
  
  theme(panel.grid = element_blank()) +
  
  xlab("Coefficient of variation")
