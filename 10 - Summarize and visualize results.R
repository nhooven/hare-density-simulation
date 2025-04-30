# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 10 - Summarize and visualize results
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 29 Apr 2025
# Date completed: 30 Apr 2025
# Date last modified: 29 Apr 2025
# R version: 4.4.3

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation

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
    
    # density as factor
    true.D = factor(true.D),
    
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
                      levels = c("4", "2")))
    
rem.2 <- rem.2 %>%
  
  mutate(
    
    # density as factor
    true.D = factor(true.D),
    
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
                      levels = c("4", "2")))

#_______________________________________________________________________
# 5. Calculate mean performance by scenario ----

# Here we'll use bootstrapping to estimate CIs on the mean bias, precision, and coverage
# across each scenario, and visualize them in informative plots

# define function
boot_perform <- function (df,
                          metric,
                          n.iter = 1000) {
  
  all.combos <- expand.grid(true.D = c("0.4", "0.8", "1.6", "3.2"),
                            fix.rate = c("0.5", "1", "4"),
                            fix.success = c("100", "60"),
                            duration = c("4", "2"))
  
  if (metric == "bias") {
    
    all.return.df <- data.frame()
    
    for (i in 1:nrow(all.combos)) {
      
      focal.combo <- all.combos[i, ]
      
      # subset df
      focal.df <- df %>% filter(true.D == focal.combo$true.D,
                                fix.rate == focal.combo$fix.rate,
                                fix.success == focal.combo$fix.success,
                                duration == focal.combo$duration)
      
      # mean
      mean.metric <- mean(focal.df$abs.perc.bias)
      
      # bootstrapping CIs
      matrix.j <- matrix(data = NA,
                         nrow = n.iter,
                         ncol = 1)
      
      for (j in 1:n.iter) {
        
        matrix.j[j, ] <- mean(sample(focal.df$abs.perc.bias, size = nrow(focal.df), replace = T))
        
      }
      
      # confidence intervals
      l95.metric <- quantile(matrix.j[ , 1], probs = 0.025)
      u95.metric <- quantile(matrix.j[ , 1], probs = 0.975)
      
      # bind into df
      metric.df <- cbind(focal.combo, mean.metric, l95.metric, u95.metric)
      
      all.return.df <- rbind(all.return.df, metric.df)
      
    }
    
  }
  
  if (metric == "precision") {
    
    all.return.df <- data.frame()
    
    for (i in 1:nrow(all.combos)) {
      
      focal.combo <- all.combos[i, ]
      
      # subset df
      focal.df <- df %>% filter(true.D == focal.combo$true.D,
                                fix.rate == focal.combo$fix.rate,
                                fix.success == focal.combo$fix.success,
                                duration == focal.combo$duration)
      
      # mean
      mean.metric <- mean(focal.df$cv.REM.D)
      
      # bootstrapping CIs
      matrix.j <- matrix(data = NA,
                         nrow = n.iter,
                         ncol = 1)
      
      for (j in 1:n.iter) {
        
        matrix.j[j, ] <- mean(sample(focal.df$cv.REM.D, size = nrow(focal.df), replace = T))
        
      }
      
      # confidence intervals
      l95.metric <- quantile(matrix.j[ , 1], probs = 0.025)
      u95.metric <- quantile(matrix.j[ , 1], probs = 0.975)
      
      # bind into df
      metric.df <- cbind(focal.combo, mean.metric, l95.metric, u95.metric)
      
      all.return.df <- rbind(all.return.df, metric.df)
      
    }
    
  }
  
  if (metric == "coverage") {
    
    all.return.df <- data.frame()
    
    for (i in 1:nrow(all.combos)) {
      
      focal.combo <- all.combos[i, ]
      
      # subset df
      focal.df <- df %>% filter(true.D == focal.combo$true.D,
                                fix.rate == focal.combo$fix.rate,
                                fix.success == focal.combo$fix.success,
                                duration == focal.combo$duration)
      
      # mean
      mean.metric <- sum(focal.df$ci.cov) / nrow(focal.df)
      
      # bootstrapping CIs
      matrix.j <- matrix(data = NA,
                         nrow = n.iter,
                         ncol = 1)
      
      for (j in 1:n.iter) {
        
        matrix.j[j, ] <- sum(sample(focal.df$ci.cov, size = nrow(focal.df), replace = T)) / nrow(focal.df)
        
      }
      
      # confidence intervals
      l95.metric <- quantile(matrix.j[ , 1], probs = 0.025)
      u95.metric <- quantile(matrix.j[ , 1], probs = 0.975)
      
      # bind into df
      metric.df <- cbind(focal.combo, mean.metric, l95.metric, u95.metric)
      
      all.return.df <- rbind(all.return.df, metric.df)
      
    }
    
  }
  
  # return
  return(all.return.df)
  
}

#_______________________________________________________________________
# 5a. Bias (absolute percent bias) ----
#_______________________________________________________________________

boot.bias.1 <- boot_perform(rem.1, "bias")
boot.bias.2 <- boot_perform(rem.2, "bias")

#_______________________________________________________________________
# 5b. Precision (coefficient of variation) ----
#_______________________________________________________________________

boot.precis.1 <- boot_perform(rem.1, "precision")
boot.precis.2 <- boot_perform(rem.2, "precision")

#_______________________________________________________________________
# 5c. Coverage (95% confidence interval coverage) ----
#_______________________________________________________________________

boot.cov.1 <- boot_perform(rem.1, "coverage")
boot.cov.2 <- boot_perform(rem.2, "coverage")

#_______________________________________________________________________
# 6. Visualize performance ----

# I originally included all the data, but it's not possible to see the CIs

#_______________________________________________________________________
# 6a. Bias ----
#_______________________________________________________________________

# Q1
ggplot(data = rem.1) +
  
  theme_bw() +
  
  # each density as columns, duration as rows
  # the nested function is neat!
  ggh4x::facet_nested(duration + fix.success ~ true.D,
                      scales = "free_y",
                      labeller = labeller(true.D = as_labeller(c("0.4" = "D = 0.4",
                                                                 "0.8" = "D = 0.8",
                                                                 "1.6" = "D = 1.6",
                                                                 "3.2" = "D = 3.2")),
                                          duration = as_labeller(c("4" = "4 weeks",
                                                                   "2" = "2 weeks")),
                                          fix.success = as_labeller(c("100" = "100%",
                                                                      "60" = "60%")))) +
  
  # CIs
  geom_errorbarh(data = boot.bias.1,
                 aes(y = interaction(fix.rate, fix.success),
                     xmin = l95.metric,
                     xmax = u95.metric,
                     color = fix.rate),
                 linewidth = 1.25,
                 height = 0) +
  
  # means
  geom_point(data = boot.bias.1,
             aes(x = mean.metric,
                 y = interaction(fix.rate, fix.success),
                 shape = fix.success,
                 size = fix.success,
                 fill = fix.rate)) +
  
  # theme arguments
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
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
        axis.title = element_text(size = 10)) +
  
  # labels
  xlab("|Percent bias|") +
  ylab("Fix rate (h)") +
  
  # shapes, sizes, and colors
  scale_shape_manual(values = c(21, 22)) +
  scale_size_manual(values = c(1.35, 1.35)) +
  scale_color_manual(values = c("palegreen4", "palegreen3", "palegreen2")) +
  scale_fill_manual(values = c("palegreen4", "palegreen3", "palegreen2")) +
  
  # reverse the y axis
  scale_y_discrete(limits = rev,
                   labels = c(rep(c("4", "1", "0.5"), 
                                  times = 4))) +
  
  # sensible breaks
  scale_x_continuous(breaks = c(20, 30, 40))

# 632 x 244

# Q2
ggplot(data = rem.2) +
  
  theme_bw() +
  
  # each density as columns, duration as rows
  # the nested function is neat!
  ggh4x::facet_nested(duration + fix.success ~ true.D,
                      scales = "free_y",
                      labeller = labeller(true.D = as_labeller(c("0.4" = "D = 0.4",
                                                                 "0.8" = "D = 0.8",
                                                                 "1.6" = "D = 1.6",
                                                                 "3.2" = "D = 3.2")),
                                          duration = as_labeller(c("4" = "4 weeks",
                                                                   "2" = "2 weeks")),
                                          fix.success = as_labeller(c("100" = "100%",
                                                                      "60" = "60%")))) +
  
  # CIs
  geom_errorbarh(data = boot.bias.2,
                 aes(y = interaction(fix.rate, fix.success),
                     xmin = l95.metric,
                     xmax = u95.metric,
                     color = fix.rate),
                 linewidth = 1.25,
                 height = 0) +
  
  # means
  geom_point(data = boot.bias.2,
             aes(x = mean.metric,
                 y = interaction(fix.rate, fix.success),
                 shape = fix.success,
                 size = fix.success,
                 fill = fix.rate)) +
  
  # theme arguments
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
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
        axis.title = element_text(size = 10)) +
  
  # labels
  xlab("|Percent bias|") +
  ylab("Fix rate (h)") +
  
  # shapes, sizes, and colors
  scale_shape_manual(values = c(21, 22)) +
  scale_size_manual(values = c(1.35, 1.35)) +
  scale_color_manual(values = c("salmon4", "salmon3", "salmon2")) +
  scale_fill_manual(values = c("salmon4", "salmon3", "salmon2")) +
  
  # reverse the y axis
  scale_y_discrete(limits = rev,
                   labels = c(rep(c("4", "1", "0.5"), 
                                  times = 4))) +
  
  # sensible breaks
  scale_x_continuous(breaks = c(40, 50))

# 632 x 244

#_______________________________________________________________________
# 6b. Precision ----
#_______________________________________________________________________

# Q1
ggplot(data = rem.1) +
  
  theme_bw() +
  
  # each density as columns, duration as rows
  # the nested function is neat!
  ggh4x::facet_nested(duration + fix.success ~ true.D,
                      scales = "free",
                      labeller = labeller(true.D = as_labeller(c("0.4" = "D = 0.4",
                                                                 "0.8" = "D = 0.8",
                                                                 "1.6" = "D = 1.6",
                                                                 "3.2" = "D = 3.2")),
                                          duration = as_labeller(c("4" = "4 weeks",
                                                                   "2" = "2 weeks")),
                                          fix.success = as_labeller(c("100" = "100%",
                                                                      "60" = "60%")))) +
  
  # CIs
  geom_errorbarh(data = boot.precis.1,
                 aes(y = interaction(fix.rate, fix.success),
                     xmin = l95.metric,
                     xmax = u95.metric,
                     color = fix.rate),
                 linewidth = 1.25,
                 height = 0) +
  
  # means
  geom_point(data = boot.precis.1,
             aes(x = mean.metric,
                 y = interaction(fix.rate, fix.success),
                 shape = fix.success,
                 size = fix.success,
                 fill = fix.rate)) +
  
  # theme arguments
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
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
        axis.title = element_text(size = 10)) +
  
  # labels
  xlab("Coefficient of variation") +
  ylab("Fix rate (h)") +
  
  # shapes, sizes, and colors
  scale_shape_manual(values = c(21, 22)) +
  scale_size_manual(values = c(1.35, 1.35)) +
  scale_color_manual(values = c("palegreen4", "palegreen3", "palegreen2")) +
  scale_fill_manual(values = c("palegreen4", "palegreen3", "palegreen2")) +
  
  # reverse the y axis
  scale_y_discrete(limits = rev,
                   labels = c(rep(c("4", "1", "0.5"), 
                                  times = 4))) +
  
  # reasonable breaks
  scale_x_continuous(breaks = seq(0, 0.35, 0.01)[-20])

# 632 x 244

# Q2
ggplot(data = rem.2) +
  
  theme_bw() +
  
  # each density as columns, duration as rows
  # the nested function is neat!
  ggh4x::facet_nested(duration + fix.success ~ true.D,
                      scales = "free",    # free scales here
                      labeller = labeller(true.D = as_labeller(c("0.4" = "D = 0.4",
                                                                 "0.8" = "D = 0.8",
                                                                 "1.6" = "D = 1.6",
                                                                 "3.2" = "D = 3.2")),
                                          duration = as_labeller(c("4" = "4 weeks",
                                                                   "2" = "2 weeks")),
                                          fix.success = as_labeller(c("100" = "100%",
                                                                      "60" = "60%")))) +
  
  # CIs
  geom_errorbarh(data = boot.precis.2,
                 aes(y = interaction(fix.rate, fix.success),
                     xmin = l95.metric,
                     xmax = u95.metric,
                     color = fix.rate),
                 linewidth = 1.25,
                 height = 0) +
  
  # means
  geom_point(data = boot.precis.2,
             aes(x = mean.metric,
                 y = interaction(fix.rate, fix.success),
                 shape = fix.success,
                 size = fix.success,
                 fill = fix.rate)) +
  
  # theme arguments
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
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
        axis.title = element_text(size = 10)) +
  
  # labels
  xlab("Coefficient of variation") +
  ylab("Fix rate (h)") +
  
  # shapes, sizes, and colors
  scale_shape_manual(values = c(21, 22)) +
  scale_size_manual(values = c(1.35, 1.35)) +
  scale_color_manual(values = c("salmon4", "salmon3", "salmon2")) +
  scale_fill_manual(values = c("salmon4", "salmon3", "salmon2")) +
  
  # reverse the y axis
  scale_y_discrete(limits = rev,
                   labels = c(rep(c("4", "1", "0.5"), 
                                  times = 4))) +
  
  scale_x_continuous(breaks = seq(0, 0.335, 0.005))


# 632 x 244

#_______________________________________________________________________
# 6c. Coverage ----

# this is the percent of replicates whose 95% CIs include the true value

#_______________________________________________________________________

# Q1
ggplot(data = rem.1) +
  
  theme_bw() +
  
  # each density as columns, duration as rows
  # the nested function is neat!
  ggh4x::facet_nested(duration + fix.success ~ true.D,
                      scales = "free_y",
                      labeller = labeller(true.D = as_labeller(c("0.4" = "D = 0.4",
                                                                 "0.8" = "D = 0.8",
                                                                 "1.6" = "D = 1.6",
                                                                 "3.2" = "D = 3.2")),
                                          duration = as_labeller(c("4" = "4 weeks",
                                                                   "2" = "2 weeks")),
                                          fix.success = as_labeller(c("100" = "100%",
                                                                      "60" = "60%")))) +
  
  # CIs
  geom_errorbarh(data = boot.cov.1,
                 aes(y = interaction(fix.rate, fix.success),
                     xmin = l95.metric,
                     xmax = u95.metric,
                     color = fix.rate),
                 linewidth = 1.25,
                 height = 0) +
  
  # means
  geom_point(data = boot.cov.1,
             aes(x = mean.metric,
                 y = interaction(fix.rate, fix.success),
                 shape = fix.success,
                 size = fix.success,
                 fill = fix.rate)) +
  
  # theme arguments
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
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
        axis.title = element_text(size = 10)) +
  
  # labels
  xlab("Proportion of 95% CIs that include true value") +
  ylab("Fix rate (h)") +
  
  # shapes, sizes, and colors
  scale_shape_manual(values = c(21, 22)) +
  scale_size_manual(values = c(1.35, 1.35)) +
  scale_color_manual(values = c("palegreen4", "palegreen3", "palegreen2")) +
  scale_fill_manual(values = c("palegreen4", "palegreen3", "palegreen2")) +
  
  # reverse the y axis
  scale_y_discrete(limits = rev,
                   labels = c(rep(c("4", "1", "0.5"), 
                                  times = 4)))

# 632 x 244

# Q2
ggplot(data = rem.2) +
  
  theme_bw() +
  
  # each density as columns, duration as rows
  # the nested function is neat!
  ggh4x::facet_nested(duration + fix.success ~ true.D,
                      scales = "free_y",
                      labeller = labeller(true.D = as_labeller(c("0.4" = "D = 0.4",
                                                                 "0.8" = "D = 0.8",
                                                                 "1.6" = "D = 1.6",
                                                                 "3.2" = "D = 3.2")),
                                          duration = as_labeller(c("4" = "4 weeks",
                                                                   "2" = "2 weeks")),
                                          fix.success = as_labeller(c("100" = "100%",
                                                                      "60" = "60%")))) +
  
  # CIs
  geom_errorbarh(data = boot.cov.2,
                 aes(y = interaction(fix.rate, fix.success),
                     xmin = l95.metric,
                     xmax = u95.metric,
                     color = fix.rate),
                 linewidth = 1.25,
                 height = 0) +
  
  # means
  geom_point(data = boot.cov.2,
             aes(x = mean.metric,
                 y = interaction(fix.rate, fix.success),
                 shape = fix.success,
                 size = fix.success,
                 fill = fix.rate)) +
  
  # theme arguments
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
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
        axis.title = element_text(size = 10)) +
  
  # labels
  xlab("Proportion of 95% CIs that include true value") +
  ylab("Fix rate (h)") +
  
  # shapes, sizes, and colors
  scale_shape_manual(values = c(21, 22)) +
  scale_size_manual(values = c(1.35, 1.35)) +
  scale_color_manual(values = c("salmon4", "salmon3", "salmon2")) +
  scale_fill_manual(values = c("salmon4", "salmon3", "salmon2")) +
  
  # reverse the y axis
  scale_y_discrete(limits = rev,
                   labels = c(rep(c("4", "1", "0.5"), 
                                  times = 4))) +
  
  scale_x_continuous(breaks = c(0.45, 0.55, 0.65, 0.75, 0.85))

#_______________________________________________________________________
# 7. Write tables ----
#_______________________________________________________________________

# if necessary for publication

#_______________________________________________________________________
# 8. Save RData ----
#_______________________________________________________________________

save.image(file = paste0(getwd(), "/Derived data/Plot files/04_30_2025.RData"))
