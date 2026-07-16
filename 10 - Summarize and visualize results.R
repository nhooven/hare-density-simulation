# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 10 - Summarize and visualize results
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 29 Apr 2025
# Date completed: 30 Apr 2025
# Date last modified: 15 Jul 2026
# R version: 4.4.3

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation
library(cowplot)

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

rem.1.TV1.mean <- readRDS("data_derived/REM_results/1_TV1_mean.rds")
rem.1.TV1.rms <- readRDS("data_derived/REM_results/1_TV1_rms.rds")

rem.1.TV2.mean <- readRDS("data_derived/REM_results/1_TV2_mean.rds")
rem.1.TV2.rms <- readRDS("data_derived/REM_results/1_TV2_rms.rds")

rem.1.TV3.mean <- readRDS("data_derived/REM_results/1_TV3_mean.rds")
rem.1.TV3.rms <- readRDS("data_derived/REM_results/1_TV3_rms.rds")

rem.2.mean <- readRDS("data_derived/REM_results/2_mean.rds")
rem.2.rms <- readRDS("data_derived/REM_results/2_rms.rds")

rem.3.mean <- readRDS("data_derived/REM_results/3_mean.rds")
rem.3.rms <- readRDS("data_derived/REM_results/3_rms.rds")

#_______________________________________________________________________
# 3. Clean data ----

# we'll add bias metrics, and add scenario variables

  # here we'll calculate:
  # percent bias: ((estimated - true) / true) * 100
  # absolute percent bias: abs(((estimated - true) / true) * 100)
  # CI coverage: is the true value within the 95% CI?

  # We'll separate out each level for comparison later
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

# function
add_bias_scen <- function (.rem) {
  
  rem.1 <- .rem |>
    
    mutate(
      
      # bias metrics
      perc.bias = ((mean.REM.D - true.D) / true.D * 100),
      abs.perc.bias = abs((mean.REM.D - true.D) / true.D * 100),
      ci.cov = ifelse(true.D >= l95.REM.D & true.D <= u95.REM.D,
                      1,
                      0),
      
      # scenarios
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
                        levels = c("4", "2"))
      
      )
  
  return(rem.1)
  
}

# use function
rem.1.TV1.mean.1 <- add_bias_scen(rem.1.TV1.mean)
rem.1.TV1.rms.1 <- add_bias_scen(rem.1.TV1.rms)

rem.1.TV2.mean.1 <- add_bias_scen(rem.1.TV2.mean)
rem.1.TV2.rms.1 <- add_bias_scen(rem.1.TV2.rms)

rem.1.TV3.mean.1 <- add_bias_scen(rem.1.TV3.mean)
rem.1.TV3.rms.1 <- add_bias_scen(rem.1.TV3.rms)

rem.2.mean.1 <- add_bias_scen(rem.2.mean)
rem.2.rms.1 <- add_bias_scen(rem.2.rms)

rem.3.mean.1 <- add_bias_scen(rem.3.mean)
rem.3.rms.1 <- add_bias_scen(rem.3.rms)

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
                                duration == focal.combo$duration) |>
        
        # drop NAs
        drop_na()
      
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
      l95.metric <- quantile(matrix.j[ , 1], probs = 0.025, na.rm = T)
      u95.metric <- quantile(matrix.j[ , 1], probs = 0.975, na.rm = T)
      
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
                                duration == focal.combo$duration) |>
        
        # drop NAs
        drop_na()
      
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
      l95.metric <- quantile(matrix.j[ , 1], probs = 0.025, na.rm = T)
      u95.metric <- quantile(matrix.j[ , 1], probs = 0.975, na.rm = T)
      
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
                                duration == focal.combo$duration) |>
        
        # drop NAs
        drop_na()
      
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
      l95.metric <- quantile(matrix.j[ , 1], probs = 0.025, na.rm = T)
      u95.metric <- quantile(matrix.j[ , 1], probs = 0.975, na.rm = T)
      
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

rem.1.TV1.mean.bias <- boot_perform(rem.1.TV1.mean.1, "bias")
rem.1.TV1.rms.bias <- boot_perform(rem.1.TV1.rms.1, "bias")

rem.1.TV2.mean.bias <- boot_perform(rem.1.TV2.mean.1, "bias")
rem.1.TV2.rms.bias <- boot_perform(rem.1.TV2.rms.1, "bias")

rem.1.TV3.mean.bias <- boot_perform(rem.1.TV3.mean.1, "bias")
rem.1.TV3.rms.bias <- boot_perform(rem.1.TV3.rms.1, "bias")

# bind together after adding identifiers
rem.1.mean.bias <- bind_rows(
  
  rem.1.TV1.mean.bias |> mutate(TV = "0.5", speed = "CVM"),
  rem.1.TV2.mean.bias |> mutate(TV = "30", speed = "CVM"),
  rem.1.TV3.mean.bias |> mutate(TV = "120", speed = "CVM")
  
) |>
  
  # factor levels
  mutate(TV = factor(TV, 
                     levels = c("0.5", "30", "120"),
                     labels = c("0.5 mins", "30 mins", "120 mins")))

rem.1.rms.bias <- bind_rows(
  
  rem.1.TV1.rms.bias |> mutate(TV = "0.5", speed = "RMS"),
  rem.1.TV2.rms.bias |> mutate(TV = "30", speed = "RMS"),
  rem.1.TV3.rms.bias |> mutate(TV = "120", speed = "RMS")
  
) |>
  
  # factor levels
  mutate(TV = factor(TV, 
                     levels = c("0.5", "30", "120"),
                     labels = c("0.5 mins", "30 mins", "120 mins")))

rem.2.mean.bias <- boot_perform(rem.2.mean.1, "bias")
rem.2.rms.bias <- boot_perform(rem.2.rms.1, "bias")

rem.3.mean.bias <- boot_perform(rem.3.mean.1, "bias")
rem.3.rms.bias <- boot_perform(rem.3.rms.1, "bias")

#_______________________________________________________________________
# 5b. Precision (coefficient of variation) ----
#_______________________________________________________________________

rem.1.TV1.mean.precis <- boot_perform(rem.1.TV1.mean.1, "precision")
rem.1.TV1.rms.precis <- boot_perform(rem.1.TV1.rms.1, "precision")

rem.1.TV2.mean.precis <- boot_perform(rem.1.TV2.mean.1, "precision")
rem.1.TV2.rms.precis <- boot_perform(rem.1.TV2.rms.1, "precision")

rem.1.TV3.mean.precis <- boot_perform(rem.1.TV3.mean.1, "precision")
rem.1.TV3.rms.precis <- boot_perform(rem.1.TV3.rms.1, "precision")

# bind together after adding identifiers
rem.1.mean.precis <- bind_rows(
  
  rem.1.TV1.mean.precis |> mutate(TV = "0.5", speed = "CVM"),
  rem.1.TV2.mean.precis |> mutate(TV = "30", speed = "CVM"),
  rem.1.TV3.mean.precis |> mutate(TV = "120", speed = "CVM")
  
) |>
  
  # factor levels
  mutate(TV = factor(TV, 
                     levels = c("0.5", "30", "120"),
                     labels = c("0.5 mins", "30 mins", "120 mins")))

rem.1.rms.precis <- bind_rows(
  
  rem.1.TV1.rms.precis |> mutate(TV = "0.5", speed = "RMS"),
  rem.1.TV2.rms.precis |> mutate(TV = "30", speed = "RMS"),
  rem.1.TV3.rms.precis |> mutate(TV = "120", speed = "RMS")
  
) |>
  
  # factor levels
  mutate(TV = factor(TV, 
                     levels = c("0.5", "30", "120"),
                     labels = c("0.5 mins", "30 mins", "120 mins")))

rem.2.mean.precis <- boot_perform(rem.2.mean.1, "precision")
rem.2.rms.precis <- boot_perform(rem.2.rms.1, "precision")

rem.3.mean.precis <- boot_perform(rem.3.mean.1, "precision")
rem.3.rms.precis <- boot_perform(rem.3.rms.1, "precision")

#_______________________________________________________________________
# 5c. Coverage (95% confidence interval coverage) ----
#_______________________________________________________________________

rem.1.TV1.mean.cov <- boot_perform(rem.1.TV1.mean.1, "coverage")
rem.1.TV1.rms.cov <- boot_perform(rem.1.TV1.rms.1, "coverage")

rem.1.TV2.mean.cov <- boot_perform(rem.1.TV2.mean.1, "coverage")
rem.1.TV2.rms.cov <- boot_perform(rem.1.TV2.rms.1, "coverage")

rem.1.TV3.mean.cov <- boot_perform(rem.1.TV3.mean.1, "coverage")
rem.1.TV3.rms.cov <- boot_perform(rem.1.TV3.rms.1, "coverage")

# bind together after adding identifiers
rem.1.mean.cov <- bind_rows(
  
  rem.1.TV1.mean.cov |> mutate(TV = "0.5", speed = "CVM"),
  rem.1.TV2.mean.cov |> mutate(TV = "30", speed = "CVM"),
  rem.1.TV3.mean.cov |> mutate(TV = "120", speed = "CVM")
  
) |>
  
  # factor levels
  mutate(TV = factor(TV, 
                     levels = c("0.5", "30", "120"),
                     labels = c("0.5 mins", "30 mins", "120 mins")))

rem.1.rms.cov <- bind_rows(
  
  rem.1.TV1.rms.cov |> mutate(TV = "0.5", speed = "RMS"),
  rem.1.TV2.rms.cov |> mutate(TV = "30", speed = "RMS"),
  rem.1.TV3.rms.cov |> mutate(TV = "120", speed = "RMS")
  
) |>
  
  # factor levels
  mutate(TV = factor(TV, 
                     levels = c("0.5", "30", "120"),
                     labels = c("0.5 mins", "30 mins", "120 mins")))

rem.2.mean.cov <- boot_perform(rem.2.mean.1, "coverage")
rem.2.rms.cov <- boot_perform(rem.2.rms.1, "coverage")

rem.3.mean.cov <- boot_perform(rem.3.mean.1, "coverage")
rem.3.rms.cov <- boot_perform(rem.3.rms.1, "coverage")

#_______________________________________________________________________
# 6. Visualize performance ----

# I originally included all the data, but it's not possible to see the CIs

#_______________________________________________________________________
# 6a. Bias ----
#_______________________________________________________________________

# Q1
rem.1.bias <- rbind(rem.1.mean.bias, rem.1.rms.bias)

ggplot(data = rem.1.bias) +
  
  theme_bw() +
  
  # each density as columns, duration as rows
  # the nested function is neat!
  ggh4x::facet_nested(TV + duration + fix.success ~ true.D,
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
  geom_errorbarh(aes(y = interaction(fix.rate, fix.success),
                     xmin = l95.metric,
                     xmax = u95.metric),
                 color = "palegreen3",
                 alpha = 0.5,
                 linewidth = 1.25,
                 height = 0) +
  
  # means
  geom_point(aes(x = mean.metric,
                 y = interaction(fix.rate, fix.success),
                 shape = speed,
                 size = speed,
                 fill = speed)) +
  
  # theme arguments
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_text(hjust = 0.01,
                                    face = "bold",
                                    size = 8),
        strip.text.y = element_text(hjust = 0.01,
                                    size = 7),
        strip.background = element_rect(fill = "white"),
        axis.text.y = element_text(size = 7,
                                   color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "pt")) +
  
  # shapes, sizes, and colors
  scale_shape_manual(values = c(21, 22)) +
  scale_size_manual(values = c(1.35, 1.25)) +
  scale_fill_manual(values = c("palegreen4", "palegreen1")) +
  
  # reverse the y axis
  scale_y_discrete(limits = rev,
                   labels = c(rep(c("4", "1", "0.5"), 
                                  times = 4))) +
  
  coord_cartesian(xlim = c(0, 275)) +
  scale_x_continuous(breaks = c(50, 150, 250)) -> bias.Q1.plot

# Q2
rem.2.bias <- rbind(rem.2.mean.bias |> mutate(speed = "CVM"),
                    rem.2.rms.bias |> mutate(speed = "RMS"))


ggplot(data = rem.2.bias) +
  
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
  geom_errorbarh(aes(y = interaction(fix.rate, fix.success),
                     xmin = l95.metric,
                     xmax = u95.metric),
                 color = "salmon3",
                 alpha = 0.5,
                 linewidth = 1.25,
                 height = 0) +
  
  # means
  geom_point(aes(x = mean.metric,
                 y = interaction(fix.rate, fix.success),
                 shape = speed,
                 size = speed,
                 fill = speed)) +
  
  # theme arguments
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_blank(),
        strip.text.y = element_text(hjust = 0.01,
                                    size = 7),
        strip.background.x = element_blank(),
        strip.background.y = element_rect(fill = "white"),
        axis.text.y = element_text(size = 7,
                                   color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0, 15.5, 0, 0), "pt")) +
  
  # shapes, sizes, and colors
  scale_shape_manual(values = c(21, 22)) +
  scale_size_manual(values = c(1.35, 1.25)) +
  scale_fill_manual(values = c("salmon4", "salmon1")) +
  
  # reverse the y axis
  scale_y_discrete(limits = rev,
                   labels = c(rep(c("4", "1", "0.5"), 
                                  times = 4))) +
  
  coord_cartesian(xlim = c(0, 275)) +
  scale_x_continuous(breaks = c(50, 150, 250)) -> bias.Q2.plot

# Q3
rem.3.bias <- rbind(rem.3.mean.bias |> mutate(speed = "CVM"),
                    rem.3.rms.bias |> mutate(speed = "RMS"))

ggplot(data = rem.3.bias) +
  
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
  geom_errorbarh(aes(y = interaction(fix.rate, fix.success),
                     xmin = l95.metric,
                     xmax = u95.metric),
                 color = "gold3",
                 alpha = 0.5,
                 linewidth = 1.25,
                 height = 0) +
  
  # means
  geom_point(aes(x = mean.metric,
                 y = interaction(fix.rate, fix.success),
                 shape = speed,
                 size = speed,
                 fill = speed)) +
  
  # theme arguments
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_blank(),
        strip.text.y = element_text(hjust = 0.01,
                                    size = 7),
        strip.background.x = element_blank(),
        strip.background.y = element_rect(fill = "white"),
        axis.text.y = element_text(size = 7,
                                   color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0, 15.5, 0, 0), "pt")) +
  
  # shapes, sizes, and colors
  scale_shape_manual(values = c(21, 22)) +
  scale_size_manual(values = c(1.35, 1.25)) +
  scale_fill_manual(values = c("gold4", "gold1")) +
  
  # reverse the y axis
  scale_y_discrete(limits = rev,
                   labels = c(rep(c("4", "1", "0.5"), 
                                  times = 4))) +
  
  coord_cartesian(xlim = c(0, 275)) +
  scale_x_continuous(breaks = c(50, 150, 250)) +
  
  xlab("|Percent bias|") -> bias.Q3.plot

# plot together
plot_grid(bias.Q1.plot, bias.Q2.plot, bias.Q3.plot,
          nrow = 3,
          rel_heights = c(3, 1, 1.3))

# 600 x 700

#_______________________________________________________________________
# 6b. Precision ----

# we originally had these as bar charts, but that won't work with Q1

#_______________________________________________________________________

# Q1
rem.1.precis <- rbind(rem.1.mean.precis, rem.1.rms.precis)

ggplot(data = rem.1.precis) +
  
  theme_bw() +
  
  # each density as columns, duration as rows
  # the nested function is neat!
  ggh4x::facet_nested(TV + duration + fix.success ~ true.D,
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
  geom_errorbarh(aes(y = interaction(fix.rate, fix.success),
                     xmin = l95.metric,
                     xmax = u95.metric),
                 color = "palegreen3",
                 alpha = 0.5,
                 linewidth = 1.25,
                 height = 0) +
  
  # means
  geom_point(aes(x = mean.metric,
                 y = interaction(fix.rate, fix.success),
                 shape = speed,
                 size = speed,
                 fill = speed)) +
  
  # theme arguments
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_text(hjust = 0.01,
                                    face = "bold",
                                    size = 8),
        strip.text.y = element_text(hjust = 0.01,
                                    size = 7),
        strip.background = element_rect(fill = "white"),
        axis.text.y = element_text(size = 7,
                                   color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "pt")) +
  
  # shapes, sizes, and colors
  scale_shape_manual(values = c(21, 22)) +
  scale_size_manual(values = c(1.35, 1.25)) +
  scale_fill_manual(values = c("palegreen4", "palegreen1")) +
  
  # reverse the y axis
  scale_y_discrete(limits = rev,
                   labels = c(rep(c("4", "1", "0.5"), 
                                  times = 4))) +
  
  coord_cartesian(xlim = c(0.1, 0.5)) +
  scale_x_continuous(breaks = c(0.2, 0.3, 0.4)) -> precis.Q1.plot

# Q2
rem.2.precis <- rbind(rem.2.mean.precis |> mutate(speed = "CVM"),
                      rem.2.rms.precis |> mutate(speed = "RMS"))

ggplot(data = rem.2.precis) +
  
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
  geom_errorbarh(aes(y = interaction(fix.rate, fix.success),
                     xmin = l95.metric,
                     xmax = u95.metric),
                 color = "salmon3",
                 alpha = 0.5,
                 linewidth = 1.25,
                 height = 0) +
  
  # means
  geom_point(aes(x = mean.metric,
                 y = interaction(fix.rate, fix.success),
                 shape = speed,
                 size = speed,
                 fill = speed)) +
  
  # theme arguments
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_blank(),
        strip.text.y = element_text(hjust = 0.01,
                                    size = 7),
        strip.background.x = element_blank(),
        strip.background.y = element_rect(fill = "white"),
        axis.text.y = element_text(size = 7,
                                   color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0, 15.5, 0, 0), "pt")) +
  
  # shapes, sizes, and colors
  scale_shape_manual(values = c(21, 22)) +
  scale_size_manual(values = c(1.35, 1.25)) +
  scale_fill_manual(values = c("salmon4", "salmon1")) +
  
  # reverse the y axis
  scale_y_discrete(limits = rev,
                   labels = c(rep(c("4", "1", "0.5"), 
                                  times = 4))) +
  
  coord_cartesian(xlim = c(0.1, 0.5)) +
  scale_x_continuous(breaks = c(0.2, 0.3, 0.4)) -> precis.Q2.plot

# Q3
rem.3.precis <- rbind(rem.3.mean.precis |> mutate(speed = "CVM"),
                      rem.3.rms.precis |> mutate(speed = "RMS"))

ggplot(data = rem.3.precis) +
  
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
  geom_errorbarh(aes(y = interaction(fix.rate, fix.success),
                     xmin = l95.metric,
                     xmax = u95.metric),
                 color = "gold3",
                 alpha = 0.5,
                 linewidth = 1.25,
                 height = 0) +
  
  # means
  geom_point(aes(x = mean.metric,
                 y = interaction(fix.rate, fix.success),
                 shape = speed,
                 size = speed,
                 fill = speed)) +
  
  # theme arguments
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_blank(),
        strip.text.y = element_text(hjust = 0.01,
                                    size = 7),
        strip.background.x = element_blank(),
        strip.background.y = element_rect(fill = "white"),
        axis.text.y = element_text(size = 7,
                                   color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0, 15.5, 0, 0), "pt")) +
  
  # shapes, sizes, and colors
  scale_shape_manual(values = c(21, 22)) +
  scale_size_manual(values = c(1.35, 1.25)) +
  scale_fill_manual(values = c("gold4", "gold1")) +
  
  # reverse the y axis
  scale_y_discrete(limits = rev,
                   labels = c(rep(c("4", "1", "0.5"), 
                                  times = 4))) +
  
  coord_cartesian(xlim = c(0.1, 0.5)) +
  scale_x_continuous(breaks = c(0.2, 0.3, 0.4)) +
  
  xlab("Coefficient of variation") -> precis.Q3.plot

# plot together
plot_grid(precis.Q1.plot, precis.Q2.plot, precis.Q3.plot,
          nrow = 3,
          rel_heights = c(3, 1, 1.3))

# 600 x 700

#_______________________________________________________________________
# 6c. Coverage ----

# this is the percent of replicates whose 95% CIs include the true value

#_______________________________________________________________________

# Q1
rem.1.cov <- rbind(rem.1.mean.cov, rem.1.rms.cov)

ggplot(data = rem.1.cov) +
  
  theme_bw() +
  
  # each density as columns, duration as rows
  # the nested function is neat!
  ggh4x::facet_nested(TV + duration + fix.success ~ true.D,
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
  geom_errorbarh(aes(y = interaction(fix.rate, fix.success),
                     xmin = l95.metric,
                     xmax = u95.metric),
                 color = "palegreen3",
                 alpha = 0.5,
                 linewidth = 1.25,
                 height = 0) +
  
  # means
  geom_point(aes(x = mean.metric,
                 y = interaction(fix.rate, fix.success),
                 shape = speed,
                 size = speed,
                 fill = speed)) +
  
  # theme arguments
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_text(hjust = 0.01,
                                    face = "bold",
                                    size = 8),
        strip.text.y = element_text(hjust = 0.01,
                                    size = 7),
        strip.background = element_rect(fill = "white"),
        axis.text.y = element_text(size = 7,
                                   color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "pt")) +
  
  # shapes, sizes, and colors
  scale_shape_manual(values = c(21, 22)) +
  scale_size_manual(values = c(1.35, 1.25)) +
  scale_fill_manual(values = c("palegreen4", "palegreen1")) +
  
  # reverse the y axis
  scale_y_discrete(limits = rev,
                   labels = c(rep(c("4", "1", "0.5"), 
                                  times = 4))) +
  
  coord_cartesian(xlim = c(0, 1.0)) +
  scale_x_continuous(breaks = c(0.1, 0.5, 0.9)) -> cov.Q1.plot

# Q2
rem.2.cov <- rbind(rem.2.mean.cov |> mutate(speed = "CVM"),
                   rem.2.rms.cov |> mutate(speed = "RMS"))

ggplot(data = rem.2.cov) +
  
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
  geom_errorbarh(aes(y = interaction(fix.rate, fix.success),
                     xmin = l95.metric,
                     xmax = u95.metric),
                 color = "salmon3",
                 alpha = 0.5,
                 linewidth = 1.25,
                 height = 0) +
  
  # means
  geom_point(aes(x = mean.metric,
                 y = interaction(fix.rate, fix.success),
                 shape = speed,
                 size = speed,
                 fill = speed)) +
  
  # theme arguments
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_blank(),
        strip.text.y = element_text(hjust = 0.01,
                                    size = 7),
        strip.background.x = element_blank(),
        strip.background.y = element_rect(fill = "white"),
        axis.text.y = element_text(size = 7,
                                   color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0, 15.5, 0, 0), "pt")) +
  
  # shapes, sizes, and colors
  scale_shape_manual(values = c(21, 22)) +
  scale_size_manual(values = c(1.35, 1.25)) +
  scale_fill_manual(values = c("salmon4", "salmon1")) +
  
  # reverse the y axis
  scale_y_discrete(limits = rev,
                   labels = c(rep(c("4", "1", "0.5"), 
                                  times = 4))) +
  
  coord_cartesian(xlim = c(0, 1.0)) +
  scale_x_continuous(breaks = c(0.1, 0.5, 0.9)) -> cov.Q2.plot

# Q3
rem.3.cov <- rbind(rem.3.mean.cov |> mutate(speed = "CVM"),
                   rem.3.rms.cov |> mutate(speed = "RMS"))

ggplot(data = rem.3.cov) +
  
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
  geom_errorbarh(aes(y = interaction(fix.rate, fix.success),
                     xmin = l95.metric,
                     xmax = u95.metric),
                 color = "gold3",
                 alpha = 0.5,
                 linewidth = 1.25,
                 height = 0) +
  
  # means
  geom_point(aes(x = mean.metric,
                 y = interaction(fix.rate, fix.success),
                 shape = speed,
                 size = speed,
                 fill = speed)) +
  
  # theme arguments
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_blank(),
        strip.text.y = element_text(hjust = 0.01,
                                    size = 7),
        strip.background.x = element_blank(),
        strip.background.y = element_rect(fill = "white"),
        axis.text.y = element_text(size = 7,
                                   color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0, 15.5, 0, 0), "pt")) +
  
  # shapes, sizes, and colors
  scale_shape_manual(values = c(21, 22)) +
  scale_size_manual(values = c(1.35, 1.25)) +
  scale_fill_manual(values = c("gold4", "gold1")) +
  
  # reverse the y axis
  scale_y_discrete(limits = rev,
                   labels = c(rep(c("4", "1", "0.5"), 
                                  times = 4))) +
  
  coord_cartesian(xlim = c(0, 1.0)) +
  scale_x_continuous(breaks = c(0.1, 0.5, 0.9)) +
  
  xlab("Coverage") -> cov.Q3.plot

# plot together
plot_grid(cov.Q1.plot, cov.Q2.plot, cov.Q3.plot,
          nrow = 3,
          rel_heights = c(3, 1, 1.3))

# 600 x 700

#_______________________________________________________________________
# 7. Write tables ----
#_______________________________________________________________________

write.table(rem.1.bias, "clipboard", sep = "\t")
write.table(rem.2.bias, "clipboard", sep = "\t")
write.table(rem.3.bias, "clipboard", sep = "\t")

write.table(rem.1.precis, "clipboard", sep = "\t")
write.table(rem.2.precis, "clipboard", sep = "\t")
write.table(rem.3.precis, "clipboard", sep = "\t")

write.table(rem.1.cov, "clipboard", sep = "\t")
write.table(rem.2.cov, "clipboard", sep = "\t")
write.table(rem.3.cov, "clipboard", sep = "\t")
