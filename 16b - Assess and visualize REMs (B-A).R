# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 16b - Assess and visualize REMs (B-A)
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 13 Dec 2024
# Date completed: 
# Date last modified: 19 Feb 2025
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation
library(ggridges)        # ridgeline plots

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

all.rem.boot <- read.csv(paste0(getwd(), "/Derived_data/REM results/all_remBA_boot.csv"))

#_______________________________________________________________________
# 3. Plot ----
#_______________________________________________________________________
# 3a. Reorder and rename factors ----
#_______________________________________________________________________

all.rem.boot.1 <- all.rem.boot %>%
  
  mutate(rep = factor(rep,
                      levels = c(1, 2, 3),
                      labels = c("Landscape 1", "Landscape 2", "Landscape 3")),
         method = factor(method,
                         labels = c("M1",
                                    "M2",
                                    "M3",
                                    "M4",
                                    "M5")))

#_______________________________________________________________________
# 3a. Direct comparison plot ----

# this is the original plot design I had -
# these show 1:1 lines, with estimated density on the x and true density on the y

#_______________________________________________________________________

# before
ggplot(data = all.rem.boot.1) +
  
  theme_bw() +
  
  facet_grid(n.cams ~ method) +
  
  coord_cartesian(xlim = c(0, 20),
                  ylim = c(0, 20)) +
  
  # 1:1 lines
  geom_abline(slope = 1, 
              intercept = 0,
              linetype = "dashed") +
  
  geom_errorbarh(aes(y = n.indiv / 10,
                     xmin = l95,
                     xmax = u95),
                 height = 0,
                 linewidth = 1.5,
                 alpha = 0.15) +
  
  geom_point(aes(x = mean,
                 y = n.indiv / 10,
                 fill = as.factor(n.indiv),
                 shape = rep),
             size = 1.5) +
  
  scale_shape_manual(values = c(21, 22, 23)) +
  
  scale_fill_viridis_d(option = "viridis") +
  scale_color_viridis_d(option = "viridis") +
  
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  
  xlab("Estimated density (individuals/ha)") +
  ylab("True density (individuals/ha)")

#_______________________________________________________________________
# 4. Calculate bias ----
#_______________________________________________________________________

# calculate density and bias
all.rem.boot.2 <- all.rem.boot.1 %>%
  
  # density
  mutate(true.density = n.indiv / 10,
         true.density.fac = factor(n.indiv / 10)) %>%
  
  # bias
  mutate(bias = mean - true.density) %>%
  
  # percent bias
  mutate(percent.bias = (bias / true.density) * 100) %>%
  
  # absolute percent bias
  mutate(abs.percent.bias = abs(percent.bias)) %>%
  
  # does the 95% confidence interval contain the true value?
  mutate(conf.cov = ifelse(true.density >= l95 &
                             true.density <= u95,
                           1,
                           0))

#_______________________________________________________________________
# 5. Bias, CV, and coverage plots ----
#_______________________________________________________________________
# 5a. Absolute bias boxplots ----
#_______________________________________________________________________

ggplot(data = all.rem.boot.2) +
  
  theme_bw() +
  
  geom_boxplot(aes(x = abs(percent.bias),
                   y = method,
                   color = method),
               outliers = FALSE) +
  
  geom_jitter(aes(x = abs(percent.bias),
                  y = method,
                  color = method),
              alpha = 0.5) +
  
  facet_grid(rep ~ n.cams) +
  
  scale_color_viridis_d(end = 0.9) +
  
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank()) +
  
  # axis labels
  xlab("Absolute percent bias")

#_______________________________________________________________________
# 5b. How does bias change with density? ----
#_______________________________________________________________________

ggplot(data = all.rem.boot.2) +
  
  theme_bw() +
  
  facet_grid(n.cams ~ rep) + 
  
  geom_line(aes(x = true.density,
                y = abs.percent.bias,
                color = method,
                group = method),
            linewidth = 1.05,
            alpha = 0.65) +
  
  geom_point(aes(x = true.density,
                 y = abs.percent.bias,
                 fill = method,
                 group = method,
                 shape = method),
             size = 1.5) + 
  
  scale_color_viridis_d(end = 0.9) +
  scale_fill_viridis_d(end = 0.9) +
  
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
  
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   size = 5)) +
  
  scale_x_log10(breaks = unique(all.rem.boot.2$true.density)) +
  
  xlab("True density (individuals/ha)") +
  ylab("Absolute percent bias")

#_______________________________________________________________________
# 5c. Precision ----

# here we'll depict how the CVs change

all.rem.boot.cv.summary <- all.rem.boot.2 %>%
  
  group_by(n.cams, n.indiv, method) %>%
  
  summarize(mean.cv = mean(cv),
            sd.cv = sd(cv)) %>%
  
  mutate(true.density = factor(n.indiv / 10))

all.rem.boot.cv.summary

#_______________________________________________________________________

# grouped column plot
ggplot(data = all.rem.boot.cv.summary) +
  
  theme_bw() +
  
  facet_grid(~ n.cams) +
  
  geom_col(aes(x = true.density,
               y = mean.cv,
               fill = method,
               group = method),
           position = position_dodge(width = 0.75)) +
  
  scale_fill_viridis_d(end = 0.9) +
  
  geom_errorbar(aes(x = true.density,
                    y = mean.cv,
                    group = method,
                    ymin = mean.cv - sd.cv,
                    ymax = mean.cv + sd.cv),
                width = 0,
                color = "gray",
                position = position_dodge(width = 0.75)) +
  
  theme(panel.grid = element_blank()) +
  
  xlab("True density (individual/ha)") +
  ylab("Density coefficient of variation")


#_______________________________________________________________________
# 9c. Coverage of confidence intervals ----

# what percentage of confidence intervals overlap wih the true value?

#_______________________________________________________________________

all.rem.cov <- all.rem.boot.2 %>%
  
  group_by(method, n.cams, true.density) %>%
  
  summarize(sum.rep = sum(conf.cov))

all.rem.cov


# column plot
ggplot(data = all.rem.cov) +
  
  theme_bw() +
  
  facet_grid(method ~ n.cams) +
  
  geom_line(aes(x = true.density,
                y = sum.rep,
                color = method),
            linewidth = 1.05) +
  
  scale_color_viridis_d(end = 0.9) +
  
  theme(panel.grid = element_blank()) +
  
  scale_x_log10() +
  
  xlab("True density (individual/ha)") +
  ylab("Number of replicates with coverage")




