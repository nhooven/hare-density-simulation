# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 11 - Figures for workflow diagram
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 06 May 2025
# Date completed: 06 May 2025
# Date last modified: 08 Jun 2025
# R version: 4.4.3

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation

#_______________________________________________________________________
# 2. Model parameter distributions ----
#_______________________________________________________________________

# "typical" individual
typical.df <- data.frame(x = seq(0, 40, 0.01),
                         y = dlnorm(seq(0, 40, 0.01), 
                                    meanlog = 2.25, 
                                    sdlog = 0.75))

ggplot(data = typical.df,
       aes(x = x,
           y = y)) +
  
  theme_bw() +
  
  geom_line(linewidth = 3,
            color = "purple3") +
  
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

# full population
full.df <- data.frame()

# loop
for (i in 1:9) {
  
  focal.df <- data.frame(x = seq(0, 40, 0.01),
                         y = dlnorm(seq(0, 40, 0.01), 
                                    meanlog = rnorm(1, 2.25, 0.65),
                                    sdlog = rnorm(1, 0.75, 0.5)),
                         i = i)
  
  full.df <- rbind(full.df, focal.df)
  
}

ggplot(data = full.df,
       aes(x = x,
           y = y,
           color = as.factor(i))) +
  
  theme_bw() +
  
  geom_line(linewidth = 3) +
  
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") +
  
  scale_color_brewer(palette = "YlOrRd")

#_______________________________________________________________________
# 3. Evaluation metrics ----

library(cowplot)
library(shades)

#_______________________________________________________________________

# bias
ggplot(data = data.frame(y = 1:3,
                         x = c(0.1, 0.5, -0.6)),
       aes(x = x,
           y = as.factor(y))) +
  
  theme_bw() +
  
  geom_vline(linetype = "dashed",
             xintercept = 0,
             color = "darkgray") +
  
  geom_point(aes(fill = as.factor(y)),
             shape = 21,
             size = 4,
             stroke = 1.4) +
  
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = "black",
                                   size = 14)) +
  
  coord_cartesian(xlim = c(-0.7, 0.7)) +
  
  scale_x_continuous(breaks = c(-0.5, 
                                0, 
                                0.5),
                     labels = c("underestimated", 
                                expression(paste("True ", italic("D"))), 
                                "overestimated")) +
  
  scale_fill_manual(values = saturation("purple3",
                                        c(0.4, 
                                          0.6, 
                                          1.0))) -> bias.plot

# precision
ggplot(data = data.frame(y = 1:3,
                         x = c(0.2, 0.8, 0.3)),
       aes(x = x,
           y = as.factor(y))) +
  
  theme_bw() +
  
  geom_col(aes(fill = as.factor(y)),
           color = "black",
           width = 0.5) +
  
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = "black",
                                   size = 14)) +
  
  coord_cartesian(xlim = c(0.1, 0.9)) +
  
  scale_x_continuous(breaks = c(0.2, 0.8),
                     labels = c("less uncertainty", "more uncertainty")) +
  
  scale_fill_manual(values = saturation("purple3",
                                        c(0.4, 
                                          0.6, 
                                          1.0))) -> precis.plot

# coverage
ggplot(data = data.frame(y = 1:3,
                         x = c(0.1, 0.5, -0.6),
                         xmin = c(-0.1, -0.2, -0.9),
                         xmax = c(0.3, 1.1, -0.3)),
       aes(x = x,
           y = as.factor(y))) +
  
  theme_bw() +
  
  geom_vline(linetype = "dashed",
             xintercept = 0,
             color = "darkgray") +
  
  geom_errorbarh(aes(y = as.factor(y),
                     xmin = xmin,
                     xmax = xmax,
                     color = as.factor(y),
                     linetype = as.factor(y),
                     linewidth = as.factor(y),
                     alpha = as.factor(y)),
                 height = 0) +
  
  geom_point(aes(fill = as.factor(y),
                 alpha = as.factor(y)),
             color = "black",
             shape = 21,
             size = 4,
             stroke = 1.4) +
  
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = "black",
                                   size = 14)) +
  
  coord_cartesian(xlim = c(-1.1, 1.1)) +
  
  scale_x_continuous(breaks = 0,
                     labels = expression(paste("True ", italic("D")))) +
  
  scale_linetype_manual(values = c("solid", 
                                   "solid", 
                                   "dashed")) +
  
  scale_linewidth_manual(values = c(4, 
                                    4, 
                                    1)) +
  
  scale_color_manual(values = saturation("purple3",
                                         c(0.4, 
                                           0.6, 
                                           1.0))) +
  
  scale_fill_manual(values = saturation("purple3",
                                        c(0.4, 
                                          0.6, 
                                          1.0))) +
  
  scale_alpha_manual(values = c(1.0, 1.0, 0.25)) -> cov.plot

# plot together
plot_grid(bias.plot, precis.plot, cov.plot,
          nrow = 3)

# 441 x 457

#_______________________________________________________________________
# 4. Track resampling ----

library(adehabitatLT)

#_______________________________________________________________________

# simulate from a CRW
base.CRW <- simm.crw(date = 1:40,
                     proj4string = CRS("epsg:32611"))

plot(base.CRW)

# dfs
base.crw.df <- base.CRW[[1]]

rate.crw.df <- base.crw.df[seq(1, 40, 2), ]    # fix rate
succ.crw.df <- base.crw.df[sample(1:40, size = 40 * 0.6), ]    # fix success
dur.crw.df <- base.crw.df[1:20, ]

# plots
# rate
ggplot() +
  
  theme_minimal() +
  
  geom_path(data = base.crw.df,
            aes(x = x,
                y = y),
            color = "darkgray",
            linewidth = 1.0) + 
  
  geom_path(data = rate.crw.df,
            aes(x = x,
                y = y),
            color = "purple3",
            linewidth = 1.2) + 
  
  geom_point(data = rate.crw.df,
             aes(x = x,
                 y = y),
             shape = 21,
             fill = "purple3",
             size = 3) +
  
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title.y = element_text(color = "black",
                                    size = 24),
        axis.title.x = element_blank()) +
  
  ylab("Fix rate") -> plot.rate

# succ
ggplot() +
  
  theme_minimal() +
  
  geom_point(data = base.crw.df,
             aes(x = x,
                 y = y),
             color = "darkgray",
             size = 2) +
  
  geom_point(data = succ.crw.df,
             aes(x = x,
                 y = y),
             shape = 21,
             fill = "purple3",
             size = 3) +
  
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title.y = element_text(color = "black",
                                    size = 24),
        axis.title.x = element_blank()) +
  
  ylab("Fix success") -> plot.succ

# duration
ggplot() +
  
  theme_minimal() +
  
  geom_path(data = base.crw.df,
            aes(x = x,
                y = y),
            color = "darkgray",
            linewidth = 1.0) + 
  
  geom_path(data = dur.crw.df,
            aes(x = x,
                y = y),
            color = "purple3",
            linewidth = 1.1) + 
  
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title.y = element_text(color = "black",
                                    size = 24),
        axis.title.x = element_blank()) +
  
  ylab("Track duration") -> plot.dur

# plot together
plot_grid(plot.rate, plot.succ, plot.dur,
          ncol = 3)

#_______________________________________________________________________
# 5. Reference vs. estimated ----
#_______________________________________________________________________

ref.est.df <- data.frame(est = c(0.6, 1.0, 1.24, 2.1, 2.5, 2.34, 3.6, 4.5),
                         ref = c(0.3, 1.0, 1.93, 1.42, 3.1, 1.87, 4.0, 2.19))

ggplot(data = ref.est.df,
       aes(x = est,
           y = ref)) +
  
  theme_bw() +
  
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dashed") +
  
  geom_point(aes(fill = as.factor(ref)),
             size = 3,
             shape = 21) +
  
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_blank()) +
  
  coord_cartesian(xlim = c(0, 5),
                  ylim = c(0, 5)) +
  
  labs(x = expression(paste("Estimated ", italic("D"))),
       y = expression(paste("Reference ", italic("D")))) +
  
  scale_fill_brewer(palette = "PuBu")
