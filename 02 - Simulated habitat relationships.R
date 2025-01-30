# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 02b - Simulated habitat relationships
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 26 Nov 2024
# Date completed: 02 Dec 2024
# Date last modified: 30 Jan 2025
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation
library(terra)           # rasters
library(amt)             # distributions

#_______________________________________________________________________
# 2. Read in rasters ----
#_______________________________________________________________________

S1 <- rast(paste0(getwd(), "/Rasters/S1.tif"))
S2 <- rast(paste0(getwd(), "/Rasters/S2.tif"))
S3 <- rast(paste0(getwd(), "/Rasters/S3.tif"))

C1 <- rast(paste0(getwd(), "/Rasters/C1.tif"))
C2 <- rast(paste0(getwd(), "/Rasters/C2.tif"))
C3 <- rast(paste0(getwd(), "/Rasters/C3.tif"))

#_______________________________________________________________________
# 3. Movement parameter distributions ----

# here we'll just approximate hare movement from preliminary data

#_______________________________________________________________________

# step lengths (gamma)
# make distribution
sl.dist <- make_gamma_distr(shape = 1.2,
                            scale = 70)

# turning angles (uniform)
# make distribution
ta.dist <- make_unif_distr(min = -pi,
                           max = pi)

# log(sl) quantiles
lsl.med <- log(qgamma(p = 0.50, shape = sl.dist$params$shape, scale = sl.dist$params$scale))
lsl.low <- log(qgamma(p = 0.025, shape = sl.dist$params$shape, scale = sl.dist$params$scale))
lsl.hig <- log(qgamma(p = 0.975, shape = sl.dist$params$shape, scale = sl.dist$params$scale))

#_______________________________________________________________________
# 4. Base iSSF coefficients ----
#_______________________________________________________________________

# importantly, these are *standardized* coefficients

coef.forage <- 1.5          # selection for forage
coef.forage.sl <- -0.3      # shorter sl with higher forage       
coef.edge <- -0.5           # avoidance of edge distance
coef.open <- -1.5           # base avoidance of open (start of step)
coef.open.sl <- 0.5         # interaction with log(sl) (longer movements when starting in open)

#_______________________________________________________________________
# 5. Random slope standard deviations  ----
#_______________________________________________________________________

# here we'll draw iSSF slopes (habitat selection only) from normal distributions
sd.low <- 0.10
sd.high <- 0.75

#_______________________________________________________________________
# 5. Parameter expectation and variation ----
#_______________________________________________________________________
# 5a. Quantiles for plotting ----
#_______________________________________________________________________

# selection for stem
iv.forage <- data.frame(coef = "forage",
                        estimate = coef.forage,
                        q.low.l = qnorm(p = 0.025, mean = coef.forage, sd = sd.low),
                        q.low.u = qnorm(p = 0.975, mean = coef.forage, sd = sd.low),
                        q.high.l = qnorm(p = 0.025, mean = coef.forage, sd = sd.high),
                        q.high.u = qnorm(p = 0.975, mean = coef.forage, sd = sd.high))

# stem step length interaction
iv.forage.sl <- data.frame(coef = "forage:log(sl)",
                           estimate = coef.forage.sl,
                           q.low.l = NA,
                           q.low.u = NA,
                           q.high.l = NA,
                           q.high.u = NA)

# avoidance of edge
iv.edge <- data.frame(coef = "edge",
                      estimate = coef.edge,
                      q.low.l = qnorm(p = 0.025, mean = coef.edge, sd = sd.low),
                      q.low.u = qnorm(p = 0.975, mean = coef.edge, sd = sd.low),
                      q.high.l = qnorm(p = 0.025, mean = coef.edge, sd = sd.high),
                      q.high.u = qnorm(p = 0.975, mean = coef.edge, sd = sd.high))

# base avoidance of open"
iv.open <- data.frame(coef = "open",
                      estimate = coef.open,
                      q.low.l = qnorm(p = 0.025, mean = coef.open, sd = sd.low),
                      q.low.u = qnorm(p = 0.975, mean = coef.open, sd = sd.low),
                      q.high.l = qnorm(p = 0.025, mean = coef.open, sd = sd.high),
                      q.high.u = qnorm(p = 0.975, mean = coef.open, sd = sd.high))

# greater steps starting in open
iv.open.sl <- data.frame(coef = "open:log(sl)",
                         estimate = coef.open.sl,
                         q.low.l = NA,
                         q.low.u = NA,
                         q.high.l = NA,
                         q.high.u = NA)

# bind together
iv.all <- rbind(iv.forage, iv.forage.sl, iv.edge, iv.open, iv.open.sl)

# reorder factor levels
iv.all$coef <- factor(iv.all$coef,
                      levels = rev(c("forage", "forage:log(sl)", "edge", "open", "open:log(sl)")))

#_______________________________________________________________________
# 5b. Plot ----
#_______________________________________________________________________

ggplot(data = iv.all,
       aes(x = estimate,
           y = coef,
           group = coef)) +
  
  theme_bw() +
  
  geom_vline(xintercept = 0,
             linetype = "dashed") +
  
  # high variability
  geom_errorbarh(aes(xmin = q.high.l,
                     xmax = q.high.u),
                 height = 0,
                 linewidth = 2.5,
                 color = "darkblue",
                 alpha = 0.25) +
  
  # low variability
  geom_errorbarh(aes(xmin = q.low.l,
                     xmax = q.low.u),
                 height = 0,
                 linewidth = 2.5,
                 color = "darkblue",
                 alpha = 0.75) +
  
  geom_point(size = 2,
             shape = 21,
             fill = "white",
             color = "black") +
  
  # labels
  ylab("") +
  xlab("Standardized coefficient") +
  
  # theme
  theme(panel.grid = element_blank())

#_______________________________________________________________________
# 6. Density plots ----

# unusued for now - in the future it might be nice to show implied relationships
# and density plots for the step length adjustments

#_______________________________________________________________________
# 6a. Stem ----
#_______________________________________________________________________

ggplot(iv.forage,
       aes(x = x,
           y = d)) +
  
  theme_bw() +
  
  geom_vline(xintercept = 0,
             linetype = "dashed") +
  
  geom_vline(xintercept = coef.forage) +
  
  geom_line(color = "darkblue",
            linewidth = 1.5) +
  
  theme(panel.grid = element_blank()) +
  
  xlab("Stem selection coefficient") +
  ylab("")

#_______________________________________________________________________
# 6c. Mature ----
#_______________________________________________________________________

# adjustment to the sl distribution
sl.dist.mature <- update_gamma(sl.dist, beta_sl = 0, beta_log_sl = 0.5)

sl.dist.df <- rbind(data.frame(x = seq(0, 350, length.out = 1000),
                               y = dgamma(seq(0, 350, length.out = 1000), 
                                          shape = sl.dist$params$shape,
                                          scale = sl.dist$params$scale),
                               type = "other"),
                    data.frame(x = seq(0, 350, length.out = 1000),
                               y = dgamma(seq(0, 350, length.out = 1000), 
                                          shape = sl.dist.mature$params$shape,
                                          scale = sl.dist.mature$params$scale),
                               type = "mature"))

# plot
ggplot(sl.dist.df,
       aes(x = x,
           y = y,
           color = type)) +
  
  theme_bw() +
  
  geom_line(linewidth = 1.15) +
  
  theme(panel.grid = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.75, 0.75)) +
  
  xlab("Step length (m)") +
  ylab("")

#_______________________________________________________________________
# 7. Response curves ----
#_______________________________________________________________________
# 7a. Stem ----
#_______________________________________________________________________

# make dfs
stem.seq <- seq(stem.range[1], stem.range[2], length.out = 100)

stem.response.med <- data.frame(x = stem.seq,
                                y = stem.seq * coef.stem + stem.seq * lsl.med * coef.stem.sl,
                                sl = "med")

stem.response.low <- data.frame(x = stem.seq,
                                y = stem.seq * coef.stem + stem.seq * lsl.low * coef.stem.sl,
                                sl = "low")

stem.response.hig <- data.frame(x = stem.seq,
                                y = stem.seq * coef.stem + stem.seq * lsl.hig * coef.stem.sl,
                                sl = "hig")

# bind together
stem.response <- rbind(stem.response.med, stem.response.low, stem.response.hig)

# plot
ggplot(data = stem.response) +
  
  theme_bw() +
  
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  
  geom_line(aes(x = x,
                y = y,
                color = sl),
            linewidth = 1.5) +
  
  theme(panel.grid = element_blank()) +
  
  xlab("Stem (standardized)") +
  ylab("log-RSS")

#_______________________________________________________________________
# 7b. Edge ----
#_______________________________________________________________________

# make df
edge.response <- data.frame(x = seq(edge.range[1], edge.range[2], length.out = 100),
                            y = seq(edge.range[1], edge.range[2], length.out = 100) * coef.edge,
                            ylow = seq(edge.range[1], edge.range[2], length.out = 100) * 
                              qnorm(p = 0.025, mean = coef.edge, sd = 1.5),
                            yhigh = seq(edge.range[1], edge.range[2], length.out = 100) * 
                              qnorm(p = 0.975, mean = coef.edge, sd = 1.5))

# plot
ggplot(edge.response,
       aes(x = x,
           y = y)) +
  
  theme_bw() +
  
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  
  geom_ribbon(aes(x = x,
                  y = y,
                  ymin = ylow,
                  ymax = yhigh),
              color = NA,
              fill = "darkgreen",
              alpha = 0.15) +
  
  geom_line(color = "darkgreen",
            linewidth = 1.5) +
  
  theme(panel.grid = element_blank()) +
  
  xlab("Distance to edge (standardized)") +
  ylab("log-RSS")

#_______________________________________________________________________
# 7c. Mature ----
#_______________________________________________________________________

# make dfs
mature.response <- data.frame(x = c("low", "med", "hig"),
                              y = c(coef.mature + coef.mature.sl * lsl.low,
                                    coef.mature + coef.mature.sl * lsl.med,
                                    coef.mature + coef.mature.sl * lsl.hig))

# plot
ggplot(data = mature.response) +
  
  theme_bw() +
  
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  
  geom_point(aes(x = x,
                 y = y,
                color = x),
             size = 4) +
  
  theme(panel.grid = element_blank()) +
  
  xlab("Step length") +
  ylab("log-RSS")

