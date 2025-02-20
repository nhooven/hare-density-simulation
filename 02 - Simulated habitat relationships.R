# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 02b - Simulated habitat relationships
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 26 Nov 2024
# Date completed: 02 Dec 2024
# Date last modified: 20 Feb 2025
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation
library(terra)           # rasters
library(amt)             # distributions
library(circular)

#_______________________________________________________________________
# 2. Read in rasters ----
#_______________________________________________________________________

B1 <- rast(paste0(getwd(), "/Rasters/B1.tif"))
B2 <- rast(paste0(getwd(), "/Rasters/B2.tif"))
B3 <- rast(paste0(getwd(), "/Rasters/B3.tif"))

A1 <- rast(paste0(getwd(), "/Rasters/A1.tif"))
A2 <- rast(paste0(getwd(), "/Rasters/A2.tif"))
A3 <- rast(paste0(getwd(), "/Rasters/A3.tif"))

#_______________________________________________________________________
# 3. Movement parameter distributions ----
#_______________________________________________________________________
# 3a. Step lengths ----
#_______________________________________________________________________

# step lengths (gamma) - mean will be 50 m in 2 hours
# the mean of a gamma is shape * scale
# make distribution
sl.dist <- make_gamma_distr(shape = 1,
                            scale = 50)

# plot
sl.dist.df <- data.frame(x = seq(0, 70, length.out = 1000),
                         y = dgamma(x = seq(0, 200, length.out = 1000),
                                    shape = sl.dist$params$shape,
                                    scale = sl.dist$params$scale))

ggplot(data = sl.dist.df,
       aes(x = x,
           y = y)) +
  
  theme_bw() +
  
  geom_line(linewidth = 1,
            color = "blue") +
  
  theme(panel.grid = element_blank()) +
  
  xlab("Step length (m)") +
  
  ylab("")

#_______________________________________________________________________
# 3b. Turning angles ----
#_______________________________________________________________________

# turning angles (von Mises) - small concentration parameter, vague directional persistence
# make distribution
ta.dist <- make_vonmises_distr(kappa = 1.5)

# plot
ta.dist.df <- data.frame(x = circular(seq(-pi, pi, length.out = 1000)),
                         y = dvonmises(x = circular(seq(-pi, pi, length.out = 1000)),
                                       kappa = ta.dist$params$kappa,
                                       mu = circular(0)))

ggplot() +
  
  theme_bw() +
  
  geom_line(data = ta.dist.df ,
            aes(x = x,
                y = y),
            linewidth = 1,
            color = "red") +
  
  xlab("Turn angle (rad)") +
  
  ylab("") +
  
  theme(panel.grid = element_blank()) +
  
  scale_x_continuous(breaks = c(-pi, -pi/2, 0, pi/2, pi),
                     labels = c(expression(-pi, -pi/2, 0, pi/2, pi)))

#_______________________________________________________________________
# 4. Base iSSF coefficients ----
#_______________________________________________________________________

# importantly, these are *standardized* coefficients

coef.fora.sl <- -0.05       # β1 - (start) fora and speed interaction = negative
coef.fora <- 1.5            # β2 = (end) fora selection = positive
coef.fora.ta <- -0.5       # β3 = (end) fora and concentration interaction = negative       
coef.elev <- 0.75           # β4 = (end) linear elev selection = positive
coef.elev2 <- -0.75         # β5 = (end) squared elev selection = negative
coef.open.sl <- 0.25        # β6 = (start) open and speed interaction = positive
coef.open <- -1.0           # β7 = (end) open selection = negative

#_______________________________________________________________________
# 5. Parameter expectation plot ----
#_______________________________________________________________________
# 5a. Data.frame ----
#_______________________________________________________________________

# bind together
coef.all <- data.frame(coef = c("fora:log(sl)",
                                "fora",
                                "fora:cos(ta)",
                                "elev",
                                "elev2",
                                "open:log(sl)",
                                "open"),
                     estimate = c(coef.fora.sl,
                                  coef.fora,
                                  coef.fora.ta,
                                  coef.elev,
                                  coef.elev2,
                                  coef.open.sl,
                                  coef.open))

# reorder factor levels
coef.all$coef <- factor(coef.all$coef,
                        levels = rev(c("fora:log(sl)",
                                       "fora",
                                       "fora:cos(ta)",
                                       "elev",
                                       "elev2",
                                       "open:log(sl)",
                                       "open")))

# add indicator for base vs. movement interaction
coef.all$type <- c("movement", "base", "movement", "base", "base", "movement", "base")

#_______________________________________________________________________
# 5b. Plot ----
#_______________________________________________________________________

ggplot(data = coef.all,
       aes(x = estimate,
           y = coef,
           fill = type,
           shape = type)) +
  
  theme_bw() +
  
  geom_vline(xintercept = 0,
             linetype = "dashed") +
  
  geom_point(size = 3,
             color = "black") +
  
  scale_fill_manual(values = c("purple", "orange")) +
  scale_shape_manual(values = c(21, 24)) +
  
  # labels
  ylab("") +
  xlab("Standardized coefficient") +
  
  # theme
  theme(panel.grid = element_blank(),
        legend.position = "top")

#_______________________________________________________________________
# 6. Implied selection relationships ----

# these are all the calculations at the "end" of the step

# movement parameter quantiles
# step length
sl.high <- qgamma(p = 0.90, shape = sl.dist$params$shape, scale = sl.dist$params$scale)
sl.low <- qgamma(p = 0.10, shape = sl.dist$params$shape, scale = sl.dist$params$scale)

# turning angle
# this is trickier - we actually just want to see what this would be at pi
ta.max <- pi

#_______________________________________________________________________
# 6a. Forage ----
#_______________________________________________________________________

# df
fora.df <- tibble(x = seq(-2, 2, length.out = 1000)) %>%
  
  mutate(y.mean = x * coef.fora,      # at the mean ta
         y.high = x * coef.fora + x * coef.fora.ta * cos(ta.max)) %>%
  
  pivot_longer(cols = c(y.mean, y.high)) %>%
  
  mutate(ta = factor(name,
                     levels = c("y.mean",
                                "y.high"),
                     labels = c("zero",
                                "maximum")))

ggplot(fora.df,
       aes(x = x,
           y = value,
           color = ta)) +
  
  theme_bw() +
  
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  
  geom_line(linewidth = 1.5) +
  
  scale_color_manual(values = c("gray", "orange")) +
  
  theme(panel.grid = element_blank(),
        legend.position = c(0.2, 0.8)) +
  
  xlab("Forage (standardized)") +
  ylab("ln(RSS)")

#_______________________________________________________________________
# 6b. Elevation ----
#_______________________________________________________________________

# df
elev.df <- tibble(x = seq(-2, 2, length.out = 1000)) %>%
  
  mutate(y.mean = x * coef.elev + x^2 * coef.elev2)

ggplot(elev.df,
       aes(x = x,
           y = y.mean)) +
  
  theme_bw() +
  
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  
  geom_line(linewidth = 1.5,
            color = "gray") +
  
  theme(panel.grid = element_blank(),
        legend.position = c(0.2, 0.8)) +
  
  xlab("Elevation (standardized)") +
  ylab("ln(RSS)")

#_______________________________________________________________________
# 6c. Openness ----
#_______________________________________________________________________

# df
open.df <- tibble(x = seq(-2, 2, length.out = 1000)) %>%
  
  mutate(y.mean = x * coef.open)

ggplot(open.df,
       aes(x = x,
           y = y.mean)) +
  
  theme_bw() +
  
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  
  geom_line(linewidth = 1.5,
            color = "gray") +
  
  theme(panel.grid = element_blank(),
        legend.position = c(0.2, 0.8)) +
  
  xlab("Openness (standardized)") +
  ylab("ln(RSS)")

#_______________________________________________________________________
# 7. Implied movement parameter distribution adjustments ----

# https://conservancy.umn.edu/server/api/core/bitstreams/63727072-87b1-4b35-b81c-8fd31b8f1e57/content

# these are all the calculations at the "start" of the step

# habitat covariate quantiles
cov.high <- quantile(seq(-2, 2, length.out = 1000), prob = 0.90)
cov.low <- quantile(seq(-2, 2, length.out = 1000), prob = 0.10)

#_______________________________________________________________________
# 7a. Forage:log(sl) ----
#_______________________________________________________________________

# new sl distributions
sl.dist.fora.high <- update_gamma(dist = sl.dist,
                                  beta_sl = 0,
                                  beta_log_sl = cov.high * coef.fora.sl)

sl.dist.fora.low <- update_gamma(dist = sl.dist,
                                 beta_sl = 0,
                                 beta_log_sl = cov.low * coef.fora.sl)

# df
fora.sl.df <- tibble(x = seq(0.1, 70, length.out = 1000)) %>%   # normal step length
  
  mutate(y.mean = dgamma(x, 
                         shape = sl.dist$params$shape,
                         scale = sl.dist$params$scale),
         y.high = dgamma(x, 
                         shape = sl.dist.fora.high$params$shape,
                         scale = sl.dist.fora.high$params$scale),
         y.low = dgamma(x, 
                        shape = sl.dist.fora.low$params$shape,
                        scale = sl.dist.fora.low$params$scale)) %>%
  
  pivot_longer(cols = c(y.mean, y.high, y.low)) %>%
  
  mutate(fora = factor(name,
                     levels = c("y.low",
                                "y.mean",
                                "y.high"),
                     labels = c("low",
                                "mean",
                                "high")))

ggplot(fora.sl.df,
       aes(x = x,
           y = value,
           color = fora)) +
  
  theme_bw() +
  
  geom_line(linewidth = 1.5) +
  
  scale_color_manual(values = c("lightblue", "gray", "orange")) +
  
  theme(panel.grid = element_blank(),
        legend.position = c(0.75, 0.75)) +
  
  xlab("Step length (m)") +
  ylab("")

#_______________________________________________________________________
# 7b. Open:log(sl) ----
#_______________________________________________________________________

# new sl distributions
sl.dist.open.high <- update_gamma(dist = sl.dist,
                                  beta_sl = 0,
                                  beta_log_sl = cov.high * coef.open.sl)

sl.dist.open.low <- update_gamma(dist = sl.dist,
                                 beta_sl = 0,
                                 beta_log_sl = cov.low * coef.open.sl)

# df
open.sl.df <- tibble(x = seq(0.1, 70, length.out = 1000)) %>%   # normal step length
  
  mutate(y.mean = dgamma(x, 
                         shape = sl.dist$params$shape,
                         scale = sl.dist$params$scale),
         y.high = dgamma(x, 
                         shape = sl.dist.open.high$params$shape,
                         scale = sl.dist.open.high$params$scale),
         y.low = dgamma(x, 
                        shape = sl.dist.open.low$params$shape,
                        scale = sl.dist.open.low$params$scale)) %>%
  
  pivot_longer(cols = c(y.mean, y.high, y.low)) %>%
  
  mutate(open = factor(name,
                       levels = c("y.low",
                                  "y.mean",
                                  "y.high"),
                       labels = c("low",
                                  "mean",
                                  "high")))

ggplot(open.sl.df,
       aes(x = x,
           y = value,
           color = open)) +
  
  theme_bw() +
  
  geom_line(linewidth = 1.5) +
  
  scale_color_manual(values = c("lightblue", "gray", "orange")) +
  
  theme(panel.grid = element_blank(),
        legend.position = c(0.75, 0.75)) +
  
  xlab("Step length (m)") +
  ylab("")

#_______________________________________________________________________
# 8. Base RSFs for generating home range centers ----

# the approach here will be to use all the base (i.e., non-interactive coefs)
# defined above to create a naive RSF surface so we aren't placing virtual HRs
# randomly - this is akin to the movement-free habitat kernel

#_______________________________________________________________________
# 8a. Scale rasters within the buffered unit ----

# (we'll say all of this is "available")
# we only need to do this for the BEFORE landscapes since we'll use the same HRCs
# for the same individuals for the AFTER landscapes

#_______________________________________________________________________

# define extent (buffered unit)
unit.buff <- st_buffer(st_read(paste0(getwd(), "/Derived_data/Shapefiles/unit_bound.shp")), 
                       dist = 100)
# crop
B1.crop <- crop(B1, unit.buff)
B2.crop <- crop(B2, unit.buff)
B3.crop <- crop(B3, unit.buff)

# scale
B1.scale <- scale(B1.crop)
B2.scale <- scale(B2.crop)
B3.scale <- scale(B3.crop)

#_______________________________________________________________________
# 8b. Calculate "naive" RSF surfaces ----
#_______________________________________________________________________

B1.rsf <- exp(coef.fora * B1.scale$fora +
              coef.elev * B1.scale$elev +
              coef.elev2 * B1.scale$elev^2 +
              coef.open * B1.scale$open)

B2.rsf <- exp(coef.fora * B2.scale$fora +
                coef.elev * B2.scale$elev +
                coef.elev2 * B2.scale$elev^2 +
                coef.open * B2.scale$open)

B3.rsf <- exp(coef.fora * B3.scale$fora +
                coef.elev * B3.scale$elev +
                coef.elev2 * B3.scale$elev^2 +
                coef.open * B3.scale$open)

# bind together and rename
B.rsf <- c(B1.rsf, B2.rsf, B3.rsf)

names(B.rsf) <- c("B1", "B2", "B3")

# plot
plot(B.rsf)

#_______________________________________________________________________
# 8c. Write raster ----
#_______________________________________________________________________

writeRaster(B.rsf, filename = "Rasters/B_rsf.tif", overwrite = T)
