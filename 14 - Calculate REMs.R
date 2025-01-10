# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 14 - Calculate REMs
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 13 Dec 2024
# Date completed: 
# Date last modified: 09 Jan 2025
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation

#_______________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________

passes <- read.csv(paste0(getwd(), "/Derived_data/Passes/final_passes.csv"))

#_______________________________________________________________________
# 3. Calculate REM density by camera ----
#_______________________________________________________________________

# define function
calc_REM_bycam <- function(x) {
  
  # point estimates
  # full naive REM
  REM.1 <- (x$total.passes / x$days) *
           (pi / ((x$dr.c.mean * 1000) * 3.5 * (2.0 + x$lens))) *
           10000
  
  # REM corrected with naive SSF prediction
  REM.2 <- REM.1 * x$ssf.mean
  
  # REM corrected with iSSF UD prediction
  REM.3 <- REM.1 * x$issf.mean
  
  # total variance with the delta method
  var.REM.1 <- x$dr.c.sd^2
  var.REM.2 <- x$dr.c.sd^2 + x$ssf.se^2
  var.REM.3 <- x$dr.se^2 + x$issf.se^2
  
  # bind onto df
  x.1 <- cbind(x, REM.1, REM.2, REM.3,
               var.REM.1, var.REM.2, var.REM.3)
  
  # return
  return(x.1)
  
}

# apply function
passes.1 <- calc_REM_bycam(passes)

#_______________________________________________________________________
# 4. Bootstrap within each landscape-var-rep-n.cams-n.indiv combinations ----
#_______________________________________________________________________

# expand grid
all.combos <- expand.grid(landscape = c("simple", "complex"),
                          variability = c("low", "high"),
                          rep = 1:3,
                          n.cams = c(4, 9, 16),
                          n.indiv = unique(passes.1$n.indiv))

# loop through all combos
all.REM <- data.frame()

# start time
start.time <- Sys.time()

for (i in 1:nrow(all.combos)) {
  
  focal.combo <- all.combos[i, ]
  
  # subset passes.1
  focal.passes <- passes.1 %>%
    
    filter(landscape == focal.combo$landscape,
           variability == focal.combo$variability,
           rep == focal.combo$rep,
           n.cams == focal.combo$n.cams,
           n.indiv == focal.combo$n.indiv)
  
  # initialize a new data.frame to hold everything, starting with weighted means
  focal.REM <- data.frame(REM.1.mean = weighted.mean(focal.passes$REM.1, 1 / focal.passes$var.REM.1),
                          REM.2.mean = weighted.mean(focal.passes$REM.2, 1 / focal.passes$var.REM.2),
                          REM.3.mean = weighted.mean(focal.passes$REM.3, 1 / focal.passes$var.REM.3))
  
  # weighted bootstrap for each REM estimate
  boot.means <- data.frame(matrix(data = NA,
                                  nrow = 5000,
                                  ncol = 3))
  
  # loop
  for (j in 1:5000) {
    
    # naive REM
    REM.1.sample <- sample(1:nrow(focal.passes), size = nrow(focal.passes), replace = TRUE, prob = 1 / focal.passes$var.REM.1)
    
    REM.1.sampled.rows <- focal.passes[REM.1.sample, ]
    
    boot.means[j, 1] <- mean(REM.1.sampled.rows$REM.1)
    
    # REM with SSF CF
    REM.2.sample <- sample(1:nrow(focal.passes), size = nrow(focal.passes), replace = TRUE, prob = 1 / focal.passes$var.REM.2)
    
    REM.2.sampled.rows <- focal.passes[REM.2.sample, ]
    
    boot.means[j, 2] <- mean(REM.2.sampled.rows$REM.2)
    
    # REM with iSSF CF and adjusted day range
    REM.3.sample <- sample(1:nrow(focal.passes), size = nrow(focal.passes), replace = TRUE, prob = 1 / focal.passes$var.REM.3)
    
    REM.3.sampled.rows <- focal.passes[REM.3.sample, ]
    
    boot.means[j, 3] <- mean(REM.3.sampled.rows$REM.3)
    
  }
  
  # add variance measures to df
  # sampling error
  focal.REM$REM.1.se <- sd(boot.means$X1)
  focal.REM$REM.2.se <- sd(boot.means$X2)
  focal.REM$REM.3.se <- sd(boot.means$X3)
  
  # coefficients of variation
  focal.REM$REM.1.cv <- focal.REM$REM.1.se / focal.REM$REM.1.mean
  focal.REM$REM.2.cv <- focal.REM$REM.2.se / focal.REM$REM.2.mean
  focal.REM$REM.3.cv <- focal.REM$REM.3.se / focal.REM$REM.3.mean
  
  # percentile confidence intervals
  focal.REM$REM.1.l95 <- quantile(boot.means$X1, prob = 0.025)
  focal.REM$REM.2.l95 <- quantile(boot.means$X2, prob = 0.025)
  focal.REM$REM.3.l95 <- quantile(boot.means$X3, prob = 0.025)
  
  focal.REM$REM.1.u95 <- quantile(boot.means$X1, prob = 0.975)
  focal.REM$REM.2.u95 <- quantile(boot.means$X2, prob = 0.975)
  focal.REM$REM.3.u95 <- quantile(boot.means$X3, prob = 0.975)
  
  # bind in identifier information
  focal.passes.info <- focal.passes %>%
    
    dplyr::select(landscape, variability, rep, n.cams, n.indiv)
  
  focal.all <- cbind(focal.passes.info, focal.REM)
  
  # bind into master df
  all.REM <- rbind(all.REM, focal.all)
  
  # status message
  time.elapsed <- round(as.numeric(difftime(Sys.time(), start.time, units = "mins")), digits = 2)
  
  print(paste0("Completed combo ", i, " of ", nrow(all.combos), " - ", time.elapsed, " mins"))
  
}

# keep distinct rows only
all.REM.1 <- all.REM %>% dplyr::distinct()

#_______________________________________________________________________
# 5. Pivot for plotting ----
#_______________________________________________________________________

# mean
all.REM.longer.mean <- all.REM %>%
  
  # select ID and value columns
  dplyr::select(c(1:5, 6:8)) %>%
  
  # pivot longer and drop the REM. prefix
  pivot_longer(cols = 6:8,
               names_prefix = "REM.") %>%
  
  # rename "value" and add "method" columns
  rename(mean = value) %>%
  mutate(method = substring(name, 1 , 1)) %>%
  
  # drop "name"
  dplyr::select(-name)

# CV
all.REM.longer.cv <- all.REM %>%
  
  # select ID and value columns
  dplyr::select(c(1:5, 12:14)) %>%
  
  # pivot longer and drop the REM. prefix
  pivot_longer(cols = 6:8,
               names_prefix = "REM.") %>%
  
  # rename "value" and add "method" columns
  rename(cv = value) %>%
  mutate(method = substring(name, 1 , 1)) %>%
  
  # drop unused columns
  dplyr::select(-c(1:6, 8))

# l95
all.REM.longer.l95 <- all.REM %>%
  
  # select ID and value columns
  dplyr::select(c(1:5, 15:17)) %>%
  
  # pivot longer and drop the REM. prefix
  pivot_longer(cols = 6:8,
               names_prefix = "REM.") %>%
  
  # rename "value" and add "method" columns
  rename(l95 = value) %>%
  mutate(method = substring(name, 1 , 1)) %>%
  
  # drop unused columns
  dplyr::select(-c(1:6, 8))

# u95
all.REM.longer.u95 <- all.REM %>%
  
  # select ID and value columns
  dplyr::select(c(1:5, 18:20)) %>%
  
  # pivot longer and drop the REM. prefix
  pivot_longer(cols = 6:8,
               names_prefix = "REM.") %>%
  
  # rename "value" and add "method" columns
  rename(u95 = value) %>%
  mutate(method = substring(name, 1 , 1)) %>%
  
  # drop unused columns
  dplyr::select(-c(1:6, 8))

# bind together
all.REM.longer <- cbind(all.REM.longer.mean,
                        all.REM.longer.cv,
                        all.REM.longer.l95,
                        all.REM.longer.u95)

# keep distinct rows only
all.REM.longer <- all.REM.longer %>% dplyr::distinct()

#_______________________________________________________________________
# 6. Plot ----
#_______________________________________________________________________
# 6a. Reorder and rename factors ----
#_______________________________________________________________________

all.REM.longer.1 <- all.REM.longer %>%
  
  mutate(method = factor(method,
                         levels = c(1, 2, 3),
                         labels = c("naive", "HS", "HS + movement")),
         landscape = factor(landscape,
                            levels = c("simple", "complex")),
         variability = factor(variability,
                              levels = c("low", "high")),
         rep = factor(rep,
                      levels = c(1, 2, 3),
                      labels = c("Landscape 1", "Landscape 2", "Landscape 3")))

#_______________________________________________________________________
# 6b. Split out by landscape and variability ----
#_______________________________________________________________________

# make sure to only keep distinct entries
all.REM.longer.SL <- all.REM.longer.1 %>% filter(landscape == "simple" & variability == "low")
all.REM.longer.SH <- all.REM.longer.1 %>% filter(landscape == "simple" & variability == "high")
all.REM.longer.CL <- all.REM.longer.1 %>% filter(landscape == "complex" & variability == "low")
all.REM.longer.CH <- all.REM.longer.1 %>% filter(landscape == "complex" & variability == "high")

#_______________________________________________________________________
# 6c. Plots ----
#_______________________________________________________________________

# SL
ggplot(data = all.REM.longer.SL) +
  
  theme_bw() +
  
  facet_grid(method ~ n.cams) +
  
  coord_cartesian(xlim = c(0, 8),
                  ylim = c(0, 8)) +
  
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

# 540 x 540

# SH
ggplot(data = all.REM.longer.SH) +
  
  theme_bw() +
  
  facet_grid(method ~ n.cams) +
  
  coord_cartesian(xlim = c(0, 7),
                  ylim = c(0, 7)) +
  
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

# CL
ggplot(data = all.REM.longer.CL) +
  
  theme_bw() +
  
  facet_grid(method ~ n.cams) +
  
  coord_cartesian(xlim = c(0, 9),
                  ylim = c(0, 9)) +
  
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

# CH
ggplot(data = all.REM.longer.CH) +
  
  theme_bw() +
  
  facet_grid(method ~ n.cams) +
  
  coord_cartesian(xlim = c(0, 8),
                  ylim = c(0, 8)) +
  
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
# 7. Calculate bias ----
#_______________________________________________________________________

# calculate density and bias
all.REM.longer <- all.REM.longer %>%
  
  # density
  mutate(true.density = n.indiv / 10) %>%
  
  # bias
  mutate(bias = mean - true.density) %>%
  
  # percent bias
  mutate(percent.bias = (bias / true.density) * 100) %>%

  # does the 95% confidence interval contain the true value?
  mutate(conf.cov = ifelse(true.density >= l95 &
                           true.density <= u95,
                           1,
                           0))

#_______________________________________________________________________
# 8. Bias, CV, and coverage plots ----
#_______________________________________________________________________
# 8a. Factor levels ----
#_______________________________________________________________________

all.REM.longer.2 <- all.REM.longer %>%
  
  mutate(method = factor(method,
                         levels = c(1, 2, 3),
                         labels = c("naive", "HS", "HS + movement")),
         landscape = factor(landscape,
                            levels = c("simple", "complex")),
         variability = factor(variability,
                              levels = c("low", "high")),
         rep = factor(rep,
                      levels = c(1, 2, 3),
                      labels = c("Landscape 1", "Landscape 2", "Landscape 3")))

# pivot wider for direct comparison
all.REM.wider.bias <- all.REM.longer.2 %>%
  
  dplyr::select(c(1:5, 7, 12:14)) %>%
  
  pivot_wider(names_from = method,
              values_from = c(bias, percent.bias, conf.cov))

#_______________________________________________________________________
# 8b. Bias ----

# here we want to show how each approach deals with bias, stratified by
# all of the other variables

#_______________________________________________________________________

# HS vs naive
ggplot(data = all.REM.wider.bias,
       aes(x = percent.bias_HS,
           y = percent.bias_naive,
           fill = as.factor(n.indiv),
           shape = rep)) +
  
  theme_bw() +
  
  coord_cartesian(xlim = c(-150, 150),
                  ylim = c(-150, 150)) +
  
  facet_grid(landscape * variability ~ n.cams) +
  
  geom_hline(yintercept = 0) +
  
  geom_vline(xintercept = 0) +
  
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  
  geom_point(size = 0.75) +
  
  scale_shape_manual(values = c(21, 22, 23)) +
  
  scale_fill_viridis_d(option = "magma") +
  
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  
  xlab("Percent bias (HS method)") +
  
  ylab("Percent bias (naive method)")

# I don't know if I like this. We'll decide what's best later

#_______________________________________________________________________
# 8c. Precision ----

# here we'll depict how the CVs change

# we'll calculate mean CV and SD from each replicate

all.REM.longer.CV.group <- all.REM.longer.2 %>%
  
  group_by(landscape, variability, n.cams, n.indiv, method) %>%
  
  summarize(mean.CV = mean(cv),
            sd.CV = sd(cv))

#_______________________________________________________________________

# grouped column plot
ggplot(data = all.REM.longer.CV.group) +
  
  theme_bw() +
  
  facet_grid(landscape * variability ~ n.cams) +
  
  geom_col(aes(x = method,
               y = mean.CV,
               fill = as.factor(n.indiv),
               group = as.factor(n.indiv)),
           color = "black",
           position = "dodge") +
  
  geom_errorbar(aes(group = as.factor(n.indiv),
                    x = method,
                    y = mean.CV,
                    ymin = mean.CV - sd.CV,
                    ymax = mean.CV + sd.CV),
                position = position_dodge()) +
  
  scale_y_continuous(breaks = c(0, 0.5, 1.0)) +
  
  scale_fill_viridis_d(option = "magma", direction = -1) +
  
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank()) +
  
  ylab("Density coefficient of variation")

#_______________________________________________________________________
# 9. Summaries ----
#_______________________________________________________________________
# 9a. Percent bias ----
#_______________________________________________________________________

all.REM.bias <- all.REM.longer.2 %>%
  
  dplyr::select(c(1:5, 7, 11, 13)) %>%
  
  group_by(method, n.cams, true.density) %>%
  
  summarize(mean.percent.bias = mean(percent.bias),
            sd.percent.bias = sd(percent.bias),
            max.percent.bias = max(percent.bias),
            min.percent.bias = min(percent.bias))

all.REM.bias

# circle plot
ggplot(all.REM.bias,
       aes(x = true.density,
           y = abs(mean.percent.bias))) +
  
  theme_bw() +
  
  facet_grid(method ~ n.cams) +
  
  coord_cartesian(ylim = c(0, 55)) +
  
  geom_point(aes(size = sd.percent.bias,
                 fill = mean.percent.bias),
             shape = 21) +
  
  scale_fill_viridis_c(direction = -1) +
  
  theme(panel.grid = element_blank()) +
  
  xlab("True density (individuals/ha)") +
  ylab("|Mean percent bias|")

# joined circle plot
ggplot(all.REM.bias) +
  
  theme_bw() +
  
  facet_grid(method ~ n.cams) +
  
  geom_hline(yintercept = 0) +
  
  geom_segment(aes(x = true.density,
                   xend = true.density,
                   y = min.percent.bias,
                   yend = max.percent.bias),
               linewidth = 1.5,
               color = "gray") +
  
  geom_point(aes(x = true.density,
                 y = max.percent.bias,
                 fill = max.percent.bias),
             shape = 21) +
    
  geom_point(aes(x = true.density,
                 y = min.percent.bias,
                 fill = min.percent.bias),
             shape = 21) +
  
  scale_fill_viridis_c(direction = -1) +
  
  scale_x_continuous(breaks = unique(all.REM.bias$true.density)) +
  
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  
  xlab("True density (individuals/ha)") +
  ylab("Mean percent bias")

# overall method performance
all.REM.bias.1 <- all.REM.longer.2 %>%
  
  dplyr::select(c(1:5, 7, 11, 13)) %>%
  
  group_by(method, n.indiv) %>%
  
  summarize(mean.percent.bias = mean(abs(percent.bias)),
            sd.percent.bias = sd(abs(percent.bias)),
            max.percent.bias = max(abs(percent.bias)),
            min.percent.bias = min(abs(percent.bias)))

all.REM.bias.1

# plot distributions of percent bias
ggplot(data = all.REM.longer.2,
       aes(x = abs(percent.bias))) +
  
  theme_bw() +
  
  facet_wrap(variability ~ method) +
  
  geom_density()


#_______________________________________________________________________
# 9c. Coverage of confidence intervals ----

# do confidence intervals overlap with the true value?

#_______________________________________________________________________

all.REM.coverage <- all.REM.longer.2 %>%
  
  dplyr::select(c(1:5, 7, 11, 14)) %>%
  
  group_by(method, n.cams, true.density) %>%
  
  summarize(percent = sum(conf.cov) / n())

all.REM.coverage

#_______________________________________________________________________
# 10. Let's model it - see which is best! ----
#_______________________________________________________________________

# the response variable will be the abs(percent bias)

# dataset for modeling
all.REM.longer.3 <- all.REM.longer.2 %>%
  
  mutate(abs.percent.bias = abs(percent.bias)) %>%
  
  # keep covariates
  dplyr::select(landscape,
                variability,
                rep,
                n.cams,
                true.density,
                method,
                abs.percent.bias,
                cv,
                conf.cov) %>%
  
  # unique "rep" covariate
  mutate(ls = as.integer(as.factor(paste0(landscape, rep)))) %>%
  
  # drop original rep
  dplyr::select(-rep)

# use a GLMM
library(glmmTMB)

model.bias <- glmmTMB(data = all.REM.longer.3,
                 formula = abs.percent.bias ~ method +
                                              landscape +
                                              variability +
                                              n.cams +
                                              true.density +
                                              (1 | ls),
                 family = gaussian())

summary(model.bias)

# looks like the HS + movement scheme doesn't help us at all with bias

# how about precision?
model.cv <- glmmTMB(data = all.REM.longer.3,
                      formula = cv ~ method +
                                     landscape +
                                     variability +
                                     n.cams +
                                     true.density +
                                     (1 | ls),
                      family = gaussian())

summary(model.cv)

# nope!

# probability of coverage?
model.conf <- glmmTMB(data = all.REM.longer.3,
                      formula = conf.cov ~ method +
                                           landscape +
                                           variability +
                                           n.cams +
                                           true.density +
                                           (1 | ls),
                    family = binomial())

summary(model.conf)

# well, no

# THE PROBLEM

# the problem is that passes are in no way related to the probability of use in these
# simulations - that's a major issue!
plot(passes$total.passes, 1 / passes$issf.mean)

# I'm guessing the UD simulations mess up my RSS calculation because of what happens out in the "fringes"
# most of the cameras are sample at RSS's OVER 1 - that really doesn't make much sense
# in general, the correction factor "smoothed" out REM estimates so that they were 
# indistinguishable from the naive REM

# none of the variation in predicted use is being captured by those simulations

# 09 Jan 2025 - what I need to go back and fix:

# - Use typical RSFs instead of SSFs to measure HS for the second method
# - Think about how the day range adjustment is calculated with end steps and so on
# - Fix those horrible UD simulations
# - Figure out the best way to propagate uncertainty (via delta and/or bootstrapping)
