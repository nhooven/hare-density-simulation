# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 15 - Calculate REMs
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 13 Dec 2024
# Date completed: 
# Date last modified: 18 Jan 2025
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation

#_______________________________________________________________________
# 2. Read in and clean data ----
#_______________________________________________________________________

passes <- read.csv(paste0(getwd(), "/Derived_data/For REM/final_passes.csv"))

# remove first column and duplicates - for some reason there are duplicates!
passes.1 <- passes %>%
  
  dplyr::select(-X) %>%
  
  dplyr::distinct()
  
#_______________________________________________________________________
# 3. Calculate REM density ----
#_______________________________________________________________________

# here we'll propagate uncertainty with a Monte Carlo resampling approach -
# all REM inputs (except constants) will be a drawn from a Normal(mean, sd)
# of that input

# we'll first calculate by-camera REM estimates, then use each draw (n = n.samp)
# to calculate a mean, SD, CV, and confidence intervals for each "combo"

#_______________________________________________________________________
# 3a. Define function ----
#_______________________________________________________________________

calc_REM_bycam <- function(x,
                           n.samp) {
  
  # loop through cameras
  all.cams <- data.frame()
  
  for (i in 1:nrow(x)) {
    
    # subset
    focal.cam <- x[i, ]
    
    # determine if the count is zero
    if (focal.cam$total.passes > 0) {
      
     # sample n.samp draws for each input value
     focal.static.dr <- rnorm(n.samp, focal.cam$static.dr.mean, focal.cam$static.dr.se)
     focal.dr <- rnorm(n.samp, focal.cam$dr.mean, focal.cam$dr.se)
     focal.rsf <- rnorm(n.samp, focal.cam$rsf.mean, focal.cam$rsf.se)
     focal.issf <- rnorm(n.samp, focal.cam$issf.mean, focal.cam$issf.se)
     
     # calculate mean, SD, and CV for each approach
     # loop through all draws
     focal.rem.calc <- data.frame()
     
     for (j in 1:n.samp) {
       
       # df of all inputs
       focal.input <- data.frame(total.passes = focal.cam$total.passes,
                                 static.dr = focal.static.dr[j],
                                 dr = focal.dr[j],
                                 rsf = focal.rsf[j],
                                 issf = focal.issf[j],
                                 lens = focal.cam$lens,
                                 days = focal.cam$days)
       
       # and df to store outputs
       focal.output <- data.frame(M1 = NA,
                                  M2 = NA,
                                  M3 = NA,
                                  M4 = NA,
                                  M5 = NA)
       
       # M1 - full naive
       focal.output$M1 <- (focal.input$total.passes / focal.input$days) *
                          (pi / ((focal.input$static.dr * 1000) * 3.5 * (2.0 + focal.input$lens))) *
                          10000
       
       # M2 - RSF correction
       focal.output$M2 <- ((focal.input$total.passes / focal.input$days) *
                          (pi / ((focal.input$static.dr * 1000) * 3.5 * (2.0 + focal.input$lens))) *
                          10000) * focal.input$rsf
       
       # M3 - iSSF correction, no movement
       focal.output$M3 <- ((focal.input$total.passes / focal.input$days) *
                          (pi / ((focal.input$static.dr * 1000) * 3.5 * (2.0 + focal.input$lens))) *
                          10000) * focal.input$issf
       
       # M4 - No correction, movement
       focal.output$M4 <- ((focal.input$total.passes / focal.input$days) *
                          (pi / ((focal.input$dr * 1000) * 3.5 * (2.0 + focal.input$lens))) *
                          10000)
       
       # M5 - iSSF correction, movement
       focal.output$M5 <- ((focal.input$total.passes / focal.input$days) *
                          (pi / ((focal.input$dr * 1000) * 3.5 * (2.0 + focal.input$lens))) *
                          10000) * focal.input$issf
       
       # bind into df
       focal.rem.calc <- rbind(focal.rem.calc, focal.output)
      
    }
    
    # calculate summary statistics
    focal.rem.calc.summary <- data.frame(M1.mean = mean(focal.rem.calc$M1),
                                         M1.sd = sd(focal.rem.calc$M1),
                                         M1.cv = sd(focal.rem.calc$M1) / mean(focal.rem.calc$M1),
                                         M2.mean = mean(focal.rem.calc$M2),
                                         M2.sd = sd(focal.rem.calc$M2),
                                         M2.cv = sd(focal.rem.calc$M2) / mean(focal.rem.calc$M2),
                                         M3.mean = mean(focal.rem.calc$M3),
                                         M3.sd = sd(focal.rem.calc$M3),
                                         M3.cv = sd(focal.rem.calc$M3) / mean(focal.rem.calc$M3),
                                         M4.mean = mean(focal.rem.calc$M4),
                                         M4.sd = sd(focal.rem.calc$M4),
                                         M4.cv = sd(focal.rem.calc$M4) / mean(focal.rem.calc$M4),
                                         M5.mean = mean(focal.rem.calc$M5),
                                         M5.sd = sd(focal.rem.calc$M5),
                                         M5.cv = sd(focal.rem.calc$M5) / mean(focal.rem.calc$M5))
    
    # bind into original df
    focal.cam <- cbind(focal.cam, focal.rem.calc.summary)
    
    # and bind into all.cams
    all.cams <- rbind(all.cams, focal.cam)
      
    } else {     # condition if n.passes = 0
      
      focal.rem.calc.summary <- data.frame(M1.mean = 0,
                                           M1.sd = 0,
                                           M1.cv = NaN,
                                           M2.mean = 0,
                                           M2.sd = 0,
                                           M2.cv = NaN,
                                           M3.mean = 0,
                                           M3.sd = 0,
                                           M3.cv = NaN,
                                           M4.mean = 0,
                                           M4.sd = 0,
                                           M4.cv = NaN,
                                           M5.mean = 0,
                                           M5.sd = 0,
                                           M5.cv = NaN)
      
      # bind into original df
      focal.cam <- cbind(focal.cam, focal.rem.calc.summary)
      
      # and bind into all.cams
      all.cams <- rbind(all.cams, focal.cam)
      
    }
    
  }
  
  # return
  return(all.cams)
  
}
  
#_______________________________________________________________________
# 3b. Use function ----
#_______________________________________________________________________

all.rem <- calc_REM_bycam(passes.1, n.samp = 500)

#_______________________________________________________________________
# 4. Calculate mean, CV, and confidence intervals by combination ----
#_______________________________________________________________________

# expand grid
all.combos <- expand.grid(trt = c("before", "after"),
                          rep = 1:3,
                          n.cams = c(4, 9, 16),
                          n.indiv = unique(passes$n.indiv))

#_______________________________________________________________________
# 4a. Define function ----
#_______________________________________________________________________

mean_rem <- function(x,
                     combos = all.combos,
                     n.samp) {
  
  # loop through all combos
  # start time
  start.time <- Sys.time()
  
  for (i in 1:nrow(combos)) {
    
    focal.combo <- combos[i, ]
    
    # subset passes
    focal.passes <- x %>%
      
      filter(trt == focal.combo$trt,
             rep == focal.combo$rep,
             n.cams == focal.combo$n.cams,
             n.indiv == focal.combo$n.indiv)
    
    # take draws
    # sample n.samp draws for each input value
    focal.M1 <- rnorm(n.samp, focal.passes$M1.mean, focal.passes$M1.sd)
    focal.M2 <- rnorm(n.samp, focal.passes$M2.mean, focal.passes$M2.sd)
    focal.M3 <- rnorm(n.samp, focal.passes$M3.mean, focal.passes$M3.sd)
    focal.M4 <- rnorm(n.samp, focal.passes$M4.mean, focal.passes$M4.sd)
    focal.M5 <- rnorm(n.samp, focal.passes$M5.mean, focal.passes$M5.sd)
    
    for (j in 1:n.samp) {
      
      # df of all REM draws
      focal.draws <- data.frame(total.passes = focal.cam$total.passes,
                                static.dr = focal.static.dr[j],
                                dr = focal.dr[j],
                                rsf = focal.rsf[j],
                                issf = focal.issf[j],
                                lens = focal.cam$lens,
                                days = focal.cam$days)
      
      # and df to store outputs
      focal.output <- data.frame(M1 = NA,
                                 M2 = NA,
                                 M3 = NA,
                                 M4 = NA,
                                 M5 = NA)
      
    }
    
  }
  
  
  
}


  
  
  
  
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
# 4. Bootstrap within each trt-rep-n.cams-n.indiv combinations ----
#_______________________________________________________________________


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
