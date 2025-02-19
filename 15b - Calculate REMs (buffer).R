# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 15b - Calculate REMs (buffer)
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 13 Dec 2024
# Date completed: 19 Feb 2025
# Date last modified: 19 Feb 2025
# R version: 4.2.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)       # tidy data cleaning and manipulation

#_______________________________________________________________________
# 2. Read in and clean data ----
#_______________________________________________________________________

passes <- read.csv(paste0(getwd(), "/Derived_data/For REM/final_passes_buff.csv"))

# remove first column and duplicates - for some reason there are duplicates!
passes.1 <- passes %>%
  
  dplyr::select(-X) %>%
  
  dplyr::distinct()

#_______________________________________________________________________
# 3. Calculate REM density by camera ----
#_______________________________________________________________________

# here we'll propagate uncertainty with a Monte Carlo resampling approach -
# all REM inputs (except constants) will be a drawn from a Normal(mean, sd)
# of that input

# we'll first calculate by-camera REM estimates, then use weighted bootstrapping
# to generate CVs and confidence intervals

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
                                           M2.mean = mean(focal.rem.calc$M2),
                                           M2.sd = sd(focal.rem.calc$M2),
                                           M3.mean = mean(focal.rem.calc$M3),
                                           M3.sd = sd(focal.rem.calc$M3),
                                           M4.mean = mean(focal.rem.calc$M4),
                                           M4.sd = sd(focal.rem.calc$M4),
                                           M5.mean = mean(focal.rem.calc$M5),
                                           M5.sd = sd(focal.rem.calc$M5))
      
      # bind into original df
      focal.cam <- cbind(focal.cam, focal.rem.calc.summary)
      
      # and bind into all.cams
      all.cams <- rbind(all.cams, focal.cam)
      
    } else {     # condition if n.passes = 0
      
      focal.rem.calc.summary <- data.frame(M1.mean = 0,
                                           M1.sd = 0,
                                           M2.mean = 0,
                                           M2.sd = 0,
                                           M3.mean = 0,
                                           M3.sd = 0,
                                           M4.mean = 0,
                                           M4.sd = 0,
                                           M5.mean = 0,
                                           M5.sd = 0)
      
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

boot_rem <- function(x,
                     combos = all.combos,
                     n.samp) {
  
  # loop through all combos
  all.rem.boot <- data.frame()
  
  # start time
  start.time <- Sys.time()
  
  for (i in 1:nrow(combos)) {
    
    focal.combo <- combos[i, ]
    
    # subset passes
    # make sure that we can also use this function for the after-only data
    if ("before" %in% unique(combos$trt)) {
      
      focal.passes <- x %>%
        
        filter(trt == focal.combo$trt,
               rep == focal.combo$rep,
               n.cams == focal.combo$n.cams,
               n.indiv == focal.combo$n.indiv)
      
    } else {
      
      focal.passes <- x %>%
        
        filter(rep == focal.combo$rep,
               n.cams == focal.combo$n.cams,
               n.indiv == focal.combo$n.indiv)  
      
    }
    
    
    
    # new df with by-camera probability weights
    focal.weights <- focal.passes %>%
      
      dplyr::select(cam.id,
                    M1.mean:M5.sd) %>%
      
      # allow zero SD having the greatest weight (i.e., not down-weighted)
      # inverse of the SD
      mutate(M1.w.raw = ifelse(M1.sd == 0,
                               NA,
                               1 / M1.sd),
             M2.w.raw = ifelse(M2.sd == 0,
                               NA,
                               1 / M2.sd),
             M3.w.raw = ifelse(M3.sd == 0,
                               NA,
                               1 / M3.sd),
             M4.w.raw = ifelse(M4.sd == 0,
                               NA,
                               1 / M4.sd),
             M5.w.raw = ifelse(M5.sd == 0,
                               NA,
                               1 / M5.sd)) %>%
      
      # summarize inverse SD weights to calculate probabilities
      mutate(M1.w.p = M1.w.raw / sum(M1.w.raw, na.rm = TRUE),
             M2.w.p = M2.w.raw / sum(M2.w.raw, na.rm = TRUE),
             M3.w.p = M3.w.raw / sum(M3.w.raw, na.rm = TRUE),
             M4.w.p = M4.w.raw / sum(M4.w.raw, na.rm = TRUE),
             M5.w.p = M5.w.raw / sum(M5.w.raw, na.rm = TRUE)) %>%
      
      # replace NAs with 1 to ensure that 0s aren't down-weighted
      mutate(M1.w.p = replace_na(M1.w.p, 1),
             M2.w.p = replace_na(M2.w.p, 1),
             M3.w.p = replace_na(M3.w.p, 1),
             M4.w.p = replace_na(M4.w.p, 1),
             M5.w.p = replace_na(M5.w.p, 1)) %>%
      
      # remove "raw" columns
      dplyr::select(-c(M1.w.raw:M5.w.raw))
    
    # conduct weighted bootstrapping
    boot.means <- data.frame(matrix(data = NA,
                                    nrow = n.samp,
                                    ncol = 5))
    
    # loop
    for (j in 1:n.samp) {
      
      # M1
      boot.means[j, 1] <- mean(focal.weights[sample(1:nrow(focal.weights), 
                                                    size = nrow(focal.weights), 
                                                    replace = TRUE, 
                                                    prob = focal.weights$M1.w.p), ]$M1.mean)
      
      # M2
      boot.means[j, 2] <- mean(focal.weights[sample(1:nrow(focal.weights), 
                                                    size = nrow(focal.weights), 
                                                    replace = TRUE, 
                                                    prob = focal.weights$M2.w.p), ]$M2.mean)
      
      # M3
      boot.means[j, 3] <- mean(focal.weights[sample(1:nrow(focal.weights), 
                                                    size = nrow(focal.weights), 
                                                    replace = TRUE, 
                                                    prob = focal.weights$M3.w.p), ]$M3.mean)
      
      # M4
      boot.means[j, 4] <- mean(focal.weights[sample(1:nrow(focal.weights), 
                                                    size = nrow(focal.weights), 
                                                    replace = TRUE, 
                                                    prob = focal.weights$M4.w.p), ]$M4.mean)
      
      # M5
      boot.means[j, 5] <- mean(focal.weights[sample(1:nrow(focal.weights), 
                                                    size = nrow(focal.weights), 
                                                    replace = TRUE, 
                                                    prob = focal.weights$M5.w.p), ]$M5.mean)
      
      
      
    }
    
    # add summaries to df
    # sampling error
    focal.combo$M1.mean <- mean(boot.means$X1)
    focal.combo$M2.mean <- mean(boot.means$X2)
    focal.combo$M3.mean <- mean(boot.means$X3)
    focal.combo$M4.mean <- mean(boot.means$X4)
    focal.combo$M5.mean <- mean(boot.means$X5)
    
    # sampling error
    focal.combo$M1.sd <- sd(boot.means$X1)
    focal.combo$M2.sd <- sd(boot.means$X2)
    focal.combo$M3.sd <- sd(boot.means$X3)
    focal.combo$M4.sd <- sd(boot.means$X4)
    focal.combo$M5.sd <- sd(boot.means$X5)
    
    # coefficients of variation
    focal.combo$M1.cv <- focal.combo$M1.sd / focal.combo$M1.mean
    focal.combo$M2.cv <- focal.combo$M2.sd / focal.combo$M2.mean
    focal.combo$M3.cv <- focal.combo$M3.sd / focal.combo$M3.mean
    focal.combo$M4.cv <- focal.combo$M4.sd / focal.combo$M4.mean
    focal.combo$M5.cv <- focal.combo$M5.sd / focal.combo$M5.mean
    
    # percentile confidence intervals
    # lower
    focal.combo$M1.l95 <- quantile(boot.means$X1, prob = 0.025)
    focal.combo$M2.l95 <- quantile(boot.means$X2, prob = 0.025)
    focal.combo$M3.l95 <- quantile(boot.means$X3, prob = 0.025)
    focal.combo$M4.l95 <- quantile(boot.means$X4, prob = 0.025)
    focal.combo$M5.l95 <- quantile(boot.means$X5, prob = 0.025)
    
    focal.combo$M1.u95 <- quantile(boot.means$X1, prob = 0.975)
    focal.combo$M2.u95 <- quantile(boot.means$X2, prob = 0.975)
    focal.combo$M3.u95 <- quantile(boot.means$X3, prob = 0.975)
    focal.combo$M4.u95 <- quantile(boot.means$X4, prob = 0.975)
    focal.combo$M5.u95 <- quantile(boot.means$X5, prob = 0.975)
    
    # bind into master df
    all.rem.boot <- rbind(all.rem.boot, focal.combo)
    
    # status message
    time.elapsed <- round(as.numeric(difftime(Sys.time(), start.time, units = "mins")), digits = 2)
    
    print(paste0("Completed combo ", i, " of ", nrow(combos), " - ", time.elapsed, " mins"))
    
  }
  
  # return
  return(all.rem.boot)
  
}

#_______________________________________________________________________
# 4b. Use function ----
#_______________________________________________________________________

all.rem.boot <- boot_rem(all.rem, combos = all.combos, n.samp = 5000)

#_______________________________________________________________________
# 5. Examine data ----
#_______________________________________________________________________

# examine observations with mean == 0
all.rem.boot %>% filter(M1.mean == 0)

all.rem %>% filter(trt == "before",
                   rep %in% c(2, 3),
                   n.indiv == 2,
                   n.cams == 4)

# no cameras detected anything! let's remove for the rest of the analysis
all.rem.boot.1 <- all.rem.boot %>%
  
  filter(M1.mean != 0)

#_______________________________________________________________________
# 6. Prepare data for plotting ----
#_______________________________________________________________________

all.rem.boot.l <- all.rem.boot.1 %>%
  
  # pivot longer and drop the M prefix
  pivot_longer(cols = 5:29,
               names_prefix = "M") %>%
  
  # add "method" and "stat" columns
  mutate(method = substring(name, 1 , 1),
         stat = substring(name, 3, 7)) %>%
  
  # drop "name"
  dplyr::select(-name) %>%
  
  # pivot_wider
  pivot_wider(names_from = stat)

#_______________________________________________________________________
# 7. Write to.csv ----
#_______________________________________________________________________

write.csv(all.rem.boot.l, paste0(getwd(), "/Derived_data/REM results/all_rem_boot_buff.csv"))

#_______________________________________________________________________
# 8. Calculate REM density by camera BA ----
#_______________________________________________________________________

# here we'll use all inputs from the before landscapes to compute the after

#_______________________________________________________________________
# 8a. Define function ----
#_______________________________________________________________________

calc_REM_BA <- function(x,
                        n.samp) {
  
  # split before and after
  xB <- x %>% filter(trt == "before")
  xA <- x %>% filter(trt == "after")
  
  # loop through cameras
  all.cams <- data.frame()
  
  for (i in 1:nrow(xB)) {
    
    # subset (assume both have the same attributes)
    focal.camB <- xB[i, ]
    focal.camA <- xA[i, ]
    
    # determine if the count in "after" is zero
    if (focal.camA$total.passes > 0) {
      
      # sample n.samp draws for each input value (from "before)
      focal.static.dr <- rnorm(n.samp, focal.camB$static.dr.mean, focal.camB$static.dr.se)
      focal.dr <- rnorm(n.samp, focal.camB$dr.mean, focal.camB$dr.se)
      focal.rsf <- rnorm(n.samp, focal.camB$rsf.mean, focal.camB$rsf.se)
      focal.issf <- rnorm(n.samp, focal.camB$issf.mean, focal.camB$issf.se)
      
      # calculate mean, SD, and CV for each approach
      # loop through all draws
      focal.rem.calc <- data.frame()
      
      for (j in 1:n.samp) {
        
        # df of all inputs (use correct ones)
        focal.input <- data.frame(total.passes = focal.camA$total.passes,
                                  static.dr = focal.static.dr[j],
                                  dr = focal.dr[j],
                                  rsf = focal.rsf[j],
                                  issf = focal.issf[j],
                                  lens = focal.camB$lens,
                                  days = focal.camB$days)
        
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
                                           M2.mean = mean(focal.rem.calc$M2),
                                           M2.sd = sd(focal.rem.calc$M2),
                                           M3.mean = mean(focal.rem.calc$M3),
                                           M3.sd = sd(focal.rem.calc$M3),
                                           M4.mean = mean(focal.rem.calc$M4),
                                           M4.sd = sd(focal.rem.calc$M4),
                                           M5.mean = mean(focal.rem.calc$M5),
                                           M5.sd = sd(focal.rem.calc$M5))
      
      # bind into original df (AFTER)
      focal.camA <- cbind(focal.camA, focal.rem.calc.summary)
      
      # and bind into all.cams
      all.cams <- rbind(all.cams, focal.camA)
      
    } else {     # condition if n.passes = 0
      
      focal.rem.calc.summary <- data.frame(M1.mean = 0,
                                           M1.sd = 0,
                                           M2.mean = 0,
                                           M2.sd = 0,
                                           M3.mean = 0,
                                           M3.sd = 0,
                                           M4.mean = 0,
                                           M4.sd = 0,
                                           M5.mean = 0,
                                           M5.sd = 0)
      
      # bind into original df (AFTER)
      focal.camA <- cbind(focal.camA, focal.rem.calc.summary)
      
      # and bind into all.cams
      all.cams <- rbind(all.cams, focal.camA)
      
    }
    
  }
  
  # return
  return(all.cams)
  
}

#_______________________________________________________________________
# 8a. Use functions ----
#_______________________________________________________________________

all.remBA <- calc_REM_BA(passes.1, n.samp = 500)

# bootstrap
a.combos <- expand.grid(rep = 1:3,
                        n.cams = c(4, 9, 16),
                        n.indiv = unique(passes$n.indiv))

all.remBA.boot <- boot_rem(all.remBA, combos = a.combos, n.samp = 5000)

#_______________________________________________________________________
# 8b. Examine data ----
#_______________________________________________________________________

# examine observations with mean == 0
all.remBA.boot %>% filter(M1.mean == 0)

#_______________________________________________________________________
# 8c. Prepare data for plotting ----
#_______________________________________________________________________

all.remBA.boot.l <- all.remBA.boot %>%
  
  # pivot longer and drop the M prefix
  pivot_longer(cols = 4:28,
               names_prefix = "M") %>%
  
  # add "method" and "stat" columns
  mutate(method = substring(name, 1 , 1),
         stat = substring(name, 3, 7)) %>%
  
  # drop "name"
  dplyr::select(-name) %>%
  
  # pivot_wider
  pivot_wider(names_from = stat)

#_______________________________________________________________________
# 8d. Write to.csv ----
#_______________________________________________________________________

write.csv(all.remBA.boot.l, paste0(getwd(), "/Derived_data/REM results/all_remBA_boot.csv"))
