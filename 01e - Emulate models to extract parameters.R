# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 01e - Emulate models to extract parameters
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 08 Apr 2025
# Date completed: 08 Apr 2025
# Date last modified: 08 Apr 2025
# R version: 4.4.3

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation
library(lubridate)            # work with dates
library(ctmm)                 # CTSP movement modeling
library(sf)                   # spatial operations
library(cowplot)              # multiple plots
library(glmmTMB)              # regression

#_______________________________________________________________________
# 2. Read in and bind data ----
#_______________________________________________________________________

all.files <- list.files("Derived data/Hares - Fitted models")

all.models <- list()

for (i in 1:length(all.files)) {
  
  load(paste0("Derived data/Hares - Fitted models/", all.files[i]))
  
  # bind into list
  all.models[[i]] <- top.model
  
  rm(top.model)
  
}

names(all.models) <- all.files

#_______________________________________________________________________
# 3. Calculate population-level mean parameters  ----
#_______________________________________________________________________

# tau position
meta.tau.p <- meta(x = all.models,
                   variable = "tau position",
                   level = 0.95,
                   sort = TRUE)

# tau velocity
meta.tau.v <- meta(x = all.models,
                   variable = "tau velocity",
                   level = 0.95,
                   sort = TRUE)

# speed
meta.speed <- meta(x = all.models,
                   variable = "speed",
                   level = 0.95,
                   sort = TRUE)

#_______________________________________________________________________
# 4. Sample possible model fits   ----

# here we'll use the "emulate" function in ctmm to make sure that models
# properly account for variance/covariance of parameters

# we'll take 500 draws from each model, making them available for
# simulation in the next step

n.draws <- 500

#_______________________________________________________________________

# list to hold all models
sampled.model.fits <- list()

# df to hold parameter values
sampled.model.params <- data.frame()

for (i in 1:length(all.models)) {
  
  # focal model
  focal.model <- all.models[[i]]
  
  focal.model.fits <- list()
  
  # sample all
  for (j in 1:n.draws) {
    
    focal.model.fits[[j]] <- emulate(focal.model,
                                     fast = TRUE)
    
    # focal.params
    focal.params <- data.frame(model = i,
                               draw = j,
                               tau.p = focal.model.fits[[j]]$tau[1],
                               tau.v = focal.model.fits[[j]]$tau[2],
                               sigma.maj = focal.model.fits[[j]]$sigma@par[1],
                               sigma.min = focal.model.fits[[j]]$sigma@par[2],
                               angle = focal.model.fits[[j]]$sigma@par[3],
                               speed.low = speed(focal.model.fits[[j]])$CI[1],
                               speed.est = speed(focal.model.fits[[j]])$CI[2],
                               speed.hig = speed(focal.model.fits[[j]])$CI[3])
    
    # bind in 
    sampled.model.params <- rbind(sampled.model.params, focal.params)
    
  }
  
  # bind in
  sampled.model.fits[[i]] <- focal.model.fits
  
}

#_______________________________________________________________________
# 5. Examine parameter distributions ----
#_______________________________________________________________________
# 5a. Tau position ----
#_______________________________________________________________________

ggplot(data = sampled.model.params,
       aes(x = tau.p)) +
  
  theme_bw() +
  
  geom_density(color = "darkblue",
               fill = "lightblue",
               alpha = 0.35,
               linewidth = 1.05) +
  
  theme(panel.grid = element_blank()) +
  
  coord_cartesian(xlim = c(0, 60000)) +
  
  scale_x_continuous(labels = round(c(0, 
                                      20000 / 3600, 
                                      40000 / 3600, 
                                      60000 / 3600))) +
  
  xlab("Tau position - hr")

#_______________________________________________________________________
# 5b. Tau velocity ----
#_______________________________________________________________________

ggplot(data = sampled.model.params,
       aes(x = tau.v)) +
  
  theme_bw() +
  
  geom_density(color = "darkblue",
               fill = "lightblue",
               alpha = 0.35,
               linewidth = 1.05) +
  
  theme(panel.grid = element_blank()) +
  
  coord_cartesian(xlim = c(0, 15000)) +
  
  scale_x_continuous(labels = round(c(0, 
                                      5000 / 3600, 
                                      10000 / 3600, 
                                      15000 / 3600))) +
  
  xlab("Tau velocity - hr")

#_______________________________________________________________________
# 5c. Speed ----
#_______________________________________________________________________

ggplot(data = sampled.model.params %>% filter(speed.est < 15),
       aes(x = speed.est)) +
  
  theme_bw() +
  
  geom_density(color = "darkblue",
               fill = "lightblue",
               alpha = 0.35,
               linewidth = 1.05) +
  
  theme(panel.grid = element_blank()) +
  
  xlab("Speed - km/day")

#_______________________________________________________________________
# 6. Examine parameter relationships ----
#_______________________________________________________________________
# 6a. Tau position and sigma major ----
#_______________________________________________________________________

ggplot(data = sampled.model.params,
       aes(x = tau.p,
           y = sigma.maj)) +
  
  theme_bw() +
  
  geom_point(alpha = 0.05,
             color = "darkblue") +
  
  theme(panel.grid = element_blank()) +
  
  xlab("Tau position - sec") +
  ylab("Sigma (major axis; m)")

#_______________________________________________________________________
# 6b. Tau velocity and sigma major ----
#_______________________________________________________________________

ggplot(data = sampled.model.params,
       aes(x = tau.v,
           y = sigma.maj)) +
  
  theme_bw() +
  
  geom_point(alpha = 0.05,
             color = "darkblue") +
  
  theme(panel.grid = element_blank()) +
  
  xlab("Tau velocity - sec") +
  ylab("Sigma (major axis; m)")

#_______________________________________________________________________
# 6c. Sigma major and sigma minor ----
#_______________________________________________________________________

ggplot(data = sampled.model.params,
       aes(x = sigma.maj,
           y = sigma.min)) +
  
  theme_bw() +
  
  geom_point(alpha = 0.05,
             color = "darkblue") +
  
  theme(panel.grid = element_blank()) +
  
  xlab("Sigma (major axis; m)") +
  ylab("Sigma (minor axis; m)")

#_______________________________________________________________________
# 6d. Tau position and speed ----
#_______________________________________________________________________

ggplot(data = sampled.model.params %>% filter(speed.est < 15),
       aes(x = tau.p,
           y = speed.est)) +
  
  theme_bw() +
  
  geom_point(alpha = 0.05,
             color = "darkblue") +
  
  theme(panel.grid = element_blank()) +
  
  xlab("Tau position - sec") +
  ylab("Speed (km/day)")

#_______________________________________________________________________
# 6e. Tau velocity and speed ----
#_______________________________________________________________________

ggplot(data = sampled.model.params %>% filter(speed.est < 15),
       aes(x = tau.v,
           y = speed.est)) +
  
  theme_bw() +
  
  geom_point(alpha = 0.05,
             color = "darkblue") +
  
  theme(panel.grid = element_blank()) +
  
  xlab("Tau velocity - sec") +
  ylab("Speed (km/day)")

#_______________________________________________________________________
# 6f. log(Tau velocity) and log(speed) ----
#_______________________________________________________________________

ggplot(data = sampled.model.params %>% filter(speed.est < 15),
       aes(x = log(tau.v),
           y = log(speed.est))) +
  
  theme_bw() +
  
  geom_point(alpha = 0.05,
             color = "darkblue") +
  
  theme(panel.grid = element_blank()) +
  
  xlab("ln(Tau velocity)") +
  ylab("ln(Speed)")

#_______________________________________________________________________
# 7. Regression ----
#_______________________________________________________________________

speed.reg <- glmmTMB(log(speed.est) ~ 
                       scale(log(tau.v)) + 
                       scale(tau.p) + 
                       scale(sigma.maj) +
                       (1 | model),
                     family = gaussian(),
                data = sampled.model.params %>% filter(speed.est < 15))

summary(speed.reg)

performance::performance(speed.reg)

# marg R2 of 0.565
# cond R2 of 0.966

#_______________________________________________________________________
# 8. Write to files ----
#_______________________________________________________________________

save(sampled.model.fits, file = "Derived data/Hares - Emulated models/sampled_fits.RData")

write.csv(sampled.model.params, "Derived data/Hares - Emulated models/sampled_params.csv")
