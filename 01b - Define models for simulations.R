# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Density - movement simulation
# Script: 01b - Define models for simulations
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 08 Apr 2025
# Date completed: 08 Apr 2025
# Date last modified: 02 Jul 2026
# R version: 4.5.2

#_______________________________________________________________________
# 1. Load in packages ----
#_______________________________________________________________________

library(tidyverse)            # data cleaning and manipulation
library(lubridate)            # work with dates
library(ctmm)                 # CTSP movement modeling
library(glmmTMB)              # regression

#_______________________________________________________________________
# 2. Read in models ----
#_______________________________________________________________________

# directory
dir.model <- "D:/hare_project/data_analysis/General/hare-gps-processing-new/data_cleaned/"

# read in
# top models
top.models <- readRDS(paste0(dir.model, "top_models.rds"))

# model selection results
model.select <- readRDS(paste0(dir.model, "all_model_select.rds"))

# shift rownames to a column
model.select$model <- rownames(model.select) |>
  
  # remove digits
  gsub('[[:digit:]]+', '', x = _)

#_______________________________________________________________________
# 3. Distributions of top models  ----
#_______________________________________________________________________

model.select.top <- model.select |> filter(mod == 1) |>
  
  mutate(model = case_when(
    
    model %in% c("IID", "IID anisotropic") ~ "IID",
    model %in% c("OU", "OU anisotropic") ~ "OU",
    model %in% c("OUF", "OUF anisotropic") ~ "OUF",
    model %in% c("OUf", "OUf anisotropic") ~ "OUf"
    
  )
  
  )

ggplot(model.select.top) +
  
  theme_bw() +
  
  geom_bar(aes(x = model)) +
  
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black")) +
  
  xlab("Top CTSP model") +
  ylab("Total tracks")

#_______________________________________________________________________
# 4. Subset tau_p and tau_v models  ----
#_______________________________________________________________________

# top models
model.select.tau.p <- model.select.top |> filter(model %in% c("OU", "OUF", "OUf"))
model.select.tau.v <- model.select.top |> filter(model %in% c("OUF", "OUf"))

# indices
indices.tau.p <- model.select.tau.p$i
indices.tau.v <- model.select.tau.v$i

# subset top model lists
models.tau.p <- top.models[indices.tau.p]
models.tau.v <- top.models[indices.tau.v]

#_______________________________________________________________________
# 5. Calculate population-level mean parameters  ----
#_______________________________________________________________________

# tau position
meta.tau.p <- meta(x = models.tau.p,
                   variable = "tau position",
                   level = 0.95,
                   sort = TRUE)

# diffusion
meta.sigma <- meta(x = models.tau.p,
                   variable = "diffusion",
                   level = 0.95,
                   sort = TRUE)

# tau velocity
meta.tau.v <- meta(x = models.tau.v,
                   variable = "tau velocity",
                   level = 0.95,
                   sort = TRUE)

# speed
meta.speed <- meta(x = models.tau.v,
                   variable = "speed",
                   level = 0.95,
                   sort = TRUE)

# sigma range
calc_sigma_range <- function (.model) {
  
  sigma.df <- data.frame(
    
    major = .model$sigma@par[1],
    minor = .model$sigma@par[2],
    cov = .model$sigma[2]
    
  )
  
  return(sigma.df)
  
}

all.sigma <- do.call(rbind, lapply(models.tau.p, calc_sigma_range))

# medians are more informative here
median(all.sigma$major)
median(all.sigma$minor)
median(all.sigma$cov) # zero

ggplot(all.sigma) +
  
  theme_bw() +
  
  geom_point(aes(x = major,
                 y = minor),
             shape = 21) +
  
  geom_point(x = median(all.sigma$major),
             y = median(all.sigma$minor),
             color = "red",
             size = 3) +
  
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black")) +
  
  xlab("Major axis variance") +
  ylab("Minor axis variance")

#_______________________________________________________________________
# 6. Q1 - Construct "empirically-inspired" models ----

# we'll pick 3 levels of tau_v (seconds up to 2 hours)
# TV1: 0.5 mins
# TV2: 30 mins
# TV3: 120 mins

# our "hypothetical" fix rate starts at 1 hr

# importantly, means are in hours but should be in seconds

# for each simulation, we'll add the stationary mean and an a random orientation

#_______________________________________________________________________

# base models
model.base.tv1 <- ctmm(
  
  tau = c(meta.tau.p[1, 2] %#% "hour", 30),
  isotropic = FALSE,
  range = T,
  error = FALSE,
  sigma = matrix(c(median(all.sigma$major),
                   0,
                   0,
                   median(all.sigma$minor)),
                 nrow = 2)
  
)

model.base.tv2 <- ctmm(
  
  tau = c(meta.tau.p[1, 2] %#% "hour", 30 %#% "minute"),
  isotropic = FALSE,
  range = T,
  error = FALSE,
  sigma = matrix(c(median(all.sigma$major),
                   0,
                   0,
                   median(all.sigma$minor)),
                 nrow = 2)
  
)

model.base.tv3 <- ctmm(
  
  tau = c(meta.tau.p[1, 2] %#% "hour", 120 %#% "minute"),
  isotropic = FALSE,
  range = T,
  error = FALSE,
  sigma = matrix(c(median(all.sigma$major),
                   0,
                   0,
                   median(all.sigma$minor)),
                 nrow = 2)
  
)

# save base models
saveRDS(model.base.tv1, "data_derived/Q1_models/tv1.rds")
saveRDS(model.base.tv2, "data_derived/Q1_models/tv2.rds")
saveRDS(model.base.tv3, "data_derived/Q1_models/tv3.rds")

#_______________________________________________________________________
# 7. Q2 - Reasonable range of individual variation ----

# generate parameters for 1,000 hypothetical movement models
# we can add these into CTMMs for sampling later
# we'll use the 95% quantiles of the empirical means

extract_tau_p <- function (.model) {
  
  return(.model$tau[1])

}

all.tau.p <- do.call(rbind, lapply(models.tau.p, extract_tau_p))

#_______________________________________________________________________

Q2.params <- data.frame(
  
  tau.p = runif(1000,
                quantile(all.tau.p, prob = 0.025),
                quantile(all.tau.p, prob = 0.975)),
  
  tau.v = runif(1000,
                30,
                120 %#% "minute"),
  
  sigma.major = runif(1000,
                      quantile(all.sigma$major, prob = 0.025),
                      quantile(all.sigma$major, prob = 0.975)),
  
  sigma.minor = runif(1000,
                      quantile(all.sigma$minor, prob = 0.025),
                      quantile(all.sigma$minor, prob = 0.975))
  
)

# ensure that tau v is never > tau p
sum(Q2.params$tau.v > Q2.params$tau.p)

# replace with the tau p (these would be "OUf" processes)
Q2.params$tau.v[which(Q2.params$tau.v > Q2.params$tau.p)] <- 
  Q2.params$tau.p[which(Q2.params$tau.v > Q2.params$tau.p)]

# save to file
saveRDS(Q2.params, "data_derived/Q2_models/Q2_params.rds")
