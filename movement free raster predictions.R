#_______________________________________________________________________
# 8. Movement-naive raster predictions ----

# here we'll calculate log-RSS (the movement-free selection kernel) from these
# movement-naive SSFs across the study area and then extract these values at 
# each camera as the "correction factor"

#_______________________________________________________________________
# 8a. Simple - low ----
#_______________________________________________________________________

SSF.params.1 <- SSF.params %>% 
  
  filter(landscape == "simple",
         variability == "low")

pred.raster.SL <- landscape.covs.simple$stem * SSF.params.1$betas[1] +
  landscape.covs.simple$edge * SSF.params.1$betas[2] +
  landscape.covs.simple$mature * SSF.params.1$betas[3]

terra::plot(pred.raster.SL)

#_______________________________________________________________________
# 8b. Complex - low ----
#_______________________________________________________________________

SSF.params.2 <- SSF.params %>% 
  
  filter(landscape == "complex",
         variability == "low")

pred.raster.CL <- landscape.covs.complex$stem * SSF.params.2$betas[1] +
  landscape.covs.complex$edge * SSF.params.2$betas[2] +
  landscape.covs.complex$mature * SSF.params.2$betas[3]

terra::plot(pred.raster.CL)

#_______________________________________________________________________
# 8c. Simple - high ----
#_______________________________________________________________________

SSF.params.3 <- SSF.params %>% 
  
  filter(landscape == "simple",
         variability == "high")

pred.raster.SH <- landscape.covs.simple$stem * SSF.params.3$betas[1] +
  landscape.covs.simple$edge * SSF.params.3$betas[2] +
  landscape.covs.simple$mature * SSF.params.3$betas[3]

terra::plot(pred.raster.SH)

#_______________________________________________________________________
# 8d. Complex - high ----
#_______________________________________________________________________

SSF.params.4 <- SSF.params %>% 
  
  filter(landscape == "complex",
         variability == "high")

pred.raster.CH <- landscape.covs.complex$stem * SSF.params.4$betas[1] +
  landscape.covs.complex$edge * SSF.params.4$betas[2] +
  landscape.covs.complex$mature * SSF.params.4$betas[3]

terra::plot(pred.raster.CH)

#_______________________________________________________________________
# 8e. Bind rasters together, rename, and save ----
#_______________________________________________________________________

SSF.pred.rasters <- c(pred.raster.SL, 
                      pred.raster.SH,
                      pred.raster.CL,
                      pred.raster.CH)

names(SSF.pred.rasters) <- c("SL", "SH", "CL", "CH")

writeRaster(SSF.pred.rasters, filename = "Rasters/SSF_pred.tif", overwrite = T)