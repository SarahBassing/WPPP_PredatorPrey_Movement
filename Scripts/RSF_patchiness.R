  #'  ====================================
  #'  Quantify patchiness of RSF surfaces
  #'  March 2023
  #'  ====================================
  #'  Script to quantify how patchy areas of low, medium, and high relative 
  #'  probability of selection are for each RSF. To be compared to average step 
  #'  lengths to determine if spatial scale of analysis and inference are compatible. 
  #'  ====================================
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(sf)
  library(stars)
  library(rgeos)
  library(terra)
  library(tidyverse)
  
  #'  RSFs: stack of 3 rasters, 1 per year
  md_smr_rsf <- rast("./Shapefiles/Predicted_RSFs/md_smr_RSFstack_global.tif") 
  md_wtr_rsf <- rast("./Shapefiles/Predicted_RSFs/md_wtr_RSFstack_global.tif")
  elk_smr_rsf <- rast("./Shapefiles/Predicted_RSFs/elk_smr_RSFstack_global.tif")
  elk_wtr_rsf <- rast("./Shapefiles/Predicted_RSFs/elk_wtr_RSFstack_global.tif")
  wtd_smr_rsf <- rast("./Shapefiles/Predicted_RSFs/wtd_smr_RSFstack_global.tif")
  wtd_wtr_rsf <- rast("./Shapefiles/Predicted_RSFs/wtd_wtr_RSFstack_global.tif")
  coug_smr_rsf <- rast("./Shapefiles/Predicted_RSFs/coug_smr_RSFstack_global.tif")
  coug_wtr_rsf <- rast("./Shapefiles/Predicted_RSFs/coug_wtr_RSFstack_global.tif")
  wolf_smr_rsf <- rast("./Shapefiles/Predicted_RSFs/wolf_smr_RSFstack_global.tif")
  wolf_wtr_rsf <- rast("./Shapefiles/Predicted_RSFs/wolf_wtr_RSFstack_global.tif")
  bob_smr_rsf <- rast("./Shapefiles/Predicted_RSFs/bob_smr_RSFstack_global.tif")
  bob_wtr_rsf <- rast("./Shapefiles/Predicted_RSFs/bob_wtr_RSFstack_global.tif")
  coy_smr_rsf <- rast("./Shapefiles/Predicted_RSFs/coy_smr_RSFstack_global.tif")
  coy_wtr_rsf <- rast("./Shapefiles/Predicted_RSFs/coy_wtr_RSFstack_global.tif")
  
  crs(md_smr_rsf, describe=TRUE, proj=TRUE)
  res(md_smr_rsf)
  
  
  
  #'  Function to reclassify raster based on high, medium, and low relative
  #'  probability of selection, detect patches of pixels, and measure average
  #'  size of patches
  patch_size <- function(rsf) {
    #'  Matrices to reclassify pixel values
    #'  High relative pr(select) = values 0.8 - 1.0
    hi_select <- matrix(c(0,0.8,0, 
                        0.8,1.0,1, 
                        1.0,2.0,0), ncol = 3, byrow = TRUE)
    #'  Medium relative pr(select) = values 0.4 - 0.8
    med_select <- matrix(c(0,0.4,0,
                         0.4,0.8,1,
                         0.8,2.0,0), ncol = 3, byrow = TRUE)
    #'  Low relative pr(select) = values 0.0 - 0.4
    low_select <- matrix(c(-1.0,0,0,
                         0,0.4,1,
                         0.4,2.0,0), ncol = 3, byrow= TRUE)
    
    #'  Reclassify raster using new classification
    hi_reclass <- classify(rsf, hi_select, include.lowest = FALSE)
    med_reclass <- classify(rsf, med_select, include.lowest = FALSE)
    low_reclass <- classify(rsf, low_select, include.lowest = FALSE)
    
    #'  Detect groups of reclassified cells surrounded by 0s
    hi_patch <- terra::patches(hi_reclass, 8, zeroAsNA = TRUE)
    med_patch <- terra::patches(med_reclass, 8, zeroAsNA = TRUE)
    low_patch <- terra::patches(hi_reclass, 8, zeroAsNA = TRUE)
    
    #'  Visualize
    plot(hi_patch)
    
    #'  Calculate area covered by cells (patch size)
    hi_area <- cellSize(hi_patch, unit="km") |> zonal(hi_patch, sum)
    med_area <- cellSize(med_patch, unit="km") |> zonal(med_patch, sum)
    low_area <- cellSize(low_patch, unit="km") |> zonal(low_patch, sum)
    
    #'  Calculate mean patch size
    hi_mean <- mean(hi_area); hi_sd <- sd(hi_area)
    med_mean <- mean(med_area); med_sd <- sd(med_area)
    low_mean <- mean(low_area); low_sd <- sd(low_area)
    
    #'  Save as table
    mean_patch_size_km <- c(hi_mean, med_mean, low_mean)
    sd_patch_size_km <- c(hi_sd, med_sd, low_sd)
    patch_value <- c("High selection", "Medium selection", "Low selection")
    patch_size_km <- cbind(patch_value, mean_patch_size_km, sd_patch_size_km)
    
    return(patch_size_km)
    
  }
  coug_smr_patch_size <- patch_size(coug_smr_rsf[[1]])
  wolf_smr_patch_size <- patch_size(wolf_smr_rsf[[1]])
  elk_smr_patch_size <- patch_size(elk_smr_rsf[[1]])
  md_smr_patch_size <- patch_size(md_smr_rsf[[1]])
  wtd_smr_patch_size <- patch_size(wtd_smr_rsf[[1]])
  
  coug_wtr_patch_size <- patch_size(coug_wtr_rsf[[1]])
  wolf_wtr_patch_size <- patch_size(wolf_wtr_rsf[[1]])
  elk_wtr_patch_size <- patch_size(elk_wtr_rsf[[1]])
  md_wtr_patch_size <- patch_size(md_wtr_rsf[[1]])
  wtd_wtr_patch_size <- patch_size(wtd_wtr_rsf[[1]])
  
  #'  Create result table and save
  smr_patch_size <- rbind(coug_smr_patch_size, wolf_smr_patch_size, elk_smr_patch_size, md_smr_patch_size, wtd_smr_patch_size)
  wtr_patch_size <- rbind(coug_wtr_patch_size, wolf_wtr_patch_size, elk_wtr_patch_size, md_wtr_patch_size, wtd_wtr_patch_size)
  patch_size <- rbind(smr_patch_size, wtr_patch_size)
  Species <- c("Cougar", "Wolf", "Elk", "Mule deer", "White-tailed deer",
               "Cougar", "Wolf", "Elk", "Mule deer", "White-tailed deer")
  Season <- c("Summer", "Summer", "Summer", "Summer", "Summer",
              "Winter", "Winter", "Winter", "Winter", "Winter")
  patch_summary_table <- cbind(Species, Season, patch_size)
  write.csv(patch_summary_table, "./Outputs/RSF_output/RSF_patch_size_summary.csv")
  
  
  
  tst_rast <- coug_smr_rsf[[1]]
  
  hi_risk <- matrix(c(0,0.8,0, 
                     0.8,1.0,1, 
                     1.0,2.0,0), ncol = 3, byrow = TRUE)
  med_risk <- matrix(c(0,0.4,0,
                       0.4,0.8,1,
                       0.8,2.0,0), ncol = 3, byrow = TRUE)
  low_risk <- matrix(c(-1.0,0,0,
                       0,0.4,1,
                       0.4,2.0,0), ncol = 3, byrow= TRUE)
  
  hi_risk_reclass <- classify(tst_rast, hi_risk, include.lowest = FALSE)
  med_risk_reclass <- classify(tst_rast, med_risk, include.lowest = FALSE)
  low_risk_reclass <- classify(tst_rast, low_risk, include.lowest = FALSE)
  
  
  hi_patch <- terra::patches(hi_risk_reclass, 8, zeroAsNA = TRUE)
  med_patch <- terra::patches(med_risk_reclass, 8, zeroAsNA = TRUE)
  low_patch <- terra::patches(hi_risk_reclass, 8, zeroAsNA = TRUE)
  
  hi_area <- cellSize(hi_patch, unit="km") |> zonal(hi_patch, sum)
  hi_mean <- mean(hi_area); hi_sd <- sd(hi_area)
  
  