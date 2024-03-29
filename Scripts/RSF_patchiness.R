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
    hi_select <- matrix(c(0,0.7,0, 
                        0.7,1.0,1, 
                        1.0,2.0,0), ncol = 3, byrow = TRUE)
    #'  Medium relative pr(select) = values 0.4 - 0.8
    med_select <- matrix(c(0,0.3,0,
                         0.3,0.7,1,
                         0.7,2.0,0), ncol = 3, byrow = TRUE)
    #'  Low relative pr(select) = values 0.0 - 0.4
    low_select <- matrix(c(-1.0,0,0,
                         0,0.3,1,
                         0.3,2.0,0), ncol = 3, byrow= TRUE)
    
    #'  Reclassify raster using new classification
    hi_reclass <- classify(rsf, hi_select, include.lowest = FALSE)
    med_reclass <- classify(rsf, med_select, include.lowest = FALSE)
    low_reclass <- classify(rsf, low_select, include.lowest = FALSE)
    
    #'  Detect groups of reclassified cells surrounded by 0s
    hi_patch <- terra::patches(hi_reclass, 8, zeroAsNA = TRUE)
    # plot(hi_patch, main = "high value patchs")
    
    med_patch <- terra::patches(med_reclass, 8, zeroAsNA = TRUE)
    # plot(med_patch, main = "medium value patches")
    
    low_patch <- terra::patches(low_reclass, 8, zeroAsNA = TRUE)
    # plot(low_patch, main = "low value patches")
    
    #'  Calculate area covered by cells (patch size)
    hi_area <- cellSize(hi_patch, unit="m") |> zonal(hi_patch, sum)
    med_area <- cellSize(med_patch, unit="m") |> zonal(med_patch, sum)
    low_area <- cellSize(low_patch, unit="m") |> zonal(low_patch, sum)
    
    # area_list <- list(hi_area, med_area, low_area)
    # return(area_list)
    
    #'  Calculate mean patch size
    hi_mean <- mean(hi_area$area); hi_sd <- sd(hi_area$area); hi_avg_length <- sqrt(hi_mean)
    med_mean <- mean(med_area$area); med_sd <- sd(med_area$area); med_avg_length <- sqrt(med_mean)
    low_mean <- mean(low_area$area); low_sd <- sd(low_area$area); low_avg_length <- sqrt(low_mean)
    
    all_patches <- rbind(hi_area, med_area, low_area)
    patch_mean <- mean(all_patches$area); patch_sd <- sd(all_patches$area); patch_length <- sqrt(patch_mean)

    #'  Save as table
    mean_patch_size_m <- c(hi_mean, med_mean, low_mean, patch_mean)
    mean_patch_size_m <- as.data.frame(mean_patch_size_m)
    sd_patch_size_m <- c(hi_sd, med_sd, low_sd, patch_sd)
    length_patch_size_m <- c(hi_avg_length, med_avg_length, low_avg_length, patch_length)
    patch_value <- c("High selection", "Medium selection", "Low selection", "All patches")
    patch_size_m <- cbind(patch_value, mean_patch_size_m, sd_patch_size_m, length_patch_size_m)
    names(patch_size_m) <- c("patch value", "mean_patch_area_m", "sd_area", "length_m")

    return(patch_size_m)
    
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
  Species <- rep(c("Cougar", "Wolf", "Elk", "Mule deer", "White-tailed deer"), times = 1, each = 4)
  Season <- rep(c("Summer", "Summer", "Summer", "Summer", "Summer"), times = 1, each = 4)
  smr_patches <- cbind(Species, Season, smr_patch_size)
  # save(smr_patches, file = "./Outputs/RSF_output/smr_patch_size.RData")
  
  wtr_patch_size <- rbind(coug_wtr_patch_size, wolf_wtr_patch_size, elk_wtr_patch_size, md_wtr_patch_size, wtd_wtr_patch_size)
  Species <- rep(c("Cougar", "Wolf", "Elk", "Mule deer", "White-tailed deer"), times = 1, each = 4)
  Season <- rep(c("Winter", "Winter", "Winter", "Winter", "Winter"), times = 1, each = 4)
  wtr_patches <- cbind(Species, Season, wtr_patch_size)
  save(wtr_patch_size, file = "./Outputs/RSF_output/wtr_patch_size.RData")
  
  patch_size <- rbind(smr_patches, wtr_patches)
  # write.csv(patch_size, "./Outputs/RSF_output/RSF_patch_size_summary.csv")
  
  predator_patch_size <- rbind(coug_smr_patch_size, wolf_smr_patch_size, coug_wtr_patch_size, wolf_wtr_patch_size) %>%
    filter(`patch value` == "All patches") %>%
    summarize(avg_patch_area_m = mean(mean_patch_area_m),
              sd_patch_area = sd(mean_patch_area_m),
              avg_patch_length_m = mean(length_m),
              sd_patch_length = sd(length_m)) %>%
    mutate(spp_group = "Predator") %>%
    relocate(spp_group, .before = avg_patch_area_m)
  prey_patch_size <- rbind(elk_smr_patch_size, md_smr_patch_size, wtd_smr_patch_size, elk_wtr_patch_size, md_wtr_patch_size, wtd_wtr_patch_size) %>%
    filter(`patch value` == "All patches") %>%
    summarize(avg_patch_area_m = mean(mean_patch_area_m),
              sd_patch_area = sd(mean_patch_area_m),
              avg_patch_length_m = mean(length_m),
              sd_patch_length = sd(length_m)) %>%
    mutate(spp_group = "Prey") %>%
    relocate(spp_group, .before = avg_patch_area_m)
  
  avg_pred_prey_patch_size <- rbind(predator_patch_size, prey_patch_size)
  write.csv(avg_pred_prey_patch_size, "./Outputs/RSF_output/avg_RSF_patch_size.csv")
  
  
  
  