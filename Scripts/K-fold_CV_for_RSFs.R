  #'  ================================
  #'  K-fold Cross-Validation for RSFs
  #'  WA Predator-Prey Project
  #'  Sarah B. Bassing
  #'  January 2022
  #'  ================================
  #'  Script to run K-fold cross-validation for RSFs following methods described 
  #'  in Boyce et al. (2002). Script pulls in used/available locations and matched 
  #'  covariates (prepared in the Collar_RSF_DataPrep.R script), then runs data 
  #'  through a series of functions that: 
  #'    1) reclassifies landcover data due to convergence issues for some 
  #'       categories identified in preliminary runs of RSFs; 
  #'    2) standardizes continuous covariates for each species, season, and 
  #'       year-specific data set & randomly assigns observations into K folds;
  #'    3) partitions K folds into training and testing data sets;
  #'    4) defines species, season, and year-specific GLMM
  #'    5) trains models on training data sets
  #'  Next it cross-validates k-fold models by:
  #'    1) Predict trained model results across study area, normalize values & plot
  #'    2) Bin normalized predictions
  #'    3) Overlay testing locations (used pts only) across binned RSF map
  #'    4) area-adjust frequency of testing pts per bin
  #'    5) Spearman's rank correlation between bin rank and frequency of test pts
  #'    ==============================

  #'  Clear memory
  rm(list=ls())
  
  #'  Load packages
  library(tidyverse)
  library(car)
  library(cvms)
  library(groupdata2)
  library(knitr)
  library(lme4)
  library(raster)
  library(sf)
  # library(doParallel)
  # library(parallel)
  # library(future.apply)

  #'  Load used and available locations, and covariate data
  load("./Outputs/RSF_pts/md_dat_all_2022-01-06.RData")
  load("./Outputs/RSF_pts/elk_dat_all_2022-01-06.RData")
  load("./Outputs/RSF_pts/wtd_dat_all_2022-01-06.RData")
  load("./Outputs/RSF_pts/coug_dat_all_2022-01-06.RData")
  load("./Outputs/RSF_pts/wolf_dat_all_2022-01-06.RData")
  load("./Outputs/RSF_pts/bob_dat_all_2022-01-06.RData")
  load("./Outputs/RSF_pts/coy_dat_all_2022-01-06.RData")
  
  #'  Read in study area grids (1km^2)
  NE_1km <- raster("./Shapefiles/NE_1km_grid.tif") #NE_1km_grid
  OK_1km <- raster("./Shapefiles/OK_1km_grid.tif") #OK_1km_grid
  
  #'  Convert rasters to pixels and extract coordinates (centroid of each cell)
  raster_dat <- function(r) {
    dots <- as(r, "SpatialPixelsDataFrame")
    gridID <- dots@grid.index
    coords <- coordinates(dots)
    pts <- as.data.frame(cbind(gridID, coords))
    pts$ID <- seq(1:nrow(pts))
    return(pts)
  }
  NE_pts <- raster_dat(NE_1km)
  OK_pts <- raster_dat(OK_1km)
  
  #'  Read in covariates extracted across each study area (1km resolution)
  load("./Outputs/Telemetry_covs/NE_covs_1km_2022-01-07.RData") 
  load("./Outputs/Telemetry_covs/OK_covs_1km_2022-01-07.RData")
  
  #'  Format study area-wide covariate data to include annually relevant data only   ##### JOINING INTRODUCES NAs
  NE.covs <- NE.covs.1km %>%      # NE.covs.30m
    mutate(StudyArea = "NE") %>%
    full_join(NE_pts, by = "ID") %>%
    dplyr::select(-gridID) %>%
    #'  Covariate values were extracted at masked locations for some reason 
    #'  need to exclude these because missing coordinate data when joined
    filter(!is.na(x))
  OK.covs <- OK.covs.1km %>%      # OK.covs.30m
    mutate(StudyArea = "OK") %>%
    full_join(OK_pts, by = "ID") %>%
    dplyr::select(-gridID) %>%
    #'  Covariate values were extracted at masked locations for some reason 
    #'  need to exclude these because missing coordinate data when joined
    filter(!is.na(x))
  SA.covs <- rbind(NE.covs, OK.covs)
  SA.covs.Year1 <- dplyr::select(SA.covs, -c(CanopyCover19, CanopyCover20, Dist2Edge19, Landcover_type19))
  names(SA.covs.Year1) <- c("ID", "Elev", "Slope", "RoadDen", "Dist2Water",
                            "HumanMod", "CanopyCover", "Dist2Edge", 
                            "Landcover_type", "StudyArea", "x", "y")
  SA.covs.Year2 <- dplyr::select(SA.covs, -c(CanopyCover18, CanopyCover20, Dist2Edge18, Landcover_type18))
  names(SA.covs.Year2) <- c("ID", "Elev", "Slope", "RoadDen", "Dist2Water",
                            "HumanMod", "CanopyCover", "Dist2Edge", 
                            "Landcover_type", "StudyArea", "x", "y")
  #'  Note: applying 2019 Dist2Edge and Landcover_type to Year3 data due to lack
  #'  of 2020 landcover data
  SA.covs.Year3 <- dplyr::select(SA.covs, -c(CanopyCover18, CanopyCover19, Dist2Edge18, Landcover_type18))
  names(SA.covs.Year3) <- c("ID", "Elev", "Slope", "RoadDen", "Dist2Water",
                            "HumanMod", "CanopyCover", "Dist2Edge", 
                            "Landcover_type", "StudyArea", "x", "y")
  #'  List study area covariates by year to mirror rest of data structure
  SA.covs_list <- list(SA.covs.Year1, SA.covs.Year2, SA.covs.Year3)
  NE.covs_list <- list(SA.covs.Year1[SA.covs.Year1$StudyArea == "NE",], SA.covs.Year2[SA.covs.Year2$StudyArea == "NE",], SA.covs.Year3[SA.covs.Year3$StudyArea == "NE",])
  OK.covs_list <- list(SA.covs.Year1[SA.covs.Year1$StudyArea == "OK",], SA.covs.Year2[SA.covs.Year2$StudyArea == "OK",], SA.covs.Year3[SA.covs.Year3$StudyArea == "OK",])
  
  
  ####  Data formatting  ####
  #'  =======================
  #'  Function to reclassify landcover_type into fewer categories
  class_landcov <- function(locs) {
    locs <- locs %>%
      mutate(
        Landcover_type = ifelse(Landcover_type == "Water", "Other", Landcover_type),
        Landcover_type = ifelse(Landcover_type == "Glacier", "Other", Landcover_type),
        Landcover_type = ifelse(Landcover_type == "Barren", "Other", Landcover_type),
        Landcover_type = ifelse(Landcover_type == "Wetland", "Wetland", Landcover_type),
        Landcover_type = ifelse(Landcover_type == "Wetland", "Wetland", Landcover_type),
        Landcover_type = ifelse(Landcover_type == "Woody Wetland", "Open Grass", Landcover_type),
        Landcover_type = ifelse(Landcover_type == "Mesic Grass", "Open Grass", Landcover_type),
        Landcover_type = ifelse(Landcover_type == "Xeric Grass", "Open Grass", Landcover_type),
        Landcover_type = ifelse(Landcover_type == "Mesic Shrub", "Shrub Mix", Landcover_type),
        Landcover_type = ifelse(Landcover_type == "Xeric Shrub", "Shrub Mix", Landcover_type),
        Landcover_type = ifelse(Landcover_type == "Forest", "Forest", Landcover_type),
        Landcover_type = ifelse(Landcover_type == "Agriculture", "Developed", Landcover_type),
        Landcover_type = ifelse(Landcover_type == "Commercial", "Developed", Landcover_type),
        Landcover_type = ifelse(Landcover_type == "Developed", "Developed", Landcover_type),
        Landcover_type = ifelse(Landcover_type == "310", "Developed", Landcover_type)
      )
    
    return(locs)
  }
  #'  Reclassify landcover data for each species
  md_dat_all <- class_landcov(md_dat_all)
  elk_dat_all <- class_landcov(elk_dat_all)
  wtd_dat_all <- class_landcov(wtd_dat_all)
  coug_dat_all <- class_landcov(coug_dat_all)
  wolf_dat_all <- class_landcov(wolf_dat_all)
  bob_dat_all <- class_landcov(bob_dat_all)
  coy_dat_all <- class_landcov(coy_dat_all)
  #'  Reclassify landcover data for study area wide data sets
  SA.covs_list <- lapply(SA.covs_list, class_landcov)
  NE.covs_list <- lapply(NE.covs_list, class_landcov)
  OK.covs_list <- lapply(OK.covs_list, class_landcov)

  
  #'  Landcover_type categories causing convergence issues for some species due 
  #'  to too few observations in some categories reducing categories further
  reclass_landcov <- function(locs) {
    locs <- locs %>%
      mutate(
        Landcover_type = as.character(as.factor(Landcover_type)),
        Landcover_type = ifelse(Landcover_type == "Developed", "Other", Landcover_type)
      )
    locs$Landcover_type <- droplevels(as.factor(locs$Landcover_type))
    locs$Landcover_type <- relevel(locs$Landcover_type, ref = "Forest")
    
    return(locs)
  }
  #'  Reclassify landcover data for each species even further
  md_dat_all_reclass <- reclass_landcov(md_dat_all)
  elk_dat_all_reclass <- reclass_landcov(elk_dat_all)
  wtd_dat_all_reclass <- reclass_landcov(wtd_dat_all)
  coug_dat_all_reclass <- reclass_landcov(coug_dat_all)
  wolf_dat_all_reclass <- reclass_landcov(wolf_dat_all)
  bob_dat_all_reclass <- reclass_landcov(bob_dat_all)
  coy_dat_all_reclass <- reclass_landcov(coy_dat_all)
  SA.covs_list_reclass <- lapply(SA.covs_list, reclass_landcov)
  NE.covs_list_reclass <- lapply(NE.covs_list, reclass_landcov)
  OK.covs_list_reclass <- lapply(OK.covs_list, reclass_landcov)
  
  #'  More reclassification required for all wolf models-
  #'  "Other", "Developed", & "Wetland" landcover types causing issues with model
  #'  convergence so lumping all together as one class
  reclass_wolf <- function(locs) {
    locs <- locs %>%
      mutate(
        Landcover_type = as.character(as.factor(Landcover_type)),
        Landcover_type = ifelse(Landcover_type == "Wetland", "Other", Landcover_type)
      )
    locs$Landcover_type <- droplevels(as.factor(locs$Landcover_type))
    locs$Landcover_type <- relevel(locs$Landcover_type, ref = "Forest")
    
    return(locs)
  }
  wolf_dat_all_reclass2 <- reclass_wolf(wolf_dat_all_reclass)
  SA.covs_list_wolf_reclass <- lapply(SA.covs_list_reclass, reclass_wolf)

  
  #'  List species-specific data by season
  #'  Using different landcover_type classifications for some species & seasons
  #'  due to convergence issues identified in initial runs of RSFs
  mdData_smr <- md_dat_all[md_dat_all$Season == "Summer18" | md_dat_all$Season == "Summer19" | md_dat_all$Season == "Summer20",]
  mdData_wtr <- md_dat_all[md_dat_all$Season == "Winter1819" | md_dat_all$Season == "Winter1920" | md_dat_all$Season == "Winter2021",]
  #'  Note the reclassified landcover_type data for elkData_winter
  elkData_smr <- elk_dat_all[elk_dat_all$Season == "Summer18" | elk_dat_all$Season == "Summer19" | elk_dat_all$Season == "Summer20",]
  elkData_wtr <- elk_dat_all_reclass[elk_dat_all_reclass$Season == "Winter1819" | elk_dat_all_reclass$Season == "Winter1920" | elk_dat_all_reclass$Season == "Winter2021",]
  #'  Note the reclassified landcover_type data for wtdData_winter 
  wtdData_smr <- wtd_dat_all[wtd_dat_all$Season == "Summer18" | wtd_dat_all$Season == "Summer19" | wtd_dat_all$Season == "Summer20",]
  wtdData_wtr <- wtd_dat_all_reclass[wtd_dat_all_reclass$Season == "Winter1819" | wtd_dat_all_reclass$Season == "Winter1920" | wtd_dat_all_reclass$Season == "Winter2021",]
  #'  Note the reclassified landcover_type data for cougData_winter
  cougData_smr <- coug_dat_all[coug_dat_all$Season == "Summer18" | coug_dat_all$Season == "Summer19" | coug_dat_all$Season == "Summer20",]
  cougData_wtr <- coug_dat_all_reclass[coug_dat_all_reclass$Season == "Winter1819" | coug_dat_all_reclass$Season == "Winter1920" | coug_dat_all_reclass$Season == "Winter2021",]
  #'  Note the double reclassified landcover_type data for wolfData
  wolfData_smr <- wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Summer18" | wolf_dat_all_reclass2$Season == "Summer19" | wolf_dat_all_reclass2$Season == "Summer20",]
  wolfData_wtr <- wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Winter1819" | wolf_dat_all_reclass2$Season == "Winter1920" | wolf_dat_all_reclass2$Season == "Winter2021",]
  #'  Note the reclassified landcover_type data for bobData
  bobData_smr <- bob_dat_all_reclass[bob_dat_all_reclass$Season == "Summer18" | bob_dat_all_reclass$Season == "Summer19" | bob_dat_all_reclass$Season == "Summer20",]
  bobData_wtr <- bob_dat_all_reclass[bob_dat_all_reclass$Season == "Winter1819" | bob_dat_all_reclass$Season == "Winter1920" | bob_dat_all_reclass$Season == "Winter2021",]
  #'  Note the reclassified landcover_type data for coyData
  coyData_smr <- coy_dat_all_reclass[coy_dat_all_reclass$Season == "Summer18" | coy_dat_all_reclass$Season == "Summer19" | coy_dat_all_reclass$Season == "Summer20",]
  coyData_wtr <- coy_dat_all_reclass[coy_dat_all_reclass$Season == "Winter1819" | coy_dat_all_reclass$Season == "Winter1920" | coy_dat_all_reclass$Season == "Winter2021",]

  
  ####  Data partitioning  ####
  #'  =========================
  #'  Define number of folds
  #'  Keep this in mind for training & testing functions below
  K <- 5
  
  #'  Function to standardize covariates & fold data sets
  #'  Note: center & scaling across all IDs but separately by spp & season
  Ztrans_Kfold <- function(dat, K){
    
    #'  Make categorical variables factors
    dat$ID <- as.factor(dat$ID)
    dat$Used <- as.factor(dat$Used)
    dat$Area <- as.factor(dat$Area)
    dat$Year <- as.factor(dat$Year)
    dat$Season <- as.factor(dat$Season)
    dat$Landcover <- as.factor(dat$Landcover)
    dat$Landcover_type <- droplevels(as.factor(dat$Landcover_type))
    dat$Landcover_type <- relevel(dat$Landcover_type, ref = "Forest")
    #'  Standardize continuous variables
    dat$Elev <- scale(dat$Elev)
    dat$Slope <- scale(dat$Slope)
    dat$TPI <- scale(dat$TPI)
    dat$RoadDen <- scale(dat$RoadDen)
    dat$Dist2Water <- scale(dat$Dist2Water)
    dat$HumanMod <- scale(dat$HumanMod)
    dat$CanopyCover <- scale(dat$CanopyCover)
    dat$Dist2Edge <- scale(dat$Dist2Edge)
    dat$PercForMix <- scale(dat$PercForMix)
    dat$PercXGrass <- scale(dat$PercXGrass)
    dat$PercXShrub <- scale(dat$PercXShrub)

    #'  Partition each data set into K folds
    #'  Use groupdata2 package with cat_col = "Used" to balance folds proportional
    #'  to 0's and 1's in each data set
    set.seed(2022)
    fold_df <- fold(dat, k = K, cat_col = "Used")
    fold_df <- as.data.frame(fold_df) 
    
    return(fold_df)
  }
  #'  Z-transform covariates and fold data for each spp & season
  mdDataz_smr <- Ztrans_Kfold(mdData_smr, K = K)
  mdDataz_wtr <- Ztrans_Kfold(mdData_wtr, K = K)
  elkDataz_smr <- Ztrans_Kfold(elkData_smr, K = K)
  elkDataz_wtr <- Ztrans_Kfold(elkData_wtr, K = K)
  wtdDataz_smr <- Ztrans_Kfold(wtdData_smr, K = K)
  wtdDataz_wtr <- Ztrans_Kfold(wtdData_wtr, K = K)
  cougDataz_smr <- Ztrans_Kfold(cougData_smr, K = K)
  cougDataz_wtr <- Ztrans_Kfold(cougData_wtr, K = K)
  wolfDataz_smr <- Ztrans_Kfold(wolfData_smr, K = K)
  wolfDataz_wtr <- Ztrans_Kfold(wolfData_wtr, K = K)
  bobDataz_smr <- Ztrans_Kfold(bobData_smr, K = K)
  bobDataz_wtr <- Ztrans_Kfold(bobData_wtr, K = K)
  coyDataz_smr <- Ztrans_Kfold(coyData_smr, K = K)
  coyDataz_wtr <- Ztrans_Kfold(coyData_wtr, K = K)
 
  
  #'  Function to partition standardized TRAINING data based on folds
  #'  First TRAINING data set excludes observations from 1st fold
  #'  Second TRAINING data set excludes observations from 2nd fold, etc.
  training_dat <- function(dat) {
    train1 <- dat[dat$.folds != 1,]
    train2 <- dat[dat$.folds != 2,] 
    train3 <- dat[dat$.folds != 3,] 
    train4 <- dat[dat$.folds != 4,] 
    train5 <- dat[dat$.folds != 5,] 
    
    training_list <- list(train1, train2, train3, train4, train5)
    
    return(training_list)
  }
  #'  Create a list of partitioned data sets using standardized data-
  #'  5 training data sets for each species and season
  mdData_smr_train <- training_dat(mdDataz_smr)
  mdData_wtr_train <- training_dat(mdDataz_wtr)
  elkData_smr_train <- training_dat(elkDataz_smr)
  elkData_wtr_train <- training_dat(elkDataz_wtr)
  wtdData_smr_train <- training_dat(wtdDataz_smr)
  wtdData_wtr_train <- training_dat(wtdDataz_wtr)
  cougData_smr_train <- training_dat(cougDataz_smr)
  cougData_wtr_train <- training_dat(cougDataz_wtr)
  wolfData_smr_train <- training_dat(wolfDataz_smr)
  wolfData_wtr_train <- training_dat(wolfDataz_wtr)
  bobData_smr_train <- training_dat(bobDataz_smr)
  bobData_wtr_train <- training_dat(bobDataz_wtr)
  coyData_smr_train <- training_dat(coyDataz_smr)
  coyData_wtr_train <- training_dat(coyDataz_wtr)

  #'  Double checked this split correctly
  unique(mdData_smr_train[[5]][".folds"]) # should have folds 1 - 4, not 5
  unique(mdData_smr_train[[2]][".folds"]) # should have folds 1 & 3-5, not 2
  
  
  #'  Function to withhold standardized TESTING data based on folds & training data
  #'  First TESTING data set includes only observations from 1st fold
  #'  Second TESTING data set includes only observations from 2nd fold, etc.
  testing_dat <- function(dat) {
    test1 <- dat[dat$.folds == 1,]
    test2 <- dat[dat$.folds == 2,]
    test3 <- dat[dat$.folds == 3,]
    test4 <- dat[dat$.folds == 4,]
    test5 <- dat[dat$.folds == 5,]
    
    testing_list <- list(test1, test2, test3, test4, test5)
    
    return(testing_list)
  }
  #'  Create a list of the withheld fold of data- 
  #'  5 testing data sets for each species and season
  mdData_smr_test <- testing_dat(mdDataz_smr)
  mdData_wtr_test <- testing_dat(mdDataz_wtr)
  elkData_smr_test <- testing_dat(elkDataz_smr)
  elkData_wtr_test <- testing_dat(elkDataz_wtr)
  wtdData_smr_test <- testing_dat(wtdDataz_smr)
  wtdData_wtr_test <- testing_dat(wtdDataz_wtr)
  cougData_smr_test <- testing_dat(cougDataz_smr)
  cougData_wtr_test <- testing_dat(cougDataz_wtr)
  wolfData_smr_test <- testing_dat(wolfDataz_smr)
  wolfData_wtr_test <- testing_dat(wolfDataz_wtr)
  bobData_smr_test <- testing_dat(bobDataz_smr)
  bobData_wtr_test <- testing_dat(bobDataz_wtr)
  coyData_smr_test <- testing_dat(coyDataz_smr)
  coyData_wtr_test <- testing_dat(coyDataz_wtr)

  #'  Double checked this split correctly
  unique(mdData_smr_test[[5]][".folds"]) # should only be fold 5
  unique(mdData_smr_test[[2]][".folds"]) # should only be fold 2
  
  
  #'  Save for later use
  save(mdData_smr_test, file = "./Outputs/RSF_output/Kfold_CV/mdData_smr_TestData.RData")
  save(mdData_wtr_test, file = "./Outputs/RSF_output/Kfold_CV/mdData_wtr_TestData.RData")
  save(elkData_smr_test, file = "./Outputs/RSF_output/Kfold_CV/elkData_smr_TestData.RData")
  save(elkData_wtr_test, file = "./Outputs/RSF_output/Kfold_CV/elkData_wtr_TestData.RData")
  save(wtdData_smr_test, file = "./Outputs/RSF_output/Kfold_CV/wtdData_smr_TestData.RData")
  save(wtdData_wtr_test, file = "./Outputs/RSF_output/Kfold_CV/wtdData_wtr_TestData.RData")
  save(cougData_smr_test, file = "./Outputs/RSF_output/Kfold_CV/cougData_smr_TestData.RData")
  save(cougData_wtr_test, file = "./Outputs/RSF_output/Kfold_CV/cougData_wtr_TestData.RData")
  save(wolfData_smr_test, file = "./Outputs/RSF_output/Kfold_CV/wolfData_smr_TestData.RData")
  save(wolfData_wtr_test, file = "./Outputs/RSF_output/Kfold_CV/wolfData_wtr_TestData.RData")
  save(bobData_smr_test, file = "./Outputs/RSF_output/Kfold_CV/bobData_smr_TestData.RData")
  save(bobData_wtr_test, file = "./Outputs/RSF_output/Kfold_CV/bobData_wtr_TestData.RData")
  save(coyData_smr_test, file = "./Outputs/RSF_output/Kfold_CV/coyData_smr_TestData.RData")
  save(coyData_wtr_test, file = "./Outputs/RSF_output/Kfold_CV/coyData_wtr_TestData.RData")
  
  
  ####  RSF models for K-fold CV  ####
  #'  Define species- and season-specific RSFs to cross-validate
  #'  Used backwards step selection to select best model for each species & season
  #'  in Resource_Selection_Function_Models.R script
  #'  Mule Deer models
  md_smr_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + Dist2Edge + Landcover_type + (1|ID)"
  md_wtr_mod <- "Used ~ 1 + Elev + Slope + RoadDen + Dist2Water + CanopyCover + Landcover_type + (1|ID)"
  #'  Elk models
  elk_smr_mod <- "Used ~ 1 + Elev + I(Elev2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  elk_wtr_mod <- "Used ~ 1 + Elev + I(Elev2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  #'  White-tailed Deer models
  wtd_smr_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + Dist2Edge + Landcover_type + (1|ID)"
  wtd_wtr_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  #'  Cougar models
  coug_smr_mod <- "Used ~ 1 + Elev + I(Elev2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  coug_wtr_mod <- "Used ~ 1 + Elev + I(Elev2) + Slope + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  #'  Wolf models
  wolf_smr_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + Dist2Water + HumanMod + Dist2Edge + Landcover_type + (1|ID)"
  wolf_wtr_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)" 
  #'  Bobcat models
  bob_smr_mod <- "Used ~ 1 + Elev + Slope + RoadDen + Dist2Water + HumanMod + Dist2Edge + Landcover_type + (1|ID)"
  bob_wtr_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + HumanMod + CanopyCover + Landcover_type + (1|ID)"
  #'  Coyote models
  coy_smr_mod <- "Used ~ 1 + Slope + RoadDen + Dist2Water + CanopyCover + Landcover_type + (1|ID)"  
  coy_wtr_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + Dist2Edge + Landcover_type + (1|ID)"
 
  
  ####  K-fold model training  ####
  #'  =============================
  #'  Function to run GLMM on K-folded data 
  glmm_fn <- function(dat, mod) {
    glmm_mod <- glmer(formula = mod, data = dat, family = binomial(link = "logit"))
    print(summary(glmm_mod))
    print(car::vif(glmm_mod))
      
    return(glmm_mod)
  }
  #'  Apply each list of species & season-specific k-folded data to glmm function
  #'  Each list contains 5 training data sets per species*season combination
  #'  This takes awhile!
  #'  Mule Deer K-fold model training
  md_kfold_smr <- lapply(mdData_smr_train, glmm_fn, mod = md_smr_mod)
  md_kfold_wtr <- lapply(mdData_wtr_train, glmm_fn, mod = md_wtr_mod)
  save(md_kfold_smr, file = paste0("./Outputs/RSF_output/Kfold_CV/md_kfold_smr_", Sys.Date(), ".RData"))
  save(md_kfold_wtr, file = paste0("./Outputs/RSF_output/Kfold_CV/md_kfold_wtr_", Sys.Date(), ".RData"))
  
  #'  Elk K-fold model training
  elk_kfold_smr <- lapply(elkData_smr_train, glmm_fn, mod = elk_smr_mod)
  elk_kfold_wtr <- lapply(elkData_wtr_train, glmm_fn, mod = elk_wtr_mod)
  save(elk_kfold_smr, file = paste0("./Outputs/RSF_output/Kfold_CV/elk_kfold_smr_", Sys.Date(), ".RData"))
  save(elk_kfold_wtr, file = paste0("./Outputs/RSF_output/Kfold_CV/elk_kfold_wtr_", Sys.Date(), ".RData"))
  
  #'  White-tailed Deer K-fold model training
  wtd_kfold_smr <- lapply(wtdData_smr_train, glmm_fn, mod = wtd_smr_mod)
  wtd_kfold_wtr <- lapply(wtdData_wtr_train, glmm_fn, mod = wtd_wtr_mod)
  save(wtd_kfold_smr, file = paste0("./Outputs/RSF_output/Kfold_CV/wtd_kfold_smr_", Sys.Date(), ".RData"))
  save(wtd_kfold_wtr, file = paste0("./Outputs/RSF_output/Kfold_CV/wtd_kfold_wtr_", Sys.Date(), ".RData"))
  
  #'  Cougar K-fold model training
  coug_kfold_smr <- lapply(cougData_smr_train, glmm_fn, mod = coug_smr_mod)
  coug_kfold_wtr <- lapply(cougData_wtr_train, glmm_fn, mod = coug_wtr_mod)
  save(coug_kfold_smr, file = paste0("./Outputs/RSF_output/Kfold_CV/coug_kfold_smr_", Sys.Date(), ".RData"))
  save(coug_kfold_wtr, file = paste0("./Outputs/RSF_output/Kfold_CV/coug_kfold_wtr_", Sys.Date(), ".RData"))
  
  #'  Wolf K-fold model training
  wolf_kfold_smr <- lapply(wolfData_smr_train, glmm_fn, mod = wolf_smr_mod)
  wolf_kfold_wtr <- lapply(wolfData_wtr_train, glmm_fn, mod = wolf_wtr_mod)
  save(wolf_kfold_smr, file = paste0("./Outputs/RSF_output/Kfold_CV/wolf_kfold_smr_", Sys.Date(), ".RData"))
  save(wolf_kfold_wtr, file = paste0("./Outputs/RSF_output/Kfold_CV/wolf_kfold_wtr_", Sys.Date(), ".RData"))

  #'  Bobcat K-fold model training
  bob_kfold_smr <- lapply(bobData_smr_train, glmm_fn, mod = bob_smr_mod)
  bob_kfold_wtr <- lapply(bobData_wtr_train, glmm_fn, mod = bob_wtr_mod)
  save(bob_kfold_smr, file = paste0("./Outputs/RSF_output/Kfold_CV/bob_kfold_smr_", Sys.Date(), ".RData"))
  save(bob_kfold_wtr, file = paste0("./Outputs/RSF_output/Kfold_CV/bob_kfold_wtr_", Sys.Date(), ".RData"))
  
  #'  Coyote K-fold model training
  coy_kfold_smr <- lapply(coyData_smr_train, glmm_fn, mod = coy_smr_mod)
  coy_kfold_wtr <- lapply(coyData_wtr_train, glmm_fn, mod = coy_wtr_mod)
  save(coy_kfold_smr, file = paste0("./Outputs/RSF_output/Kfold_CV/coy_kfold_smr_", Sys.Date(), ".RData"))
  save(coy_kfold_wtr, file = paste0("./Outputs/RSF_output/Kfold_CV/coy_kfold_wtr_", Sys.Date(), ".RData"))

  
  
  ####  Projecting trained models  ####
  #'  =================================
  #'  Function to find mean & standard deviation for original (raw) covariates
  #'  Necessary for standardizing study area-wide covs based on original models
  #'  Note: summarizes data by spp & season, same as how the data were separated 
  #'  for k-fold cv 
  cov_summary <- function(covs) {
    mu.cov <- covs %>% 
      summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)))
    sd.cov <- covs %>% 
      summarise(across(where(is.numeric), ~sd(.x, na.rm = TRUE)))
    mu.sd.cov <- rbind(mu.cov, sd.cov)
    parameter <- as.data.frame(c("Mean", "SD"))
    colnames(parameter) <- "Parameter"
    cov_summary <- cbind(parameter, mu.sd.cov)
    return(cov_summary)
  }
  #'  Summarize raw spp & season-specific covariate values 
  #'  Requires the untransformed covariates for each species & year
  mdCov_smr_summary <- cov_summary(mdData_smr)
  mdCov_wtr_summary <- cov_summary(mdData_wtr)
  elkCov_smr_summary <- cov_summary(elkData_smr)
  elkCov_wtr_summary <- cov_summary(elkData_wtr)
  wtdCov_smr_summary <- cov_summary(wtdData_smr)
  wtdCov_wtr_summary <- cov_summary(wtdData_wtr)
  cougCov_smr_summary <- cov_summary(cougData_smr)
  cougCov_wtr_summary <- cov_summary(cougData_wtr)
  wolfCov_smr_summary <- cov_summary(wolfData_smr)
  wolfCov_wtr_summary <- cov_summary(wolfData_wtr)
  bobCov_smr_summary <- cov_summary(bobData_smr)
  bobCov_wtr_summary <- cov_summary(bobData_wtr)
  coyCov_smr_summary <- cov_summary(coyData_smr)
  coyCov_wtr_summary <- cov_summary(coyData_wtr)

  #'  Standardize study area-wide covariates based on the mean [1] and SD [2] of 
  #'  the original covariate values that went into the training models
  scaling_covs <- function(mu.sd, covs) {
    scaling_covs <- covs %>%
      transmute(
        ID = ID,
        Elev = (Elev - mu.sd$Elev[1]) / mu.sd$Elev[2],
        Slope = (Slope - mu.sd$Slope[1]) / mu.sd$Slope[2],
        RoadDen = (RoadDen - mu.sd$RoadDen[1]) / mu.sd$RoadDen[2],
        Dist2Water = (Dist2Water - mu.sd$Dist2Water[1]) / mu.sd$Dist2Water[2],
        HumanMod = (HumanMod - mu.sd$HumanMod[1]) / mu.sd$HumanMod[2],
        CanopyCover = (CanopyCover - mu.sd$CanopyCover[1]) / mu.sd$CanopyCover[2],
        Dist2Edge = (Dist2Edge - mu.sd$Dist2Edge[1]) / mu.sd$Dist2Edge[2],
        Landcover = as.factor(Landcover_type),
        Landcover_Developed = as.numeric(ifelse(Landcover_type == "Developed", 1, 0)),
        Landcover_Grass = as.numeric(ifelse(Landcover_type == "Open Grass", 1, 0)),
        Landcover_Other = as.numeric(ifelse(Landcover_type == "Other", 1, 0)),
        Landcover_Shrub = as.numeric(ifelse(Landcover_type == "Shrub Mix", 1, 0)),
        Landcover_Wetland = as.numeric(ifelse(Landcover_type == "Wetland", 1, 0)),
        StudyArea = as.factor(StudyArea),
        x = as.numeric(x),
        y = as.numeric(y))
    return(scaling_covs)
  }
  #'  ATTENTION: Be sure to use the correct classification/reclassification of 
  #'  the landcover_type variables for each species and season!
  md_smr_zcovs <- lapply(OK.covs_list, scaling_covs, mu.sd = mdCov_smr_summary)
  md_wtr_zcovs <- lapply(OK.covs_list, scaling_covs, mu.sd = mdCov_wtr_summary)
  elk_smr_zcovs <- lapply(NE.covs_list, scaling_covs, mu.sd = elkCov_smr_summary)
  elk_wtr_zcovs <- lapply(NE.covs_list_reclass, scaling_covs, mu.sd = elkCov_wtr_summary)
  wtd_smr_zcovs <- lapply(NE.covs_list, scaling_covs, mu.sd = wtdCov_smr_summary)
  wtd_wtr_zcovs <- lapply(NE.covs_list_reclass, scaling_covs, mu.sd = wtdCov_wtr_summary)
  coug_smr_zcovs <- lapply(SA.covs_list, scaling_covs, mu.sd = cougCov_smr_summary)
  coug_wtr_zcovs <- lapply(SA.covs_list_reclass, scaling_covs, mu.sd = cougCov_wtr_summary)
  wolf_smr_zcovs <- lapply(SA.covs_list_wolf_reclass, scaling_covs, mu.sd = wolfCov_smr_summary)
  wolf_wtr_zcovs <- lapply(SA.covs_list_wolf_reclass, scaling_covs, mu.sd = wolfCov_wtr_summary)
  bob_smr_zcovs <- lapply(SA.covs_list_reclass, scaling_covs, mu.sd = bobCov_smr_summary)
  bob_wtr_zcovs <- lapply(SA.covs_list_reclass, scaling_covs, mu.sd = bobCov_wtr_summary)
  coy_smr_zcovs <- lapply(SA.covs_list_reclass, scaling_covs, mu.sd = coyCov_smr_summary)
  coy_wtr_zcovs <- lapply(SA.covs_list_reclass, scaling_covs, mu.sd = coyCov_wtr_summary)
  #' #'  fyi: mapply allows each call of the scaling_covs function to get the 1st,
  #' #'  2nd, 3rd, etc. element of BOTH lists, SIMPLIFY keeps output in list format
  #' coy_smr_zcovs <- mapply(scaling_covs, mu.sd = coyCov_smr_summary, covs = SA.covs_list_reclass, SIMPLIFY = FALSE)
  #' coy_wtr_zcovs <- mapply(scaling_covs, mu.sd = coyCov_wtr_summary, covs = SA.covs_list_reclass, SIMPLIFY = FALSE)
  
  
  #'  Read in saved k-fold trained model results
  load("./Outputs/RSF_output/Kfold_CV/coy_kfold_smr_2022-01-13.RData")
  load("./Outputs/RSF_output/Kfold_CV/coy_kfold_wtr_2022-01-13.RData")
  
  #'  Function to save parameter estimates from each trained model
  #'  Use coef(mod) to look at random effects estimates
  rounddig <- 2
  
  rsf_out <- function(mod, spp, season){
    betas <- mod@beta
    se <- sqrt(diag(vcov(mod)))
    z <- summary(mod)$coef[,3]
    pval <- summary(mod)$coef[,4]
    out <- as.data.frame(cbind(betas, se, pval)) %>%
      transmute(
        Species = rep(spp, nrow(.)),
        Season = rep(season, nrow(.)),
        Parameter = row.names(.),
        Estimate = round(betas, rounddig),
        SE = round(se, rounddig),
        Z = round(z, rounddig),
        Pval = round(pval, rounddig)) %>%
      dplyr::select(c(Species, Season, Parameter, Estimate)) %>%
      mutate(Parameter = ifelse(Parameter == "(Intercept)", "alpha", Parameter),
             Parameter = ifelse(Parameter == "Elev", "b.elev", Parameter),
             Parameter = ifelse(Parameter == "I(Elev^2)", "b.elev2", Parameter),
             Parameter = ifelse(Parameter == "Slope", "b.slope", Parameter),
             Parameter = ifelse(Parameter == "RoadDen", "b.road", Parameter),
             Parameter = ifelse(Parameter == "Dist2Water", "b.water", Parameter),
             Parameter = ifelse(Parameter == "HumanMod", "b.hm", Parameter),
             Parameter = ifelse(Parameter == "CanopyCover", "b.canopy", Parameter),
             Parameter = ifelse(Parameter == "Dist2Edge", "b.edge", Parameter),
             Parameter = ifelse(Parameter == "Landcover_typeDeveloped", "b.developed", Parameter),
             Parameter = ifelse(Parameter == "Landcover_typeOpen Grass", "b.grass", Parameter),
             Parameter = ifelse(Parameter == "Landcover_typeOther", "b.other", Parameter),
             Parameter = ifelse(Parameter == "Landcover_typeShrub Mix", "b.shrub", Parameter),
             Parameter = ifelse(Parameter == "Landcover_typeWetland", "b.wetland", Parameter)) %>%
      #'  Spread data so each row represents model coefficients for a single season, single species model
      pivot_wider(names_from = Parameter, values_from = Estimate) 
    
    #'  Covariates excluded from species-specific models not included in this data 
    #'  frame but necessary for predicting function to work below
    #'  Vector of columns names that need to be included in this data frame
    nms <- c("Species", "Season", "alpha", "b.elev", "b.elev2", "b.slope", "b.road", "b.water", "b.hm", "b.canopy", "b.edge", "b.grass", "b.other", "b.shrub", "b.wetland")
    #'  Identify if there are any missing column names in the data frame
    Missing <- setdiff(nms, names(out))
    #'  Add missing columns and fill with 0's
    out[Missing] <- 0
    #'  Return data frame based on full list of column names
    out <- out[nms]
    
    return(out)
  }
  #'  Extract coefficient estimates for each trained model
  coy_smr_trainout <- lapply(coy_kfold_smr, rsf_out, spp = "Coyote", season = "Summer")
  coy_wtr_trainout <- lapply(coy_kfold_wtr, rsf_out, spp = "Coyote", season = "Winter")
  
 
  #'  Function to predict across all grid cells based on RSF results
  #'  Should end up with 1 predicted value per grid cell
  #'  NOTE: I want the predict relative probability of selection from RSF so not 
  #'  using a logit transformation like I normally would with logistic regression.
  #'  Instead, dropping the intercept from the model and just exponentiating the 
  #'  coeffs*covs (Fieberg et al. 2020)
  predict_rsf <- function(coef, cov) {
    predict_rsf <- c()
    #'  Predict across each grid cell
    for(i in 1:nrow(cov)) {
      predict_rsf[i] <- exp(coef$b.elev*cov$Elev[i] + coef$b.elev2*I(cov$Elev[i]^2) + 
                            coef$b.slope*cov$Slope[i] + coef$b.road*cov$RoadDen[i] +
                            coef$b.water*cov$Dist2Water[i] + coef$b.hm*cov$HumanMod[i] +
                            coef$b.canopy*cov$CanopyCover[i] + coef$b.edge*cov$Dist2Edge[i] +
                              coef$b.developed*cov$Landcover_Developed[i] +
                            coef$b.grass*cov$Landcover_Grass[i] + coef$b.other*cov$Landcover_Other[i] +
                            coef$b.shrub*cov$Landcover_Shrub[i] + coef$b.wetland*cov$Landcover_Wetland[i])
    }
    predict_rsf <- as.data.frame(predict_rsf)
    predict_rsf <- cbind(cov$ID, cov$x, cov$y, cov$StudyArea, predict_rsf)
    colnames(predict_rsf) <- c("ID", "x", "y", "StudyArea", "predict_rsf")
    
    return(predict_rsf)
  }
  #'  Run estimated coefficients from trained models through function to predict 
  #'  relative probability of selection for each species, season, and year -- 
  #'  Prediction should be year-specific b/c some habitat variables change from
  #'  year to year so need to allow predictions to vary accordingly
  #'  Results in 5 predicted RSFs from the K-fold training models for each year 
  #'  per season and species
  coy_smr18_Kpredict <- lapply(coy_smr_trainout, predict_rsf, cov = coy_smr_zcovs[[1]])
  coy_smr19_Kpredict <- lapply(coy_smr_trainout, predict_rsf, cov = coy_smr_zcovs[[2]])
  coy_smr20_Kpredict <- lapply(coy_smr_trainout, predict_rsf, cov = coy_smr_zcovs[[3]])
  coy_smr_Kpredict <- list(coy_smr18_Kpredict, coy_smr19_Kpredict, coy_smr20_Kpredict)
  # save(coy_smr_Kpredict, file = paste0("./Outputs/RSF_output/Kfold_CV/coy_smr_Kpredict_", Sys.Date(), ".RData"))
  coy_wtr1819_Kpredict <- lapply(coy_wtr_trainout, predict_rsf, cov = coy_wtr_zcovs[[1]])
  coy_wtr1920_Kpredict <- lapply(coy_wtr_trainout, predict_rsf, cov = coy_wtr_zcovs[[2]])
  coy_wtr2021_Kpredict <- lapply(coy_wtr_trainout, predict_rsf, cov = coy_wtr_zcovs[[3]])
  coy_wtr_Kpredict <- list(coy_wtr1819_Kpredict, coy_wtr1920_Kpredict, coy_wtr2021_Kpredict)
  # save(coy_wtr_Kpredict, file = paste0("./Outputs/RSF_output/Kfold_CV/coy_wtr_Kpredict_", Sys.Date(), ".RData"))

  
  
  
  
  ####
  ####  NEED TO DEAL WITH THESE OUTLIERS BEFORE RE-SCALING!!!
  ####
  
  
  
  #'  Re-scale predicted RSF values between 0 & 1 for plotting
  RSF_rescale <- function(out) {
    rescale_val <- out %>%
      mutate(
        rescale_rsf = round(predict_rsf/(max(predict_rsf, na.rm = T)), digits = 4)
      ) %>%
      dplyr::select(-c(ID, StudyArea, predict_rsf))
    return(rescale_val)
  }
  #'  Rescale predicted RSF values within each list of lists
  coy_smr_Krsf <- lapply(coy_smr_Kpredict, lapply, RSF_rescale) # definitely some extreme outliers
  c1.1 <- coy_smr_Krsf[[1]][[1]]
  c2.1 <- coy_smr_Krsf[[1]][[1]]
  c3.1 <- coy_smr_Krsf[[1]][[1]]
  coy_wtr_Krsf <- lapply(coy_wtr_Kpredict, lapply, RSF_rescale)
   
  
  #'  Define desired projections
  sa_proj <- projection("EPSG:2855")  # NAD83(HARN) / Washington North
  
  #'  Load study area shapefiles
  OK.SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "METHOW_SA") %>%
    st_transform(crs = sa_proj)
  OK.SA <- as(OK.SA, "Spatial")
  NE.SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "NE_SA") %>%
    st_transform(crs = sa_proj)
  NE.SA <- as(NE.SA, "Spatial")
  
  #'  Rasterize predicted RSF values
  rasterize_rsf <- function(rsf_list) {
    df <- rsf_list
    #'  Identify coordinates of rsf predictions
    coordinates(df) <- ~ x + y
    #'  Coerce predictions to SpatialPixelsDataFrame
    gridded(df) <- TRUE
    #'  Coerce to raster
    rasterRSF <- raster(df)
    #'  Define projection
    crs(rasterRSF) <- sa_proj
    plot(rasterRSF)
    plot(NE.SA, add = T)
    plot(OK.SA, add = T)
    
    return(rasterRSF)
  }
  coy_smr_Kfoldraster <- lapply(coy_smr_Krsf, lapply, rasterize_rsf)
  coy_wtr_Kfoldraster <- lapply(coy_wtr_Krsf, lapply, rasterize_rsf)
  
  #'  Bin RSF predictions (equal-sample sized bin or quantiles)
  #'  https://stackoverflow.com/questions/57922248/r-version-of-esri-slice-tool
  bin_rsf <- function(rast, season, species) {
    trained_rsf <- rast
    #'  Create 10 breaks for re-scaled RSF values ranging 0 to 1
    breaks <- seq(0, 1, 1/10)
    #'  Group re-scaled RSF values into bins based on cutoffs
    quants <- quantile(sampleRegular(trained_rsf, ncell(trained_rsf)), breaks, na.rm = TRUE)
    #'  Create new raster of binned RSF values
    binned_rsf <- cut(trained_rsf, quants)
    plot(binned_rsf, legend = F, main = paste(season, species, "Predicted RSF Bins"))
    legend("topright", legend = c(seq(1,length(breaks)-1,1)),
           fill = terrain.colors(n, alpha = 1, rev = T), cex=0.85, bty = "n")
    plot(NE.SA, add = T)
    plot(OK.SA, add = T)

    return(binned_rsf)
  }
  #'  Bin k-fold RSF predictions
  #'  End up with 5 binned RSF rasters for each K-fold trained model per year
  coy_smr_Kbins <- lapply(coy_smr_Kfoldraster, lapply, bin_rsf, season = "Summer", species = "Coyote")
  coy_wtr_Kbins <- lapply(coy_wtr_Kfoldraster, lapply, bin_rsf, season = "Winter", species = "Coyote")
  
  
  
  
  
  #'  Calculate area of each bin
  
  
  
  
  
  
  
  ####  Cross-validate RSFs  ####
  #'  ===========================
  #'  Read in saved testing data
  load("./Outputs/RSF_output/Kfold_CV/mdData_smr_TestData.RData")
  load("./Outputs/RSF_output/Kfold_CV/mdData_wtr_TestData.RData")
  load("./Outputs/RSF_output/Kfold_CV/elkData_smr_TestData.RData")
  load("./Outputs/RSF_output/Kfold_CV/elkData_wtr_TestData.RData")
  load("./Outputs/RSF_output/Kfold_CV/wtdData_smr_TestData.RData")
  load("./Outputs/RSF_output/Kfold_CV/wtdData_wtr_TestData.RData")
  load("./Outputs/RSF_output/Kfold_CV/cougData_smr_TestData.RData")
  load("./Outputs/RSF_output/Kfold_CV/cougData_wtr_TestData.RData")
  load("./Outputs/RSF_output/Kfold_CV/wolfData_smr_TestData.RData")
  load("./Outputs/RSF_output/Kfold_CV/wolfData_wtr_TestData.RData")
  load("./Outputs/RSF_output/Kfold_CV/bobData_smr_TestData.RData")
  load("./Outputs/RSF_output/Kfold_CV/bobData_wtr_TestData.RData")
  load("./Outputs/RSF_output/Kfold_CV/coyData_smr_TestData.RData")
  load("./Outputs/RSF_output/Kfold_CV/coyData_wtr_TestData.RData")
  
  #'  Retain just the used locations from each testing data set and make spatial
  testing_pts <- function(test_dat) {
    used_locs <- test_dat %>%
      filter(Used == 1) %>%
      dplyr::select(c(x, y, ID, Season, Year, Area, Used, .folds)) %>%
      #'  Convert to an sf object
      st_as_sf(coords = c("x", "y"), crs = sa_proj)
    #' #'  Create a SpatialPointsDataFrame of the used locations
    #' coords <- used_locs[, c("x", "y")]
    #' dat <- used_locs[, c(3:8)]
    #' sa_proj <- CRS("EPSG:2855")
    #' spdf <- SpatialPointsDataFrame(coords = coords, data = dat, proj4string = sa_proj)
    #' 
    #' return(spdf)
    return(used_locs)
  }
  md_smr_test_used <- lapply(mdData_smr_test, testing_pts)
  md_wtr_test_used <- lapply(mdData_wtr_test, testing_pts)
  elk_smr_test_used <- lapply(elkData_smr_test, testing_pts)
  elk_wtr_test_used <- lapply(elkData_wtr_test, testing_pts)
  wtd_smr_test_used <- lapply(wtdData_smr_test, testing_pts)
  wtd_wtr_test_used <- lapply(wtdData_wtr_test, testing_pts)
  coug_smr_test_used <- lapply(cougData_smr_test, testing_pts)
  coug_wtr_test_used <- lapply(cougData_wtr_test, testing_pts)
  wolf_smr_test_used <- lapply(wolfData_smr_test, testing_pts)
  wolf_wtr_test_used <- lapply(wolfData_wtr_test, testing_pts)
  bob_smr_test_used <- lapply(bobData_smr_test, testing_pts)
  bob_wtr_test_used <- lapply(bobData_wtr_test, testing_pts)
  coy_smr_test_used <- lapply(coyData_smr_test, testing_pts)
  coy_wtr_test_used <- lapply(coyData_wtr_test, testing_pts)
  
  
  #'  Plot used testing observations across binned RSF predictions for respective
  #'  trained models (e.g., test1 against train1, test2 against train2, etc.)
  #'  extract raster value at each location
  used_bins <- function(Kbins, test_used){
    bin_freq <- c()
    for(i in seq_along(coy_smr_test_used)) {
      used_bins <- raster::extract(Kbins[[i]], test_used[[i]], df = TRUE)
      bin_freq <- rbind(bin_freq, used_bins)
    }
    hist(bin_freq$layer)
    return(bin_freq)
  }
  #'  Extract annual bin values at used locations from k-fold testing data
  #'  Remember: test_used is a list of 5 df & Kbins is 3 lists of 5 rasters
  #'  where Kbins = [[1]] for Year1, [[2]] for Year2, & [[3]] for Year3
  coy_smr_binfreq <- lapply(coy_smr_Kbins, used_bins, test_used = coy_smr_test_used)
  coy_wtr_binfreq <- lapply(coy_wtr_Kbins, used_bins, test_used = coy_wtr_test_used)
  
  # bin_freq <- c()
  # for(i in seq_along(coy_smr_test_used)) {
  #   used_bins <- raster::extract(coy_smr_Kbins[[1]][[i]], coy_smr_test_used[[i]], df = TRUE)
  #   bin_freq <- rbind(bin_freq, used_bins)
  # }
  # hist(bin_freq$layer)
  
  #'  Count frequency of used bins and area-weight frequency
  #'  Spearman's rank correlation for each k-fold model
  #'  Estimate mean rank correlation as measure of predictive ability for each 
  #'  spp/season/year-specific model
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  