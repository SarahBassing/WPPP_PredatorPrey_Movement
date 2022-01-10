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
  #'    6) 
  #'   

  #'  Clear memory
  rm(list=ls())
  
  #'  Load packages
  library(tidyverse)
  library(car)
  library(cvms)
  library(groupdata2)
  library(knitr)
  library(lme4)
  library(doParallel)
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
  
  ####  Data formatting  ####
  #'  =======================
  #'  Function to reclassify landcover_type categories
  class_landcov <- function(locs) {
    #'  Re-classify landcover into fewer categories
    #'  Based on Taylor's feedback
    #'  Open grass: mesic grass, xeric grass, wetland woody
    #'  Shrubby mix: mesic shrub, xeric shrub
    #'  Forest
    #'  Wetland
    #'  Other: water, barren, glacier
    #'  Developed: agriculture, commercial, developed
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
  #'  Reclassify landcover data for each species
  md_dat_all <- class_landcov(md_dat_all)
  elk_dat_all <- class_landcov(elk_dat_all)
  wtd_dat_all <- class_landcov(wtd_dat_all)
  coug_dat_all <- class_landcov(coug_dat_all)
  wolf_dat_all <- class_landcov(wolf_dat_all)
  bob_dat_all <- class_landcov(bob_dat_all)
  coy_dat_all <- class_landcov(coy_dat_all)
  
  #'  Landcover_type categories causing convergence issues for some species due to
  #'  too few observations in some categories (e.g., "Other", "Wetland") so
  #'  reclassifying into fewer categories
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

  
  #'  List species-specific data by season
  #'  Using different landcover_type classifications for some species & seasons
  #'  due to convergence issues identified in initial runs of RSFs
  mdData_smr <- list(md_dat_all[md_dat_all$Season == "Summer18",], md_dat_all[md_dat_all$Season == "Summer19",], md_dat_all[md_dat_all$Season == "Summer20",])
  mdData_wtr <- list(md_dat_all[md_dat_all$Season == "Winter1819",], md_dat_all[md_dat_all$Season == "Winter1920",], md_dat_all[md_dat_all$Season == "Winter2021",])
  #'  Note the reclassified landcover_type data for elkData_winter
  elkData_smr <- list(elk_dat_all[elk_dat_all$Season == "Summer18",], elk_dat_all[elk_dat_all$Season == "Summer19",], elk_dat_all[elk_dat_all$Season == "Summer20",])
  elkData_wtr <- list(elk_dat_all_reclass[elk_dat_all_reclass$Season == "Winter1819",], elk_dat_all_reclass[elk_dat_all_reclass$Season == "Winter1920",], elk_dat_all_reclass[elk_dat_all_reclass$Season == "Winter2021",])
  #'  Note the reclassified landcover_type data for wtdData_winter 
  wtdData_smr <- list(wtd_dat_all[wtd_dat_all$Season == "Summer18",], wtd_dat_all[wtd_dat_all$Season == "Summer19",], wtd_dat_all[wtd_dat_all$Season == "Summer20",])
  wtdData_wtr <- list(wtd_dat_all_reclass[wtd_dat_all_reclass$Season == "Winter1819",], wtd_dat_all_reclass[wtd_dat_all_reclass$Season == "Winter1920",], wtd_dat_all_reclass[wtd_dat_all_reclass$Season == "Winter2021",])
  #'  Note the reclassified landcover_type data for cougData_winter
  cougData_smr <- list(coug_dat_all[coug_dat_all$Season == "Summer18",], coug_dat_all[coug_dat_all$Season == "Summer19",], coug_dat_all[coug_dat_all$Season == "Summer20",])
  cougData_wtr <- list(coug_dat_all_reclass[coug_dat_all_reclass$Season == "Winter1819",], coug_dat_all_reclass[coug_dat_all_reclass$Season == "Winter1920",], coug_dat_all_reclass[coug_dat_all_reclass$Season == "Winter2021",])
  #'  Note the double reclassified landcover_type data for wolfData
  wolfData_smr <- list(wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Summer18",], wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Summer19",], wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Summer20",])
  wolfData_wtr <- list(wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Winter1819",], wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Winter1920",], wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Winter2021",])
  #'  Note the reclassified landcover_type data for bobData
  bobData_smr <- list(bob_dat_all_reclass[bob_dat_all_reclass$Season == "Summer18",], bob_dat_all_reclass[bob_dat_all_reclass$Season == "Summer19",], bob_dat_all_reclass[bob_dat_all_reclass$Season == "Summer20",])
  bobData_wtr <- list(bob_dat_all_reclass[bob_dat_all_reclass$Season == "Winter1819",], bob_dat_all_reclass[bob_dat_all_reclass$Season == "Winter1920",], bob_dat_all_reclass[bob_dat_all_reclass$Season == "Winter2021",])
  #'  Note the reclassified landcover_type data for coyData
  coyData_smr <- list(coy_dat_all_reclass[coy_dat_all_reclass$Season == "Summer18",], coy_dat_all_reclass[coy_dat_all_reclass$Season == "Summer19",], coy_dat_all_reclass[coy_dat_all_reclass$Season == "Summer20",])
  coyData_wtr <- list(coy_dat_all_reclass[coy_dat_all_reclass$Season == "Winter1819",], coy_dat_all_reclass[coy_dat_all_reclass$Season == "Winter1920",], coy_dat_all_reclass[coy_dat_all_reclass$Season == "Winter2021",])

  
  ####  Data partitioning  ####
  #'  =========================
  #'  Define number of folds
  #'  Keep this in mind for training & testing functions below
  K <- 5
  
  #'  Function to standardize covariates & fold data sets
  #'  Note: center & scaling across all IDs but separately by spp, season, & year
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
  #'  Z-transform covariates and fold data for each spp, season, & year
  mdDataz_smr <- lapply(mdData_smr, Ztrans_Kfold, K = K)
  mdDataz_wtr <- lapply(mdData_wtr, Ztrans_Kfold, K = K)
  elkDataz_smr <- lapply(elkData_smr, Ztrans_Kfold, K = K)
  elkDataz_wtr <- lapply(elkData_wtr, Ztrans_Kfold, K = K)
  wtdDataz_smr <- lapply(wtdData_smr, Ztrans_Kfold, K = K)
  wtdDataz_wtr <- lapply(wtdData_wtr, Ztrans_Kfold, K = K)
  cougDataz_smr <- lapply(cougData_smr, Ztrans_Kfold, K = K)
  cougDataz_wtr <- lapply(cougData_wtr, Ztrans_Kfold, K = K)
  wolfDataz_smr <- lapply(wolfData_smr, Ztrans_Kfold, K = K)
  wolfDataz_wtr <- lapply(wolfData_wtr, Ztrans_Kfold, K = K)
  bobDataz_smr <- lapply(bobData_smr, Ztrans_Kfold, K = K)
  bobDataz_wtr <- lapply(bobData_wtr, Ztrans_Kfold, K = K)
  coyDataz_smr <- lapply(coyData_smr, Ztrans_Kfold, K = K)
  coyDataz_wtr <- lapply(coyData_wtr, Ztrans_Kfold, K = K)

  
  #'  Function to partition standardized TRAINING data based on folds
  #'  First TRAINING data set excludes observations from 1st fold
  #'  Second TRAINING data set excludtes observations from 2nd fold, etc.
  training_dat <- function(dat) {
    train1 <- dat[dat$.folds != 1,]
    train2 <- dat[dat$.folds != 2,] 
    train3 <- dat[dat$.folds != 3,] 
    train4 <- dat[dat$.folds != 4,] 
    train5 <- dat[dat$.folds != 5,] 
    
    training_list <- list(train1, train2, train3, train4, train5)
    
    return(training_list)
  }
  #'  Create a list of lists of partitioned data sets using standardized data
  #'  5 training data sets for each of 3 season-specific data sets per species
  mdData_smr_train <- lapply(mdDataz_smr, training_dat)
  mdData_wtr_train <- lapply(mdDataz_wtr, training_dat)
  elkData_smr_train <- lapply(elkDataz_smr, training_dat)
  elkData_wtr_train <- lapply(elkDataz_wtr, training_dat)
  wtdData_smr_train <- lapply(wtdDataz_smr, training_dat)
  wtdData_wtr_train <- lapply(wtdDataz_wtr, training_dat)
  cougData_smr_train <- lapply(cougDataz_smr, training_dat)
  cougData_wtr_train <- lapply(cougDataz_wtr, training_dat)
  wolfData_smr_train <- lapply(wolfDataz_smr, training_dat)
  wolfData_wtr_train <- lapply(wolfDataz_wtr, training_dat)
  bobData_smr_train <- lapply(bobDataz_smr, training_dat)
  bobData_wtr_train <- lapply(bobDataz_wtr, training_dat)
  coyData_smr_train <- lapply(coyDataz_smr, training_dat)
  coyData_wtr_train <- lapply(coyDataz_wtr, training_dat)
  
  #'  Double checked this split correctly
  unique(mdData_smr_train[[1]][[5]][".folds"]) # should have folds 1 - 4, not 5
  unique(mdData_smr_train[[1]][[2]][".folds"]) # should have folds 1 & 3-5, not 2
  
  
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
  #'  Create a list of lists of the withheld fold of data 
  #'  5 testing data sets for each of 3 season-specific data sets per species
  mdData_smr_test <- lapply(mdDataz_smr, testing_dat)
  mdData_wtr_test <- lapply(mdDataz_wtr, testing_dat)
  elkData_smr_test <- lapply(elkDataz_smr, testing_dat)
  elkData_wtr_test <- lapply(elkDataz_wtr, testing_dat)
  wtdData_smr_test <- lapply(wtdDataz_smr, testing_dat)
  wtdData_wtr_test <- lapply(wtdDataz_wtr, testing_dat)
  cougData_smr_test <- lapply(cougDataz_smr, testing_dat)
  cougData_wtr_test <- lapply(cougDataz_wtr, testing_dat)
  wolfData_smr_test <- lapply(wolfDataz_smr, testing_dat)
  wolfData_wtr_test <- lapply(wolfDataz_wtr, testing_dat)
  bobData_smr_test <- lapply(bobDataz_smr, testing_dat)
  bobData_wtr_test <- lapply(bobDataz_wtr, testing_dat)
  coyData_smr_test <- lapply(coyDataz_smr, testing_dat)
  coyData_wtr_test <- lapply(coyDataz_wtr, testing_dat)

  #'  Double checked this split correctly
  unique(mdData_smr_test[[1]][[5]][".folds"]) # should only be fold 5
  unique(mdData_smr_test[[1]][[2]][".folds"]) # should only be fold 2
  
  #'  Save for later use
  save(mdData_smr_test, "./Outputs/RSF_output/Kfold_CV/mdData_smr_TestData.RData")
  save(mdData_wtr_test, "./Outputs/RSF_output/Kfold_CV/mdData_wtr_TestData.RData")
  save(elkData_smr_test, "./Outputs/RSF_output/Kfold_CV/elkData_smr_TestData.RData")
  save(elkData_wtr_test, "./Outputs/RSF_output/Kfold_CV/elkData_wtr_TestData.RData")
  save(wtdData_smr_test, "./Outputs/RSF_output/Kfold_CV/wtdData_smr_TestData.RData")
  save(wtdData_wtr_test, "./Outputs/RSF_output/Kfold_CV/wtdData_wtr_TestData.RData")
  save(cougData_smr_test, "./Outputs/RSF_output/Kfold_CV/cougData_smr_TestData.RData")
  save(cougData_wtr_test, "./Outputs/RSF_output/Kfold_CV/cougData_wtr_TestData.RData")
  save(wolfData_smr_test, "./Outputs/RSF_output/Kfold_CV/wolfData_smr_TestData.RData")
  save(wolfData_wtr_test, "./Outputs/RSF_output/Kfold_CV/wolfData_wtr_TestData.RData")
  save(bobData_smr_test, "./Outputs/RSF_output/Kfold_CV/bobData_smr_TestData.RData")
  save(bobData_wtr_test, "./Outputs/RSF_output/Kfold_CV/bobData_wtr_TestData.RData")
  save(coyData_smr_test, "./Outputs/RSF_output/Kfold_CV/coyData_smr_TestData.RData")
  save(coyData_wtr_test, "./Outputs/RSF_output/Kfold_CV/coyData_wtr_TestData.RData")
  
  
  #'  Define species- and season-specific models for RSFs to cross-validate
  #'  Used backwards step selection to select covariates for each model
  ####  Mule Deer models  ####
  md_smr18_mod <- "Used ~ 1 + Elev + I(Elev^2) + RoadDen + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  md_smr19_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Landcover_type + (1|ID)"
  md_smr20_mod <- "Used ~ 1 + Elev + I(Elev^2) + RoadDen + Dist2Water + CanopyCover + Landcover_type + (1|ID)"
  md_wtr1819_mod <- "Used ~ 1 + Elev + Slope + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  md_wtr1920_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  md_wtr2021_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  ####  Elk models  ####
  elk_smr18_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  elk_smr19_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  elk_smr20_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  elk_wtr1819_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  elk_wtr1920_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  elk_wtr2021_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  ####  White-tailed Deer models  ####
  wtd_smr18_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Edge + Landcover_type + (1|ID)"
  wtd_smr19_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  wtd_smr20_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  wtd_wtr1819_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + Dist2Edge + Landcover_type + (1|ID)"
  wtd_wtr1920_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  wtd_wtr2021_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  ####  Cougar models  ####
  coug_smr18_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)" 
  coug_smr19_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  coug_smr20_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  coug_wtr1819_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  coug_wtr1920_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  coug_wtr2021_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  ####  Wolf models  ####
  wolf_smr18_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + Dist2Water + Dist2Edge + Landcover_type + (1|ID)"
  wolf_smr19_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  wolf_smr20_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  wolf_wtr1819_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)" 
  wolf_wtr1920_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + HumanMod + Dist2Edge + Landcover_type + (1|ID)"
  wolf_wtr2021_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  ####  Bobcat models  ####
  #'  Only data for MVBOB90M in smr18--- not enough data to make inference about bobcat resource selection across 2 study areas
  # bob_smr18_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  bob_smr19_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + HumanMod + Dist2Edge + Landcover_type + (1|ID)"
  bob_smr20_mod <- "Used ~ 1 + RoadDen + Dist2Water + HumanMod + Landcover_type + (1|ID)"
  #'  Only data for MVBOB88M & MVBOB90M wtr1819--- not enough data to make inference about bobcat resource selection across 2 study areas
  # bob_wtr1819_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  bob_wtr1920_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  bob_wtr2021_mod <- "Used ~ 1 + Elev + Slope + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  ####  Coyote models  ####
  #'  Data from only MVCOY68F, NECOY1F, NECOY2M, & NECOY3F in snmr18--- hesitant to extrapolate selection across study areas
  coy_smr18_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"  
  coy_smr19_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Landcover_type + (1|ID)"
  coy_smr20_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
  #'  Data from only MVCOY68F, NECOY1F, NECOY2M, NECOY3F & NECOY4M in wtr1819--- hesitant to extrapolate selection across study areas
  coy_wtr1819_mod <- "Used ~ 1 + Elev + Dist2Water + Dist2Edge + Landcover_type + (1|ID)"
  coy_wtr1920_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + Landcover_type + (1|ID)"
  coy_wtr2021_mod <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)"
 
  
  ####  K-fold model training  ####
  #'  =============================
  #'  Function to run GLMM on K folded data 
  glmm_fn <- function(dat, mod) {
    glmm_mod <- glmer(formula = mod, data = dat, family = binomial(link = "logit"))
    print(summary(glmm_mod))
    print(car::vif(glmm_mod))
      
    return(glmm_mod)
  }
  #'  Apply spp/season/year-specific glmm to list of k-folded data
  #'  Each list contains 5 training data sets per year [[1]], year [[2]], year [[3]]
  ####  Mule Deer K-fold model training  ####
  md_smr18 <- lapply(mdData_smr_train[[1]], glmm_fn, mod = md_smr18_mod)
  md_smr19 <- lapply(mdData_smr_train[[2]], glmm_fn, mod = md_smr19_mod)
  md_smr20 <- lapply(mdData_smr_train[[3]], glmm_fn, mod = md_smr20_mod)
  md_wtr1819 <- lapply(mdData_wtr_train[[1]], glmm_fn, mod = md_wtr1819_mod)
  md_wtr1920 <- lapply(mdData_wtr_train[[2]], glmm_fn, mod = md_wtr1920_mod)
  md_wtr2021 <- lapply(mdData_wtr_train[[3]], glmm_fn, mod = md_wtr2021_mod)
  #'  Save trained models
  md_kfold_smr <- list(md_smr18, md_smr19, md_smr20)
  md_kfold_wtr <- list(md_wtr1819, md_wtr1920, md_wtr2021)
  save(md_kfold_smr, file = paste0("./Outputs/RSF_output/Kfold_CV/md_kfold_smr_", Sys.Date(), ".RData"))
  save(md_kfold_wtr, file = paste0("./Outputs/RSF_output/Kfold_CV/md_kfold_wtr_", Sys.Date(), ".RData"))
  
  ####  Elk K-fold model training  ####
  elk_smr18 <- lapply(elkData_smr_train[[1]], glmm_fn, mod = elk_smr18_mod)
  elk_smr19 <- lapply(elkData_smr_train[[2]], glmm_fn, mod = elk_smr19_mod)
  elk_smr20 <- lapply(elkData_smr_train[[3]], glmm_fn, mod = elk_smr20_mod)
  elk_wtr1819 <- lapply(elkData_wtr_train[[1]], glmm_fn, mod = elk_wtr1819_mod)
  elk_wtr1920 <- lapply(elkData_wtr_train[[2]], glmm_fn, mod = elk_wtr1920_mod)
  elk_wtr2021 <- lapply(elkData_wtr_train[[3]], glmm_fn, mod = elk_wtr2021_mod)
  #'  Save trained models
  elk_kfold_smr <- list(elk_smr18, elk_smr19, elk_smr20)
  elk_kfold_wtr <- list(elk_wtr1819, elk_wtr1920, elk_wtr2021)
  save(elk_kfold_smr, file = paste0("./Outputs/RSF_output/Kfold_CV/elk_kfold_smr_", Sys.Date(), ".RData"))
  save(elk_kfold_wtr, file = paste0("./Outputs/RSF_output/Kfold_CV/elk_kfold_wtr_", Sys.Date(), ".RData"))
  
  ####  White-tailed Deer K-fold model training  ####
  wtd_smr18 <- lapply(wtdData_smr_train[[1]], glmm_fn, mod = wtd_smr18_mod)
  wtd_smr19 <- lapply(wtdData_smr_train[[2]], glmm_fn, mod = wtd_smr19_mod)
  wtd_smr20 <- lapply(wtdData_smr_train[[3]], glmm_fn, mod = wtd_smr20_mod)
  wtd_wtr1819 <- lapply(wtdData_wtr_train[[1]], glmm_fn, mod = wtd_wtr1819_mod)
  wtd_wtr1920 <- lapply(wtdData_wtr_train[[2]], glmm_fn, mod = wtd_wtr1920_mod)
  wtd_wtr2021 <- lapply(wtdData_wtr_train[[3]], glmm_fn, mod = wtd_wtr2021_mod)
  #'  Save trained models
  wtd_kfold_smr <- list(wtd_smr18, wtd_smr19, wtd_smr20)
  wtd_kfold_wtr <- list(wtd_wtr1819, wtd_wtr1920, wtd_wtr2021)
  save(wtd_kfold_smr, file = paste0("./Outputs/RSF_output/Kfold_CV/wtd_kfold_smr_", Sys.Date(), ".RData"))
  save(wtd_kfold_wtr, file = paste0("./Outputs/RSF_output/Kfold_CV/wtd_kfold_wtr_", Sys.Date(), ".RData"))
  
  ####  Cougar K-fold model training  ####
  coug_smr18 <- lapply(cougData_smr_train[[1]], glmm_fn, mod = coug_smr18_mod)
  coug_smr19 <- lapply(cougData_smr_train[[2]], glmm_fn, mod = coug_smr19_mod)
  coug_smr20 <- lapply(cougData_smr_train[[3]], glmm_fn, mod = coug_smr20_mod)
  coug_wtr1819 <- lapply(cougData_wtr_train[[1]], glmm_fn, mod = coug_wtr1819_mod)
  coug_wtr1920 <- lapply(cougData_wtr_train[[2]], glmm_fn, mod = coug_wtr1920_mod)
  coug_wtr2021 <- lapply(cougData_wtr_train[[3]], glmm_fn, mod = coug_wtr2021_mod)
  #'  Save trained models
  coug_kfold_smr <- list(coug_smr18, coug_smr19, coug_smr20)
  coug_kfold_wtr <- list(coug_wtr1819, coug_wtr1920, coug_wtr2021)
  save(coug_kfold_smr, file = paste0("./Outputs/RSF_output/Kfold_CV/coug_kfold_smr_", Sys.Date(), ".RData"))
  save(coug_kfold_wtr, file = paste0("./Outputs/RSF_output/Kfold_CV/coug_kfold_wtr_", Sys.Date(), ".RData"))
  
  ####  Wolf K-fold model training  ####
  wolf_smr18 <- lapply(wolfData_smr_train[[1]], glmm_fn, mod = wolf_smr18_mod)
  wolf_smr19 <- lapply(wolfData_smr_train[[2]], glmm_fn, mod = wolf_smr19_mod)
  wolf_smr20 <- lapply(wolfData_smr_train[[3]], glmm_fn, mod = wolf_smr20_mod)
  wolf_wtr1819 <- lapply(wolfData_wtr_train[[1]], glmm_fn, mod = wolf_wtr1819_mod)
  wolf_wtr1920 <- lapply(wolfData_wtr_train[[2]], glmm_fn, mod = wolf_wtr1920_mod)
  wolf_wtr2021 <- lapply(wolfData_wtr_train[[3]], glmm_fn, mod = wolf_wtr2021_mod)
  #'  Save trained models
  wolf_kfold_smr <- list(wolf_smr18, wolf_smr19, wolf_smr20)
  wolf_kfold_wtr <- list(wolf_wtr1819, wolf_wtr1920, wolf_wtr2021)
  save(wolf_kfold_smr, file = paste0("./Outputs/RSF_output/Kfold_CV/wolf_kfold_smr_", Sys.Date(), ".RData"))
  save(wolf_kfold_wtr, file = paste0("./Outputs/RSF_output/Kfold_CV/wolf_kfold_wtr_", Sys.Date(), ".RData"))
  
  ####  Bobcat K-fold model training  ####
  #'  Note: smr18 & wtr1819 not run due to too few collars in year 1
  # bob_smr18 <- lapply(bobData_smr_train[[1]], glmm_fn, mod = bob_smr18_mod)
  bob_smr19 <- lapply(bobData_smr_train[[2]], glmm_fn, mod = bob_smr19_mod)
  bob_smr20 <- lapply(bobData_smr_train[[3]], glmm_fn, mod = bob_smr20_mod)
  # bob_wtr1819 <- lapply(bobData_wtr_train[[1]], glmm_fn, mod = bob_wtr1819_mod)
  bob_wtr1920 <- lapply(bobData_wtr_train[[2]], glmm_fn, mod = bob_wtr1920_mod)
  bob_wtr2021 <- lapply(bobData_wtr_train[[3]], glmm_fn, mod = bob_wtr2021_mod)
  #'  Save trained models
  bob_kfold_smr <- list(bob_smr19, bob_smr20)
  bob_kfold_wtr <- list(bob_wtr1920, bob_wtr2021)
  save(bob_kfold_smr, file = paste0("./Outputs/RSF_output/Kfold_CV/bob_kfold_smr_", Sys.Date(), ".RData"))
  save(bob_kfold_wtr, file = paste0("./Outputs/RSF_output/Kfold_CV/bob_kfold_wtr_", Sys.Date(), ".RData"))
  
  ####  Coyote K-fold model training  ####
  coy_smr18 <- lapply(coyData_smr_train[[1]], glmm_fn, mod = coy_smr18_mod)
  coy_smr19 <- lapply(coyData_smr_train[[2]], glmm_fn, mod = coy_smr19_mod)
  coy_smr20 <- lapply(coyData_smr_train[[3]], glmm_fn, mod = coy_smr20_mod)
  coy_wtr1819 <- lapply(coyData_wtr_train[[1]], glmm_fn, mod = coy_wtr1819_mod)
  coy_wtr1920 <- lapply(coyData_wtr_train[[2]], glmm_fn, mod = coy_wtr1920_mod)
  coy_wtr2021 <- lapply(coyData_wtr_train[[3]], glmm_fn, mod = coy_wtr2021_mod)
  #'  Save trained models
  coy_kfold_smr <- list(coy_smr18, coy_smr19, coy_smr20)
  coy_kfold_wtr <- list(coy_wtr1819, coy_wtr1920, coy_wtr2021)
  save(coy_kfold_smr, file = paste0("./Outputs/RSF_output/Kfold_CV/coy_kfold_smr_", Sys.Date(), ".RData"))
  save(coy_kfold_wtr, file = paste0("./Outputs/RSF_output/Kfold_CV/coy_kfold_wtr_", Sys.Date(), ".RData"))
  
  
  ####  Projecting trained models  ####
  #'  =================================
  #'  Read in saved k-fold trained model results
  #'  Extract betas for each trained model
  #'  Fit study area-wide covs to GLMM with estimated betas
  #'  Plot predicted results across study areas (using 1km grid)
  #'  Bin RSF predictions (area-weighted or quantiles)
  #'  Save binned predictions as rasters
  
  
  ####  Cross-validate RSFs  ####
  #'  ===========================
  #'  Read in saved testing data
  #'  Plot used testing observations across binned RSF predictions for respective
  #'  trained models (e.g., test1 against train1, test2 against train2, etc.)
  #'  Count frequency of used locations in each bin
  #'  Spearman's rank correlation for each k-fold model
  #'  Estimate mean rank correlation as measure of predictive ability for each 
  #'  spp/season/year-specific model
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  