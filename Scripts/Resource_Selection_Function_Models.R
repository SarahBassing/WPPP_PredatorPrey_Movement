  #'  ============================================
  #'  3rd Order Resource Selection Functions (RSFs)
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  December 2021
  #'  ============================================
  #'  Script to run 3rd order resource selection functions (RSFs) and predict
  #'  relative probability of selection across study areas for each species and
  #'  season. These predictive surfaces will be used to represent probability of
  #'  predator/prey presence as a covariate in HMMs.
  #'  
  #'  Using backwards step selection for initial model selection, then K-fold
  #'  cross-validation to evaluate each model's predictive performance 
  #'  (see K-fold_CV_for_RSF.R for K-fold CV analysis).
  #'  
  #'  
  #'  Cleaned telemetry and covariate data were prepared with
  #'  Collar_RSF_DataPrep.R script. 
  #'  ============================================
  
  #'  Clear memory
  rm(list=ls())

  #'  Load packages
  library(tidyverse)
  library(car)
  library(cvms)
  library(groupdata2)
  library(knitr)
  library(lme4)
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
  
  
  #'  Function to re-classify landcover into fewer categories
  #'  Based on T. Ganz's input:
  #'  Open grass: mesic grass, xeric grass, wetland woody
  #'  Shrubby mix: mesic shrub, xeric shrub
  #'  Other: water, barren, glacier
  #'  Developed: agriculture, commercial, developed
  #'  Forest
  #'  Wetland
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
  
  #'  Function to reclassify land cover into fewer categories
  #'  Landcover_type categories causing convergence issues for some species due to
  #'  too few observations in some categories (e.g., "Other", "Developed") 
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
  
  
  #'  Center & scale covariates 
  #'  Note: standardizing across all IDs but separately by species & season
  standardize_covs <- function(locs){
    #'  Make categorical variables factors
    locs$ID <- as.factor(locs$ID)
    locs$Used <- as.factor(locs$Used)
    locs$Area <- as.factor(locs$Area)
    locs$Year <- as.factor(locs$Year)
    locs$Season <- as.factor(locs$Season)
    locs$Landcover <- as.factor(locs$Landcover)
    locs$Landcover_type <- droplevels(as.factor(locs$Landcover_type))
    locs$Landcover_type <- relevel(locs$Landcover_type, ref = "Forest")
    #'  Standardize continuous variables
    locs$Elev <- scale(locs$Elev)
    locs$Slope <- scale(locs$Slope)
    locs$TPI <- scale(locs$TPI)
    locs$RoadDen <- scale(locs$RoadDen)
    locs$Dist2Water <- scale(locs$Dist2Water)
    locs$HumanMod <- scale(locs$HumanMod)
    locs$CanopyCover <- scale(locs$CanopyCover)
    locs$Dist2Edge <- scale(locs$Dist2Edge)
    locs$PercForMix <- scale(locs$PercForMix)
    locs$PercXGrass <- scale(locs$PercXGrass)
    locs$PercXShrub <- scale(locs$PercXShrub)
    
    locs <- as.data.frame(locs)
    
    return(locs)
  }
  #'  List datasets by season & standardize covariates
  mdData_smr <- md_dat_all[md_dat_all$Season == "Summer18" | md_dat_all$Season == "Summer19" | md_dat_all$Season == "Summer20",]
  mdData_wtr <- md_dat_all[md_dat_all$Season == "Winter1819" | md_dat_all$Season == "Winter1920" | md_dat_all$Season == "Winter2021",]
  mdDataz_smr <- standardize_covs(mdData_smr)
  mdDataz_wtr <- standardize_covs(mdData_wtr)
  #'  Note the reclassified landcover_type data for elkData_winter
  elkData_smr <- elk_dat_all[elk_dat_all$Season == "Summer18" | elk_dat_all$Season == "Summer19" | elk_dat_all$Season == "Summer20",]
  elkData_wtr <- elk_dat_all_reclass[elk_dat_all_reclass$Season == "Winter1819" | elk_dat_all_reclass$Season == "Winter1920" | elk_dat_all_reclass$Season == "Winter2021",]
  elkDataz_smr <- standardize_covs(elkData_smr)
  elkDataz_wtr <- standardize_covs(elkData_wtr)
  #'  Note the reclassified landcover_type data for wtdData_winter 
  wtdData_smr <- wtd_dat_all[wtd_dat_all$Season == "Summer18" | wtd_dat_all$Season == "Summer19" | wtd_dat_all$Season == "Summer20",]
  wtdData_wtr <- wtd_dat_all_reclass[wtd_dat_all_reclass$Season == "Winter1819" | wtd_dat_all_reclass$Season == "Winter1920" | wtd_dat_all_reclass$Season == "Winter2021",]
  wtdDataz_smr <- standardize_covs(wtdData_smr)
  wtdDataz_wtr <- standardize_covs(wtdData_wtr)
  #'  Note the reclassified landcover_type data for cougData_winter
  cougData_smr <- coug_dat_all[coug_dat_all$Season == "Summer18" | coug_dat_all$Season == "Summer19" | coug_dat_all$Season == "Summer20",]
  cougData_wtr <- coug_dat_all_reclass[coug_dat_all_reclass$Season == "Winter1819" | coug_dat_all_reclass$Season == "Winter1920" | coug_dat_all_reclass$Season == "Winter2021",]
  cougDataz_smr <- standardize_covs(cougData_smr)
  cougDataz_wtr <- standardize_covs(cougData_wtr)
  #'  Note the double reclassified landcover_type data for wolfData
  wolfData_smr <- wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Summer18" | wolf_dat_all_reclass2$Season == "Summer19" | wolf_dat_all_reclass2$Season == "Summer20",]
  wolfData_wtr <- wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Winter1819" | wolf_dat_all_reclass2$Season == "Winter1920" | wolf_dat_all_reclass2$Season == "Winter2021",]
  wolfDataz_smr <- standardize_covs(wolfData_smr)
  wolfDataz_wtr <- standardize_covs(wolfData_wtr)
  #'  Note the reclassified landcover_type data for bobData
  bobData_smr <- bob_dat_all_reclass[bob_dat_all_reclass$Season == "Summer18" | bob_dat_all_reclass$Season == "Summer19" | bob_dat_all_reclass$Season == "Summer20",]
  bobData_wtr <- bob_dat_all_reclass[bob_dat_all_reclass$Season == "Winter1819" | bob_dat_all_reclass$Season == "Winter1920" | bob_dat_all_reclass$Season == "Winter2021",]
  bobDataz_smr <- standardize_covs(bobData_smr)
  bobDataz_wtr <- standardize_covs(bobData_wtr)
  #'  Note the reclassified landcover_type data for coyData
  coyData_smr <- coy_dat_all_reclass[coy_dat_all_reclass$Season == "Summer18" | coy_dat_all_reclass$Season == "Summer19" | coy_dat_all_reclass$Season == "Summer20",]
  coyData_wtr <- coy_dat_all_reclass[coy_dat_all_reclass$Season == "Winter1819" | coy_dat_all_reclass$Season == "Winter1920" | coy_dat_all_reclass$Season == "Winter2021",]
  coyDataz_smr <- standardize_covs(coyData_smr)
  coyDataz_wtr <- standardize_covs(coyData_wtr)

  #'  Correlation Matrix
  #'  ==================
  #'  Function to create correlation matrix for all continuous covariates at once
  cov_correlation <- function(dat) {
    used <- dat[dat$Used == 1,]
    covs <- used[,c("Elev", "Slope", "TPI", "RoadDen",
                    "Dist2Water", "HumanMod", "CanopyCover", "Dist2Edge",
                    "PercForMix", "PercXGrass", "PercXShrub")]
    cor_matrix <- cor(covs, use = "complete.obs")
    return(cor_matrix)
  }
  #'  Generate correlation matrix for each species and season
  (md_smr_corr <- cov_correlation(mdData_smr))
  (md_wtr_corr <- cov_correlation(mdData_wtr))
  (elk_smr_corr <- cov_correlation(elkData_smr))
  (elk_wtr_corr <- cov_correlation(elkData_wtr))
  (wtd_smr_corr <- cov_correlation(wtdData_smr))
  (wtd_wtr_corr <- cov_correlation(wtdData_wtr))
  (coug_smr_corr <- cov_correlation(cougData_smr))
  (coug_wtr_corr <- cov_correlation(cougData_wtr))
  (wolf_smr_corr <- cov_correlation(wolfData_smr))
  (wolf_wtr_corr <- cov_correlation(wolfData_wtr))
  (bob_smr_corr <- cov_correlation(bobData_smr))
  (bob_wtr_corr <- cov_correlation(bobData_wtr))
  (coy_smr_corr <- cov_correlation(coyData_smr))
  (coy_wtr_corr <- cov_correlation(coyData_wtr))

  #'  Elevation & TPI are highly correlated (almost 100%) for all datasets so nixing TPI entirely
  #'  Elevation & Human Modified correlated in MD smr/wtr, ELK smr, COUG smr, & COY smr/wtr- nixing HumanMod for those models
  #'  Using Landcover_type instead of % Forest, % Grass, & % Shrub because...
  #'  % Shrub correlated with Elevation, TPI and Human Modified in MD smr
  #'  Dist2Edg correlated with % Forest & % Grass in ELK wtr
  #'  % Forest & Grass correlated in COUG and WOLF wtr
  #'  % Grass & Shrub correlated for BOb smr
  #'  % Forest & Canopy Cover, % Forest & Shrub correlated for BOB wtr
  
  
  
  #'  Resource Selection Function Models
  #'  ==================================
  #'  Functions to run logistic mixed effects models that include random effect 
  #'  for individual. Habitat covariates excluded based on species (see notes 
  #'  above) and convergence issues. Seasonal models run separately so RSFs are 
  #'  specific to the species & season but with all years combined. Some covariates
  #'  vary annually however so each species & season-specific RSF is predicted 
  #'  across an annual study area map (under prediction section below). 
  #'  
  #'  Using backwards step selection to identify top model with all significant
  #'  covariates. Covariates excluded through this process are commented out next
  #'  to model in order that the variable was removed.
  #'  ==================================
  
  glmm_fn <- function(mod, dat) {
    glmm_mod <- glmer(formula = mod, data = dat, family = binomial(link = "logit"))
    print(summary(glmm_mod))
    print(car::vif(glmm_mod))
    
    return(glmm_mod)
  }
  #'  Run species, season, and year specific models through glmm function
  
  ####  Mule Deer RSFs  ####
  #'  Dropping HumanMod in mulie models due to high correlation with other covariates
  md_smr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + Dist2Edge + Landcover_type + (1|ID)", dat = mdDataz_smr) # + CanopyCover
  md_wtr <- glmm_fn(mod = "Used ~ 1 + Elev + Slope + RoadDen + Dist2Water + CanopyCover + Landcover_type + (1|ID)", dat = mdDataz_wtr) # + I(Elev^2) + Dist2Edge

  ####  Elk RSFs  ####
  #'  Dropping HumanMod in elk summer models due to high correlation with other covariates
  elk_smr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = elkDataz_smr)
  #'  Note: using reclassified version of landcover for winter elk models
  elk_wtr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = elkDataz_wtr)

  ####  White-tailed Deer RSFs  ####
  wtd_smr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + Dist2Edge + Landcover_type + (1|ID)",  dat = wtdDataz_smr) # + CanopyCover
  #'  Note: using reclassified version of landcover for winter WTD models
  wtd_wtr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wtdDataz_wtr)

  ####  Cougar RSFs  ####
  #'  Dropping HumanMod in cougar summer models due to high correlation with other covariates
  coug_smr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = cougDataz_smr) 
  #'  Note: using reclassified version of landcover for winter cougar models
  coug_wtr <- glmm_fn(mod = "Used ~ 1 + Elev + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Landcover_type + (1|ID)",  dat = cougDataz_wtr) # + I(Elev^2) + Dist2Edge

  ####  Wolf RSFs  ####
  #'  NOTE: using 2nd reclassified version of landcover categories for wolf models
  wolf_smr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + Dist2Water + HumanMod + Dist2Edge + Landcover_type + (1|ID)",  dat = wolfDataz_smr) # + RoadDen + CanopyCover
  wolf_wtr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wolfDataz_wtr)

  ####  Bobcat RSFs  ####
  #'  Note: using reclassified version of landcover for summer bobcat models
  bob_smr <- glmm_fn(mod = "Used ~ 1 + Elev + Slope + RoadDen + Dist2Water + HumanMod + Dist2Edge + Landcover_type + (1|ID)",  dat = bobDataz_smr) # + I(Elev^2) + CanopyCover
  bob_wtr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + HumanMod + CanopyCover + Landcover_type + (1|ID)",  dat = bobDataz_wtr) # + Dist2Water + Dist2Edge

  ####  Coyote RSFs  ####
  #'  Note: using reclassified version of landcover for all coyote models
  #'  Dropping HumanMod from all coyote models due to high correlation with other covariates
  coy_smr <- glmm_fn(mod = "Used ~ 1 + Slope + RoadDen + Dist2Water + CanopyCover + Landcover_type + (1|ID)", dat = coyDataz_smr)  # + I(Elev^2) + Dist2Edge + Elev
  coy_wtr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + Dist2Edge + Landcover_type + (1|ID)", dat = coyDataz_wtr) # + CanopyCover

  #'  Group species-specific models
  RSF_MD_list <- list(md_smr, md_wtr)
  RSF_ELK_list <- list(elk_smr, elk_wtr)
  RSF_WTD_list <- list(wtd_smr, wtd_wtr)
  RSF_COUG_list <- list(coug_smr, coug_wtr)
  RSF_WOLF_list <- list(wolf_smr, wolf_wtr)
  RSF_BOB_list <- list(bob_smr, bob_wtr)
  RSF_COY_list <- list(coy_smr, coy_wtr)
  
  #'  Save
  save(RSF_MD_list, file = paste0("./Outputs/RSF_output/RSF_MD_list_", Sys.Date(), ".RData"))
  save(RSF_ELK_list, file = paste0("./Outputs/RSF_output/RSF_ELK_list_", Sys.Date(), ".RData"))
  save(RSF_WTD_list, file = paste0("./Outputs/RSF_output/RSF_WTD_list_", Sys.Date(), ".RData"))
  save(RSF_COUG_list, file = paste0("./Outputs/RSF_output/RSF_COUG_list_", Sys.Date(), ".RData"))
  save(RSF_WOLF_list, file = paste0("./Outputs/RSF_output/RSF_WOLF_list_", Sys.Date(), ".RData"))
  save(RSF_BOB_list, file = paste0("./Outputs/RSF_output/RSF_BOB_list_", Sys.Date(), ".RData"))
  save(RSF_COY_list, file = paste0("./Outputs/RSF_output/RSF_COY_list_", Sys.Date(), ".RData"))
  
  

  ####  Project RSF results across study areas  ####
  #'  ==============================================
  
  #'  Load RSFs
  load("./Outputs/RSF_output/RSF_MD_list_2022-01-13.RData")
  load("./Outputs/RSF_output/RSF_ELK_list_2022-01-13.RData")
  load("./Outputs/RSF_output/RSF_WTD_list_2022-01-13.RData")
  load("./Outputs/RSF_output/RSF_COUG_list_2022-01-13.RData")
  load("./Outputs/RSF_output/RSF_WOLF_list_2022-01-13.RData")
  load("./Outputs/RSF_output/RSF_BOB_list_2022-01-13.RData")
  load("./Outputs/RSF_output/RSF_COY_list_2022-01-13.RData")
  
  #'  Load spatial libraries
  library(sf)
  library(raster)
  
  #'  Define desired projections
  sa_proj <- projection("EPSG:2855")  # NAD83(HARN) / Washington North
  
  #'  Read in study area grids
  NE_1km <- raster("./Shapefiles/NE_1km_grid_mask.tif")
  OK_1km <- raster("./Shapefiles/OK_1km_grid_mask.tif")
  # NE_30m <- raster("./Shapefiles/NE_30m_grid_mask.tif")
  # OK_30m <- raster("./Shapefiles/OK_30m_grid_mask.tif")
  
  #'  Load study area shapefiles
  OK.SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "METHOW_SA") %>%
    st_transform(crs = sa_proj)
  OK.SA <- as(OK.SA, "Spatial")
  NE.SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "NE_SA") %>%
    st_transform(crs = sa_proj)
  NE.SA <- as(NE.SA, "Spatial")
  
  #'  Convert rasters to pixels and extract coordinates (centroid of each cell)
  #'  FYI: grid.index is the original cell number but b/c some cells were masked 
  #'  out for large water bodies, the gridID does not match the extracted study 
  #'  area-wide covariate df so need to create a new ID specific to masked grid
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
  # NE_pts <- raster_dat(NE_30m)
  # OK_pts <- raster_dat(OK_30m)
  
  #'  Read in covariates extracted across each study area 
  load("./Outputs/Telemetry_covs/NE_covs_1km_2022-01-13.RData") 
  load("./Outputs/Telemetry_covs/OK_covs_1km_2022-01-13.RData")
  # load("./Outputs/Telemetry_covs/NE_covs_30m_2022-01-13.RData") 
  # load("./Outputs/Telemetry_covs/OK_covs_30m_2022-01-13.RData")
  
  #'  Format study area-wide covariate data to include annually relevant data only
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
  
  
  #'  Call landcover and scaling functions from above to re-format covariates
  SA.covs_list <- lapply(SA.covs_list, class_landcov)
  SA.covs_list_reclass <- lapply(SA.covs_list, reclass_landcov)
  SA.covs_list_wolf_reclass <- lapply(SA.covs_list_reclass, reclass_wolf)
  #'  Reformat study area specific covariate dfs for ungulate models
  NE.covs_list <- lapply(NE.covs_list, class_landcov)
  OK.covs_list <- lapply(OK.covs_list, class_landcov)
  NE.covs_list_reclass <- lapply(NE.covs_list, reclass_landcov)
  OK.covs_list_reclass <- lapply(OK.covs_list, reclass_landcov)
  

  #'  Function to find mean & standard deviation for raw covariates in RSFs
  #'  Necessary for standardizing study area-wide covs based on original models
  #'  Note: summarizes data by spp & season, same as data structure in RSF
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
  #'  the original covariate values that went into the RSFs
  scaling_covs <- function(covs, mu.sd) {
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
  #'  Standardize study area-wide covariates based on species & season-specific
  #'  model covariate means & SDs
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

  
  #'  Function to save parameter estimates from each RSF
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
             Parameter = ifelse(Parameter == "Landcover_typeOpen Grass", "b.grass", Parameter),
             Parameter = ifelse(Parameter == "Landcover_typeOther", "b.other", Parameter),
             Parameter = ifelse(Parameter == "Landcover_typeShrub Mix", "b.shrub", Parameter),
             Parameter = ifelse(Parameter == "Landcover_typeWetland", "b.wetland", Parameter)) %>%
      #'  Spread data so each row represents model coefficients for a single season, single species model
      pivot_wider(names_from = Parameter, values_from = Estimate) 
    
    #'  Covariates excluded from species-specific models not included in the data
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
  #'  Extract coefficient estimates for each model (list: summer rsf [[1]], winter rsf [[2]])
  md_smr_rsfout <- rsf_out(RSF_MD_list[[1]], spp = "Mule Deer", season = "Summer")
  md_wtr_rsfout <- rsf_out(RSF_MD_list[[2]], spp = "Mule Deer", season = "Winter")
  elk_smr_rsfout <- rsf_out(RSF_ELK_list[[1]], spp = "Elk", season = "Summer")
  elk_wtr_rsfout <- rsf_out(RSF_ELK_list[[2]], spp = "Elk", season = "Winter")
  wtd_smr_rsfout <- rsf_out(RSF_WTD_list[[1]], spp = "White-tailed Deer", season = "Summer")
  wtd_wtr_rsfout <- rsf_out(RSF_WTD_list[[2]], spp = "White-tailed Deer", season = "Winter")
  coug_smr_rsfout <- rsf_out(RSF_COUG_list[[1]], spp = "Cougar", season = "Summer")
  coug_wtr_rsfout <- rsf_out(RSF_COUG_list[[2]], spp = "Cougar", season = "Winter")
  wolf_smr_rsfout <- rsf_out(RSF_WOLF_list[[1]], spp = "Wolf", season = "Summer")
  wolf_wtr_rsfout <- rsf_out(RSF_WOLF_list[[2]], spp = "Wolf", season = "Winter")
  bob_smr_rsfout <- rsf_out(RSF_BOB_list[[1]], spp = "Bobcat", season = "Summer")
  bob_wtr_rsfout <- rsf_out(RSF_BOB_list[[2]], spp = "Bobcat", season = "Winter")
  coy_smr_rsfout <- rsf_out(RSF_COY_list[[1]], spp = "Coyote", season = "Summer")
  coy_wtr_rsfout <- rsf_out(RSF_COY_list[[2]], spp = "Coyote", season = "Winter")
  
  
  #'  Function to predict across all grid cells based on RSF results
  #'  Should end up with 1 predicted value per grid cell
  #'  NOTE: I want the predict relative probability of selection from RSF dropping 
  #'  the intercept from the model and just exponentiating the coeffs*covs (Fieberg et al. 2020)
  predict_rsf <- function(cov, coef) {
    predict_rsf <- c()
    #'  Predict across each grid cell
    for(i in 1:nrow(cov)) {
      predict_rsf[i] <- exp(coef$b.elev*cov$Elev[i] + coef$b.elev2*I(cov$Elev[i]^2) + 
                              coef$b.slope*cov$Slope[i] + coef$b.road*cov$RoadDen[i] +
                              coef$b.water*cov$Dist2Water[i] + coef$b.hm*cov$HumanMod[i] +
                              coef$b.canopy*cov$CanopyCover[i] + coef$b.edge*cov$Dist2Edge[i] +
                              coef$b.grass*cov$Landcover_Grass[i] + coef$b.other*cov$Landcover_Other[i] +
                              coef$b.shrub*cov$Landcover_Shrub[i] + coef$b.wetland*cov$Landcover_Wetland[i])
    }
    predict_rsf <- as.data.frame(predict_rsf)
    predict_rsf <- cbind(cov$ID, cov$x, cov$y, cov$StudyArea, predict_rsf)
    colnames(predict_rsf) <- c("ID", "x", "y", "StudyArea", "predict_rsf")
    
    return(predict_rsf)
  }
  #'  Predict species & season-specific RSFs for each year across the study areas
  #'  NOTE: Applying annually varying covariate data to the same species & season-
  #'  specific RSF model because I expect the general relationships between a 
  #'  species and the covariates to be the same across years but that the spatial
  #'  distribution of those resource units may change annually.
  #'  Remember- *zcovs is a list of standardized covariates, 1 data frame per year
  md_smr_rsf_sa <- lapply(md_smr_zcovs, predict_rsf, coef = md_smr_rsfout)
  md_wtr_rsf_sa <- lapply(md_wtr_zcovs, predict_rsf, coef = md_wtr_rsfout)
  elk_smr_rsf_sa <- lapply(elk_smr_zcovs, predict_rsf, coef = elk_smr_rsfout)
  elk_wtr_rsf_sa <- lapply(elk_wtr_zcovs, predict_rsf, coef = elk_wtr_rsfout)
  wtd_smr_rsf_sa <- lapply(wtd_smr_zcovs, predict_rsf, coef = wtd_smr_rsfout)
  wtd_wtr_rsf_sa <- lapply(wtd_wtr_zcovs, predict_rsf, coef = wtd_wtr_rsfout)
  coug_smr_rsf_sa <- lapply(coug_smr_zcovs, predict_rsf, coef = coug_smr_rsfout)
  coug_wtr_rsf_sa <- lapply(coug_wtr_zcovs, predict_rsf, coef = coug_wtr_rsfout)
  wolf_smr_rsf_sa <- lapply(wolf_smr_zcovs, predict_rsf, coef = wolf_smr_rsfout)
  wolf_wtr_rsf_sa <- lapply(wolf_wtr_zcovs, predict_rsf, coef = wolf_wtr_rsfout)
  bob_smr_rsf_sa <- lapply(bob_smr_zcovs, predict_rsf, coef = bob_smr_rsfout)
  bob_wtr_rsf_sa <- lapply(bob_wtr_zcovs, predict_rsf, coef = bob_wtr_rsfout)
  coy_smr_rsf_sa <- lapply(coy_smr_zcovs, predict_rsf, coef = coy_smr_rsfout)
  coy_wtr_rsf_sa <- lapply(coy_wtr_zcovs, predict_rsf, coef = coy_wtr_rsfout)
  
  chk <- coy_smr_rsf_sa[[1]]
  
  #'  List and save
  all_spp_RSF_predicted <- list(md_smr_rsf_sa, md_wtr_rsf_sa, elk_smr_rsf_sa, 
                                elk_wtr_rsf_sa, wtd_smr_rsf_sa, wtd_wtr_rsf_sa, 
                                coug_smr_rsf_sa, coug_wtr_rsf_sa, wolf_smr_rsf_sa, 
                                wolf_wtr_rsf_sa, bob_smr_rsf_sa, bob_wtr_rsf_sa, 
                                coy_smr_rsf_sa, coy_wtr_rsf_sa)
  # save(all_spp_RSF_predicted, file = paste0("./Outputs/RSF_output/all_spp_RSF_predicted_", Sys.Date(), ".RData"))
  coy_rsf_predicted <- list(coy_smr_rsf_sa, coy_wtr_rsf_sa)

  
  #'  Identify any outliers
  outliers <- function(predicted, title, covs_list) {
    #'  Summarize predicted values
    hist(predicted$predict_rsf, breaks = 100, main = title)
    boxplot(predicted$predict_rsf, main = title)
    #'  What value represents the 99th percentile in the predicted RSF values
    quant <- quantile(predicted$predict_rsf, c(0.99), na.rm = TRUE)
    #'  Print that value and maximum prediction
    print(quant); print(max(predicted$predict_rsf, na.rm = TRUE))
    #'  Identify the 1% most extreme values and set to 99th percentile value
    predicted <- predicted %>%
      mutate(outlier = ifelse(predict_rsf > quant, "outlier", "not_outlier"),
             adjusted_rsf = ifelse(outlier == "outlier", quant, predict_rsf))
    #'  How many predicted values are above the 99th percentile?
    outlier <- predicted[predicted$outlier == "outlier",]
    outlier <- filter(outlier, !is.na(outlier))
    print(nrow(outlier))
    
    #' Summarize covariates associated with extreme values
    bigvalues <- full_join(predicted, covs_list, by = c("ID", "StudyArea", "x", "y"))
    print(summary(bigvalues)) 
    hist(bigvalues$Elev, breaks = 25, main = "Frequency of Elevation values", xlab = "Standardize Elevation")
    hist(bigvalues$Elev[bigvalues$outlier == "outlier"], breaks = 25, main = "Frequency of Outlier Elevation values", xlab = "Standardize Elevation")
    hist(bigvalues$Slope, breaks = 25, main = "Frequency of Slope values", xlab = "Standardize Slope")
    hist(bigvalues$Slope[bigvalues$outlier == "outlier"], breaks = 25, main = "Frequency of Outlier Slope values", xlab = "Standardize Slope")
    hist(bigvalues$RoadDen, breaks = 25, main = "Frequency of Road Density values", xlab = "Standardize Road Density")
    hist(bigvalues$RoadDen[bigvalues$outlier == "outlier"], breaks = 25, main = "Frequency of Outlier Road Density values", xlab = "Standardize Road Density")
    hist(bigvalues$Dist2Water, breaks = 25, main = "Frequency of Distance to Water values", xlab = "Standardize Dist. to Water")
    hist(bigvalues$Dist2Water[bigvalues$outlier == "outlier"], breaks = 25, main = "Frequency of Outlier Distance to Water values", xlab = "Standardize Dist. to Water")
    hist(bigvalues$HumanMod, breaks = 25, main = "Frequency of Human Modified Landscape values", xlab = "Standardize Human Mod")
    hist(bigvalues$HumanMod[bigvalues$outlier == "outlier"], breaks = 25, main = "Frequency of Outlier Human Modified Landscape values", xlab = "Standardize Human Mod")
    hist(bigvalues$CanopyCover, breaks = 25, main = "Frequency of Canopy Cover values", xlab = "Standardize Canopy Cover")
    hist(bigvalues$CanopyCover[bigvalues$outlier == "outlier"], breaks = 25, main = "Frequency of Outlier Canopy Cover values", xlab = "Standardize Canopy Cover")
    hist(bigvalues$Dist2Edge, breaks = 25, main = "Frequency of Distance to Edge values", xlab = "Standardize Dist. to Edge")
    hist(bigvalues$Dist2Edge[bigvalues$outlier == "outlier"], breaks = 25, main = "Frequency of Outlier Distance to Edge values", xlab = "Standardize Dist. to Edge")
    
    #' plot(bigvalues$RoadDen, bigvalues$predict_rsf, main = "Correlation between Predicted RSF and Road Density", xlab = "Road Density", ylab = "Predicted RSF Values")
    #' #'  Simple regression between covaraites
    #' abline(lm(bigvalues$predict_rsf ~ bigvalues$RoadDen), col = "red")
    #' #'  Mark the 99th percentil value
    #' abline(h = quant, col = "blue")
    #' #'  Frequency of covariates associated with non-extreme & extreme RSF values
    #' hist(bigvalues$RoadDen, breaks = 10, main = "Frequency of Road Density values", xlab = "Road Density")
    #' hist(bigvalues$RoadDen[bigvalues$outlier == "not_outlier"], breaks = 10, main = "Frequency of Road Density values with NO RSF outliers", xlab = "Road Density")
    #' hist(bigvalues$RoadDen[bigvalues$outlier == "outlier"], breaks = 10, main = "Frequency of Road Density values for RSF outliers", xlab = "Road Density")
    #' 
    #' bigvalues <- cbind(predicted, covs_list)
    #' plot(bigvalues$Dist2Water, bigvalues$predict_rsf, main = "Correlation between Predicted RSF and Distance to Water", xlab = "Distance to Water", ylab = "Predicted RSF Values")
    #' #'  Simple regression between covaraites
    #' abline(lm(bigvalues$predict_rsf ~ bigvalues$Dist2Water), col = "red")
    #' #'  Mark the 99th percentil value
    #' abline(h = quant, col = "blue")
    #' #'  Frequency of covariates associated with non-extreme & extreme RSF values
    #' hist(bigvalues$Dist2Water, breaks = 10, main = "Frequency of Distance to Water values", xlab = "Distance to Water")
    #' hist(bigvalues$Dist2Water[bigvalues$outlier == "not_outlier"], breaks = 10, main = "Frequency of Distance to Water values with NO RSF outliers", xlab = "Distance to Water")
    #' hist(bigvalues$Dist2Water[bigvalues$outlier == "outlier"], breaks = 10, main = "Frequency of Distance to Water values for RSF outliers", xlab = "Distance to Water")

    return(bigvalues)
  }
  #'  Identify outlier predictions and possible covariates associated with those
  #'  Be sure to used standardized covariates for evaluation
  md_smr_outliers <- lapply(md_smr_rsf_sa, outliers, title = "Mule Deer Summer RSF Predictions", covs_list = md_smr_zcovs[[1]]) 
  md_wtr_outliers <- lapply(md_wtr_rsf_sa, outliers, title = "Mule Deer Winter RSF Predictions", covs_list = md_wtr_zcovs[[1]]) #' one extreme slope value associated w/ one big outlier (mask out pixels Slope >6)
  elk_smr_outliers <- lapply(elk_smr_rsf_sa, outliers, title = "Elk Summer RSF Predictions", covs_list = elk_smr_zcovs[[1]]) #' one extreme dist2water value associated w/ one big outlier (mask out pixels Dist2Water >9)
  elk_wtr_outliers <- lapply(elk_wtr_rsf_sa, outliers, title = "Elk Winter RSF Predictions", covs_list = elk_wtr_zcovs[[1]]) #' extreme dist2water values associated w/ outliers (mask out pixels Dist2Water >9)
  wtd_smr_outliers <- lapply(wtd_smr_rsf_sa, outliers, title = "White-tailed Deer Summer RSF Predictions", covs_list = wtd_smr_zcovs[[1]]) #' extreme road density values associated w/ outliers (mask out pixels RoadDen >13)
  wtd_wtr_outliers <- lapply(wtd_wtr_rsf_sa, outliers, title = "White-tailed Deer Winter RSF Predictions", covs_list = wtd_wtr_zcovs[[1]]) #' extreme road density & Dist2Water values associated w/ outliers (mask out pixels RoadDen >13, Dist2Water >9)
  coug_smr_outliers <- lapply(coug_smr_rsf_sa, outliers, title = "Cougar Summer RSF Predictions", covs_list = coug_smr_zcovs[[1]])
  coug_wtr_outliers <- lapply(coug_wtr_rsf_sa, outliers, title = "Cougar Winter RSF Predictions", covs_list = coug_wtr_zcovs[[1]])
  wolf_smr_outliers <- lapply(wolf_smr_rsf_sa, outliers, title = "Wolf Summer RSF Predictions", covs_list = wolf_smr_zcovs[[1]])
  wolf_wtr_outliers <- lapply(wolf_wtr_rsf_sa, outliers, title = "Wolf Winter RSF Predictions", covs_list = wolf_wtr_zcovs[[1]])
  bob_smr_outliers <- lapply(bob_smr_rsf_sa, outliers, title = "Bobcat Summer RSF Predictions", covs_list = bob_smr_zcovs[[1]])
  bob_wtr_outliers <- lapply(bob_wtr_rsf_sa, outliers, title = "Bobcat Winter RSF Predictions", covs_list = bob_wtr_zcovs[[1]])
  coy_smr_outliers <- lapply(coy_smr_rsf_sa, outliers, title = "Coyote Summer RSF Predictions", covs_list = bob_smr_zcovs[[1]]) #' extreme road density, Dist2Water, & HumanMod values for some outliers but also for non-outliers
  coy_wtr_outliers <- lapply(coy_wtr_rsf_sa, outliers, title = "Coyote Winter RSF Predictions", covs_list = coy_wtr_zcovs[[1]]) #' extreme road density values associated w/ outliers (mask out pixels RoadDen >13)
  #'  After reviewing extreme RSF values and covaraite values associated with 
  #'  those locations I am masking pixels for:
  #'   1. Winter mule deer where SLOPE > 6
  #'   2. Summer and Winter elk where DIST2WATER > 9
  #'   3. Summer white-tailed deer where ROADDEN > 13
  #'   4. Winter white-tailed deer where ROADDEN > 13 & DIST2WATER > 9
  #'   5. Winter coyote where ROADDEN > 13
  
  #'  Exclude extreme outliers as identified above
  #'  Mule Deer
  mask_md <- function(predicted) {
    predicted <- predicted %>%
      mutate(masked_rsf = ifelse(Slope > 6, "NA", predict_rsf),
             masked_rsf = as.numeric(masked_rsf))
  }
  md_wtr_outliers <- lapply(md_wtr_outliers, mask_md)
  #'  Elk
  mask_elk <- function(predicted) {
    predicted <- predicted %>%
      mutate(masked_rsf = ifelse(Dist2Water > 9, "NA", predict_rsf),
             masked_rsf = as.numeric(masked_rsf))
  }
  elk_smr_outliers <- lapply(elk_smr_outliers, mask_elk)
  elk_wtr_outliers <- lapply(elk_wtr_outliers, mask_elk)
  #'  Winter White-tailed Deer
  mask_wtd <- function(predicted) {
    predicted <- predicted %>%
      mutate(masked_rsf = ifelse(RoadDen > 13, "NA", predict_rsf),
             masked_rsf = ifelse(Dist2Water > 9, "NA", masked_rsf),
             masked_rsf = as.numeric(masked_rsf))
  }
  wtd_wtr_outliers <- lapply(wtd_wtr_outliers, mask_wtd)
  #'  Summer White-tailed Deer, Winter Coyote
  mask_rd <- function(predicted) {
    predicted <- predicted %>%
      mutate(masked_rsf = ifelse(RoadDen > 13, "NA", predict_rsf),
             masked_rsf = as.numeric(masked_rsf))
  }
  wtd_smr_outliers <- lapply(wtd_smr_outliers, mask_rd)
  coy_wtr_outliers <- lapply(coy_wtr_outliers, mask_rd)
  
  ####  Even after masking out extreme values related to high covariate values
  ####  there are still extreme RSF values that skew the rescaling and lead to
  ####  misleading plots. Using adjusted RSF values instead where extreme 1% are
  ####  set to 99th percentile value.
  

  #'  Re-scale predicted RSF values between 0 & 1 for plotting
  RSF_rescale <- function(out) {
    rescale_val <- out %>%
      mutate(
        rescale_rsf = round(adjusted_rsf/(max(adjusted_rsf, na.rm = T)), digits = 4)
        # rescale_rsf = round(masked_rsf/(max(masked_rsf, na.rm = T)), digits = 4)
        # rescale_rsf = round(predict_rsf/(max(predict_rsf, na.rm = T)), digits = 4)
      ) %>%
      dplyr::select(c(x, y, rescale_rsf))
    return(rescale_val)
  }
  #'  Rescale predicted RSF values within each list of lists
  md_smr_rescale_sa <- lapply(md_smr_outliers, RSF_rescale)  #md_smr_rsf_sa
  md_wtr_rescale_sa <- lapply(md_wtr_outliers, RSF_rescale) 
  elk_smr_rescale_sa <- lapply(elk_smr_outliers, RSF_rescale) 
  elk_wtr_rescale_sa <- lapply(elk_wtr_outliers, RSF_rescale) 
  wtd_smr_rescale_sa <- lapply(wtd_smr_outliers, RSF_rescale)  
  wtd_wtr_rescale_sa <- lapply(wtd_wtr_outliers, RSF_rescale) 
  coug_smr_rescale_sa <- lapply(coug_smr_outliers, RSF_rescale)  
  coug_wtr_rescale_sa <- lapply(coug_wtr_outliers, RSF_rescale)
  wolf_smr_rescale_sa <- lapply(wolf_smr_outliers, RSF_rescale)  
  wolf_wtr_rescale_sa <- lapply(wolf_wtr_outliers, RSF_rescale)
  bob_smr_rescale_sa <- lapply(bob_smr_outliers, RSF_rescale)  
  bob_wtr_rescale_sa <- lapply(bob_wtr_outliers, RSF_rescale)
  coy_smr_rescale_sa <- lapply(coy_smr_outliers, RSF_rescale)  
  coy_wtr_rescale_sa <- lapply(coy_wtr_outliers, RSF_rescale) 

  
  chk <- md_smr_rescale_sa[[1]]
  
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
  md_smr_RSFraster <- lapply(md_smr_rescale_sa, rasterize_rsf)
  md_wtr_RSFraster <- lapply(md_wtr_rescale_sa, rasterize_rsf)
  elk_smr_RSFraster <- lapply(elk_smr_rescale_sa, rasterize_rsf)
  elk_wtr_RSFraster <- lapply(elk_wtr_rescale_sa, rasterize_rsf)
  wtd_smr_RSFraster <- lapply(wtd_smr_rescale_sa, rasterize_rsf) 
  wtd_wtr_RSFraster <- lapply(wtd_wtr_rescale_sa, rasterize_rsf) 
  coug_smr_RSFraster <- lapply(coug_smr_rescale_sa, rasterize_rsf)
  coug_wtr_RSFraster <- lapply(coug_wtr_rescale_sa, rasterize_rsf)
  wolf_smr_RSFraster <- lapply(wolf_smr_rescale_sa, rasterize_rsf)
  wolf_wtr_RSFraster <- lapply(wolf_wtr_rescale_sa, rasterize_rsf)
  bob_smr_RSFraster <- lapply(bob_smr_rescale_sa, rasterize_rsf)
  bob_wtr_RSFraster <- lapply(bob_wtr_rescale_sa, rasterize_rsf)
  coy_smr_RSFraster <- lapply(coy_smr_rescale_sa, rasterize_rsf)
  coy_wtr_RSFraster <- lapply(coy_wtr_rescale_sa, rasterize_rsf)
  
  #'  Rename rasters
  rename_raster <- function(raster_list) {
    L <- setNames(raster_list, c("Year1", "Year2", "Year3"))
    S <- stack(L)
    return(S)
  }
  md_smr_RSFstack <- rename_raster(md_smr_RSFraster)
  md_wtr_RSFstack <- rename_raster(md_wtr_RSFraster)
  elk_smr_RSFstack <- rename_raster(elk_smr_RSFraster)
  elk_wtr_RSFstack <- rename_raster(elk_wtr_RSFraster)
  wtd_smr_RSFstack <- rename_raster(wtd_smr_RSFraster)
  wtd_wtr_RSFstack <- rename_raster(wtd_wtr_RSFraster)
  coug_smr_RSFstack <- rename_raster(coug_smr_RSFraster)
  coug_wtr_RSFstack <- rename_raster(coug_wtr_RSFraster)
  wolf_smr_RSFstack <- rename_raster(wolf_smr_RSFraster)
  wolf_wtr_RSFstack <- rename_raster(wolf_wtr_RSFraster)
  bob_smr_RSFstack <- rename_raster(bob_smr_RSFraster)
  bob_wtr_RSFstack <- rename_raster(bob_wtr_RSFraster)
  coy_smr_RSFstack <- rename_raster(coy_smr_RSFraster)
  coy_wtr_RSFstack <- rename_raster(coy_wtr_RSFraster)

  
  #'  Plot & Save
  pdf(file = "./Outputs/RSF_Maps.pdf")
  plot(md_smr_RSFstack[[1]], main = "Summer Mule Deer Predicted RSF (2018)"); plot(OK.SA, add = T)
  plot(md_wtr_RSFstack[[1]], main = "Winter Mule Deer Predicted RSF (2018)"); plot(OK.SA, add = T)
  plot(elk_smr_RSFstack[[1]], main = "Summer Elk Predicted RSF (2018)"); plot(NE.SA, add = T)
  plot(elk_wtr_RSFstack[[1]], main = "Winter Elk Predicted RSF (2018)"); plot(NE.SA, add = T)
  plot(wtd_smr_RSFstack[[1]], main = "Summer White-tailed Deer Predicted RSF (2018)"); plot(NE.SA, add = T)
  plot(wtd_wtr_RSFstack[[1]], main = "Winter White-tailed Deer Predicted RSF (2018)"); plot(NE.SA, add = T)
  plot(coug_smr_RSFstack[[1]], main = "Summer Cougar Predicted RSF (2018)"); plot(OK.SA, add = T); plot(NE.SA, add = T)
  plot(coug_wtr_RSFstack[[1]], main = "Winter Cougar Predicted RSF (2018)"); plot(OK.SA, add = T); plot(NE.SA, add = T)
  plot(wolf_smr_RSFstack[[1]], main = "Summer Wolf Predicted RSF (2018)"); plot(OK.SA, add = T); plot(NE.SA, add = T)
  plot(wolf_wtr_RSFstack[[1]], main = "Winter Wolf Predicted RSF (2018)"); plot(OK.SA, add = T); plot(NE.SA, add = T)
  plot(bob_smr_RSFstack[[1]], main = "Summer Bobcat Predicted RSF (2018)"); plot(OK.SA, add = T); plot(NE.SA, add = T)
  plot(bob_wtr_RSFstack[[1]], main = "Winter Bobcat Predicted RSF (2018)"); plot(OK.SA, add = T); plot(NE.SA, add = T)
  plot(coy_smr_RSFstack[[1]], main = "Summer Coyote Predicted RSF (2018)"); plot(OK.SA, add = T); plot(NE.SA, add = T)
  plot(coy_wtr_RSFstack[[1]], main = "Winter Coyote Predicted RSF (2018)"); plot(OK.SA, add = T); plot(NE.SA, add = T)
  dev.off()
  
  
  
  #'  SAVE!
  writeRaster(md_smr_RSFstack, filename = "./Shapefiles/Predicted_RSFs/md_smr_RSFstack.tif", bylayer = FALSE, format = 'GTiff', overwrite = TRUE)
  writeRaster(md_wtr_RSFstack, filename = "./Shapefiles/Predicted_RSFs/md_wtr_RSFstack.tif", bylayer = FALSE, format = 'GTiff', overwrite = TRUE)
  writeRaster(elk_smr_RSFstack, filename = "./Shapefiles/Predicted_RSFs/elk_smr_RSFstack.tif", bylayer = FALSE, format = 'GTiff', overwrite = TRUE)
  writeRaster(elk_wtr_RSFstack, filename = "./Shapefiles/Predicted_RSFs/elk_wtr_RSFstack.tif", bylayer = FALSE, format = 'GTiff', overwrite = TRUE)
  writeRaster(wtd_smr_RSFstack, filename = "./Shapefiles/Predicted_RSFs/wtd_smr_RSFstack.tif", bylayer = FALSE, format = 'GTiff', overwrite = TRUE)
  writeRaster(wtd_wtr_RSFstack, filename = "./Shapefiles/Predicted_RSFs/wtd_wtr_RSFstack.tif", bylayer = FALSE, format = 'GTiff', overwrite = TRUE)
  writeRaster(coug_smr_RSFstack, filename = "./Shapefiles/Predicted_RSFs/coug_smr_RSFstack.tif", bylayer = FALSE, format = 'GTiff', overwrite = TRUE)
  writeRaster(coug_wtr_RSFstack, filename = "./Shapefiles/Predicted_RSFs/coug_wtr_RSFstack.tif", bylayer = FALSE, format = 'GTiff', overwrite = TRUE)
  writeRaster(wolf_smr_RSFstack, filename = "./Shapefiles/Predicted_RSFs/wolf_smr_RSFstack.tif", bylayer = FALSE, format = 'GTiff', overwrite = TRUE)
  writeRaster(wolf_wtr_RSFstack, filename = "./Shapefiles/Predicted_RSFs/wolf_wtr_RSFstack.tif", bylayer = FALSE, format = 'GTiff', overwrite = TRUE)
  writeRaster(bob_smr_RSFstack, filename = "./Shapefiles/Predicted_RSFs/bob_smr_RSFstack.tif", bylayer = FALSE, format = 'GTiff', overwrite = TRUE)
  writeRaster(bob_wtr_RSFstack, filename = "./Shapefiles/Predicted_RSFs/bob_wtr_RSFstack.tif", bylayer = FALSE, format = 'GTiff', overwrite = TRUE)
  writeRaster(coy_smr_RSFstack, filename = "./Shapefiles/Predicted_RSFs/coy_smr_RSFstack.tif", bylayer = FALSE, format = 'GTiff', overwrite = TRUE)
  writeRaster(coy_wtr_RSFstack, filename = "./Shapefiles/Predicted_RSFs/coy_wtr_RSFstack.tif", bylayer = FALSE, format = 'GTiff', overwrite = TRUE)
  

  
  
    #' ####  Summary tables  ####
  #' #'  Save model outputs in table format 
  #' #'  Functions extract outputs for each sub-model and appends species/season info
  #' 
  #' #'  Pull out RSF results
  #' load("./Outputs/RSF_output/RSF_MD_list_2022-01-02.RData")  
  #' load("./Outputs/RSF_output/RSF_ELK_list_2022-01-02.RData")
  #' load("./Outputs/RSF_output/RSF_WTD_list_2022-01-02.RData")
  #' load("./Outputs/RSF_output/RSF_COUG_list_2022-01-02.RData") 
  #' load("./Outputs/RSF_output/RSF_WOLF_list_2022-01-02.RData")
  #' load("./Outputs/RSF_output/RSF_BOB_list_2022-01-02.RData")
  #' load("./Outputs/RSF_output/RSF_COY_list_2022-01-02.RData")
  #' 
  #' 
  #' #'  Function to save parameter estimates & p-values
  #' #'  use coef(mod) to look at random effects estimates
  #' rounddig <- 2
  #' 
  #' rsf_out <- function(mod, spp, season){
  #'   betas <- mod@beta
  #'   se <- sqrt(diag(vcov(mod)))
  #'   z <- summary(mod)$coef[,3]
  #'   pval <- summary(mod)$coef[,4]
  #'   out <- as.data.frame(cbind(betas, se, pval)) %>%
  #'     transmute(
  #'       Species = rep(spp, nrow(.)),
  #'       Season = rep(season, nrow(.)),
  #'       Parameter = row.names(.),
  #'       Estimate = round(betas, rounddig),
  #'       SE = round(se, rounddig),
  #'       Z = round(z, rounddig),
  #'       Pval = round(pval, rounddig)) 
  #'   rownames(out) <- NULL
  #'   return(out)
  #' }
  #' md_s1819_rsf <- rsf_out(md_global_smr, "Mule Deer", "Summer")
  #' md_w1820_rsf <- rsf_out(md_global_wtr, "Mule Deer", "Winter")
  #' elk_s1819_rsf <- rsf_out(elk_global_smr, "Elk", "Summer")
  #' elk_w1820_rsf <- rsf_out(elk_global_wtr, "Elk", "Winter")
  #' wtd_s1819_rsf <- rsf_out(wtd_global_smr, "White-tailed Deer", "Summer")
  #' wtd_w1820_rsf <- rsf_out(wtd_global_wtr, "White-tailed Deer", "Winter")
  #' coug_s1819_rsf <- rsf_out(coug_global_smr, "Cougar", "Summer")
  #' coug_w1820_rsf <- rsf_out(coug_global_wtr, "Cougar", "Winter")
  #' wolf_s1819_rsf <- rsf_out(wolf_global_smr, "Wolf", "Summer")
  #' wolf_w1820_rsf <- rsf_out(wolf_global_wtr, "Wolf", "Winter")
  #' bob_s1819_rsf <- rsf_out(bob_global_smr, "Bobcat", "Summer")
  #' bob_w1820_rsf <- rsf_out(bob_global_wtr, "Bobcat", "Winter")
  #' coy_s1819_rsf <- rsf_out(coy_global_smr, "Coyote", "Summer")
  #' coy_w1820_rsf <- rsf_out(coy_global_wtr, "Coyote", "Winter")
  #' 
  #' md_sSAonly_rsf <- rsf_out(md_SA_only_smr, "Mule Deer", "Summer")
  #' # write.csv(md_sSAonly_rsf, paste0("./Outputs/Tables/RSF_Results_MuleDeer_SAonly_", Sys.Date(), ".csv"))  
  #' 
  #' 
  #' 
  #' #'  Merge into larger data frames for easy comparison
  #' summer_rsf <- rbind(bob_s1819_rsf, coug_s1819_rsf, coy_s1819_rsf, wolf_s1819_rsf,
  #'                     elk_s1819_rsf, md_s1819_rsf, wtd_s1819_rsf) 
  #' winter_rsf <- rbind(bob_w1820_rsf, coug_w1820_rsf, coy_w1820_rsf, wolf_w1820_rsf,
  #'                     elk_w1820_rsf, md_w1820_rsf, wtd_w1820_rsf) 
  #' rsf_results <- rbind(summer_rsf, winter_rsf) %>%
  #'   arrange(Species)
  #' colnames(rsf_results) <- c("Species", "Season", "Parameter", "Estimate", "SE", "z", "Pval")
  #' 
  #' 
  #' #'  Spread this out so the coefficient effects are easier to compare across species
  #' rsf_results_wide <- rsf_results %>% 
  #'   dplyr::select(-z) %>%
  #'   mutate(
  #'     SE = round(SE, 2),
  #'     SE = paste0("(", SE, ")")
  #'   ) %>%
  #'   #'  Bold significant variables- doesn't work if continue manipulating data frame
  #'   # condformat(.) %>%
  #'   # rule_text_bold(c(Estimate, SE, Pval), expression = Pval <= 0.05) %>%
  #'   unite(Est_SE, Estimate, SE, sep = " ") %>%
  #'   unite(Est_SE_Pval, Est_SE, Pval, sep = "_") %>%
  #'   spread(Parameter, Est_SE_Pval) %>%
  #'   separate("(Intercept)", c("Intercept (SE)", "Intercept Pval"), sep = "_") %>%
  #'   # separate("AreaOK", c("AreaOK (SE)", "AreaOK Pval"), sep = "_") %>%
  #'   separate("Elev", c("Elev (SE)", "Elev Pval"), sep = "_") %>%
  #'   separate("Slope", c("Slope (SE)", "Slope Pval"), sep = "_") %>%
  #'   separate("PercForMix", c("PercForMix (SE)", "PercForMix Pval"), sep = "_") %>%
  #'   separate("PercXGrass", c("PercXGrass (SE)", "PercXGrass Pval"), sep = "_") %>%
  #'   separate("PercXShrub", c("PercXShrub (SE)", "PercXShrub Pval"), sep = "_") %>%
  #'   separate("RoadDen", c("Road Density (SE)", "Road Density Pval"), sep = "_") %>%
  #'   # separate("HumanMod", c("HumanMod (SE)", "HumanMod Pval"), sep = "_") %>%
  #'   arrange(match(Species, c("Bobcat", "Cougar", "Coyote", "Wolf", "Mule Deer", "Elk", "White-tailed Deer"))) %>%
  #'   arrange(match(Season, c("Summer", "Winter")))
  #' 
  #' 
  #' #'  Save!
  #' write.csv(rsf_results, paste0("./Outputs/Tables/RSF_Results_NoHM_", Sys.Date(), ".csv"))  #'  KEEP TRACK of whether Human Modified was excluded from models!
  #' write.csv(rsf_results_wide, paste0("./Outputs/Tables/RSF_Results_wide_NoHM_", Sys.Date(), ".csv"))
  #' 
  #' save.image("./Outputs/RSF_script_results.RData")
  #' # save(md_smr18, file = paste0("./Outputs/RSF_output/md_RSF_smr18_", Sys.Date(), ".RData"))
  #' # save(md_smr19, file = paste0("./Outputs/RSF_output/md_RSF_smr19_", Sys.Date(), ".RData"))
  #' # save(md_smr20, file = paste0("./Outputs/RSF_output/md_RSF_smr20_", Sys.Date(), ".RData"))
  #' # save(md_wtr1819, file = paste0("./Outputs/RSF_output/md_RSF_wtr1819_", Sys.Date(), ".RData"))
  #' # save(md_wtr1920, file = paste0("./Outputs/RSF_output/md_RSF_wtr1920_", Sys.Date(), ".RData"))
  #' # save(md_wtr2021, file = paste0("./Outputs/RSF_output/md_RSF_wtr2021_", Sys.Date(), ".RData"))
  #' # save(elk_smr18, file = paste0("./Outputs/RSF_output/elk_RSF_smr18_", Sys.Date(), ".RData"))
  #' # save(elk_smr19, file = paste0("./Outputs/RSF_output/elk_RSF_smr19_", Sys.Date(), ".RData"))
  #' # save(elk_smr20, file = paste0("./Outputs/RSF_output/elk_RSF_smr20_", Sys.Date(), ".RData"))
  #' # save(elk_wtr1819, file = paste0("./Outputs/RSF_output/elk_RSF_wtr1819_", Sys.Date(), ".RData"))
  #' # save(elk_wtr1920, file = paste0("./Outputs/RSF_output/elk_RSF_wtr1920_", Sys.Date(), ".RData"))
  #' # save(elk_wtr2021, file = paste0("./Outputs/RSF_output/elk_RSF_wtr2021_", Sys.Date(), ".RData"))
  #' # save(wtd_smr18, file = paste0("./Outputs/RSF_output/wtd_RSF_smr18_", Sys.Date(), ".RData"))
  #' # save(wtd_smr19, file = paste0("./Outputs/RSF_output/wtd_RSF_smr19_", Sys.Date(), ".RData"))
  #' # save(wtd_smr20, file = paste0("./Outputs/RSF_output/wtd_RSF_smr20_", Sys.Date(), ".RData"))
  #' # save(wtd_wtr1819, file = paste0("./Outputs/RSF_output/wtd_RSF_wtr1819_", Sys.Date(), ".RData"))
  #' # save(wtd_wtr1920, file = paste0("./Outputs/RSF_output/wtd_RSF_wtr1920_", Sys.Date(), ".RData"))
  #' # save(wtd_wtr2021, file = paste0("./Outputs/RSF_output/wtd_RSF_wtr2021_", Sys.Date(), ".RData"))
  #' # save(coug_smr18, file = paste0("./Outputs/RSF_output/coug_RSF_smr18_", Sys.Date(), ".RData"))
  #' # save(coug_smr19, file = paste0("./Outputs/RSF_output/coug_RSF_smr19_", Sys.Date(), ".RData"))
  #' # save(coug_smr20, file = paste0("./Outputs/RSF_output/coug_RSF_smr20_", Sys.Date(), ".RData"))
  #' # save(coug_wtr1819, file = paste0("./Outputs/RSF_output/coug_RSF_wtr1819_", Sys.Date(), ".RData"))
  #' # save(coug_wtr1920, file = paste0("./Outputs/RSF_output/coug_RSF_wtr1920_", Sys.Date(), ".RData"))
  #' # save(coug_wtr2021, file = paste0("./Outputs/RSF_output/coug_RSF_wtr2021_", Sys.Date(), ".RData"))
  #' # save(wolf_smr18, file = paste0("./Outputs/RSF_output/wolf_RSF_smr18_", Sys.Date(), ".RData"))
  #' # save(wolf_smr19, file = paste0("./Outputs/RSF_output/wolf_RSF_smr19_", Sys.Date(), ".RData"))
  #' # save(wolf_smr20, file = paste0("./Outputs/RSF_output/wolf_RSF_smr20_", Sys.Date(), ".RData"))
  #' # save(wolf_wtr1819, file = paste0("./Outputs/RSF_output/wolf_RSF_wtr1819_", Sys.Date(), ".RData"))
  #' # save(wolf_wtr1920, file = paste0("./Outputs/RSF_output/wolf_RSF_wtr1920_", Sys.Date(), ".RData"))
  #' # save(wolf_wtr2021, file = paste0("./Outputs/RSF_output/wolf_RSF_wtr2021_", Sys.Date(), ".RData"))
  #' # save(coy_smr18, file = paste0("./Outputs/RSF_output/coy_RSF_smr18_", Sys.Date(), ".RData"))
  #' # save(coy_smr19, file = paste0("./Outputs/RSF_output/coy_RSF_smr19_", Sys.Date(), ".RData"))
  #' # save(coy_smr20, file = paste0("./Outputs/RSF_output/coy_RSF_smr20_", Sys.Date(), ".RData"))
  #' # save(coy_wtr1819, file = paste0("./Outputs/RSF_output/coy_RSF_wtr1819_", Sys.Date(), ".RData"))
  #' # save(coy_wtr1920, file = paste0("./Outputs/RSF_output/coy_RSF_wtr1920_", Sys.Date(), ".RData"))
  #' # save(coy_wtr2021, file = paste0("./Outputs/RSF_output/coy_RSF_wtr2021_", Sys.Date(), ".RData"))
  #' # save(bob_smr18, file = paste0("./Outputs/RSF_output/bob_RSF_smr18_", Sys.Date(), ".RData"))
  #' # save(bob_smr19, file = paste0("./Outputs/RSF_output/bob_RSF_smr19_", Sys.Date(), ".RData"))
  #' # save(bob_smr20, file = paste0("./Outputs/RSF_output/bob_RSF_smr20_", Sys.Date(), ".RData"))
  #' # save(bob_wtr1819, file = paste0("./Outputs/RSF_output/bob_RSF_wtr1819_", Sys.Date(), ".RData"))
  #' # save(bob_wtr1920, file = paste0("./Outputs/RSF_output/bob_RSF_wtr1920_", Sys.Date(), ".RData"))
  #' # save(bob_wtr2021, file = paste0("./Outputs/RSF_output/bob_RSF_wtr2021_", Sys.Date(), ".RData"))

  