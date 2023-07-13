  #'  ============================================
  #'  3rd Order Resource Selection Functions (RSFs)
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  July 2023
  #'  ============================================
  #'  Script to run 3rd order resource selection functions (RSFs) and predict
  #'  relative probability of selection across study areas for each species and
  #'  season.
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
  
  #'  Load used and available locations, and covariate data
  load("./Outputs/RSF_pts/coug_dat_all_forTRG_2022-06-15.RData")
  load("./Outputs/RSF_pts/wolf_dat_all_forTRG_2022-06-15.RData")
  load("./Outputs/RSF_pts/bob_dat_all_forTRG_2022-06-15.RData")
  load("./Outputs/RSF_pts/coy_dat_all_forTRG_2022-06-15.RData")

  #'  Ensure use/available data labeled correctly by study area for wolves
  wolf_dat_all <- mutate(wolf_dat_all,
                         Area = ifelse(grepl("W61M", ID), "OK", "NE"),  
                         Area = ifelse(grepl("W88M", ID), "OK", Area),
                         Area = ifelse(grepl("W93M", ID), "OK", Area),
                         Area = ifelse(grepl("W94M", ID), "OK", Area),
                         Area = ifelse(grepl("W110M", ID), "OK", Area))
  
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
    #'  Leave weights as is
    locs$w <- locs$w
    
    locs <- as.data.frame(locs)
    
    return(locs)
  }
  #'  Subset data sets by season & standardize covariates
  #'  Note the reclassified landcover_type data for cougData_winter
  cougData_smr <- coug_dat_all_reclass[coug_dat_all_reclass$Season == "Summer17" | coug_dat_all_reclass$Season == "Summer18" | coug_dat_all_reclass$Season == "Summer19" | coug_dat_all_reclass$Season == "Summer20",]
  cougData_fll <- coug_dat_all_reclass[coug_dat_all_reclass$Season == "Fall17" | coug_dat_all_reclass$Season == "Fall18" | coug_dat_all_reclass$Season == "Fall19" | coug_dat_all_reclass$Season == "Fall20",]
  cougData_wtr <- coug_dat_all_reclass[coug_dat_all_reclass$Season == "Winter1718" | coug_dat_all_reclass$Season == "Winter1819" | coug_dat_all_reclass$Season == "Winter1920" | coug_dat_all_reclass$Season == "Winter2021",]
  cougData_sprg <- coug_dat_all_reclass[coug_dat_all_reclass$Season == "Spring18" | coug_dat_all_reclass$Season == "Spring19" | coug_dat_all_reclass$Season == "Spring20" | coug_dat_all_reclass$Season == "Spring21",]
  cougDataz_smr <- standardize_covs(cougData_smr)
  cougDataz_fll <- standardize_covs(cougData_fll)
  cougDataz_wtr <- standardize_covs(cougData_wtr)
  cougDataz_sprg <- standardize_covs(cougData_sprg)
  #'  Note the double reclassified landcover_type data for wolfData
  wolfData_smr <- wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Summer17" | wolf_dat_all_reclass2$Season == "Summer18" | wolf_dat_all_reclass2$Season == "Summer19" | wolf_dat_all_reclass2$Season == "Summer20",]
  wolfData_fll <- wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Fall17" | wolf_dat_all_reclass2$Season == "Fall18" | wolf_dat_all_reclass2$Season == "Fall19" | wolf_dat_all_reclass2$Season == "Fall20",]
  wolfData_wtr <- wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Winter1718" | wolf_dat_all_reclass2$Season == "Winter1819" | wolf_dat_all_reclass2$Season == "Winter1920" | wolf_dat_all_reclass2$Season == "Winter2021",]
  wolfData_sprg <- wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Spring18" | wolf_dat_all_reclass2$Season == "Spring19" | wolf_dat_all_reclass2$Season == "Spring20" | wolf_dat_all_reclass2$Season == "Spring21",]
  wolfDataz_smr <- standardize_covs(wolfData_smr)
  wolfDataz_fll <- standardize_covs(wolfData_fll)
  wolfDataz_wtr <- standardize_covs(wolfData_wtr)
  wolfDataz_sprg <- standardize_covs(wolfData_sprg)
  #'  Note the reclassified landcover_type data for bobData
  bobData_smr <- bob_dat_all_reclass[bob_dat_all_reclass$Season == "Summer17" | bob_dat_all_reclass$Season == "Summer18" | bob_dat_all_reclass$Season == "Summer19" | bob_dat_all_reclass$Season == "Summer20",]
  bobData_fll <- bob_dat_all_reclass[bob_dat_all_reclass$Season == "Fall17" | bob_dat_all_reclass$Season == "Fall18" | bob_dat_all_reclass$Season == "Fall19" | bob_dat_all_reclass$Season == "Fall20",]
  bobData_wtr <- bob_dat_all_reclass[bob_dat_all_reclass$Season == "Winter1718" | bob_dat_all_reclass$Season == "Winter1819" | bob_dat_all_reclass$Season == "Winter1920" | bob_dat_all_reclass$Season == "Winter2021",]
  bobData_sprg <- bob_dat_all_reclass[bob_dat_all_reclass$Season == "Summer18" | bob_dat_all_reclass$Season == "Summer19" | bob_dat_all_reclass$Season == "Summer20"| bob_dat_all_reclass$Season == "Spring21",]
  bobDataz_smr <- standardize_covs(bobData_smr)
  bobDataz_fll <- standardize_covs(bobData_fll)
  bobDataz_wtr <- standardize_covs(bobData_wtr)
  bobDataz_sprg <- standardize_covs(bobData_sprg)
  #'  Note the reclassified landcover_type data for coyData
  coyData_smr <- coy_dat_all_reclass[coy_dat_all_reclass$Season == "Summer17" | coy_dat_all_reclass$Season == "Summer18" | coy_dat_all_reclass$Season == "Summer19" | coy_dat_all_reclass$Season == "Summer20",]
  coyData_fll <- coy_dat_all_reclass[coy_dat_all_reclass$Season == "Fall17" | coy_dat_all_reclass$Season == "Fall18" | coy_dat_all_reclass$Season == "Fall19" | coy_dat_all_reclass$Season == "Fall20",]
  coyData_wtr <- coy_dat_all_reclass[coy_dat_all_reclass$Season == "Winter1718" | coy_dat_all_reclass$Season == "Winter1819" | coy_dat_all_reclass$Season == "Winter1920" | coy_dat_all_reclass$Season == "Winter2021",]
  coyData_sprg <- coy_dat_all_reclass[coy_dat_all_reclass$Season == "Spring18" | coy_dat_all_reclass$Season == "Spring19" | coy_dat_all_reclass$Season == "Spring20"| coy_dat_all_reclass$Season == "Spring21",]
  coyDataz_smr <- standardize_covs(coyData_smr)
  coyDataz_fll <- standardize_covs(coyData_fll)
  coyDataz_wtr <- standardize_covs(coyData_wtr)
  coyDataz_sprg <- standardize_covs(coyData_sprg)
  
  #'  Correlation Matrix
  #'  ==================
  #'  Function to create correlation matrix for all continuous covariates at once
  #'  Ignore PercForMix, PercXGrass, PercXShrub- they aren't being used in RSFs
  cov_correlation <- function(dat) {
    used <- dat[dat$Used == 1,]
    covs <- used[,c("Elev", "Slope", "TPI", "RoadDen",
                    "Dist2Water", "HumanMod", "CanopyCover", "Dist2Edge",
                    "PercForMix", "PercXGrass", "PercXShrub")]
    cor_matrix <- cor(covs, use = "complete.obs")
    return(cor_matrix)
  }
  #'  Generate correlation matrix for each species and season
  (coug_smr_corr <- cov_correlation(cougData_smr))
  (coug_fll_corr <- cov_correlation(cougData_fll))
  (coug_wtr_corr <- cov_correlation(cougData_wtr))
  (coug_sprg_corr <- cov_correlation(cougData_sprg))
  (wolf_smr_corr <- cov_correlation(wolfData_smr))
  (wolf_fll_corr <- cov_correlation(wolfData_fll))
  (wolf_wtr_corr <- cov_correlation(wolfData_wtr))
  (wolf_sprg_corr <- cov_correlation(wolfData_sprg))
  (bob_smr_corr <- cov_correlation(bobData_smr))
  (bob_fll_corr <- cov_correlation(bobData_fll))
  (bob_wtr_corr <- cov_correlation(bobData_wtr))
  (bob_sprg_corr <- cov_correlation(bobData_sprg))
  (coy_smr_corr <- cov_correlation(coyData_smr))
  (coy_fll_corr <- cov_correlation(coyData_fll))
  (coy_wtr_corr <- cov_correlation(coyData_wtr))
  (coy_sprg_corr <- cov_correlation(coyData_sprg))
  
  #'  Just checking but will still use all covariates in each RSF even if highly
  #'  correlated b/c trying to make as predictive as possible
  
  
  #'  Resource Selection Function Models
  #'  ==================================
  #'  Functions to run logistic mixed effects models that include random effect 
  #'  for individual. Habitat covariates excluded based on species (see notes 
  #'  above) and convergence issues. Seasonal models run separately so RSFs are 
  #'  specific to the species & season but with all years combined. Some covariates
  #'  vary annually however so each species & season-specific RSF is predicted 
  #'  across an annual study area map (see prediction section below). 
  #'  
  #'  Using backwards step selection to identify top model with all significant
  #'  covariates. Covariates excluded through this process are commented out next
  #'  to model in the order each variable was removed.
  #'  ==================================
  
  glmm_fn <- function(mod, dat) {
    glmm_mod <- glmer(formula = mod, data = dat, weights = w, family = binomial(link = "logit"))
    print(summary(glmm_mod))
    print(car::vif(glmm_mod))
    
    return(glmm_mod)
  }
  #'  Run species, season, and year specific models through glmm function
  
  ####  Cougar RSFs  ####
  coug_smr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = cougDataz_smr) 
  coug_fll <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = cougDataz_fll) 
  coug_wtr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = cougDataz_wtr) 
  coug_sprg <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = cougDataz_sprg) 
  
  ####  Wolf RSFs  ####
  #'  NOTE: using 2nd reclassified version of landcover categories for wolf models
  wolf_smr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wolfDataz_smr) 
  wolf_fll <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wolfDataz_fll) 
  wolf_wtr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wolfDataz_wtr)
  wolf_sprg <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wolfDataz_sprg) 
  
  ####  Bobcat RSFs  ####
  bob_smr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = bobDataz_smr) 
  bob_fll <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = bobDataz_fll) 
  bob_wtr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = bobDataz_wtr) 
  bob_sprg <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = bobDataz_sprg) 
  
  ####  Coyote RSFs  ####
  coy_smr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)", dat = coyDataz_smr)  
  coy_fll <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)", dat = coyDataz_fll)  
  coy_wtr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)", dat = coyDataz_wtr) 
  coy_sprg <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)", dat = coyDataz_sprg)  
  
  #'  Group species-specific models
  RSF_COUG_list <- list(coug_smr, coug_fll, coug_wtr, coug_sprg)
  RSF_WOLF_list <- list(wolf_smr, wolf_fll, wolf_wtr, wolf_sprg)
  RSF_BOB_list <- list(bob_smr, bob_fll, bob_wtr, bob_sprg)
  RSF_COY_list <- list(coy_smr, coy_fll, coy_wtr, coy_sprg)
  
  #'  Save
  save(RSF_COUG_list, file = paste0("./Outputs/RSF_output/RSF_COUG_list_forTRG_", Sys.Date(), ".RData"))
  save(RSF_WOLF_list, file = paste0("./Outputs/RSF_output/RSF_WOLF_list_forTRG_", Sys.Date(), ".RData"))
  save(RSF_BOB_list, file = paste0("./Outputs/RSF_output/RSF_BOB_list_forTRG_", Sys.Date(), ".RData"))
  save(RSF_COY_list, file = paste0("./Outputs/RSF_output/RSF_COY_list_forTRG_", Sys.Date(), ".RData"))

  
  
  ####  Project RSF results across study areas  ####
  #'  ==============================================
  
  #'  Load RSFs
  load("./Outputs/RSF_output/RSF_COUG_list_forTRG_2023-07-07.RData") # _forTRG_2022-06-15
  load("./Outputs/RSF_output/RSF_WOLF_list_forTRG_2023-07-07.RData")
  load("./Outputs/RSF_output/RSF_BOB_list_forTRG_2023-07-07.RData")
  load("./Outputs/RSF_output/RSF_COY_list_forTRG_2023-07-07.RData")
  
  #'  Load spatial libraries
  library(sf)
  library(raster)
  
  #'  Define desired projections
  sa_proj <- projection("EPSG:2855")  # NAD83(HARN) / Washington North
  # sa_proj <- "+proj=lcc +lat_0=47 +lon_0=-120.833333333333 +lat_1=48.7333333333333 +lat_2=47.5 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs"
  
  #'  Read in study area grids
  # NE_1km <- raster("./Shapefiles/NE_1km_grid_mask.tif")
  # OK_1km <- raster("./Shapefiles/OK_1km_grid_mask.tif")
  NE_30m <- raster("./Shapefiles/NE_30m_grid_mask.tif") # MASSIVE
  OK_30m <- raster("./Shapefiles/OK_30m_grid_mask.tif") # MASSIVE
  
  plot(OK_30m)
  projection(OK_30m)
  
  #'  Load study area shapefiles
  OK.SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "METHOW_SA") %>%
    st_transform(crs = sa_proj) 
  OK.SA <- as(OK.SA, "Spatial")
  NE.SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "NE_SA") %>%
    st_transform(crs = sa_proj)
  NE.SA <- as(NE.SA, "Spatial")
  
  #'  Convert rasters to pixels and extract coordinates (centroid of each cell)
  #'  FYI: "data" are the grid cell IDs from the original WPPP reference grid;
  #'  grid.index are the cell IDs based on renumbered cells in cropped rasters.
  #'  Because some cells were masked out for large water bodies in both versions 
  #'  of these rasters, the gridID does not match the extracted study area-wide 
  #'  covariate df so need to create a new ID specific to masked grid
  raster_dat <- function(r) {
    dots <- as(r, "SpatialPixelsDataFrame")
    ref_grid_ID <- dots@data
    gridID <- dots@grid.index
    coords <- coordinates(dots)
    pts <- as.data.frame(cbind(ref_grid_ID, gridID, coords))
    pts$ID <- seq(1:nrow(pts))
    names(pts) <- c("ref_gridID", "gridID", "x", "y", "ID")
    return(pts)
  }
  # NE_pts <- raster_dat(NE_1km)
  # OK_pts <- raster_dat(OK_1km)
  NE_pts <- raster_dat(NE_30m)
  OK_pts <- raster_dat(OK_30m)
  
  #'  Read in covariates extracted across each study area 
  # load("./Outputs/Telemetry_covs/NE_covs_1km_2022-05-17.RData") 
  # load("./Outputs/Telemetry_covs/OK_covs_1km_2022-05-17.RData")
  load("./Outputs/Telemetry_covs/NE_covs_30m_2022-05-17.RData")
  load("./Outputs/Telemetry_covs/OK_covs_30m_2022-05-17.RData")
  
  #'  Format study area-wide covariate data to include annually relevant data only
  NE.covs <- NE.covs.30m %>%      # NE.covs.1km
    mutate(StudyArea = "NE") %>%
    full_join(NE_pts, by = "ID") %>%
    dplyr::select(-gridID) %>%
    #'  In case covariates were extracted at masked locations- drop these because 
    #'  missing coordinate data when joined (this shouldn't actually happen though)
    filter(!is.na(x))
  OK.covs <- OK.covs.30m %>%      # OK.covs.1km
    mutate(StudyArea = "OK") %>%
    full_join(OK_pts, by = "ID") %>%
    dplyr::select(-gridID) %>%
    filter(!is.na(x))
  SA.covs <- rbind(NE.covs, OK.covs)
  SA.covs.Year0 <- dplyr::select(SA.covs, -c(CanopyCover18, CanopyCover19, CanopyCover20, 
                                             Dist2Edge18, Dist2Edge19, Dist2Edge20, 
                                             Landcover_type18, Landcover_type19, Landcover_type20))
  names(SA.covs.Year0) <- c("ID", "Elev", "Slope", "RoadDen", "Dist2Water",
                            "HumanMod", "CanopyCover", "Dist2Edge", 
                            "Landcover_type", "StudyArea", "ref_gridID", "x", "y")
  SA.covs.Year1 <- dplyr::select(SA.covs, -c(CanopyCover17, CanopyCover19, CanopyCover20, 
                                             Dist2Edge17, Dist2Edge19, Dist2Edge20, 
                                             Landcover_type17, Landcover_type19, Landcover_type20))
  names(SA.covs.Year1) <- c("ID", "Elev", "Slope", "RoadDen", "Dist2Water",
                            "HumanMod", "CanopyCover", "Dist2Edge", 
                            "Landcover_type", "StudyArea", "ref_gridID", "x", "y")
  SA.covs.Year2 <- dplyr::select(SA.covs, -c(CanopyCover17, CanopyCover18, CanopyCover20, 
                                             Dist2Edge17, Dist2Edge18, Dist2Edge20, 
                                             Landcover_type17, Landcover_type18, Landcover_type20))
  names(SA.covs.Year2) <- c("ID", "Elev", "Slope", "RoadDen", "Dist2Water",
                            "HumanMod", "CanopyCover", "Dist2Edge", 
                            "Landcover_type", "StudyArea", "ref_gridID", "x", "y")
  SA.covs.Year3 <- dplyr::select(SA.covs, -c(CanopyCover17, CanopyCover18, CanopyCover19, 
                                             Dist2Edge17, Dist2Edge18, Dist2Edge19, 
                                             Landcover_type17, Landcover_type18, Landcover_type19))
  names(SA.covs.Year3) <- c("ID", "Elev", "Slope", "RoadDen", "Dist2Water",
                            "HumanMod", "CanopyCover", "Dist2Edge", 
                            "Landcover_type", "StudyArea", "ref_gridID", "x", "y")
  #'  List study area covariates by year to mirror rest of data structure
  SA.covs_list <- list(SA.covs.Year0, SA.covs.Year1, SA.covs.Year2, SA.covs.Year3)
  NE.covs_list <- list(SA.covs.Year0[SA.covs.Year0$StudyArea == "NE",], SA.covs.Year1[SA.covs.Year1$StudyArea == "NE",], SA.covs.Year2[SA.covs.Year2$StudyArea == "NE",], SA.covs.Year3[SA.covs.Year3$StudyArea == "NE",])
  OK.covs_list <- list(SA.covs.Year0[SA.covs.Year0$StudyArea == "OK",], SA.covs.Year1[SA.covs.Year1$StudyArea == "OK",], SA.covs.Year2[SA.covs.Year2$StudyArea == "OK",], SA.covs.Year3[SA.covs.Year3$StudyArea == "OK",])
  
  
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
  cougCov_smr_summary <- cov_summary(cougData_smr)
  cougCov_fll_summary <- cov_summary(cougData_fll)
  cougCov_wtr_summary <- cov_summary(cougData_wtr)
  cougCov_sprg_summary <- cov_summary(cougData_sprg)
  wolfCov_smr_summary <- cov_summary(wolfData_smr)
  wolfCov_fll_summary <- cov_summary(wolfData_fll)
  wolfCov_wtr_summary <- cov_summary(wolfData_wtr)
  wolfCov_sprg_summary <- cov_summary(wolfData_sprg)
  bobCov_smr_summary <- cov_summary(bobData_smr)
  bobCov_fll_summary <- cov_summary(bobData_fll)
  bobCov_wtr_summary <- cov_summary(bobData_wtr)
  bobCov_sprg_summary <- cov_summary(bobData_sprg)
  coyCov_smr_summary <- cov_summary(coyData_smr)
  coyCov_fll_summary <- cov_summary(coyData_fll)
  coyCov_wtr_summary <- cov_summary(coyData_wtr)
  coyCov_sprg_summary <- cov_summary(coyData_sprg)
  
  
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
        #'  Dummy variables for Landcover_type, Forest represents the intercept
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
  coug_smr_zcovs <- lapply(SA.covs_list_reclass, scaling_covs, mu.sd = cougCov_smr_summary)
  coug_fll_zcovs <- lapply(SA.covs_list_reclass, scaling_covs, mu.sd = cougCov_fll_summary)
  coug_wtr_zcovs <- lapply(SA.covs_list_reclass, scaling_covs, mu.sd = cougCov_wtr_summary)
  coug_sprg_zcovs <- lapply(SA.covs_list_reclass, scaling_covs, mu.sd = cougCov_sprg_summary)
  wolf_smr_zcovs <- lapply(SA.covs_list_wolf_reclass, scaling_covs, mu.sd = wolfCov_smr_summary)
  wolf_fll_zcovs <- lapply(SA.covs_list_wolf_reclass, scaling_covs, mu.sd = wolfCov_fll_summary)
  wolf_wtr_zcovs <- lapply(SA.covs_list_wolf_reclass, scaling_covs, mu.sd = wolfCov_wtr_summary)
  wolf_sprg_zcovs <- lapply(SA.covs_list_wolf_reclass, scaling_covs, mu.sd = wolfCov_sprg_summary)
  bob_smr_zcovs <- lapply(SA.covs_list_reclass, scaling_covs, mu.sd = bobCov_smr_summary)
  bob_fll_zcovs <- lapply(SA.covs_list_reclass, scaling_covs, mu.sd = bobCov_fll_summary)
  bob_wtr_zcovs <- lapply(SA.covs_list_reclass, scaling_covs, mu.sd = bobCov_wtr_summary)
  bob_sprg_zcovs <- lapply(SA.covs_list_reclass, scaling_covs, mu.sd = bobCov_sprg_summary)
  coy_smr_zcovs <- lapply(SA.covs_list_reclass, scaling_covs, mu.sd = coyCov_smr_summary)
  coy_fll_zcovs <- lapply(SA.covs_list_reclass, scaling_covs, mu.sd = coyCov_fll_summary)
  coy_wtr_zcovs <- lapply(SA.covs_list_reclass, scaling_covs, mu.sd = coyCov_wtr_summary)
  coy_sprg_zcovs <- lapply(SA.covs_list_reclass, scaling_covs, mu.sd = coyCov_sprg_summary)
  
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
      #'  Use p-values to change non-significant coefficients (alpha-level > 0.05) 
      #'  to 0 so there is no effect. Only exception is when Elev^2 is significant
      #'  but Elev is not. Coefficient for Elev still needs to be in the model.
      # mutate(Estimate = ifelse(Pval > 0.05 & Estimate == 0, Estimate == 0.00, Estimate)) %>%
      # mutate(Estimate = ifelse(Pval > 0.05 & Parameter != "Elev", Estimate == 0, Estimate)) %>%
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
    
    #'  Covariates excluded from species-specific models not included in the data
    #'  frame but necessary for predicting function to work below
    #'  Vector of columns names that need to be included in this data frame
    nms <- c("Species", "Season", "alpha", "b.elev", "b.elev2", "b.slope", "b.road", 
             "b.water", "b.hm", "b.canopy", "b.edge", "b.developed", "b.grass",  
             "b.other", "b.shrub", "b.wetland")
    #'  Identify if there are any missing column names in the data frame
    Missing <- setdiff(nms, names(out))
    #'  Add missing columns and fill with 0's
    out[Missing] <- 0
    #'  Return data frame based on full list of column names
    out <- out[nms]
    
    return(out)
  }
  #'  Extract coefficient estimates for each model (list: summer rsf [[1]], summer rsf [[2]], winter rsf [[3]], summer rsf [[4]])
  coug_smr_rsfout <- rsf_out(RSF_COUG_list[[1]], spp = "Cougar", season = "Summer")
  coug_fll_rsfout <- rsf_out(RSF_COUG_list[[2]], spp = "Cougar", season = "Fall")
  coug_wtr_rsfout <- rsf_out(RSF_COUG_list[[3]], spp = "Cougar", season = "Winter")
  coug_sprg_rsfout <- rsf_out(RSF_COUG_list[[4]], spp = "Cougar", season = "Spring")
  wolf_smr_rsfout <- rsf_out(RSF_WOLF_list[[1]], spp = "Wolf", season = "Summer")
  wolf_fll_rsfout <- rsf_out(RSF_WOLF_list[[2]], spp = "Wolf", season = "Fall")
  wolf_wtr_rsfout <- rsf_out(RSF_WOLF_list[[3]], spp = "Wolf", season = "Winter")
  wolf_sprg_rsfout <- rsf_out(RSF_WOLF_list[[4]], spp = "Wolf", season = "Spring")
  bob_smr_rsfout <- rsf_out(RSF_BOB_list[[1]], spp = "Bobcat", season = "Summer")
  bob_fll_rsfout <- rsf_out(RSF_BOB_list[[2]], spp = "Bobcat", season = "Fall")
  bob_wtr_rsfout <- rsf_out(RSF_BOB_list[[3]], spp = "Bobcat", season = "Winter")
  bob_sprg_rsfout <- rsf_out(RSF_BOB_list[[4]], spp = "Bobcat", season = "Spring")
  coy_smr_rsfout <- rsf_out(RSF_COY_list[[1]], spp = "Coyote", season = "Summer")
  coy_fll_rsfout <- rsf_out(RSF_COY_list[[2]], spp = "Coyote", season = "Fall")
  coy_wtr_rsfout <- rsf_out(RSF_COY_list[[3]], spp = "Coyote", season = "Winter")
  coy_sprg_rsfout <- rsf_out(RSF_COY_list[[4]], spp = "Coyote", season = "Spring")
  
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
                              coef$b.developed*cov$Landcover_Developed[i] + 
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
  coug_smr_rsf_sa <- lapply(coug_smr_zcovs, predict_rsf, coef = coug_smr_rsfout)
  coug_fll_rsf_sa <- lapply(coug_fll_zcovs, predict_rsf, coef = coug_fll_rsfout)
  coug_wtr_rsf_sa <- lapply(coug_wtr_zcovs, predict_rsf, coef = coug_wtr_rsfout)
  coug_sprg_rsf_sa <- lapply(coug_sprg_zcovs, predict_rsf, coef = coug_sprg_rsfout)
  wolf_smr_rsf_sa <- lapply(wolf_smr_zcovs, predict_rsf, coef = wolf_smr_rsfout)
  wolf_fll_rsf_sa <- lapply(wolf_fll_zcovs, predict_rsf, coef = wolf_fll_rsfout)
  wolf_wtr_rsf_sa <- lapply(wolf_wtr_zcovs, predict_rsf, coef = wolf_wtr_rsfout)
  wolf_sprg_rsf_sa <- lapply(wolf_sprg_zcovs, predict_rsf, coef = wolf_sprg_rsfout)
  bob_smr_rsf_sa <- lapply(bob_smr_zcovs, predict_rsf, coef = bob_smr_rsfout)
  bob_fll_rsf_sa <- lapply(bob_fll_zcovs, predict_rsf, coef = bob_fll_rsfout)
  bob_wtr_rsf_sa <- lapply(bob_wtr_zcovs, predict_rsf, coef = bob_wtr_rsfout)
  bob_sprg_rsf_sa <- lapply(bob_sprg_zcovs, predict_rsf, coef = bob_sprg_rsfout)
  coy_smr_rsf_sa <- lapply(coy_smr_zcovs, predict_rsf, coef = coy_smr_rsfout)
  coy_fll_rsf_sa <- lapply(coy_fll_zcovs, predict_rsf, coef = coy_fll_rsfout)
  coy_wtr_rsf_sa <- lapply(coy_wtr_zcovs, predict_rsf, coef = coy_wtr_rsfout)
  coy_sprg_rsf_sa <- lapply(coy_sprg_zcovs, predict_rsf, coef = coy_sprg_rsfout)
  
  #'  List and save
  all_spp_RSF_predicted <- list(coug_smr_rsf_sa, coug_fll_rsf_sa, coug_wtr_rsf_sa, coug_sprg_rsf_sa,  
                                wolf_smr_rsf_sa, wolf_fll_rsf_sa, wolf_wtr_rsf_sa, wolf_sprg_rsf_sa,
                                bob_smr_rsf_sa, bob_fll_rsf_sa, bob_wtr_rsf_sa, bob_sprg_rsf_sa,
                                coy_smr_rsf_sa, coy_fll_rsf_sa, coy_wtr_rsf_sa, coy_sprg_rsf_sa)
  save(all_spp_RSF_predicted, file = paste0("./Outputs/RSF_output/all_spp_RSF_predicted_forTRG_", Sys.Date(), ".RData"))
  
  load("./Outputs/RSF_output/all_spp_RSF_predicted_forTRG_2023-07-08.RData")  # _forTRG_2022-06-17
  
  #'  Function to identify any outliers
  outliers <- function(predicted, title, covs_list) {
    #' #'  Summarize predicted values
    #' hist(predicted$predict_rsf, breaks = 100, main = title)
    #' boxplot(predicted$predict_rsf, main = title)
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
    
    return(predicted)
  }
  #'  Identify outlier predictions and possible covariates associated with those
  #'  Be sure to used standardized covariates for evaluation (using 2018 covariate data)
  coug_smr_outliers <- lapply(all_spp_RSF_predicted[[1]], outliers, title = "Cougar Summer RSF Predictions", covs_list = coug_smr_zcovs[[2]])
  coug_fll_outliers <- lapply(all_spp_RSF_predicted[[2]], outliers, title = "Cougar Fall RSF Predictions", covs_list = coug_fll_zcovs[[2]])
  coug_wtr_outliers <- lapply(all_spp_RSF_predicted[[3]], outliers, title = "Cougar Winter RSF Predictions", covs_list = coug_wtr_zcovs[[2]])
  coug_sprg_outliers <- lapply(all_spp_RSF_predicted[[4]], outliers, title = "Cougar Spring RSF Predictions", covs_list = coug_sprg_zcovs[[2]])
  wolf_smr_outliers <- lapply(all_spp_RSF_predicted[[5]], outliers, title = "Wolf Summer RSF Predictions", covs_list = wolf_smr_zcovs[[2]])
  wolf_fll_outliers <- lapply(all_spp_RSF_predicted[[6]], outliers, title = "Wolf Fall RSF Predictions", covs_list = wolf_fll_zcovs[[2]])
  wolf_wtr_outliers <- lapply(all_spp_RSF_predicted[[7]], outliers, title = "Wolf Winter RSF Predictions", covs_list = wolf_wtr_zcovs[[2]])
  wolf_sprg_outliers <- lapply(all_spp_RSF_predicted[[8]], outliers, title = "Wolf Spring RSF Predictions", covs_list = wolf_sprg_zcovs[[2]])
  bob_smr_outliers <- lapply(all_spp_RSF_predicted[[9]], outliers, title = "Bobcat Summer RSF Predictions", covs_list = bob_smr_zcovs[[2]])
  bob_fll_outliers <- lapply(all_spp_RSF_predicted[[10]], outliers, title = "Bobcat Fall RSF Predictions", covs_list = bob_fll_zcovs[[2]])
  bob_wtr_outliers <- lapply(all_spp_RSF_predicted[[11]], outliers, title = "Bobcat Winter RSF Predictions", covs_list = bob_wtr_zcovs[[2]])
  bob_sprg_outliers <- lapply(all_spp_RSF_predicted[[12]], outliers, title = "Bobcat Spring RSF Predictions", covs_list = bob_sprg_zcovs[[2]])
  coy_smr_outliers <- lapply(all_spp_RSF_predicted[[13]], outliers, title = "Coyote Summer RSF Predictions", covs_list = coy_smr_zcovs[[2]])
  coy_fll_outliers <- lapply(all_spp_RSF_predicted[[14]], outliers, title = "Coyote Fall RSF Predictions", covs_list = coy_fll_zcovs[[2]])
  coy_wtr_outliers <- lapply(all_spp_RSF_predicted[[15]], outliers, title = "Coyote Winter RSF Predictions", covs_list = coy_wtr_zcovs[[2]])
  coy_sprg_outliers <- lapply(all_spp_RSF_predicted[[16]], outliers, title = "Coyote Spring RSF Predictions", covs_list = coy_sprg_zcovs[[2]])
  
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
  coug_smr_rescale_sa <- lapply(coug_smr_outliers, RSF_rescale)
  coug_fll_rescale_sa <- lapply(coug_fll_outliers, RSF_rescale)
  coug_wtr_rescale_sa <- lapply(coug_wtr_outliers, RSF_rescale)
  coug_sprg_rescale_sa <- lapply(coug_sprg_outliers, RSF_rescale)
  wolf_smr_rescale_sa <- lapply(wolf_smr_outliers, RSF_rescale)
  wolf_fll_rescale_sa <- lapply(wolf_fll_outliers, RSF_rescale)
  wolf_wtr_rescale_sa <- lapply(wolf_wtr_outliers, RSF_rescale)
  wolf_sprg_rescale_sa <- lapply(wolf_sprg_outliers, RSF_rescale)
  bob_smr_rescale_sa <- lapply(bob_smr_outliers, RSF_rescale) 
  bob_fll_rescale_sa <- lapply(bob_fll_outliers, RSF_rescale) 
  bob_wtr_rescale_sa <- lapply(bob_wtr_outliers, RSF_rescale)
  bob_sprg_rescale_sa <- lapply(bob_sprg_outliers, RSF_rescale) 
  coy_smr_rescale_sa <- lapply(coy_smr_outliers, RSF_rescale)  
  coy_fll_rescale_sa <- lapply(coy_fll_outliers, RSF_rescale)  
  coy_wtr_rescale_sa <- lapply(coy_wtr_outliers, RSF_rescale) 
  coy_sprg_rescale_sa <- lapply(coy_sprg_outliers, RSF_rescale)  
  
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
  coug_smr_RSFraster <- lapply(coug_smr_rescale_sa, rasterize_rsf)
  coug_fll_RSFraster <- lapply(coug_fll_rescale_sa, rasterize_rsf)
  coug_wtr_RSFraster <- lapply(coug_wtr_rescale_sa, rasterize_rsf)
  coug_sprg_RSFraster <- lapply(coug_sprg_rescale_sa, rasterize_rsf) 
  wolf_smr_RSFraster <- lapply(wolf_smr_rescale_sa, rasterize_rsf)
  wolf_fll_RSFraster <- lapply(wolf_fll_rescale_sa, rasterize_rsf)
  wolf_wtr_RSFraster <- lapply(wolf_wtr_rescale_sa, rasterize_rsf)
  wolf_sprg_RSFraster <- lapply(wolf_sprg_rescale_sa, rasterize_rsf)
  bob_smr_RSFraster <- lapply(bob_smr_rescale_sa, rasterize_rsf)
  bob_fll_RSFraster <- lapply(bob_fll_rescale_sa, rasterize_rsf) 
  bob_wtr_RSFraster <- lapply(bob_wtr_rescale_sa, rasterize_rsf)
  bob_sprg_RSFraster <- lapply(bob_sprg_rescale_sa, rasterize_rsf)
  coy_smr_RSFraster <- lapply(coy_smr_rescale_sa, rasterize_rsf)
  coy_fll_RSFraster <- lapply(coy_fll_rescale_sa, rasterize_rsf)
  coy_wtr_RSFraster <- lapply(coy_wtr_rescale_sa, rasterize_rsf)
  coy_sprg_RSFraster <- lapply(coy_sprg_rescale_sa, rasterize_rsf)
  
  #'  Rename rasters
  rename_raster <- function(raster_list) {
    L <- setNames(raster_list, c("Year1", "Year2", "Year3"))
    S <- stack(L)
    return(S)
  }
  coug_smr_RSFstack <- rename_raster(coug_smr_RSFraster)
  coug_fll_RSFstack <- rename_raster(coug_fll_RSFraster)
  coug_wtr_RSFstack <- rename_raster(coug_wtr_RSFraster)
  coug_sprg_RSFstack <- rename_raster(coug_sprg_RSFraster)
  wolf_smr_RSFstack <- rename_raster(wolf_smr_RSFraster)
  wolf_fll_RSFstack <- rename_raster(wolf_fll_RSFraster)
  wolf_wtr_RSFstack <- rename_raster(wolf_wtr_RSFraster)
  wolf_sprg_RSFstack <- rename_raster(wolf_sprg_RSFraster)
  bob_smr_RSFstack <- rename_raster(bob_smr_RSFraster)
  bob_fll_RSFstack <- rename_raster(bob_fll_RSFraster)
  bob_wtr_RSFstack <- rename_raster(bob_wtr_RSFraster)
  bob_sprg_RSFstack <- rename_raster(bob_sprg_RSFraster)
  coy_smr_RSFstack <- rename_raster(coy_smr_RSFraster)
  coy_fll_RSFstack <- rename_raster(coy_fll_RSFraster)
  coy_wtr_RSFstack <- rename_raster(coy_wtr_RSFraster)
  coy_sprg_RSFstack <- rename_raster(coy_sprg_RSFraster)
  
  #'  Plot & Save
  pdf(file = "./Outputs/RSF_output/RSF_Maps_forTRG.pdf")
  plot(coug_smr_RSFstack[[2]], main = "Summer Cougar Predicted RSF (2018)"); plot(OK.SA, add = T); plot(NE.SA, add = T)
  plot(coug_fll_RSFstack[[2]], main = "Fall Cougar Predicted RSF (2018)"); plot(OK.SA, add = T); plot(NE.SA, add = T)
  plot(coug_wtr_RSFstack[[2]], main = "Winter Cougar Predicted RSF (2018)"); plot(OK.SA, add = T); plot(NE.SA, add = T)
  plot(coug_sprg_RSFstack[[2]], main = "Spring Cougar Predicted RSF (2018)"); plot(OK.SA, add = T); plot(NE.SA, add = T)
  plot(wolf_smr_RSFstack[[2]], main = "Summer Wolf Predicted RSF (2018)"); plot(OK.SA, add = T); plot(NE.SA, add = T)
  plot(wolf_fll_RSFstack[[2]], main = "Fall Wolf Predicted RSF (2018)"); plot(OK.SA, add = T); plot(NE.SA, add = T)
  plot(wolf_wtr_RSFstack[[2]], main = "Winter Wolf Predicted RSF (2018)"); plot(OK.SA, add = T); plot(NE.SA, add = T)
  plot(wolf_sprg_RSFstack[[2]], main = "Spring Wolf Predicted RSF (2018)"); plot(OK.SA, add = T); plot(NE.SA, add = T)
  plot(bob_smr_RSFstack[[2]], main = "Summer Bobcat Predicted RSF (2018)"); plot(OK.SA, add = T); plot(NE.SA, add = T)
  plot(bob_fll_RSFstack[[2]], main = "Fall Bobcat Predicted RSF (2018)"); plot(OK.SA, add = T); plot(NE.SA, add = T)
  plot(bob_wtr_RSFstack[[2]], main = "Winter Bobcat Predicted RSF (2018)"); plot(OK.SA, add = T); plot(NE.SA, add = T)
  plot(bob_sprg_RSFstack[[2]], main = "Spring Bobcat Predicted RSF (2018)"); plot(OK.SA, add = T); plot(NE.SA, add = T)
  plot(coy_smr_RSFstack[[2]], main = "Summer Coyote Predicted RSF (2018)"); plot(OK.SA, add = T); plot(NE.SA, add = T)
  plot(coy_fll_RSFstack[[2]], main = "Fall Coyote Predicted RSF (2018)"); plot(OK.SA, add = T); plot(NE.SA, add = T)
  plot(coy_wtr_RSFstack[[2]], main = "Winter Coyote Predicted RSF (2018)"); plot(OK.SA, add = T); plot(NE.SA, add = T)
  plot(coy_sprg_RSFstack[[2]], main = "Spring Coyote Predicted RSF (2018)"); plot(OK.SA, add = T); plot(NE.SA, add = T)
  dev.off()
  
  #'  Convert to rast format for terra
  library(terra)
  coug_smr_RSFstack <- rast(coug_smr_RSFstack)
  coug_fll_RSFstack <- rast(coug_fll_RSFstack)
  coug_wtr_RSFstack <- rast(coug_wtr_RSFstack)
  coug_sprg_RSFstack <- rast(coug_sprg_RSFstack)
  wolf_smr_RSFstack <- rast(wolf_smr_RSFstack)
  wolf_fll_RSFstack <- rast(wolf_fll_RSFstack)
  wolf_wtr_RSFstack <- rast(wolf_wtr_RSFstack)
  wolf_sprg_RSFstack <- rast(wolf_sprg_RSFstack)
  bob_smr_RSFstack <- rast(bob_smr_RSFstack)
  bob_fll_RSFstack <- rast(bob_fll_RSFstack)
  bob_wtr_RSFstack <- rast(bob_wtr_RSFstack)
  bob_sprg_RSFstack <- rast(bob_sprg_RSFstack)
  coy_smr_RSFstack <- rast(coy_smr_RSFstack)
  coy_fll_RSFstack <- rast(coy_fll_RSFstack)
  coy_wtr_RSFstack <- rast(coy_wtr_RSFstack)
  coy_sprg_RSFstack <- rast(coy_sprg_RSFstack)
  
  #'  SAVE!
  terra::writeRaster(coug_smr_RSFstack, filename = "./Shapefiles/Predicted_RSFs/Predicted_Predator_RSF_forTRG/coug_smr_RSFstack_global.tif", overwrite = TRUE)
  terra::writeRaster(coug_fll_RSFstack, filename = "./Shapefiles/Predicted_RSFs/Predicted_Predator_RSF_forTRG/coug_fll_RSFstack_global.tif", overwrite = TRUE)
  terra::writeRaster(coug_wtr_RSFstack, filename = "./Shapefiles/Predicted_RSFs/Predicted_Predator_RSF_forTRG/coug_wtr_RSFstack_global.tif", overwrite = TRUE)
  terra::writeRaster(coug_sprg_RSFstack, filename = "./Shapefiles/Predicted_RSFs/Predicted_Predator_RSF_forTRG/coug_sprg_RSFstack_global.tif", overwrite = TRUE)  
  terra::writeRaster(wolf_smr_RSFstack, filename = "./Shapefiles/Predicted_RSFs/Predicted_Predator_RSF_forTRG/wolf_smr_RSFstack_global.tif", overwrite = TRUE)
  terra::writeRaster(wolf_fll_RSFstack, filename = "./Shapefiles/Predicted_RSFs/Predicted_Predator_RSF_forTRG/wolf_fll_RSFstack_global.tif", overwrite = TRUE)
  terra::writeRaster(wolf_wtr_RSFstack, filename = "./Shapefiles/Predicted_RSFs/Predicted_Predator_RSF_forTRG/wolf_wtr_RSFstack_global.tif", overwrite = TRUE)
  terra::writeRaster(wolf_sprg_RSFstack, filename = "./Shapefiles/Predicted_RSFs/Predicted_Predator_RSF_forTRG/wolf_sprg_RSFstack_global.tif", overwrite = TRUE)
  terra::writeRaster(bob_smr_RSFstack, filename = "./Shapefiles/Predicted_RSFs/Predicted_Predator_RSF_forTRG/bob_smr_RSFstack_global.tif", overwrite = TRUE)
  terra::writeRaster(bob_fll_RSFstack, filename = "./Shapefiles/Predicted_RSFs/Predicted_Predator_RSF_forTRG/bob_fll_RSFstack_global.tif", overwrite = TRUE)  
  terra::writeRaster(bob_wtr_RSFstack, filename = "./Shapefiles/Predicted_RSFs/Predicted_Predator_RSF_forTRG/bob_wtr_RSFstack_global.tif", overwrite = TRUE)
  terra::writeRaster(bob_sprg_RSFstack, filename = "./Shapefiles/Predicted_RSFs/Predicted_Predator_RSF_forTRG/bob_sprg_RSFstack_global.tif", overwrite = TRUE)
  terra::writeRaster(coy_smr_RSFstack, filename = "./Shapefiles/Predicted_RSFs/Predicted_Predator_RSF_forTRG/coy_smr_RSFstack_global.tif", overwrite = TRUE)
  terra::writeRaster(coy_fll_RSFstack, filename = "./Shapefiles/Predicted_RSFs/Predicted_Predator_RSF_forTRG/coy_fll_RSFstack_global.tif", overwrite = TRUE)
  terra::writeRaster(coy_wtr_RSFstack, filename = "./Shapefiles/Predicted_RSFs/Predicted_Predator_RSF_forTRG/coy_wtr_RSFstack_global.tif", overwrite = TRUE)
  terra::writeRaster(coy_sprg_RSFstack, filename = "./Shapefiles/Predicted_RSFs/Predicted_Predator_RSF_forTRG/coy_sprg_RSFstack_global.tif", overwrite = TRUE)
  
  #' ####  Summary tables  ####
  #'  Save model outputs in table format
  #'  Function to save parameter estimates & p-values
  #'  use coef(mod) to look at random effects estimates
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
        Pval = round(pval, 3))
    rownames(out) <- NULL
    return(out)
  }
  coug_smr_rsf_out <- rsf_out(RSF_COUG_list[[1]], "Cougar", "Summer")
  coug_fll_rsf_out <- rsf_out(RSF_COUG_list[[2]], "Cougar", "Fall")
  coug_wtr_rsf_out <- rsf_out(RSF_COUG_list[[3]], "Cougar", "Winter")
  coug_sprg_rsf_out <- rsf_out(RSF_COUG_list[[4]], "Cougar", "Spring")
  wolf_smr_rsf_out <- rsf_out(RSF_WOLF_list[[1]], "Wolf", "Summer")
  wolf_fll_rsf_out <- rsf_out(RSF_WOLF_list[[2]], "Wolf", "Fall")
  wolf_wtr_rsf_out <- rsf_out(RSF_WOLF_list[[3]], "Wolf", "Winter")
  wolf_sprg_rsf_out <- rsf_out(RSF_WOLF_list[[4]], "Wolf", "Spring")
  bob_smr_rsf_out <- rsf_out(RSF_BOB_list[[1]], "Bobcat", "Summer")
  bob_fll_rsf_out <- rsf_out(RSF_BOB_list[[2]], "Bobcat", "Fall")
  bob_wtr_rsf_out <- rsf_out(RSF_BOB_list[[3]], "Bobcat", "Winter")
  bob_sprg_rsf_out <- rsf_out(RSF_BOB_list[[4]], "Bobcat", "Spring")
  coy_smr_rsf_out <- rsf_out(RSF_COY_list[[1]], "Coyote", "Summer")
  coy_fll_rsf_out <- rsf_out(RSF_COY_list[[2]], "Coyote", "Fall")
  coy_wtr_rsf_out <- rsf_out(RSF_COY_list[[3]], "Coyote", "Winter")
  coy_sprg_rsf_out <- rsf_out(RSF_COY_list[[4]], "Coyote", "Spring")
  
  #'  Merge into larger data frames for easy comparison
  summer_rsf <- rbind(coug_smr_rsf_out, wolf_smr_rsf_out, bob_smr_rsf_out, coy_smr_rsf_out)
  fall_rsf <- rbind(coug_fll_rsf_out, wolf_fll_rsf_out, bob_fll_rsf_out, coy_fll_rsf_out)
  winter_rsf <- rbind(coug_wtr_rsf_out, wolf_wtr_rsf_out, bob_wtr_rsf_out, coy_wtr_rsf_out)
  spring_rsf <- rbind(coug_sprg_rsf_out, wolf_sprg_rsf_out, bob_sprg_rsf_out, coy_sprg_rsf_out)
  rsf_results <- rbind(summer_rsf, fall_rsf, winter_rsf, spring_rsf) %>%
    arrange(Species)
  colnames(rsf_results) <- c("Species", "Season", "Parameter", "Estimate", "SE", "z", "Pval")
  
  
  #'  Spread this out so the coefficient effects are easier to compare across species
  rsf_results_wide <- rsf_results %>%
    dplyr::select(-z) %>%
    mutate(
      SE = round(SE, 2),
      SE = paste0("(", SE, ")")
    ) %>%
    #'  Bold significant variables- doesn't work if continue manipulating data frame
    # condformat(.) %>%
    # rule_text_bold(c(Estimate, SE, Pval), expression = Pval <= 0.05) %>%
    unite(Est_SE, Estimate, SE, sep = " ") %>%
    unite(Est_SE_Pval, Est_SE, Pval, sep = "_") %>%
    spread(Parameter, Est_SE_Pval) %>%
    separate("(Intercept)", c("Intercept (SE)", "Intercept Pval"), sep = "_") %>%
    # separate("AreaOK", c("AreaOK (SE)", "AreaOK Pval"), sep = "_") %>%
    separate("Elev", c("Elev (SE)", "Elev Pval"), sep = "_") %>%
    separate("I(Elev^2)", c("I(Elev^2) (SE)", "I(Elev^2) Pval"), sep = "_") %>%
    separate("Slope", c("Slope (SE)", "Slope Pval"), sep = "_") %>%
    separate("RoadDen", c("RoadDen (SE)", "RoadDen Pval"), sep = "_") %>%
    separate("Dist2Water", c("Dist2Water (SE)", "Dist2Water Pval"), sep = "_") %>%
    separate("CanopyCover", c("CanopyCover (SE)", "CanopyCover Pval"), sep = "_") %>%
    separate("Dist2Edge", c("Dist2Edge (SE)", "Dist2Edge Pval"), sep = "_") %>%
    #separate("Landcover_typeDeveloped", c("Landcover_typeDeveloped (SE)", "Landcover_typeDeveloped Pval"), sep = "_") %>%
    separate("Landcover_typeOpen Grass", c("Landcover_typeOpen Grass (SE)", "Landcover_typeOpen Grass Pval"), sep = "_") %>%
    separate("Landcover_typeOther", c("Landcover_typeOther (SE)", "Landcover_typeOther Pval"), sep = "_") %>%
    separate("Landcover_typeShrub Mix", c("Landcover_typeShrub Mix (SE)", "Landcover_typeShrub Mix Pval"), sep = "_") %>%
    separate("Landcover_typeWetland", c("Landcover_typeWetland (SE)", "Landcover_typeWetland Pval"), sep = "_") %>%
    arrange(match(Species, c("Cougar", "Wolf", "Bobcat", "Coyote"))) %>%
    arrange(match(Season, c("Summer", "Fall", "Winter", "Spring")))
  
  
  #'  Save!
  write.csv(rsf_results, paste0("./Outputs/RSF_output/RSF_Results_forTRG_", Sys.Date(), ".csv"))  
  write.csv(rsf_results_wide, paste0("./Outputs/RSF_output/RSF_Results_wide_forTRG_", Sys.Date(), ".csv"))
  
  
