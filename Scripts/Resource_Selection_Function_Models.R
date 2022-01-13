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
  #'  Based on Taylor's feedback
  #'  Open grass: mesic grass, xeric grass, wetland woody
  #'  Shrubby mix: mesic shrub, xeric shrub
  #'  Forest
  #'  Wetland
  #'  Other: water, barren, glacier
  #'  Developed: agriculture, commercial, developed
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
  #'  Note: standardizing across all IDs but separately by species, season, & year
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
  # mdData_smr <- list(md_dat_all[md_dat_all$Season == "Summer18",], md_dat_all[md_dat_all$Season == "Summer19",], md_dat_all[md_dat_all$Season == "Summer20",])
  # mdData_wtr <- list(md_dat_all[md_dat_all$Season == "Winter1819",], md_dat_all[md_dat_all$Season == "Winter1920",], md_dat_all[md_dat_all$Season == "Winter2021",])
  # mdData_smr <- lapply(mdData_smr, standardize_covs)
  # mdData_wtr <- lapply(mdData_wtr, standardize_covs)
  #'  Note the reclassified landcover_type data for elkData_winter
  elkData_smr <- elk_dat_all[elk_dat_all$Season == "Summer18" | elk_dat_all$Season == "Summer19" | elk_dat_all$Season == "Summer20",]
  elkData_wtr <- elk_dat_all_reclass[elk_dat_all_reclass$Season == "Winter1819" | elk_dat_all_reclass$Season == "Winter1920" | elk_dat_all_reclass$Season == "Winter2021",]
  elkDataz_smr <- standardize_covs(elkData_smr)
  elkDataz_wtr <- standardize_covs(elkData_wtr)
  # elkData_smr <- list(elk_dat_all[elk_dat_all$Season == "Summer18",], elk_dat_all[elk_dat_all$Season == "Summer19",], elk_dat_all[elk_dat_all$Season == "Summer20",])
  # elkData_wtr <- list(elk_dat_all_reclass[elk_dat_all_reclass$Season == "Winter1819",], elk_dat_all_reclass[elk_dat_all_reclass$Season == "Winter1920",], elk_dat_all_reclass[elk_dat_all_reclass$Season == "Winter2021",])
  # elkData_smr <- lapply(elkData_smr, standardize_covs)
  # elkData_wtr <- lapply(elkData_wtr, standardize_covs)
  #'  Note the reclassified landcover_type data for wtdData_winter 
  wtdData_smr <- wtd_dat_all[wtd_dat_all$Season == "Summer18" | wtd_dat_all$Season == "Summer19" | wtd_dat_all$Season == "Summer20",]
  wtdData_wtr <- wtd_dat_all_reclass[wtd_dat_all_reclass$Season == "Winter1819" | wtd_dat_all_reclass$Season == "Winter1920" | wtd_dat_all_reclass$Season == "Winter2021",]
  wtdDataz_smr <- standardize_covs(wtdData_smr)
  wtdDataz_wtr <- standardize_covs(wtdData_wtr)
  # wtdData_smr <- list(wtd_dat_all[wtd_dat_all$Season == "Summer18",], wtd_dat_all[wtd_dat_all$Season == "Summer19",], wtd_dat_all[wtd_dat_all$Season == "Summer20",])
  # wtdData_wtr <- list(wtd_dat_all_reclass[wtd_dat_all_reclass$Season == "Winter1819",], wtd_dat_all_reclass[wtd_dat_all_reclass$Season == "Winter1920",], wtd_dat_all_reclass[wtd_dat_all_reclass$Season == "Winter2021",])
  # wtdData_smr <- lapply(wtdData_smr, standardize_covs)
  # wtdData_wtr <- lapply(wtdData_wtr, standardize_covs)
  #'  Note the reclassified landcover_type data for cougData_winter
  cougData_smr <- coug_dat_all[coug_dat_all$Season == "Summer18" | coug_dat_all$Season == "Summer19" | coug_dat_all$Season == "Summer20",]
  cougData_wtr <- coug_dat_all_reclass[coug_dat_all_reclass$Season == "Winter1819" | coug_dat_all_reclass$Season == "Winter1920" | coug_dat_all_reclass$Season == "Winter2021",]
  cougDataz_smr <- standardize_covs(cougData_smr)
  cougDataz_wtr <- standardize_covs(cougData_wtr)
  # cougData_smr <- list(coug_dat_all[coug_dat_all$Season == "Summer18",], coug_dat_all[coug_dat_all$Season == "Summer19",], coug_dat_all[coug_dat_all$Season == "Summer20",])
  # cougData_wtr <- list(coug_dat_all_reclass[coug_dat_all_reclass$Season == "Winter1819",], coug_dat_all_reclass[coug_dat_all_reclass$Season == "Winter1920",], coug_dat_all_reclass[coug_dat_all_reclass$Season == "Winter2021",])
  # cougData_smr <- lapply(cougData_smr, standardize_covs)
  # cougData_wtr <- lapply(cougData_wtr, standardize_covs)
  #'  Note the double reclassified landcover_type data for wolfData
  wolfData_smr <- wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Summer18" | wolf_dat_all_reclass2$Season == "Summer19" | wolf_dat_all_reclass2$Season == "Summer20",]
  wolfData_wtr <- wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Winter1819" | wolf_dat_all_reclass2$Season == "Winter1920" | wolf_dat_all_reclass2$Season == "Winter2021",]
  wolfDataz_smr <- standardize_covs(wolfData_smr)
  wolfDataz_wtr <- standardize_covs(wolfData_wtr)
  # wolfData_smr <- list(wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Summer18",], wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Summer19",], wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Summer20",])
  # wolfData_wtr <- list(wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Winter1819",], wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Winter1920",], wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Winter2021",])
  # wolfData_smr <- lapply(wolfData_smr, standardize_covs)
  # wolfData_wtr <- lapply(wolfData_wtr, standardize_covs)
  #'  Note the reclassified landcover_type data for bobData
  bobData_smr <- bob_dat_all_reclass[bob_dat_all_reclass$Season == "Summer18" | bob_dat_all_reclass$Season == "Summer19" | bob_dat_all_reclass$Season == "Summer20",]
  bobData_wtr <- bob_dat_all_reclass[bob_dat_all_reclass$Season == "Winter1819" | bob_dat_all_reclass$Season == "Winter1920" | bob_dat_all_reclass$Season == "Winter2021",]
  bobDataz_smr <- standardize_covs(bobData_smr)
  bobDataz_wtr <- standardize_covs(bobData_wtr)
  # bobData_smr <- list(bob_dat_all_reclass[bob_dat_all_reclass$Season == "Summer18",], bob_dat_all_reclass[bob_dat_all_reclass$Season == "Summer19",], bob_dat_all_reclass[bob_dat_all_reclass$Season == "Summer20",])
  # bobData_wtr <- list(bob_dat_all_reclass[bob_dat_all_reclass$Season == "Winter1819",], bob_dat_all_reclass[bob_dat_all_reclass$Season == "Winter1920",], bob_dat_all_reclass[bob_dat_all_reclass$Season == "Winter2021",])
  # bobData_smr <- lapply(bobData_smr, standardize_covs)
  # bobData_wtr <- lapply(bobData_wtr, standardize_covs)
  #'  Note the reclassified landcover_type data for coyData
  coyData_smr <- coy_dat_all_reclass[coy_dat_all_reclass$Season == "Summer18" | coy_dat_all_reclass$Season == "Summer19" | coy_dat_all_reclass$Season == "Summer20",]
  coyData_wtr <- coy_dat_all_reclass[coy_dat_all_reclass$Season == "Winter1819" | coy_dat_all_reclass$Season == "Winter1920" | coy_dat_all_reclass$Season == "Winter2021",]
  coyDataz_smr <- standardize_covs(coyData_smr)
  coyDataz_wtr <- standardize_covs(coyData_wtr)
  # coyData_smr <- list(coy_dat_all_reclass[coy_dat_all_reclass$Season == "Summer18",], coy_dat_all_reclass[coy_dat_all_reclass$Season == "Summer19",], coy_dat_all_reclass[coy_dat_all_reclass$Season == "Summer20",])
  # coyData_wtr <- list(coy_dat_all_reclass[coy_dat_all_reclass$Season == "Winter1819",], coy_dat_all_reclass[coy_dat_all_reclass$Season == "Winter1920",], coy_dat_all_reclass[coy_dat_all_reclass$Season == "Winter2021",])
  # coyData_smr <- lapply(coyData_smr, standardize_covs)
  # coyData_wtr <- lapply(coyData_wtr, standardize_covs)
  
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
  # (md_smr_corr <- lapply(mdData_smr, cov_correlation))
  # (md_wtr_corr <- lapply(mdData_wtr, cov_correlation))
  # (elk_smr_corr <- lapply(elkData_smr, cov_correlation))
  # (elk_wtr_corr <- lapply(elkData_wtr, cov_correlation))
  # (wtd_smr_corr <- lapply(wtdData_smr, cov_correlation))
  # (wtd_wtr_corr <- lapply(wtdData_wtr, cov_correlation))
  # (coug_smr_corr <- lapply(cougData_smr, cov_correlation))
  # (coug_wtr_corr <- lapply(cougData_wtr, cov_correlation))
  # (wolf_smr_corr <- lapply(wolfData_smr, cov_correlation))
  # (wolf_wtr_corr <- lapply(wolfData_wtr, cov_correlation))
  # (bob_smr_corr <- lapply(bobData_smr, cov_correlation))
  # (bob_wtr_corr <- lapply(bobData_wtr, cov_correlation))
  # (coy_smr_corr <- lapply(coyData_smr, cov_correlation))
  # (coy_wtr_corr <- lapply(coyData_wtr, cov_correlation))
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
  #'  above) and convergence issues. Annual models run separately so predicted 
  #'  distributions are specific to the species, season, and year. Landcover_type
  #'  reclassified for some species where too few observations occurred in "Other"
  #'  category- "Other" and "Developed" combined in these cases (input data labeled
  #'  with _reclass to distinguish the differences in landcover classifications).
  
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
  md_wtr <- glmm_fn(mod = "Used ~ 1 + Elev + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)", dat = mdDataz_wtr) # + I(Elev^2)
  #' #'  Dropping HumanMod in mulie models due to high correlation with other covariates
  #' md_smr18 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + RoadDen + CanopyCover + Dist2Edge + Landcover_type + (1|ID)", dat = mdData_smr[[1]]) # + Slope  + Dist2Water
  #' md_smr19 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Landcover_type + (1|ID)", dat = mdData_smr[[2]]) # + Dist2Edge
  #' md_smr20 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + RoadDen + Dist2Water + CanopyCover + Landcover_type + (1|ID)", dat = mdData_smr[[3]]) # + Dist2Edge + Slope
  #' md_wtr1819 <- glmm_fn(mod = "Used ~ 1 + Elev + Slope + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)", dat = mdData_wtr[[1]]) # + I(Elev^2) + RoadDen
  #' md_wtr1920 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)", dat = mdData_wtr[[2]])
  #' md_wtr2021 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)", dat = mdData_wtr[[3]])
  
  ####  Elk RSFs  ####
  #'  Dropping HumanMod in elk summer models due to high correlation with other covariates
  elk_smr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = elkDataz_smr)
  #'  Note: using reclassified version of landcover for winter elk models
  elk_wtr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = elkDataz_wtr)
  #' #'  Dropping HumanMod in elk summer models due to high correlation with other covariates
  #' elk_smr18 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = elkData_smr[[1]])
  #' elk_smr19 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = elkData_smr[[2]]) 
  #' elk_smr20 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = elkData_smr[[3]])
  #' #'  Note: using reclassified version of landcover for winter elk models
  #' elk_wtr1819 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = elkData_wtr[[1]])
  #' elk_wtr1920 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = elkData_wtr[[2]]) # + HumanMod
  #' elk_wtr2021 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = elkData_wtr[[3]])
  
  ####  White-tailed Deer RSFs  ####
  wtd_smr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + Dist2Edge + Landcover_type + (1|ID)",  dat = wtdDataz_smr) # + CanopyCover
  #'  Note: using reclassified version of landcover for winter WTD models
  wtd_wtr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wtdDataz_wtr)
  #' #'  Dropping HumanMod in wtd summer models due to high correlation with other covariates
  #' wtd_smr18 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Edge + Landcover_type + (1|ID)",  dat = wtdData_smr[[1]]) # + Dist2Water
  #' wtd_smr19 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wtdData_smr[[2]])
  #' wtd_smr20 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wtdData_smr[[3]])
  #' #'  Note: using reclassified version of landcover for winter WTD models
  #' wtd_wtr1819 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + Dist2Edge + Landcover_type + (1|ID)",  dat = wtdData_wtr[[1]]) # + Dist2Water  + CanopyCover + RoadDen + HumanMod
  #' wtd_wtr1920 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wtdData_wtr[[2]])
  #' wtd_wtr2021 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wtdData_wtr[[3]])
  
  ####  Cougar RSFs  ####
  #'  Dropping HumanMod in cougar summer models due to high correlation with other covariates
  coug_smr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = cougDataz_smr) 
  #'  Note: using reclassified version of landcover for winter cougar models
  coug_wtr <- glmm_fn(mod = "Used ~ 1 + Elev + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Landcover_type + (1|ID)",  dat = cougDataz_wtr) # + I(Elev^2) + Dist2Edge
  #' coug_smr18 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = cougData_smr[[1]])  
  #' coug_smr19 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = cougData_smr[[2]])
  #' coug_smr20 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = cougData_smr[[3]]) 
  #' #'  Note: using reclassified version of landcover for winter cougar models
  #' coug_wtr1819 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = cougData_wtr[[1]]) # + RoadDen
  #' coug_wtr1920 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = cougData_wtr[[2]]) # + RoadDen
  #' coug_wtr2021 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = cougData_wtr[[3]]) # + RoadDen  + Dist2Water
  
  ####  Wolf RSFs  ####
  #'  NOTE: using 2nd reclassified version of landcover categories for wolf models
  wolf_smr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + Dist2Water + HumanMod + Dist2Edge + Landcover_type + (1|ID)",  dat = wolfDataz_smr) # + RoadDen + CanopyCover
  wolf_wtr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wolfDataz_wtr)
  # wolf_smr18 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + Dist2Water + Dist2Edge + Landcover_type + (1|ID)",  dat = wolfData_smr[[1]]) #+ HumanMod  + RoadDen  + CanopyCover
  # wolf_smr19 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wolfData_smr[[2]])
  # wolf_smr20 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wolfData_smr[[3]])
  # wolf_wtr1819 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wolfData_wtr[[1]])  # + Dist2Water 
  # wolf_wtr1920 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + HumanMod + Dist2Edge + Landcover_type + (1|ID)",  dat = wolfData_wtr[[2]])  # + Dist2Water  + CanopyCover
  # wolf_wtr2021 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wolfData_wtr[[3]])
  
  ####  Bobcat RSFs  ####
  #'  Note: using reclassified version of landcover for summer bobcat models
  bob_smr <- glmm_fn(mod = "Used ~ 1 + Elev + Slope + RoadDen + Dist2Water + HumanMod + Dist2Edge + Landcover_type + (1|ID)",  dat = bobDataz_smr) # + I(Elev^2) + CanopyCover
  bob_wtr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + HumanMod + CanopyCover + Landcover_type + (1|ID)",  dat = bobDataz_wtr) # + Dist2Water + Dist2Edge
  #' #'  Only data for MVBOB90M in smr18--- not enough data to make inference about bobcat resource selection across 2 study areas
  #' # bob_smr18 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = bobData_smr[[1]])
  #' bob_smr19 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + HumanMod + Dist2Edge + Landcover_type + (1|ID)",  dat = bobData_smr[[2]]) # + CanopyCover + Dist2Water
  #' bob_smr20 <- glmm_fn(mod = "Used ~ 1 + RoadDen + Dist2Water + HumanMod + Landcover_type + (1|ID)",  dat = bobData_smr[[3]])  # + Slope + CanopyCover + I(Elev^2) + Dist2Edge + Elev
  #' #'  Only data for MVBOB88M & MVBOB90M wtr1819--- not enough data to make inference about bobcat resource selection across 2 study areas
  #' # bob_wtr1819 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = bobData_wtr[[1]])
  #' bob_wtr1920 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = bobData_wtr[[2]]) #  + Dist2Water
  #' bob_wtr2021 <- glmm_fn(mod = "Used ~ 1 + Elev + Slope + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = bobData_wtr[[3]]) # + Dist2Water + RoadDen + I(Elev^2)
  
  ####  Coyote RSFs  ####
  #'  Note: using reclassified version of landcover for all coyote models
  #'  Dropping HumanMod from all coyote models due to high correlation with other covariates
  coy_smr <- glmm_fn(mod = "Used ~ 1 + Slope + RoadDen + Dist2Water + CanopyCover + Landcover_type + (1|ID)", dat = coyDataz_smr)  # + I(Elev^2) + Dist2Edge + Elev
  coy_wtr <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + Dist2Edge + Landcover_type + (1|ID)", dat = coyDataz_wtr) # + CanopyCover
  #' #'  Data from only MVCOY68F, NECOY1F, NECOY2M, & NECOY3F in snmr18--- hesitant to extrapolate selection across study areas
  #' coy_smr18 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)", dat = coyData_smr[[1]])  
  #' coy_smr19 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Landcover_type + (1|ID)", dat = coyData_smr[[2]])  #  + Dist2Edge
  #' coy_smr20 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)", dat = coyData_smr[[3]])
  #' #'  Data from only MVCOY68F, NECOY1F, NECOY2M, NECOY3F & NECOY4M in wtr1819--- hesitant to extrapolate selection across study areas
  #' coy_wtr1819 <- glmm_fn(mod = "Used ~ 1 + Elev + Dist2Water + Dist2Edge + Landcover_type + (1|ID)", dat = coyData_wtr[[1]]) # + CanopyCover  + RoadDen + I(Elev^2)  + Slope
  #' coy_wtr1920 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + Landcover_type + (1|ID)", dat = coyData_wtr[[2]])  # + CanopyCover  + Dist2Edge
  #' coy_wtr2021 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)", dat = coyData_wtr[[3]])
  
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
  
  #'  Load spatial libraries
  library(sf)
  library(raster)
  
  #'  Define desired projections
  sa_proj <- projection("EPSG:2855")  # NAD83(HARN) / Washington North
  wgs84 <- projection("+proj=longlat +datum=WGS84 +no_defs")
  
  #'  Read in study area grids (1km^2)
  NE_1km <- raster("./Shapefiles/NE_1km_grid.tif")
  OK_1km <- raster("./Shapefiles/OK_1km_grid.tif")
  # NE_30m <- raster("./Shapefiles/NE_30m_grid.tif")
  # OK_30m <- raster("./Shapefiles/OK_30m_grid.tif")
  
  #'  Load study area shapefiles
  OK.SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "METHOW_SA") %>%
    st_transform(crs = sa_proj)
  OK.SA <- as(OK.SA, "Spatial")
  NE.SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "NE_SA") %>%
    st_transform(crs = sa_proj)
  NE.SA <- as(NE.SA, "Spatial")
  
  #'  Convert rasters to pixels and extract coordinates (centroid of each cell)
  raster_dat <- function(r) {
    dots <- as(r, "SpatialPixelsDataFrame")
    ID <- dots@grid.index
    coords <- coordinates(dots)
    pts <- as.data.frame(cbind(ID, coords))
    return(pts)
  }
  NE_pts <- raster_dat(NE_1km)
  OK_pts <- raster_dat(OK_1km)
  
  #'  Read in covariates extracted across each study area (1km resolution)
  load("./Outputs/Telemetry_covs/NE_covs_1km_2022-01-07.RData") 
  load("./Outputs/Telemetry_covs/OK_covs_1km_2022-01-07.RData")
  # load("./Outputs/Telemetry_covs/NE_covs_30m_2022-01-07.RData") 
  # load("./Outputs/Telemetry_covs/OK_covs_30m_2022-01-13.RData")
  
  #'  Format study area-wide covariate data to include annually relevant data only
  NE.covs <- NE.covs.1km %>%      # NE.covs.30m
    mutate(StudyArea = "NE") %>%
    full_join(NE_pts, by = "ID")
  OK.covs <- OK.covs.1km %>%      # OK.covs.30m
    mutate(StudyArea = "OK") %>%
    full_join(OK_pts, by = "ID")
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
  
  
  #'  Call landcover and scaling functions from above to format covariates
  SA.covs_list <- lapply(SA.covs_list, class_landcov)
  SA.covs_list_reclass <- lapply(SA.covs_list, reclass_landcov)
  SA.covs_list_wolf_reclass <- lapply(SA.covs_list_reclass, reclass_wolf)
  
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
  # mdCov_smr_summary <- lapply(mdData_smr, cov_summary)
  # mdCov_wtr_summary <- lapply(mdData_wtr, cov_summary)
  # elkCov_smr_summary <- lapply(elkData_smr, cov_summary)
  # elkCov_wtr_summary <- lapply(elkData_wtr, cov_summary)
  # wtdCov_smr_summary <- lapply(wtdData_smr, cov_summary)
  # wtdCov_wtr_summary <- lapply(wtdData_wtr, cov_summary)
  # cougCov_smr_summary <- lapply(cougData_smr, cov_summary)
  # cougCov_wtr_summary <- lapply(cougData_wtr, cov_summary)
  # wolfCov_smr_summary <- lapply(wolfData_smr, cov_summary)
  # wolfCov_wtr_summary <- lapply(wolfData_wtr, cov_summary)
  # bobCov_smr_summary <- lapply(bobData_smr, cov_summary)
  # bobCov_wtr_summary <- lapply(bobData_wtr, cov_summary)
  # coyCov_smr_summary <- lapply(coyData_smr, cov_summary)
  # coyCov_wtr_summary <- lapply(coyData_wtr, cov_summary)
  
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

  #'  Function to save parameter estimates & p-values from each RSF
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
  #'  Extract coefficient estimates for each trained model
  md_smr_rsfout <- rsf_out(md_smr, spp = "Mule Deer", season = "Summer")
  md_wtr_rsfout <- rsf_out(md_wtr, spp = "Mule Deer", season = "Winter")
  elk_smr_rsfout <- rsf_out(elk_smr, spp = "Elk", season = "Summer")
  elk_wtr_rsfout <- rsf_out(elk_wtr, spp = "Elk", season = "Winter")
  wtd_smr_rsfout <- rsf_out(wtd_smr, spp = "White-tailed Deer", season = "Summer")
  wtd_wtr_rsfout <- rsf_out(wtd_wtr, spp = "White-tailed Deer", season = "Winter")
  coug_smr_rsfout <- rsf_out(coug_smr, spp = "Cougar", season = "Summer")
  coug_wtr_rsfout <- rsf_out(coug_wtr, spp = "Cougar", season = "Winter")
  wolf_smr_rsfout <- rsf_out(wolf_smr, spp = "Wolf", season = "Summer")
  wolf_wtr_rsfout <- rsf_out(wolf_wtr, spp = "Wolf", season = "Winter")
  bob_smr_rsfout <- rsf_out(bob_smr, spp = "Bobcat", season = "Summer")
  bob_wtr_rsfout <- rsf_out(bob_wtr, spp = "Bobcat", season = "Winter")
  coy_smr_rsfout <- rsf_out(coy_smr, spp = "Coyote", season = "Summer")
  coy_wtr_rsfout <- rsf_out(coy_wtr, spp = "Coyote", season = "Winter")
  
  
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
  
  #'  List and save
  # all_spp_RSF_predicted <- list(md_smr_rsf_sa, md_wtr_rsf_sa, elk_smr_rsf_sa, elk_wtr_rsf_sa, wtd_smr_rsf_sa, wtd_wtr_rsf_sa, coug_smr_rsf_sa, coug_wtr_rsf_sa, wolf_smr_rsf_sa, wolf_wtr_rsf_sa, bob_smr_rsf_sa, bob_wtr_rsf_sa, coy_smr_rsf_sa, coy_wtr_rsf_sa)
  # save(all_spp_RSF_predicted, file = "./Outputs/RSF_output/all_spp_RSF_predicted_", Sys.Date(), ".RData")
  coy_rsf_predicted <- list(coy_smr_rsf_sa, coy_wtr_rsf_sa)
  save(coy_rsf_predicted, file = "./Outputs/RSF_output/coy_rsf_predicted_", Sys.Date(), ".RData")
  
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
  md_smr_rescale_sa <- lapply(md_smr_rsf_sa, RSF_rescale)  
  md_wtr_rescale_sa <- lapply(md_wtr_rsf_sa, RSF_rescale)
  elk_smr_rescale_sa <- lapply(elk_smr_rsf_sa, RSF_rescale)  
  elk_wtr_rescale_sa <- lapply(elk_wtr_rsf_sa, RSF_rescale)
  wtd_smr_rescale_sa <- lapply(wtd_smr_rsf_sa, RSF_rescale)  
  wtd_wtr_rescale_sa <- lapply(wtd_wtr_rsf_sa, RSF_rescale)
  coug_smr_rescale_sa <- lapply(coug_smr_rsf_sa, RSF_rescale)  
  coug_wtr_rescale_sa <- lapply(coug_wtr_rsf_sa, RSF_rescale)
  wolf_smr_rescale_sa <- lapply(wolf_smr_rsf_sa, RSF_rescale)  
  wolf_wtr_rescale_sa <- lapply(wolf_wtr_rsf_sa, RSF_rescale)
  bob_smr_rescale_sa <- lapply(bob_smr_rsf_sa, RSF_rescale)  
  bob_wtr_rescale_sa <- lapply(bob_wtr_rsf_sa, RSF_rescale)
  coy_smr_rescale_sa <- lapply(coy_smr_rsf_sa, RSF_rescale)  # a few major outliers on the high end
  coy_wtr_rescale_sa <- lapply(coy_wtr_rsf_sa, RSF_rescale)
  
  ####  NEED TO DEAL WITH THESE OUTLIERS BEFORE RESCALING!!!
  ####  NEED TO DEAL WITH MASKING WATERBODIES OUT OF GRID BEFORE RESCALING TOO
  
  
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
    # plot(NE.SA, add = T)
    # plot(OK.SA, add = T)
    
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
  coy_smr_RSFstack <- rename_raster(coy_smr_RSFraster)
  coy_wtr_RSFstack <- rename_raster(coy_wtr_RSFraster)

  writeRaster(coy_smr_RSFstack, filename = "./Shapefiles/Predicted_RSFs/coy_smr_RSFstack.tif", bylayer = FALSE, format = 'GTiff', overwrite = TRUE)
  writeRaster(coy_wtr_RSFstack, filename = "./Shapefiles/Predicted_RSFs/coy_wtr_RSFstack.tif", bylayer = FALSE, format = 'GTiff', overwrite = TRUE)
  
  tst <- stack("./Shapefiles/Predicted_RSFs/coy_wtr_RSFstack.tif")
  
  
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

  