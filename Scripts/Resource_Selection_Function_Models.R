  #'  ============================================
  #'  3rd Order Resource Selection Functions (RSFs)
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  December 2021
  #'  ============================================
  #'  Script to run 3rd order resource selection functions (RSFs) and evaluate
  #'  predictive capacity with K-fold cross validation. Final RSFs will be used
  #'  as covariates to represent probability of predator/prey presence across 
  #'  study areas.
  #'  
  #'  Using backwards step selection for initial covariate selection, then K-fold
  #'  cross validation to test model performance for predicting.
  #'  
  #'  Notes on K-fold CV metrics:
  #'  From cvms vignette (pg. 6)
  #'  Confusion matrix: The Pos_ columns tells you whether a row is a True
  #'  Positive (TP), True Negative (TN), False Positive (FP), or False Negative 
  #'  (FN), depending on which level is the "positive" class. I.e. the level you 
  #'  wish to predict.
  #'  Pg. 24-26 of cvms vignette describe how metrics are calculated
  #'  Sensitivity = true positive rate
  #'  Specificity = true negative rate
  #'  
  #'  Kappa statistic is a measure of how closely the observations predicted by
  #'  the model matched the testing data, controlling for the accuracy of the 
  #'  model as measured by the expected accuracy. The kappa statistic is a 
  #'  measure of model performance and is directly comparable to other kappa 
  #'  statistics for competing models. Landis and Koch considers 0-0.20 as slight, 
  #'  0.21-0.40 as fair, 0.41-0.60 as moderate, 0.61-0.80 as substantial, and 
  #'  0.81-1 as almost perfect. Fleiss considers kappas > 0.75 as excellent, 
  #'  0.40-0.75 as fair to good, and < 0.40 as poor. Note that both scales are 
  #'  somewhat arbitrary and what's a "good" kappa value is going to depend on 
  #'  the context and objective of the model. Always consider the kappa statistic
  #'  within context of the confusion matrix. 
  #'  (https://stats.stackexchange.com/questions/82162/cohens-kappa-in-plain-english)
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
    locs$Elev2 <- scale((locs$Elev)^2)
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
    
    locs <- as.data.frame(locs) %>%
      relocate(Elev2, .after = Elev)
    
    return(locs)
  }
  #'  List datasets by season & standardize covariates
  mdData_smr <- list(md_dat_all[md_dat_all$Season == "Summer18",], md_dat_all[md_dat_all$Season == "Summer19",], md_dat_all[md_dat_all$Season == "Summer20",])
  mdData_wtr <- list(md_dat_all[md_dat_all$Season == "Winter1819",], md_dat_all[md_dat_all$Season == "Winter1920",], md_dat_all[md_dat_all$Season == "Winter2021",])
  mdData_smr <- lapply(mdData_smr, standardize_covs)
  mdData_wtr <- lapply(mdData_wtr, standardize_covs)
  #'  Note the reclassified landcover_type data for elkData_winter
  elkData_smr <- list(elk_dat_all[elk_dat_all$Season == "Summer18",], elk_dat_all[elk_dat_all$Season == "Summer19",], elk_dat_all[elk_dat_all$Season == "Summer20",])
  elkData_wtr <- list(elk_dat_all_reclass[elk_dat_all_reclass$Season == "Winter1819",], elk_dat_all_reclass[elk_dat_all_reclass$Season == "Winter1920",], elk_dat_all_reclass[elk_dat_all_reclass$Season == "Winter2021",])
  elkData_smr <- lapply(elkData_smr, standardize_covs)
  elkData_wtr <- lapply(elkData_wtr, standardize_covs)
  #'  Note the reclassified landcover_type data for wtdData_winter 
  wtdData_smr <- list(wtd_dat_all[wtd_dat_all$Season == "Summer18",], wtd_dat_all[wtd_dat_all$Season == "Summer19",], wtd_dat_all[wtd_dat_all$Season == "Summer20",])
  wtdData_wtr <- list(wtd_dat_all_reclass[wtd_dat_all_reclass$Season == "Winter1819",], wtd_dat_all_reclass[wtd_dat_all_reclass$Season == "Winter1920",], wtd_dat_all_reclass[wtd_dat_all_reclass$Season == "Winter2021",])
  wtdData_smr <- lapply(wtdData_smr, standardize_covs)
  wtdData_wtr <- lapply(wtdData_wtr, standardize_covs)
  #'  Note the reclassified landcover_type data for cougData_winter
  cougData_smr <- list(coug_dat_all[coug_dat_all$Season == "Summer18",], coug_dat_all[coug_dat_all$Season == "Summer19",], coug_dat_all[coug_dat_all$Season == "Summer20",])
  cougData_wtr <- list(coug_dat_all_reclass[coug_dat_all_reclass$Season == "Winter1819",], coug_dat_all_reclass[coug_dat_all_reclass$Season == "Winter1920",], coug_dat_all_reclass[coug_dat_all_reclass$Season == "Winter2021",])
  cougData_smr <- lapply(cougData_smr, standardize_covs)
  cougData_wtr <- lapply(cougData_wtr, standardize_covs)
  #'  Note the double reclassified landcover_type data for wolfData
  wolfData_smr <- list(wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Summer18",], wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Summer19",], wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Summer20",])
  wolfData_wtr <- list(wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Winter1819",], wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Winter1920",], wolf_dat_all_reclass2[wolf_dat_all_reclass2$Season == "Winter2021",])
  wolfData_smr <- lapply(wolfData_smr, standardize_covs)
  wolfData_wtr <- lapply(wolfData_wtr, standardize_covs)
  #'  Note the reclassified landcover_type data for bobData
  bobData_smr <- list(bob_dat_all_reclass[bob_dat_all_reclass$Season == "Summer18",], bob_dat_all_reclass[bob_dat_all_reclass$Season == "Summer19",], bob_dat_all_reclass[bob_dat_all_reclass$Season == "Summer20",])
  bobData_wtr <- list(bob_dat_all_reclass[bob_dat_all_reclass$Season == "Winter1819",], bob_dat_all_reclass[bob_dat_all_reclass$Season == "Winter1920",], bob_dat_all_reclass[bob_dat_all_reclass$Season == "Winter2021",])
  bobData_smr <- lapply(bobData_smr, standardize_covs)
  bobData_wtr <- lapply(bobData_wtr, standardize_covs)
  #'  Note the reclassified landcover_type data for coyData
  coyData_smr <- list(coy_dat_all_reclass[coy_dat_all_reclass$Season == "Summer18",], coy_dat_all_reclass[coy_dat_all_reclass$Season == "Summer19",], coy_dat_all_reclass[coy_dat_all_reclass$Season == "Summer20",])
  coyData_wtr <- list(coy_dat_all_reclass[coy_dat_all_reclass$Season == "Winter1819",], coy_dat_all_reclass[coy_dat_all_reclass$Season == "Winter1920",], coy_dat_all_reclass[coy_dat_all_reclass$Season == "Winter2021",])
  coyData_smr <- lapply(coyData_smr, standardize_covs)
  coyData_wtr <- lapply(coyData_wtr, standardize_covs)
  
  
  
  
  #' #'  Run season & species-specific data through prep function
  #' mdData_smr <- spp_dataPrep(md_dat_all[md_dat_all$Season == "Summer18" | md_dat_all$Season == "Summer19" | md_dat_all$Season == "Summer20",])
  #' mdData_wtr <- spp_dataPrep(md_dat_all[md_dat_all$Season == "Winter1819" | md_dat_all$Season == "Winter1920" | md_dat_all$Season == "Winter2021",])
  #' elkData_smr <- spp_dataPrep(elk_dat_all[elk_dat_all$Season == "Summer18" | elk_dat_all$Season == "Summer19" | elk_dat_all$Season == "Summer20",])
  #' elkData_wtr <- spp_dataPrep(elk_dat_all[elk_dat_all$Season == "Winter1819" | elk_dat_all$Season == "Winter1920" | elk_dat_all$Season == "Winter2021",])
  #' wtdData_smr <- spp_dataPrep(wtd_dat_all[wtd_dat_all$Season == "Summer18" | wtd_dat_all$Season == "Summer19" | wtd_dat_all$Season == "Summer20",])
  #' wtdData_wtr <- spp_dataPrep(wtd_dat_all[wtd_dat_all$Season == "Winter1819" | wtd_dat_all$Season == "Winter1920" | wtd_dat_all$Season == "Winter2021",])
  #' cougData_smr <- spp_dataPrep(coug_dat_all[coug_dat_all$Season == "Summer18" | coug_dat_all$Season == "Summer19" | coug_dat_all$Season == "Summer20",])
  #' cougData_wtr <- spp_dataPrep(coug_dat_all[coug_dat_all$Season == "Winter1819" | coug_dat_all$Season == "Winter1920" | coug_dat_all$Season == "Winter2021",])
  #' wolfData_smr <- spp_dataPrep(wolf_dat_all[wolf_dat_all$Season == "Summer18" | wolf_dat_all$Season == "Summer19"| wolf_dat_all$Season == "Summer20",])
  #' wolfData_wtr <- spp_dataPrep(wolf_dat_all[wolf_dat_all$Season == "Winter1819" | wolf_dat_all$Season == "Winter1920"| wolf_dat_all$Season == "Winter2021",])
  #' bobData_smr <- spp_dataPrep(bob_dat_all[bob_dat_all$Season == "Summer18" | bob_dat_all$Season == "Summer19"| bob_dat_all$Season == "Summer20",])
  #' bobData_wtr <- spp_dataPrep(bob_dat_all[bob_dat_all$Season == "Winter1819" | bob_dat_all$Season == "Winter1920"| bob_dat_all$Season == "Winter2021",])
  #' coyData_smr <- spp_dataPrep(coy_dat_all[coy_dat_all$Season == "Summer18" | coy_dat_all$Season == "Summer19" | coy_dat_all$Season == "Summer20",])
  #' coyData_wtr <- spp_dataPrep(coy_dat_all[coy_dat_all$Season == "Winter1819" | coy_dat_all$Season == "Winter1920" | coy_dat_all$Season == "Winter2021",])
  #' 
  #' #'  Landcover_type categories causing convergence issues for some species due to
  #' #'  too few observations in some categories (e.g., "Other", "Wetland") so
  #' #'  reclassifying into fewer categories
  #' reclass_landcov <- function(locs) {
  #'   locs <- locs %>%
  #'     mutate(
  #'       Landcover_type = as.character(as.factor(Landcover_type)),
  #'       Landcover_type = ifelse(Landcover_type == "Developed", "Other", Landcover_type)
  #'     )
  #'   locs$Landcover_type <- droplevels(as.factor(locs$Landcover_type))
  #'   locs$Landcover_type <- relevel(locs$Landcover_type, ref = "Forest")
  #'   
  #'   return(locs)
  #' }
  #' #'  Reclassify landcover categories
  #' mdData_smr_reclass <- reclass_landcov(mdData_smr)
  #' mdData_wtr_reclass <- reclass_landcov(mdData_wtr)
  #' elkData_smr_reclass <- reclass_landcov(elkData_smr)
  #' elkData_wtr_reclass <- reclass_landcov(elkData_wtr)
  #' wtdData_smr_reclass <- reclass_landcov(wtdData_smr)
  #' wtdData_wtr_reclass <- reclass_landcov(wtdData_wtr)
  #' cougData_smr_reclass <- reclass_landcov(cougData_smr)
  #' cougData_wtr_reclass <- reclass_landcov(cougData_wtr)
  #' wolfData_smr_reclass <- reclass_landcov(wolfData_smr)
  #' wolfData_wtr_reclass <- reclass_landcov(wolfData_wtr)
  #' bobData_smr_reclass <- reclass_landcov(bobData_smr)
  #' bobData_wtr_reclass <- reclass_landcov(bobData_wtr)
  #' coyData_smr_reclass <- reclass_landcov(coyData_smr)
  #' coyData_wtr_reclass <- reclass_landcov(coyData_wtr)
  #' 
  #' #'  More reclassification required for all wolf models-
  #' #'  "Other", "Developed", & "Wetland" landcover types causing issues with model
  #' #'  convergence so lumping all together as one class
  #' reclass_landcov <- function(locs) {
  #'   locs <- locs %>%
  #'     mutate(
  #'       Landcover_type = as.character(as.factor(Landcover_type)),
  #'       Landcover_type = ifelse(Landcover_type == "Wetland", "Other", Landcover_type)
  #'     )
  #'   locs$Landcover_type <- droplevels(as.factor(locs$Landcover_type))
  #'   locs$Landcover_type <- relevel(locs$Landcover_type, ref = "Forest")
  #'   
  #'   return(locs)
  #' }
  #' wolfData_smr_reclass2 <- reclass_landcov(wolfData_smr_reclass)
  #' wolfData_wtr_reclass2 <- reclass_landcov(wolfData_wtr_reclass)
  
  
  #'  Correlation Matrix
  #'  ==================
  #'  Function to create correlation matrix for all continuous covariates at once
  cov_correlation <- function(dat) {
    used <- dat[dat$Used == 1,]
    covs <- used[,c("Elev", "Elev2", "Slope", "TPI", "RoadDen",
                    "Dist2Water", "HumanMod", "CanopyCover", "Dist2Edge",
                    "PercForMix", "PercXGrass", "PercXShrub")]
    cor_matrix <- cor(covs, use = "complete.obs")
    return(cor_matrix)
  }
  #'  Generate correlation matrix for each species and season
  (md_smr_corr <- lapply(mdData_smr, cov_correlation))
  (md_wtr_corr <- lapply(mdData_wtr, cov_correlation))
  (elk_smr_corr <- lapply(elkData_smr, cov_correlation))
  (elk_wtr_corr <- lapply(elkData_wtr, cov_correlation))
  (wtd_smr_corr <- lapply(wtdData_smr, cov_correlation))
  (wtd_wtr_corr <- lapply(wtdData_wtr, cov_correlation))
  (coug_smr_corr <- lapply(cougData_smr, cov_correlation))
  (coug_wtr_corr <- lapply(cougData_wtr, cov_correlation))
  (wolf_smr_corr <- lapply(wolfData_smr, cov_correlation))
  (wolf_wtr_corr <- lapply(wolfData_wtr, cov_correlation))
  (bob_smr_corr <- lapply(bobData_smr, cov_correlation))
  (bob_wtr_corr <- lapply(bobData_wtr, cov_correlation))
  (coy_smr_corr <- lapply(coyData_smr, cov_correlation))
  (coy_wtr_corr <- lapply(coyData_wtr, cov_correlation))
  # (md_smr_corr <- cov_correlation(mdData_smr))
  # (md_wtr_corr <- cov_correlation(mdData_wtr))
  # (elk_smr_corr <- cov_correlation(elkData_smr))
  # (elk_wtr_corr <- cov_correlation(elkData_wtr))
  # (wtd_smr_corr <- cov_correlation(wtdData_smr))
  # (wtd_wtr_corr <- cov_correlation(wtdData_wtr))
  # (coug_smr_corr <- cov_correlation(cougData_smr))
  # (coug_wtr_corr <- cov_correlation(cougData_wtr))
  # (wolf_smr_corr <- cov_correlation(wolfData_smr))
  # (wolf_wtr_corr <- cov_correlation(wolfData_wtr))
  # (bob_smr_corr <- cov_correlation(bobData_smr))
  # (bob_wtr_corr <- cov_correlation(bobData_wtr))
  # (coy_smr_corr <- cov_correlation(coyData_smr))
  # (coy_wtr_corr <- cov_correlation(coyData_wtr))

  #'  Elevation & TPI are highly correlated (almost 100%) for all datasets so nixing TPI entirely
  #'  Elevation & Human Modified correlated in MD smr/wtr, ELK smr, WTD smr, COUG smr, & COY smr/wtr- nixing Human Mod for those models
  #'  Elevation & Human Modified correlated with WOLF smr20
  #'  Slope & Elev highly correlated in WTD wtr1819; CanopyCover & Elev correlated (0.6) in WTD smr18
  #'  Using Landcover_type instead of % Forest, % Grass, & % Shrub because...
  #'  % Shrub correlated with Elevation, RoadDen, and Human Modified in MD smr
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
  md_smr18 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + RoadDen + CanopyCover + Dist2Edge + Landcover_type + (1|ID)", dat = mdData_smr[[1]]) # + Slope  + Dist2Water
  md_smr19 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Water + CanopyCover + Landcover_type + (1|ID)", dat = mdData_smr[[2]]) # + Dist2Edge
  md_smr20 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + RoadDen + Dist2Water + CanopyCover + Landcover_type + (1|ID)", dat = mdData_smr[[3]]) # + Dist2Edge + Slope
  
  md_wtr1819 <- glmm_fn(mod = "Used ~ 1 + Elev + Slope + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)", dat = mdData_wtr[[1]]) # + I(Elev^2) + RoadDen
  md_wtr1920 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)", dat = mdData_wtr[[2]])
  md_wtr2021 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)", dat = mdData_wtr[[3]])
  
  ####  Elk RSFs  ####
  #'  Dropping HumanMod in elk summer models due to high correlation with other covariates
  elk_smr18 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = elkData_smr[[1]])
  elk_smr19 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = elkData_smr[[2]]) 
  elk_smr20 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = elkData_smr[[3]])
  
  #'  Note: using reclassified version of landcover for winter elk models
  elk_wtr1819 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = elkData_wtr[[1]])
  elk_wtr1920 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = elkData_wtr[[2]]) # + HumanMod
  elk_wtr2021 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = elkData_wtr[[3]])
  
  ####  White-tailed Deer RSFs  ####
  #'  Dropping HumanMod in wtd summer models due to high correlation with other covariates
  #'  Removed CanopyCover from summer18 model due to collinearity with Elevation
  wtd_smr18 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Edge + Landcover_type + (1|ID)",  dat = wtdData_smr[[1]]) # + Dist2Water
  wtd_smr19 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wtdData_smr[[2]])
  wtd_smr20 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wtdData_smr[[3]])
  
  #'  Note: using reclassified version of landcover for winter WTD models
  #'  Take a close look at Slope- correlated with Elev but Elev^2 might address this issue
  wtd_wtr1819 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + Dist2Edge + Landcover_type + (1|ID)",  dat = wtdData_wtr[[1]]) # + Dist2Water  + CanopyCover + RoadDen + HumanMod
  wtd_wtr1920 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wtdData_wtr[[2]])
  wtd_wtr2021 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wtdData_wtr[[3]])
  
  ####  Cougar RSFs  ####
  #'  Dropping HumanMod in cougar summer models due to high correlation with other covariates
  coug_smr18 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = cougData_smr[[1]])  
  coug_smr19 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = cougData_smr[[2]])
  coug_smr20 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = cougData_smr[[3]]) 
  
  #'  Note: using reclassified version of landcover for winter cougar models
  coug_wtr1819 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = cougData_wtr[[1]]) # + RoadDen
  coug_wtr1920 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = cougData_wtr[[2]]) # + RoadDen
  coug_wtr2021 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = cougData_wtr[[3]]) # + RoadDen  + Dist2Water
  
  ####  Wolf RSFs  ####
  #'  NOTE: using 2nd reclassified version of landcover categories for wolf models
  #'  Dropped HumanMod from summer20 model due to collinearity with Elevation
  wolf_smr18 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + Dist2Water + Dist2Edge + Landcover_type + (1|ID)",  dat = wolfData_smr[[1]]) #+ HumanMod  + RoadDen  + CanopyCover
  wolf_smr19 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wolfData_smr[[2]])
  wolf_smr20 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wolfData_smr[[3]])
  
  wolf_wtr1819 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wolfData_wtr[[1]])  # + Dist2Water 
  wolf_wtr1920 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + HumanMod + Dist2Edge + Landcover_type + (1|ID)",  dat = wolfData_wtr[[2]])  # + Dist2Water  + CanopyCover
  wolf_wtr2021 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wolfData_wtr[[3]])
  
  ####  Bobcat RSFs  ####
  #'  Note: using reclassified version of landcover for summer bobcat models
  #'  Only data for MVBOB90M in smr18--- not enough data to make inference about bobcat resource selection across 2 study areas
  # bob_smr18 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = bobData_smr[[1]])
  bob_smr19 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + HumanMod + Dist2Edge + Landcover_type + (1|ID)",  dat = bobData_smr[[2]]) # + CanopyCover + Dist2Water
  bob_smr20 <- glmm_fn(mod = "Used ~ 1 + RoadDen + Dist2Water + HumanMod + Landcover_type + (1|ID)",  dat = bobData_smr[[3]])  # + Slope + CanopyCover + I(Elev^2) + Dist2Edge + Elev
  
  #'  Only data for MVBOB88M & MVBOB90M wtr1819--- not enough data to make inference about bobcat resource selection across 2 study areas
  # bob_wtr1819 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = bobData_wtr[[1]])
  bob_wtr1920 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = bobData_wtr[[2]]) #  + Dist2Water
  bob_wtr2021 <- glmm_fn(mod = "Used ~ 1 + Elev + Slope + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = bobData_wtr[[3]]) # + Dist2Water + RoadDen + I(Elev^2)
  
  ####  Coyote RSFs  ####
  #'  Note: using reclassified version of landcover for all coyote models
  #'  Dropping HumanMod from all coyote models due to high correlation with other covariates
  #'  Data from only MVCOY68F, NECOY1F, NECOY2M, & NECOY3F in snmr18--- hesitant to extrapolate selection across study areas
  coy_smr18 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)", dat = coyData_smr[[1]])  
  coy_smr19 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Water + CanopyCover + Landcover_type + (1|ID)", dat = coyData_smr[[2]])  #  + Dist2Edge
  coy_smr20 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)", dat = coyData_smr[[3]])
  
  #'  Data from only MVCOY68F, NECOY1F, NECOY2M, NECOY3F & NECOY4M in wtr1819--- hesitant to extrapolate selection across study areas
  coy_wtr1819 <- glmm_fn(mod = "Used ~ 1 + Elev + Dist2Water + Dist2Edge + Landcover_type + (1|ID)", dat = coyData_wtr[[1]]) # + CanopyCover  + RoadDen + Elev2  + Slope
  coy_wtr1920 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Water + Landcover_type + (1|ID)", dat = coyData_wtr[[2]])  # + CanopyCover  + Dist2Edge
  coy_wtr2021 <- glmm_fn(mod = "Used ~ 1 + Elev + Elev2 + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)", dat = coyData_wtr[[3]])
  
  #'  Group species-specific models
  RSF_MD_list <- list(md_smr18, md_smr19, md_smr20, md_wtr1819, md_wtr1920, md_wtr2021)
  RSF_ELK_list <- list(elk_smr18, elk_smr19, elk_smr20, elk_wtr1819, elk_wtr1920, elk_wtr2021)
  RSF_WTD_list <- list(wtd_smr18, wtd_smr19, wtd_smr20, wtd_wtr1819, wtd_wtr1920, wtd_wtr2021)
  RSF_COUG_list <- list(coug_smr18, coug_smr19, coug_smr20, coug_wtr1819, coug_wtr1920, coug_wtr2021)
  RSF_WOLF_list <- list(wolf_smr18, wolf_smr19, wolf_smr20, wolf_wtr1819, wolf_wtr1920, wolf_wtr2021)
  RSF_BOB_list <- list(bob_smr18, bob_smr19, bob_smr20, bob_wtr1819, bob_wtr1920, bob_wtr2021)
  RSF_COY_list <- list(coy_smr18, coy_smr19, coy_smr20, coy_wtr1819, coy_wtr1920, coy_wtr2021)
  
  #'  Save
  save(RSF_MD_list, file = paste0("./Outputs/RSF_output/RSF_MD_list_", Sys.Date(), ".RData"))
  save(RSF_ELK_list, file = paste0("./Outputs/RSF_output/RSF_ELK_list_", Sys.Date(), ".RData"))
  save(RSF_WTD_list, file = paste0("./Outputs/RSF_output/RSF_WTD_list_", Sys.Date(), ".RData"))
  save(RSF_COUG_list, file = paste0("./Outputs/RSF_output/RSF_COUG_list_", Sys.Date(), ".RData"))
  save(RSF_WOLF_list, file = paste0("./Outputs/RSF_output/RSF_WOLF_list_", Sys.Date(), ".RData"))
  save(RSF_BOB_list, file = paste0("./Outputs/RSF_output/RSF_BOB_list_", Sys.Date(), ".RData"))
  save(RSF_COY_list, file = paste0("./Outputs/RSF_output/RSF_COY_list_", Sys.Date(), ".RData"))
  

  #' #'  K-fold Cross Validation
  #' #'  =======================
  #' #'  Assess predictive capacity for each RSF by splitting data into K folds and
  #' #'  leaving one fold of data out to test how predictive the model is. Running
  #' #'  in parallel to speed up the K-fold analysis. Using cvms package to handle 
  #' #'  random effects in RSFs. Code was adapted from the following websites:
  #' #'  https://github.com/LudvigOlsen/cvms#examples
  #' #'  https://cran.r-project.org/web/packages/cvms/cvms.pdf
  #' 
  #' #'  Monitor time
  #' start.time <- Sys.time()
  #' #'  Set up to run in parallel
  #' #'  Identify how many cores I want to use
  #' detectCores(logical = FALSE)
  #' #' Run in parallel on local computer with specified number of cores
  #' registerDoParallel(16) # reduce when not using lab computer
  #' 
  #' #'  Function to run k-fold cross validation on each RSF
  #' #'  Requires the data set, number of folds (K), and regression model be defined
  #' k_fold_rsf <- function(dat, K, mod) {
  #'   
  #'   #' #'  Set seed so its reproducible
  #'   #' set.seed(2021)
  #'   
  #'   #'  Partition data into folds
  #'   #'  Use groupdata2 package with cat_col = "Used" to balance folds proportional
  #'   #'  to 0's and 1's in dataset
  #'   fold_df <- fold(dat, k = K, cat_col = "Used")
  #'   
  #'   #'  Run K-fold cross-validation
  #'   #'  Positive argument indicates the level from `targets` (Used) to predict 
  #'   #'  (positive = 1 means predict the 0's, postive = 2 means predict the 1's)
  #'   #'  Preprocessing argument centers & scales UNstandardized continuous covariates 
  #'   #'  for each fold, but this always produces errors when I try to use it
  #'   CV1 <- cross_validate(fold_df, formulas = mod, family = "binomial", positive = 2, parallel = TRUE) #, parallel = TRUE, preprocessing = "standardize", REML = FALSE,  #control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  #'   
  #'   #'  View coefficients per fold
  #'   print(CV1$`Coefficients`[[1]] %>% kable())
  #'   # Metrics
  #'   print(CV1 %>% select(1:15) %>% kable(digits = 5))
  #'   #' Confusion matrix- important when considering Kappa statistic and calculating accuracy metrics
  #'   print(CV1$`Confusion Matrix`[[1]] %>% kable())
  #'   
  #'   return(CV1)
  #' }
  #' #'  Define number of folds
  #' K <- 5
  #' 
  #' #'  Run data & model for each species, season, & year in parallel
  #' ####  Mule Deer K-fold  ####
  #' md_smr18_cv <- k_fold_rsf(dat =  mdData_smr[mdData_smr$Year == "Year1",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + RoadDen + CanopyCover + Landcover_type + (1|ID)")  # + Slope  + Dist2Water  + Dist2Edge
  #' md_smr19_cv <- k_fold_rsf(dat =  mdData_smr[mdData_smr$Year == "Year2",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Landcover_type + (1|ID)") # + Dist2Edge
  #' md_smr20_cv <- k_fold_rsf(dat =  mdData_smr[mdData_smr$Year == "Year3",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + RoadDen + Dist2Water + CanopyCover + Landcover_type + (1|ID)") # + Dist2Edge  + Slope
  #' md_wtr1819_cv <- k_fold_rsf(dat =  mdData_wtr[mdData_wtr$Year == "Year1",], K = K, mod = "Used ~ 1 + Elev + Slope + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)") # + RoadDen  + I(Elev^2)
  #' md_wtr1920_cv <- k_fold_rsf(dat =  mdData_wtr[mdData_wtr$Year == "Year2",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)")
  #' md_wtr2021_cv <- k_fold_rsf(dat =  mdData_wtr[mdData_wtr$Year == "Year3",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)")
  #' 
  #' ####  Elk K-fold  ####
  #' elk_smr18_cv <- k_fold_rsf(dat = elkData_smr[elkData_smr$Year == "Year1",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)")
  #' elk_smr19_cv <- k_fold_rsf(dat = elkData_smr[elkData_smr$Year == "Year2",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)") # + Slope
  #' elk_smr20_cv <- k_fold_rsf(dat = elkData_smr[elkData_smr$Year == "Year3",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)")
  #' elk_wtr1819_cv <- k_fold_rsf(dat = elkData_wtr_reclass[elkData_wtr_reclass$Year == "Year1",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)")
  #' elk_wtr1920_cv <- k_fold_rsf(dat = elkData_wtr_reclass[elkData_wtr_reclass$Year == "Year2",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)")
  #' elk_wtr2021_cv <- k_fold_rsf(dat = elkData_wtr_reclass[elkData_wtr_reclass$Year == "Year3",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)")
  #' 
  #' ####  White-tailed Deer K-fold  ####
  #' wtd_smr18_cv <- k_fold_rsf(dat = wtdData_smr[wtdData_smr$Year == "Year1",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + HumanMod + Dist2Edge + Landcover_type + (1|ID)") # + CanopyCover + Dist2Water
  #' wtd_smr19_cv <- k_fold_rsf(dat = wtdData_smr[wtdData_smr$Year == "Year2",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)")
  #' wtd_smr20_cv <- k_fold_rsf(dat = wtdData_smr[wtdData_smr$Year == "Year3",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)")
  #' wtd_wtr1819_cv <- k_fold_rsf(dat = wtdData_wtr_reclass[wtdData_wtr_reclass$Year == "Year1",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + Dist2Edge + Landcover_type + (1|ID)") # + Dist2Water + CanopyCover + RoadDen + HumanMod
  #' wtd_wtr1920_cv <- k_fold_rsf(dat = wtdData_wtr_reclass[wtdData_wtr_reclass$Year == "Year2",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)")
  #' wtd_wtr2021_cv <- k_fold_rsf(dat = wtdData_wtr_reclass[wtdData_wtr_reclass$Year == "Year3",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)")
  #' 
  #' ####  Cougar K-fold  ####
  #' coug_smr18_cv <- k_fold_rsf(dat = cougData_smr[cougData_smr$Year == "Year1",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)")  # + HumanMod
  #' coug_smr19_cv <- k_fold_rsf(dat = cougData_smr[cougData_smr$Year == "Year2",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)")
  #' coug_smr20_cv <- k_fold_rsf(dat = cougData_smr[cougData_smr$Year == "Year3",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)") #  + HumanMod
  #' coug_wtr1819_cv <- k_fold_rsf(dat = cougData_wtr_reclass[cougData_wtr_reclass$Year == "Year1",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)") # + RoadDen
  #' coug_wtr1920_cv <- k_fold_rsf(dat = cougData_wtr_reclass[cougData_wtr_reclass$Year == "Year2",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)") # + RoadDen
  #' coug_wtr2021_cv <- k_fold_rsf(dat = cougData_wtr_reclass[cougData_wtr_reclass$Year == "Year3",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)") # + RoadDen + Dist2Water
  #' 
  #' ####  Wolf K-fold  ####
  #' #'  NOTE: using 2nd reclassified version of landcover categories for wolf models
  #' wolf_smr18_cv <- k_fold_rsf(dat = wolfData_smr_reclass2[wolfData_smr_reclass2$Year == "Year1",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)")
  #' wolf_smr19_cv <- k_fold_rsf(dat = wolfData_smr_reclass2[wolfData_smr_reclass2$Year == "Year2",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)")
  #' wolf_smr20_cv <- k_fold_rsf(dat = wolfData_smr_reclass2[wolfData_smr_reclass2$Year == "Year3",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)")
  #' wolf_wtr1819_cv <- k_fold_rsf(dat = wolfData_wtr_reclass2[wolfData_wtr_reclass2$Year == "Year1",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)")  # + Dist2Water
  #' wolf_wtr1920_cv <- k_fold_rsf(dat = wolfData_wtr_reclass2[wolfData_wtr_reclass2$Year == "Year2",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)")  # + Dist2Water
  #' wolf_wtr2021_cv <- k_fold_rsf(dat = wolfData_wtr_reclass2[wolfData_wtr_reclass2$Year == "Year3",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)")
  #' 
  #' ####  Bobcat K-fold  ####
  #' #'  Note: using reclassified version of landcover for summer bobcat models
  #' #'  Only data for MVBOB90M in smr18--- not enough data to make inference about bobcat resource selection across 2 study areas & no need for random effect
  #' # bob_smr18_cv <- k_fold_rsf(dat = bobData_smr[bobData_smr$Year == "Year1",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)")
  #' bob_smr19_cv <- k_fold_rsf(dat = bobData_smr_reclass[bobData_smr_reclass$Year == "Year2",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)")
  #' bob_smr20_cv <- k_fold_rsf(dat = bobData_smr_reclass[bobData_smr_reclass$Year == "Year3",], K = K, mod = "Used ~ 1 + Elev + Dist2Water + HumanMod + Dist2Edge + Landcover_type + (1|ID)")  # + CanopyCover + I(Elev^2) + Slope  + RoadDen
  #' #'  Only data for MVBOB88M & MVBOB90M wtr1819--- not enough data to make inference about bobcat resource selection across 2 study areas
  #' # bob_wtr1819_cv <- k_fold_rsf(dat = bobData_wtr[bobData_wtr$Year == "Year1",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)")
  #' bob_wtr1920_cv <- k_fold_rsf(dat = bobData_wtr[bobData_wtr$Year == "Year2",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + HumanMod + Dist2Edge + Landcover_type + (1|ID)") # + Dist2Water  + CanopyCover
  #' bob_wtr2021_cv <- k_fold_rsf(dat = bobData_wtr[bobData_wtr$Year == "Year3",], K = K, mod = "Used ~ 1 + Elev + Slope + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)") # + Dist2Water  + RoadDen  + I(Elev^2)
  #' 
  #' ####  Coyote K- fold  ####
  #' #'  Note: using reclassified version of landcover for all coyote models
  #' #'  Data from only MVCOY68F, NECOY1F, NECOY2M, & NECOY3F in snmr18--- hesitant to extrapolate selection across study areas
  #' coy_smr18_cv <- k_fold_rsf(dat = coyData_smr_reclass[coyData_smr_reclass$Year == "Year1",], K = K, mod = "Used ~ 1 + Elev + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)") # + I(Elev^2) 
  #' coy_smr19_cv <- k_fold_rsf(dat = coyData_smr_reclass[coyData_smr_reclass$Year == "Year2",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Landcover_type + (1|ID)") # + Dist2Edge
  #' coy_smr20_cv <- k_fold_rsf(dat = coyData_smr_reclass[coyData_smr_reclass$Year == "Year3",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)")
  #' #'  Data from only MVCOY68F, NECOY1F, NECOY2M, NECOY3F & NECOY4M in wtr1819--- hesitant to extrapolate selection across study areas
  #' coy_wtr1819_cv <- k_fold_rsf(dat = coyData_wtr_reclass[coyData_wtr_reclass$Year == "Year1",], K = K, mod = "Used ~ 1 + Elev + Dist2Water + Dist2Edge + Landcover_type + (1|ID)") # + RoadDen + I(Elev^2) + CanopyCover + Slope
  #' coy_wtr1920_cv <- k_fold_rsf(dat = coyData_wtr_reclass[coyData_wtr_reclass$Year == "Year2",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + Landcover_type + (1|ID)") # + CanopyCover  + Dist2Edge
  #' coy_wtr2021_cv <- k_fold_rsf(dat = coyData_wtr_reclass[coyData_wtr_reclass$Year == "Year3",], K = K, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)")
  #' 
  #' #'  IGNORE WARNING MESSAGE!!!
  #' #'  Get a warning message associated with running in parallel but K-fold  
  #' #'  appears to be working fine so ignoring it for now
  #' #'  Warning messages: <anonymous>: ... may be used in an incorrect context: ‘.fun(piece, ...)’
  #' 
  #' #'  End time keeping
  #' end.time <- Sys.time()
  #' #'  How long did this take?
  #' difftime(end.time, start.time, units = "hours")
  #' 
  #' #'  Compute the average of the k recorded errors. This is called the 
  #' #'  cross-validation error serving as the performance metric for the model.
  #' #'  Report Accuracy and Kappa statistics!
  #'  
  #' 
  #' #'  Group species-specific K-fold results
  #' CV_MD_list <- list(md_smr18_cv, md_smr19_cv, md_smr20_cv, md_wtr1819_cv, md_wtr1920_cv, md_wtr2021_cv)
  #' CV_ELK_list <- list(elk_smr18_cv, elk_smr19_cv, elk_smr20_cv, elk_wtr1819_cv, elk_wtr1920_cv, elk_wtr2021_cv)
  #' CV_WTD_list <- list(wtd_smr18_cv, wtd_smr19_cv, wtd_smr20_cv, wtd_wtr1819_cv, wtd_wtr1920_cv, wtd_wtr2021_cv)
  #' CV_COUG_list <- list(coug_smr18_cv, coug_smr19_cv, coug_smr20_cv, coug_wtr1819_cv, coug_wtr1920_cv, coug_wtr2021_cv)
  #' CV_WOLF_list <- list(wolf_smr18_cv, wolf_smr19_cv, wolf_smr20_cv, wolf_wtr1819_cv, wolf_wtr1920_cv, wolf_wtr2021_cv)
  #' CV_BOB_list <- list(bob_smr18_cv, bob_smr19_cv, bob_smr20_cv, bob_wtr1819_cv, bob_wtr1920_cv, bob_wtr2021_cv)
  #' CV_COY_list <- list(coy_smr18_cv, coy_smr19_cv, coy_smr20_cv, coy_wtr1819_cv, coy_wtr1920_cv, coy_wtr2021_cv)
  #' 
  #' #'  Save
  #' save(CV_MD_list, file = paste0("./Outputs/RSF_output/CV_MD_list_", Sys.Date(), ".RData"))
  #' save(CV_ELK_list, file = paste0("./Outputs/RSF_output/CV_ELK_list_", Sys.Date(), ".RData"))
  #' save(CV_WTD_list, file = paste0("./Outputs/RSF_output/CV_WTD_list_", Sys.Date(), ".RData"))
  #' save(CV_COUG_list, file = paste0("./Outputs/RSF_output/CV_COUG_list_", Sys.Date(), ".RData"))
  #' save(CV_WOLF_list, file = paste0("./Outputs/RSF_output/CV_WOLF_list_", Sys.Date(), ".RData"))
  #' save(CV_BOB_list, file = paste0("./Outputs/RSF_output/CV_BOB_list_", Sys.Date(), ".RData"))
  #' save(CV_COY_list, file = paste0("./Outputs/RSF_output/CV_COY_list_", Sys.Date(), ".RData"))
  #' 
  #' 
  #' 
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

  