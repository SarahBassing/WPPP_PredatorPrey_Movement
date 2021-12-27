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
  library(parallel)
  library(doParallel)
  library(future.apply)
  
  #'  Load used and available locations, and covariate data
  load("./Outputs/RSF_pts/md_dat_all_2021-12-10.RData")
  load("./Outputs/RSF_pts/elk_dat_all_2021-12-10.RData")
  load("./Outputs/RSF_pts/wtd_dat_all_2021-12-10.RData")
  load("./Outputs/RSF_pts/coug_dat_all_2021-12-10.RData")
  load("./Outputs/RSF_pts/wolf_dat_all_2021-12-10.RData")
  load("./Outputs/RSF_pts/bob_dat_all_2021-12-10.RData")
  load("./Outputs/RSF_pts/coy_dat_all_2021-12-10.RData")
  
  
  #'  Center & scale covariates 
  #'  Note: standardizing across all IDs & years, but separately by season & spp
  spp_dataPrep <- function(locs){
    #'  Make categorical variables factors
    locs$ID <- as.factor(locs$ID)
    locs$Used <- as.factor(locs$Used)
    locs$Area <- as.factor(locs$Area)
    locs$Year <- as.factor(locs$Year)
    locs$Season <- as.factor(locs$Season)
    locs$Landcover <- as.factor(locs$Landcover)
    # locs$Landcover_type <- as.character(as.factor(locs$Landcover_type))
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
    locs$Landcover_type <- droplevels(as.factor(locs$Landcover_type))
    locs$Landcover_type <- relevel(locs$Landcover_type, ref = "Forest")
  
    locs <- as.data.frame(locs)
    
    return(locs)
  }
  #'  Run season & species-specific data through prep function
  mdData_smr <- spp_dataPrep(md_dat_all[md_dat_all$Season == "Summer18" | md_dat_all$Season == "Summer19" | md_dat_all$Season == "Summer20",])
  mdData_wtr <- spp_dataPrep(md_dat_all[md_dat_all$Season == "Winter1819" | md_dat_all$Season == "Winter1920" | md_dat_all$Season == "Winter2021",])
  elkData_smr <- spp_dataPrep(elk_dat_all[elk_dat_all$Season == "Summer18" | elk_dat_all$Season == "Summer19" | elk_dat_all$Season == "Summer20",])
  elkData_wtr <- spp_dataPrep(elk_dat_all[elk_dat_all$Season == "Winter1819" | elk_dat_all$Season == "Winter1920" | elk_dat_all$Season == "Winter2021",])
  wtdData_smr <- spp_dataPrep(wtd_dat_all[wtd_dat_all$Season == "Summer18" | wtd_dat_all$Season == "Summer19" | wtd_dat_all$Season == "Summer20",])
  wtdData_wtr <- spp_dataPrep(wtd_dat_all[wtd_dat_all$Season == "Winter1819" | wtd_dat_all$Season == "Winter1920" | wtd_dat_all$Season == "Winter2021",])
  cougData_smr <- spp_dataPrep(coug_dat_all[coug_dat_all$Season == "Summer18" | coug_dat_all$Season == "Summer19" | coug_dat_all$Season == "Summer20",])
  cougData_wtr <- spp_dataPrep(coug_dat_all[coug_dat_all$Season == "Winter1819" | coug_dat_all$Season == "Winter1920" | coug_dat_all$Season == "Winter2021",])
  wolfData_smr <- spp_dataPrep(wolf_dat_all[wolf_dat_all$Season == "Summer18" | wolf_dat_all$Season == "Summer19"| wolf_dat_all$Season == "Summer20",])
  wolfData_wtr <- spp_dataPrep(wolf_dat_all[wolf_dat_all$Season == "Winter1819" | wolf_dat_all$Season == "Winter1920"| wolf_dat_all$Season == "Winter2021",])
  bobData_smr <- spp_dataPrep(bob_dat_all[bob_dat_all$Season == "Summer18" | bob_dat_all$Season == "Summer19"| bob_dat_all$Season == "Summer20",])
  bobData_wtr <- spp_dataPrep(bob_dat_all[bob_dat_all$Season == "Winter1819" | bob_dat_all$Season == "Winter1920"| bob_dat_all$Season == "Winter2021",])
  coyData_smr <- spp_dataPrep(coy_dat_all[coy_dat_all$Season == "Summer18" | coy_dat_all$Season == "Summer19" | coy_dat_all$Season == "Summer20",])
  coyData_wtr <- spp_dataPrep(coy_dat_all[coy_dat_all$Season == "Winter1819" | coy_dat_all$Season == "Winter1920" | coy_dat_all$Season == "Winter2021",])
  
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
  #'  Reclassify landcover categories
  mdData_smr_reclass <- reclass_landcov(mdData_smr)
  mdData_wtr_reclass <- reclass_landcov(mdData_wtr)
  elkData_smr_reclass <- reclass_landcov(elkData_smr)
  elkData_wtr_reclass <- reclass_landcov(elkData_wtr)
  wtdData_smr_reclass <- reclass_landcov(wtdData_smr)
  wtdData_wtr_reclass <- reclass_landcov(wtdData_wtr)
  cougData_smr_reclass <- reclass_landcov(cougData_smr)
  cougData_wtr_reclass <- reclass_landcov(cougData_wtr)
  wolfData_smr_reclass <- reclass_landcov(wolfData_smr)
  wolfData_wtr_reclass <- reclass_landcov(wolfData_wtr)
  bobData_smr_reclass <- reclass_landcov(bobData_smr)
  bobData_wtr_reclass <- reclass_landcov(bobData_wtr)
  coyData_smr_reclass <- reclass_landcov(coyData_smr)
  coyData_wtr_reclass <- reclass_landcov(coyData_wtr)
  
  
  
  #' #'  Function to create correlation matrix for all continuous covariates at once
  #' cov_correlation <- function(dat) {
  #'   used <- dat[dat$Used == 1,]
  #'   covs <- used[,c("Elev", "Slope", "TPI", "RoadDen", 
  #'                   "Dist2Water", "HumanMod", "CanopyCover", "Dist2Edge",
  #'                   "PercForMix", "PercXGrass", "PercXShrub")]
  #'   cor_matrix <- cor(covs, use = "complete.obs")
  #'   return(cor_matrix)
  #' }
  #' #'  Generate correlation matrix for each species and season
  #' (md_smr_corr <- cov_correlation(mdData_smr)) 
  #' (md_wtr_corr <- cov_correlation(mdData_wtr)) 
  #' (elk_smr_corr <- cov_correlation(elkData_smr)) 
  #' (elk_wtr_corr <- cov_correlation(elkData_wtr)) 
  #' (wtd_smr_corr <- cov_correlation(wtdData_smr))
  #' (wtd_wtr_corr <- cov_correlation(wtdData_wtr))
  #' (coug_smr_corr <- cov_correlation(cougData_smr)) 
  #' (coug_wtr_corr <- cov_correlation(cougData_wtr)) 
  #' (wolf_smr_corr <- cov_correlation(wolfData_smr))
  #' (wolf_wtr_corr <- cov_correlation(wolfData_wtr)) 
  #' (bob_smr_corr <- cov_correlation(bobData_smr))
  #' (bob_wtr_corr <- cov_correlation(bobData_wtr)) 
  #' (coy_smr_corr <- cov_correlation(coyData_smr)) 
  #' (coy_wtr_corr <- cov_correlation(coyData_wtr)) 
  #' 
  #' #'  Elevation & TPI are highly correlated (almost 100%) for all datasets so nixing TPI entirely
  #' #'  Elevation & Human Modified correlated in MD smr/wtr, ELK smr, & COY smr- nixing Human Mod for those models
  #' #'  Using Landcover_type instead of % Forest, % Grass, & % Shrub because...
  #' #'  % Shrub correlated with Elevation, RoadDen, and Human Modified in MD smr
  #' #'  Dist2Edg correlated with % Forest & % Grass in ELK wtr
  #' #'  % Forest & Grass correlated in COUG and WOLF wtr
  #' #'  % Grass & Shrub correlated for BOb smr
  #' #'  % Forest & Canopy Cover, % Forest & Shrub correlated for BOB wtr
  

  ####  K-fold CV  ####
  #https://github.com/LudvigOlsen/cvms#examples
  #https://cran.r-project.org/web/packages/cvms/cvms.pdf
  
  #'  Monitor time
  start.time <- Sys.time()
  #'  Set up to run in parallel
  #'  Identify how many cores I want to use
  detectCores(logical = FALSE)
  cl <- parallel::makeCluster(4)  # change to 14 when working correctly on lab computer
  #'  Run in parallel on local computer with specified number of cores
  plan(cluster, workers = cl)
  # registerDoParallel(4)
  
  #'  Function to run k-fold cross validation on each RSF
  #'  Requires the data set, number of folds (K), and regression model be defined
  k_fold_rsf <- function(dat, K, mod) {
    
    #'  Set seed so its reproducible
    set.seed(2021)
    
    #'  Partition data into folds
    #'  Use groupdata2 package with cat_col = "Used" to balance folds proportional
    #'  to 0's and 1's in dataset
    fold_df <- fold(dat, k = K) #k = 2
    
    #'  Define model
    # mod <- mod #"Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)" 
    
    #'  Run K-fold cross-validation
    #'  Positive argument indicates the level from `targets` (Used) to predict 
    #'  (positive = 1 means predict the 0's, postive = 2 means predict the 1's)
    #'  Preprocessing argument centers & scales UNstandardized continuous covariates 
    #'  for each fold, but this always produces errors when I try to use it
    CV1 <- cross_validate(fold_df, formulas = mod, family = "binomial", positive = 2) #preprocessing = "standardize", REML = FALSE,  #control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    
    #'  View coefficients per fold
    print(CV1$`Coefficients`[[1]] %>% kable())
    # Metrics
    print(CV1 %>% select(1:15) %>% kable(digits = 5))
    #' Confusion matrix- important when considering Kappa statistic and calculating accuracy metrics
    print(CV1$`Confusion Matrix`[[1]] %>% kable())
    
    return(CV1)
  }
  #'  Assess predictive capacity for each RSF
  #'  Run data & model for each species, season, & year in parallel
  # bob_wtr1920_cv <- k_fold_rsf(dat = bobData_wtr[bobData_wtr$Year == "Year2",], K = 2, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + HumanMod + Dist2Edge + Landcover_type + (1|ID)")
  # coy_smr20_cv <- k_fold_rsf(dat = coyData_smr_reclass[coyData_smr_reclass$Year == "Year3",], K = 2, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)")
  coy_smr20_cv <- future_lapply(dat = coyData_smr_reclass[coyData_smr_reclass$Year == "Year3",], FUN = k_fold_rsf, K = 2, mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)")
  
  
  
  #'  End time keeping
  end.time <- Sys.time()
  #'  Stop running in parallel
  parallel::stopCluster(cl)
  #'  How long did this take?
  difftime(end.time, start.time, units = "hours")
  
  
  
  #' mddata = mdData_smr[mdData_smr$Year == "Year3",]
  #' coydata = coyData_smr_reclass[coyData_smr_reclass$Year == "Year3",]
  #' #'  Partition training data
  #' #'  Use groupdata2 package with cat_col = "Used" to balance folds proportional
  #' #'  to 0's and 1's in dataset
  #' set.seed(2021)
  #' fold_df <- fold(mddata, k = 2)
  #' 
  #' #'  Define model
  #' formulas <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)" 
  #' 
  #' #'  Run K-fold cross-validation
  #' #'  Positive argument indicates the level from `targets` (Used) to predict (positive = 1 means predict the 0's, postive = 2 means predict the 1's)
  #' #'  Preprocessing argument centers & scales UNstandardized continuous covariates for each fold, but this always produces errors when I try to use it
  #' CV1 <- cross_validate(fold_df, formulas = formulas, family = "binomial", positive = 2) #preprocessing = "standardize", REML = FALSE, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  #' 
  #' #'  View coefficients per fold
  #' CV1$`Coefficients`[[1]] %>% kable()
  #' # Metrics
  #' CV1 %>% select(1:15) %>% kable(digits = 5)
  #' #' Confusion matrix- important when considering Kappa statistic and calculating accuracy metrics
  #' CV1$`Confusion Matrix`[[1]] %>% kable()
  
  
  
  
  #'  Compute the average of the k recorded errors. This is called the cross-validation error serving as the performance metric for the model.
  #'  Report Accuracy and Kappa statistics!
   
  
  
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
  md_smr18 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + RoadDen + CanopyCover + Landcover_type + (1|ID)", dat =  mdData_smr[mdData_smr$Year == "Year1",])  # + Slope  + Dist2Water  + Dist2Edge
  md_smr19 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Landcover_type + (1|ID)", dat =  mdData_smr[mdData_smr$Year == "Year2",]) # + Dist2Edge
  md_smr20 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + RoadDen + Dist2Water + CanopyCover + Landcover_type + (1|ID)", dat =  mdData_smr[mdData_smr$Year == "Year3",]) # + Dist2Edge  + Slope
  
  md_wtr1819 <- glmm_fn(mod = "Used ~ 1 + Elev + Slope + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)", dat =  mdData_wtr[mdData_wtr$Year == "Year1",]) # + RoadDen  + I(Elev^2)
  md_wtr1920 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)", dat =  mdData_wtr[mdData_wtr$Year == "Year2",])
  md_wtr2021 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)", dat =  mdData_wtr[mdData_wtr$Year == "Year3",])
  
  ####  Elk RSFs  ####
  #'  Dropping HumanMod in elk summer models due to high correlation with other covariates
  elk_smr18 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = elkData_smr[elkData_smr$Year == "Year1",])
  elk_smr19 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = elkData_smr[elkData_smr$Year == "Year2",]) # + Slope
  elk_smr20 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = elkData_smr[elkData_smr$Year == "Year3",])
  
  elk_wtr1819 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = elkData_wtr_reclass[elkData_wtr_reclass$Year == "Year1",])
  elk_wtr1920 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = elkData_wtr_reclass[elkData_wtr_reclass$Year == "Year2",])
  elk_wtr2021 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = elkData_wtr_reclass[elkData_wtr_reclass$Year == "Year3",])
  
  ####  White-tailed Deer RSFs  ####
  wtd_smr18 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + HumanMod + Dist2Edge + Landcover_type + (1|ID)",  dat = wtdData_smr[wtdData_smr$Year == "Year1",]) # + CanopyCover + Dist2Water
  wtd_smr19 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wtdData_smr[wtdData_smr$Year == "Year2",])
  wtd_smr20 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wtdData_smr[wtdData_smr$Year == "Year3",])

  wtd_wtr1819 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + Dist2Edge + Landcover_type + (1|ID)",  dat = wtdData_wtr_reclass[wtdData_wtr_reclass$Year == "Year1",]) # + Dist2Water + CanopyCover + RoadDen + HumanMod
  wtd_wtr1920 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wtdData_wtr_reclass[wtdData_wtr_reclass$Year == "Year2",])
  wtd_wtr2021 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wtdData_wtr_reclass[wtdData_wtr_reclass$Year == "Year3",])
  
  ####  Cougar RSFs  ####
  coug_smr18 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = cougData_smr[cougData_smr$Year == "Year1",])  # + HumanMod
  coug_smr19 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = cougData_smr[cougData_smr$Year == "Year2",])
  coug_smr20 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = cougData_smr[cougData_smr$Year == "Year3",]) #  + HumanMod
  
  coug_wtr1819 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = cougData_wtr_reclass[cougData_wtr_reclass$Year == "Year1",]) # + RoadDen
  coug_wtr1920 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = cougData_wtr_reclass[cougData_wtr_reclass$Year == "Year2",]) # + RoadDen
  coug_wtr2021 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = cougData_wtr_reclass[cougData_wtr_reclass$Year == "Year3",]) # + RoadDen + Dist2Water
  
  ####  Wolf RSFs  ####
  #'  "Other", "Developed", & "Wetland" landcover types causing issues with model
  #'  convergence for all wolf models so lumping all together as one class
  reclass_landcov <- function(locs) {
    locs <- locs %>%
    mutate(
      Landcover_type = as.character(as.factor(Landcover_type)),
      Landcover_type = ifelse(Landcover_type == "Wetland", "Other", Landcover_type)
    )
  locs$Landcover_type <- droplevels(as.factor(locs$Landcover_type))
  locs$Landcover_type <- relevel(locs$Landcover_type, ref = "Forest")
  
  return(locs)
  }
  wolfData_smr_reclass2 <- reclass_landcov(wolfData_smr_reclass)
  wolfData_wtr_reclass2 <- reclass_landcov(wolfData_wtr_reclass)
  
  wolf_smr18 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wolfData_smr_reclass2[wolfData_smr_reclass2$Year == "Year1",])
  wolf_smr19 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wolfData_smr_reclass2[wolfData_smr_reclass2$Year == "Year2",])
  wolf_smr20 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wolfData_smr_reclass2[wolfData_smr_reclass2$Year == "Year3",])
  
  wolf_wtr1819 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wolfData_wtr_reclass2[wolfData_wtr_reclass2$Year == "Year1",])  # + Dist2Water
  wolf_wtr1920 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wolfData_wtr_reclass2[wolfData_wtr_reclass2$Year == "Year2",])  # + Dist2Water
  wolf_wtr2021 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = wolfData_wtr_reclass2[wolfData_wtr_reclass2$Year == "Year3",])
  
  ####  Bobcat RSFs  ####
  #'  Only data for MVBOB90M in smr18--- not enough data to make inference about bobcat resource selection across 2 study areas
  # bob_smr18 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = bobData_smr[bobData_smr$Year == "Year1",])
  bob_smr19 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = bobData_smr_reclass[bobData_smr_reclass$Year == "Year2",])
  bob_smr20 <- glmm_fn(mod = "Used ~ 1 + Elev + Dist2Water + HumanMod + Dist2Edge + Landcover_type + (1|ID)",  dat = bobData_smr_reclass[bobData_smr_reclass$Year == "Year3",])  # + CanopyCover + I(Elev^2) + Slope  + RoadDen
  
  #'  Only data for MVBOB88M & MVBOB90M wtr1819--- not enough data to make inference about bobcat resource selection across 2 study areas
  # bob_wtr1819 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = bobData_wtr[bobData_wtr$Year == "Year1",])
  bob_wtr1920 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + HumanMod + Dist2Edge + Landcover_type + (1|ID)",  dat = bobData_wtr[bobData_wtr$Year == "Year2",]) # + Dist2Water  + CanopyCover
  bob_wtr2021 <- glmm_fn(mod = "Used ~ 1 + Elev + Slope + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID)",  dat = bobData_wtr[bobData_wtr$Year == "Year3",]) # + Dist2Water  + RoadDen  + I(Elev^2)
  
  ####  Coyote RSFs  ####
  #'  Data from only MVCOY68F, NECOY1F, NECOY2M, & NECOY3F in snmr18--- hesitant to extrapolate selection across study areas
  #'  Dropping HumanMod due to high correlation with other covariates
  coy_smr18 <- glmm_fn(mod = "Used ~ 1 + Elev + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)", dat = coyData_smr_reclass[coyData_smr_reclass$Year == "Year1",]) # + I(Elev^2) 
  coy_smr19 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Landcover_type + (1|ID)", dat = coyData_smr_reclass[coyData_smr_reclass$Year == "Year2",]) # + Dist2Edge
  coy_smr20 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)", dat = coyData_smr_reclass[coyData_smr_reclass$Year == "Year3",])
  
  #'  Data from only MVCOY68F, NECOY1F, NECOY2M, NECOY3F & NECOY4M in wtr1819--- hesitant to extrapolate selection across study areas
  coy_wtr1819 <- glmm_fn(mod = "Used ~ 1 + Elev + Dist2Water + Dist2Edge + Landcover_type + (1|ID)", dat = coyData_wtr_reclass[coyData_wtr_reclass$Year == "Year1",]) # + RoadDen + I(Elev^2) + CanopyCover + Slope
  coy_wtr1920 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + Landcover_type + (1|ID)", dat = coyData_wtr_reclass[coyData_wtr_reclass$Year == "Year2",]) # + CanopyCover  + Dist2Edge
  coy_wtr2021 <- glmm_fn(mod = "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID)", dat = coyData_wtr_reclass[coyData_wtr_reclass$Year == "Year3",])

  #'  Save
  RSF_MD_list <- list(md_smr18, md_smr19, md_smr20, md_wtr1819, md_wtr1920, md_wtr2021)
  RSF_ELK_list <- list(elk_smr18, elk_smr19, elk_smr20, elk_wtr1819, elk_wtr1920, elk_wtr2021)
  RSF_WTD_list <- list(wtd_smr18, wtd_smr19, wtd_smr20, wtd_wtr1819, wtd_wtr1920, wtd_wtr2021)
  RSF_COUG_list <- list(coug_smr18, coug_smr19, coug_smr20, coug_wtr1819, coug_wtr1920, coug_wtr2021)
  RSF_WOLF_list <- list(wolf_smr18, wolf_smr19, wolf_smr20, wolf_wtr1819, wolf_wtr1920, wolf_wtr2021)
  RSF_BOB_list <- list(bob_smr18, bob_smr19, bob_smr20, bob_wtr1819, bob_wtr1920, bob_wtr2021)
  RSF_COY_list <- list(coy_smr18, coy_smr19, coy_smr20, coy_wtr1819, coy_wtr1920, coy_wtr2021)
  
  save(RSF_MD_list, file = paste0("./Outputs/RSF_output/RSF_MD_list_", Sys.Date(), ".RData"))
  save(RSF_ELK_list, file = paste0("./Outputs/RSF_output/RSF_ELK_list_", Sys.Date(), ".RData"))
  save(RSF_WTD_list, file = paste0("./Outputs/RSF_output/RSF_WTD_list_", Sys.Date(), ".RData"))
  save(RSF_COUG_list, file = paste0("./Outputs/RSF_output/RSF_COUG_list_", Sys.Date(), ".RData"))
  save(RSF_WOLF_list, file = paste0("./Outputs/RSF_output/RSF_WOLF_list_", Sys.Date(), ".RData"))
  save(RSF_BOB_list, file = paste0("./Outputs/RSF_output/RSF_BOB_list_", Sys.Date(), ".RData"))
  save(RSF_COY_list, file = paste0("./Outputs/RSF_output/RSF_COY_list_", Sys.Date(), ".RData"))
  
  
  # save(md_smr18, file = paste0("./Outputs/RSF_output/md_RSF_smr18_", Sys.Date(), ".RData"))
  # save(md_smr19, file = paste0("./Outputs/RSF_output/md_RSF_smr19_", Sys.Date(), ".RData"))
  # save(md_smr20, file = paste0("./Outputs/RSF_output/md_RSF_smr20_", Sys.Date(), ".RData"))
  # save(md_wtr1819, file = paste0("./Outputs/RSF_output/md_RSF_wtr1819_", Sys.Date(), ".RData"))
  # save(md_wtr1920, file = paste0("./Outputs/RSF_output/md_RSF_wtr1920_", Sys.Date(), ".RData"))
  # save(md_wtr2021, file = paste0("./Outputs/RSF_output/md_RSF_wtr2021_", Sys.Date(), ".RData"))
  # save(elk_smr18, file = paste0("./Outputs/RSF_output/elk_RSF_smr18_", Sys.Date(), ".RData"))
  # save(elk_smr19, file = paste0("./Outputs/RSF_output/elk_RSF_smr19_", Sys.Date(), ".RData"))
  # save(elk_smr20, file = paste0("./Outputs/RSF_output/elk_RSF_smr20_", Sys.Date(), ".RData"))
  # save(elk_wtr1819, file = paste0("./Outputs/RSF_output/elk_RSF_wtr1819_", Sys.Date(), ".RData"))
  # save(elk_wtr1920, file = paste0("./Outputs/RSF_output/elk_RSF_wtr1920_", Sys.Date(), ".RData"))
  # save(elk_wtr2021, file = paste0("./Outputs/RSF_output/elk_RSF_wtr2021_", Sys.Date(), ".RData"))
  # save(wtd_smr18, file = paste0("./Outputs/RSF_output/wtd_RSF_smr18_", Sys.Date(), ".RData"))
  # save(wtd_smr19, file = paste0("./Outputs/RSF_output/wtd_RSF_smr19_", Sys.Date(), ".RData"))
  # save(wtd_smr20, file = paste0("./Outputs/RSF_output/wtd_RSF_smr20_", Sys.Date(), ".RData"))
  # save(wtd_wtr1819, file = paste0("./Outputs/RSF_output/wtd_RSF_wtr1819_", Sys.Date(), ".RData"))
  # save(wtd_wtr1920, file = paste0("./Outputs/RSF_output/wtd_RSF_wtr1920_", Sys.Date(), ".RData"))
  # save(wtd_wtr2021, file = paste0("./Outputs/RSF_output/wtd_RSF_wtr2021_", Sys.Date(), ".RData"))
  # save(coug_smr18, file = paste0("./Outputs/RSF_output/coug_RSF_smr18_", Sys.Date(), ".RData"))
  # save(coug_smr19, file = paste0("./Outputs/RSF_output/coug_RSF_smr19_", Sys.Date(), ".RData"))
  # save(coug_smr20, file = paste0("./Outputs/RSF_output/coug_RSF_smr20_", Sys.Date(), ".RData"))
  # save(coug_wtr1819, file = paste0("./Outputs/RSF_output/coug_RSF_wtr1819_", Sys.Date(), ".RData"))
  # save(coug_wtr1920, file = paste0("./Outputs/RSF_output/coug_RSF_wtr1920_", Sys.Date(), ".RData"))
  # save(coug_wtr2021, file = paste0("./Outputs/RSF_output/coug_RSF_wtr2021_", Sys.Date(), ".RData"))
  # save(wolf_smr18, file = paste0("./Outputs/RSF_output/wolf_RSF_smr18_", Sys.Date(), ".RData"))
  # save(wolf_smr19, file = paste0("./Outputs/RSF_output/wolf_RSF_smr19_", Sys.Date(), ".RData"))
  # save(wolf_smr20, file = paste0("./Outputs/RSF_output/wolf_RSF_smr20_", Sys.Date(), ".RData"))
  # save(wolf_wtr1819, file = paste0("./Outputs/RSF_output/wolf_RSF_wtr1819_", Sys.Date(), ".RData"))
  # save(wolf_wtr1920, file = paste0("./Outputs/RSF_output/wolf_RSF_wtr1920_", Sys.Date(), ".RData"))
  # save(wolf_wtr2021, file = paste0("./Outputs/RSF_output/wolf_RSF_wtr2021_", Sys.Date(), ".RData"))
  # save(coy_smr18, file = paste0("./Outputs/RSF_output/coy_RSF_smr18_", Sys.Date(), ".RData"))
  # save(coy_smr19, file = paste0("./Outputs/RSF_output/coy_RSF_smr19_", Sys.Date(), ".RData"))
  # save(coy_smr20, file = paste0("./Outputs/RSF_output/coy_RSF_smr20_", Sys.Date(), ".RData"))
  # save(coy_wtr1819, file = paste0("./Outputs/RSF_output/coy_RSF_wtr1819_", Sys.Date(), ".RData"))
  # save(coy_wtr1920, file = paste0("./Outputs/RSF_output/coy_RSF_wtr1920_", Sys.Date(), ".RData"))
  # save(coy_wtr2021, file = paste0("./Outputs/RSF_output/coy_RSF_wtr2021_", Sys.Date(), ".RData"))
  # save(bob_smr18, file = paste0("./Outputs/RSF_output/bob_RSF_smr18_", Sys.Date(), ".RData"))
  # save(bob_smr19, file = paste0("./Outputs/RSF_output/bob_RSF_smr19_", Sys.Date(), ".RData"))
  # save(bob_smr20, file = paste0("./Outputs/RSF_output/bob_RSF_smr20_", Sys.Date(), ".RData"))
  # save(bob_wtr1819, file = paste0("./Outputs/RSF_output/bob_RSF_wtr1819_", Sys.Date(), ".RData"))
  # save(bob_wtr1920, file = paste0("./Outputs/RSF_output/bob_RSF_wtr1920_", Sys.Date(), ".RData"))
  # save(bob_wtr2021, file = paste0("./Outputs/RSF_output/bob_RSF_wtr2021_", Sys.Date(), ".RData"))
  
  
  ####  Summary tables  ####
  #'  Save model outputs in table format 
  #'  Functions extract outputs for each sub-model and appends species/season info
  
  #'  Pull out RSF results
  load("./Outputs/RSF_output/RSF_MD_list_2021-12-18.RData")
  load("./Outputs/RSF_output/RSF_ELK_list_2021-12-18.RData")
  load("./Outputs/RSF_output/RSF_WTD_list_2021-12-18.RData")
  load("./Outputs/RSF_output/RSF_COUG_list_2021-12-18.RData") 
  load("./Outputs/RSF_output/RSF_WOLF_list_2021-12-18.RData")
  load("./Outputs/RSF_output/RSF_BOB_list_2021-12-18.RData")
  load("./Outputs/RSF_output/RSF_COY_list_2021-12-18.RData")


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
        Pval = round(pval, rounddig)) 
    rownames(out) <- NULL
    return(out)
  }
  md_s1819_rsf <- rsf_out(md_global_smr, "Mule Deer", "Summer")
  md_w1820_rsf <- rsf_out(md_global_wtr, "Mule Deer", "Winter")
  elk_s1819_rsf <- rsf_out(elk_global_smr, "Elk", "Summer")
  elk_w1820_rsf <- rsf_out(elk_global_wtr, "Elk", "Winter")
  wtd_s1819_rsf <- rsf_out(wtd_global_smr, "White-tailed Deer", "Summer")
  wtd_w1820_rsf <- rsf_out(wtd_global_wtr, "White-tailed Deer", "Winter")
  coug_s1819_rsf <- rsf_out(coug_global_smr, "Cougar", "Summer")
  coug_w1820_rsf <- rsf_out(coug_global_wtr, "Cougar", "Winter")
  wolf_s1819_rsf <- rsf_out(wolf_global_smr, "Wolf", "Summer")
  wolf_w1820_rsf <- rsf_out(wolf_global_wtr, "Wolf", "Winter")
  bob_s1819_rsf <- rsf_out(bob_global_smr, "Bobcat", "Summer")
  bob_w1820_rsf <- rsf_out(bob_global_wtr, "Bobcat", "Winter")
  coy_s1819_rsf <- rsf_out(coy_global_smr, "Coyote", "Summer")
  coy_w1820_rsf <- rsf_out(coy_global_wtr, "Coyote", "Winter")
  
  md_sSAonly_rsf <- rsf_out(md_SA_only_smr, "Mule Deer", "Summer")
  # write.csv(md_sSAonly_rsf, paste0("./Outputs/Tables/RSF_Results_MuleDeer_SAonly_", Sys.Date(), ".csv"))  
  
  

  #'  Merge into larger data frames for easy comparison
  summer_rsf <- rbind(bob_s1819_rsf, coug_s1819_rsf, coy_s1819_rsf, wolf_s1819_rsf,
                      elk_s1819_rsf, md_s1819_rsf, wtd_s1819_rsf) 
  winter_rsf <- rbind(bob_w1820_rsf, coug_w1820_rsf, coy_w1820_rsf, wolf_w1820_rsf,
                      elk_w1820_rsf, md_w1820_rsf, wtd_w1820_rsf) 
  rsf_results <- rbind(summer_rsf, winter_rsf) %>%
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
    separate("Slope", c("Slope (SE)", "Slope Pval"), sep = "_") %>%
    separate("PercForMix", c("PercForMix (SE)", "PercForMix Pval"), sep = "_") %>%
    separate("PercXGrass", c("PercXGrass (SE)", "PercXGrass Pval"), sep = "_") %>%
    separate("PercXShrub", c("PercXShrub (SE)", "PercXShrub Pval"), sep = "_") %>%
    separate("RoadDen", c("Road Density (SE)", "Road Density Pval"), sep = "_") %>%
    # separate("HumanMod", c("HumanMod (SE)", "HumanMod Pval"), sep = "_") %>%
    arrange(match(Species, c("Bobcat", "Cougar", "Coyote", "Wolf", "Mule Deer", "Elk", "White-tailed Deer"))) %>%
    arrange(match(Season, c("Summer", "Winter")))
  
  
  #'  Save!
  write.csv(rsf_results, paste0("./Outputs/Tables/RSF_Results_NoHM_", Sys.Date(), ".csv"))  #'  KEEP TRACK of whether Human Modified was excluded from models!
  write.csv(rsf_results_wide, paste0("./Outputs/Tables/RSF_Results_wide_NoHM_", Sys.Date(), ".csv"))
  
  save.image("./Outputs/RSF_script_results.RData")
  
  
  
  #'  SHOULD I BE CONSIDERING A VARIENCE INFLATION FACTOR ON THESE ESTIMATES????
  