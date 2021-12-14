  #'  ============================================
  #'  Resource Selection Functions (cam vs collar analysis)
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  May 2021
  #'  ============================================
  #'  Script to build resource selection function models for comparison to 
  #'  occupancy models and HMMs.
  #'  
  #'  Cleaned telemetry and covariate data were prepared with
  #'  Collar_RSF_DataPrep.R script which take FOREVER to run. 
  #'  ============================================
  
  #'  Clear memory
  rm(list=ls())

  #'  Load packages for selecting available points
  library(tidyverse)
  library(car)
  library(caret)
  library(lme4)
  
  #'  Load used and available locations, and covariate data
  #'  3rd Order Selection
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
    locs$Landcover_type <- as.character(as.factor(locs$Landcover_type))
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
        Landcover_type = ifelse(Landcover_type == "Developed", "Developed", Landcover_type)
      )
  
    locs <- as.data.frame(locs)
    
    return(locs)
  }
  #'  Run season & species-specific data through prep function
  mdData_smr <- spp_dataPrep(md_dat_all[md_dat_all$Season == "Summer18" | md_dat_all$Season == "Summer19",])
  mdData_wtr <- spp_dataPrep(md_dat_all[md_dat_all$Season == "Winter1819" | md_dat_all$Season == "Winter1920",])
  elkData_smr <- spp_dataPrep(elk_dat_all[elk_dat_all$Season == "Summer18" | elk_dat_all$Season == "Summer19",])
  elkData_wtr <- spp_dataPrep(elk_dat_all[elk_dat_all$Season == "Winter1819" | elk_dat_all$Season == "Winter1920",])
  wtdData_smr <- spp_dataPrep(wtd_dat_all[wtd_dat_all$Season == "Summer18" | wtd_dat_all$Season == "Summer19",])
  wtdData_wtr <- spp_dataPrep(wtd_dat_all[wtd_dat_all$Season == "Winter1819" | wtd_dat_all$Season == "Winter1920",])
  cougData_smr <- spp_dataPrep(coug_dat_all[coug_dat_all$Season == "Summer18" | coug_dat_all$Season == "Summer19",])
  cougData_wtr <- spp_dataPrep(coug_dat_all[coug_dat_all$Season == "Winter1819" | coug_dat_all$Season == "Winter1920",])
  wolfData_smr <- spp_dataPrep(wolf_dat_all[wolf_dat_all$Season == "Summer18" | wolf_dat_all$Season == "Summer19",])
  wolfData_wtr <- spp_dataPrep(wolf_dat_all[wolf_dat_all$Season == "Winter1819" | wolf_dat_all$Season == "Winter1920",])
  bobData_smr <- spp_dataPrep(bob_dat_all[bob_dat_all$Season == "Summer18" | bob_dat_all$Season == "Summer19",])
  bobData_wtr <- spp_dataPrep(bob_dat_all[bob_dat_all$Season == "Winter1819" | bob_dat_all$Season == "Winter1920",])
  coyData_smr <- spp_dataPrep(coy_dat_all[coy_dat_all$Season == "Summer18" | coy_dat_all$Season == "Summer19",])
  coyData_wtr <- spp_dataPrep(coy_dat_all[coy_dat_all$Season == "Winter1819" | coy_dat_all$Season == "Winter1920",])
  
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
  #'  Elevation & Human Modified correlated in MD smr/wtr, ELK smr, & COY smr- nixing Human Mod for those models
  #'  % Shrub correlated with Elevation, RoadDen, and Human Modified in MD smr- nixing % Shrub
  #'  Dist2Edg correlated with % Forest & % Grass in ELK wtr- nixing Dist2Edge for this model
  #'  % Forest & Grass correlated in COUG and WOLF wtr- nixing % Grass for those models
  #'  % Grass & Shrub correlated for BOb smr- nixing % Grass for this model
  #'  % Forest & Canopy Cover, % Forest & Shrub correlated for BOB wtr- nixing Canopy & % Shrub for this model
  
  
  
  
  
  
  #'  Resource Selection Function Models
  #'  ==================================
  #'  Functions to run logistic mixed effects models that include random effect 
  #'  for individual. Habitat covariates excluded based on species (see notes 
  #'  above) and convergence issues. Annual models run separately so predicted 
  #'  distributions are specific to the species, season, and year.
  
  ####  Mule Deer RSF  ####
  #'  SUMMER 2018
  #'  Dropping HumanMod due to high correlation with other covariates
  md_smr18 <- glmer(Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID),  
                         data = mdData_smr[mdData_smr$Year == "Year1",], family = binomial(link = "logit"))
  summary(md_smr18)
  car::vif(md_smr18)
  #'  SUMMER 2019
  #'  Dropping HumanMod due to high correlation with other covariates
  md_smr19 <- glmer(Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID), 
                    data = mdData_smr[mdData_smr$Year == "Year2",], family = binomial(link = "logit"))
  summary(md_smr19)
  car::vif(md_smr19)
  #'  SUMMER 2020
  #'  Dropping HumanMod due to high correlation with other covariates
  md_smr20 <- glmer(Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID), 
                    data = mdData_smr[mdData_smr$Year == "Year3",], family = binomial(link = "logit"))
  summary(md_smr20)
  car::vif(md_smr20)

  #'  WINTER 2018-2019
  #'  Dropping HumanMod due to high correlation with other covariates
  md_wtr1819 <- glmer(Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID),
                          data = mdData_wtr[mdData_wtr$Year == "Year1",], family = binomial(link = "logit"))
  summary(md_wtr1819)
  car::vif(md_wtr1819)
  #'  WINTER 2019-2020
  #'  Dropping HumanMod due to high correlation with other covariates
  md_wtr1920 <- glmer(Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID),
                      data = mdData_wtr[mdData_wtr$Year == "Year2",], family = binomial(link = "logit"))
  summary(md_wtr1920)
  car::vif(md_wtr1920)
  #'  WINTER 2020-2021
  #'  Dropping HumanMod due to high correlation with other covariates
  md_wtr2021 <- glmer(Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID),
                      data = mdData_wtr[mdData_wtr$Year == "Year3",], family = binomial(link = "logit"))
  summary(md_wtr2021)
  car::vif(md_wtr2021)
  
  
  
  ####  Elk RSF  ####
  #'  SUMMER 2018
  #'  Dropping HumanMod due to high correlation with other covariates
  elk_smr18 <- glmer(Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID),  
                        data = elkData_smr[elkData_smr  == "Year1",], family = binomial(link = "logit"))
  summary(elk_smr18)
  car::vif(elk_smr18)
  #'  SUMMER 2019
  #'  Dropping HumanMod due to high correlation with other covariates
  elk_smr19 <- glmer(Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + PLandcover_type + (1|ID),  
                     data = elkData_smr[elkData_smr  == "Year2",], family = binomial(link = "logit"))
  summary(elk_smr19)
  car::vif(elk_smr19)
  #'  SUMMER 2020
  #'  Dropping HumanMod due to high correlation with other covariates
  elk_smr20 <- glmer(Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID),  
                     data = elkData_smr[elkData_smr  == "Year3",], family = binomial(link = "logit"))
  summary(elk_smr20)
  car::vif(elk_smr20)
  
  #'  WINTER 2018-2019
  elk_wtr1819 <- glmer(Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Dist2Edge + Landcover_type + (1|ID), 
                          data = elkData_wtr[elkData_wtr$Year == "Year1",], family = binomial(link = "logit"))
  summary(elk_wtr1819)
  car::vif(elk_wtr1819)
  #'  WINTER 2019-2020
  #'  Dropping Dist2Edge due to high correlation with other covariates
  elk_wtr1920 <- glmer(Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Landcover_type + (1|ID), 
                       data = elkData_wtr[elkData_wtr$Year == "Year2",], family = binomial(link = "logit"))
  summary(elk_wtr1920)
  car::vif(elk_wtr1920)
  #'  WINTER 2020-2021
  #'  Dropping Dist2Edge due to high correlation with other covariates
  elk_wtr2021 <- glmer(Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + HumanMod + CanopyCover + Landcover_type + (1|ID), 
                       data = elkData_wtr[elkData_wtr$Year == "Year3",], family = binomial(link = "logit"))
  summary(elk_wtr2021)
  car::vif(elk_wtr2021)
  
  
  
  ####  White-tailed Deer RSF  ####
  #'  SUMMER 2018
  #'  Dropping HumanMod due to high correlation with other covariates
  wtd_global_smr <- glmer(Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1|ID), 
                          data = wtdData_smr, family = binomial(link = "logit"))
  summary(wtd_global_smr)
  car::vif(wtd_global_smr)

  #'  WINTERS 2018-2019 & 2019-2020
  wtd_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID) + (1|Year), 
                          data = wtdData_wtr, family = binomial(link = "logit"))
  summary(wtd_global_wtr)
  car::vif(wtd_global_wtr)
  
  
  ####  Cougar RSF  ####
  #'  Random effect for individual & year
  #'  SUMMERS 2018 & 2019
  coug_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID) + (1|Year), 
                          data = cougData_smr, family = binomial(link = "logit"))
  summary(coug_global_smr)
  car::vif(coug_global_smr)
  
  #'  WINTERS 2018-2019 & 2019-2020
  #'  Dropping PercXGrass due to high correlation with PercForMix
  # coug_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID) + (1|Year), 
  #                         data = cougData_wtr, family = binomial(link = "logit"))
  coug_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXShrub + RoadDen + (1|ID) + (1|Year), 
                           data = cougData_wtr, family = binomial(link = "logit"))
  summary(coug_global_wtr)
  car::vif(coug_global_wtr)
  
  
  ####  Wolf RSF  ####
  #'  Random effect for individual; too few collars active both years for year effect
  #'  SUMMERS 2018 & 2019
  wolf_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID), 
                       data = wolfData_smr, family = binomial(link = "logit"))  #1/9 collars active both years
  summary(wolf_global_smr)
  car::vif(wolf_global_smr)
  
  #'  WINTERS 2018-2019 & 2019-2020
  #'  Dropping PercXGrass due to high correlation with PercForMix  
  # wolf_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID), 
  #                          data = wolfData_wtr, family = binomial(link = "logit")) #no collars active both years
  wolf_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXShrub + RoadDen + (1|ID), 
                           data = wolfData_wtr, family = binomial(link = "logit")) #no collars active both years
  summary(wolf_global_wtr)
  car::vif(wolf_global_wtr)

  
  ####  Bobcat RSF  ####
  #'  Random effect for individual; too few collars active both years for year effect
  #'  SUMMERS 2018 & 2019
  bob_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID), 
                           data = bobData_smr, family = binomial(link = "logit")) #1/10 collars active both years
  summary(bob_global_smr)
  car::vif(bob_global_smr)
  
  #'  WINTERS 2018-2019 & 2019-2020  
  #'  Dropping PercXShrub due to high correlation with PercForMix
  # bob_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID), 
  #                          data = bobData_wtr, family = binomial(link = "logit")) #1/11 collars active both years
  bob_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + RoadDen + (1|ID), 
                          data = bobData_wtr, family = binomial(link = "logit")) #1/11 collars active both years
  summary(bob_global_wtr)
  car::vif(bob_global_wtr)
  #  Winter1920 MVBOB71M disperses somewhere on the Colville reservation between 
  #  two study areas and make the OK bobcat MCP super big as a result. Cut his
  #  winter 1920 data completely (exclude from RSF and MCP).
  
  
  ####  Coy RSF  ####
  
  #https://www.youtube.com/watch?v=BQ1VAZ7jNYQ helpful video but I don't think you can use this for models with random effects
  coydata = coyData_smr[coyData_smr$Year == "Year1",]
  library(caret)
  #'  Partitioning of the data- create index matrix of selected values
  #'  Set seed to split data
  set.seed(2021)
  #'  Randomly split data once into training and test data (80:20)
  index <- createDataPartition(coydata$Used, p = 0.8, list = FALSE, times = 1)
  #'  Create training and test data frame based on indexed data
  train_df <- coydata[index,]
  test_df <- coydata[-index,]
  #'  Specify type of training method and number of folds
  ctrlspecs <- trainControl(method = "cv", number = 5, savePredictions = "all",
                            classProbs = TRUE) # save predictions & class probabilities
  #'  Set random seed split training data into folds
  set.seed(2021)
  #'  Specify the statistical model
  model1 <- train(Used ~ Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover,
                  data = train_df, method = "glm", family = binomial, trCrontrol = ctrlspecs)
  #'  Predict outcome using model from train_df applied to the test_df
  predicted <- predict(model1, newdata = test_df)
  #'  Create confusion matrix to get the classification accuracy, etc.
  confusionMatrix(data = predicted, test_df$Used)
  
  
  library(cvms)
  library(groupdata2)
  mddata = mdData_smr[mdData_smr$Year == "Year1",]
  #'  Set seed to split data
  set.seed(2021)
  #'  Randomly split data once into training and test data (80:20) using caret package
  index <- createDataPartition(mddata$Used, p = 0.8, list = FALSE, times = 1)
  #'  Create training and test data frame based on indexed data
  train_df <- mddata[index,]
  test_df <- mddata[-index,]
  #'  Set seed to create folds
  set.seed(2021)
  #'  Partition training data into balanced folds using groupdata2 package
  train_df <- fold(train_df, k = 5, cat_col = "Used")   # Can't group with id_col b/c number of observations differ with each ID
  #'  Define model
  formulas <- "Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + Landcover_type + (1 | ID)"
  #'  Run K-fold cross validation with training data and model
  coytst <- cross_validate(train_df, formulas = formulas, family = "binomial", REML = FALSE)
  coytst
  #'  Report Accuracy and Kappa statistics!
  

  
  
  
  #'  SUMMER 2018
  #'  Dropping HumanMod due to high correlation with other covariates
  coy_smr18 <- glmer(Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + PercForMix + PercXGrass + PercXShrub + (1|ID), 
                           data = coyData_smr[coyData_smr$Year == "Year1",], family = binomial(link = "logit")) 
  summary(coy_smr18)
  car::vif(coy_smr18)
  #'  SUMMER 2019
  #'  Dropping HumanMod due to high correlation with other covariates
  coy_smr19 <- glmer(Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + PercForMix + PercXGrass + PercXShrub + (1|ID), 
                     data = coyData_smr[coyData_smr$Year == "Year2",], family = binomial(link = "logit")) 
  summary(coy_smr19)
  car::vif(coy_smr19)
  #'  SUMMER 2020
  #'  Dropping HumanMod due to high correlation with other covariates
  coy_smr20 <- glmer(Used ~ 1 + Elev + I(Elev^2) + Slope + RoadDen + Dist2Water + CanopyCover + Dist2Edge + PercForMix + PercXGrass + PercXShrub + (1|ID), 
                     data = coyData_smr[coyData_smr$Year == "Year3",], family = binomial(link = "logit")) 
  summary(coy_smr20)
  car::vif(coy_smr20)
  
  #'  WINTERS 2018-2019 & 2019-2020  
  #'  Dropping PercXGrass due to high correlation with PercForMix
  coy_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXShrub + RoadDen + (1|ID) + (1|Year), 
                          data = coyData_wtr, family = binomial(link = "logit")) 
  summary(coy_global_wtr)
  car::vif(coy_global_wtr)
  
  
  #'  Save
  save(md_global_smr, file = paste0("./Outputs/RSF_output/md_RSF_smr_NoHM_", Sys.Date(), ".RData"))  #'  KEEP TRACK of whether human modified was excluded from models!
  save(md_global_wtr, file = paste0("./Outputs/RSF_output/md_RSF_wtr_NoHM_", Sys.Date(), ".RData"))
  save(elk_global_smr, file = paste0("./Outputs/RSF_output/elk_RSF_smr_NoHM_", Sys.Date(), ".RData"))
  save(elk_global_wtr, file = paste0("./Outputs/RSF_output/elk_RSF_wtr_NoHM_", Sys.Date(), ".RData"))
  save(wtd_global_smr, file = paste0("./Outputs/RSF_output/wtd_RSF_smr_NoHM_", Sys.Date(), ".RData"))
  save(wtd_global_wtr, file = paste0("./Outputs/RSF_output/wtd_RSF_wtr_NoHM_", Sys.Date(), ".RData"))
  save(coug_global_smr, file = paste0("./Outputs/RSF_output/coug_RSF_smr_NoHM_", Sys.Date(), ".RData"))
  save(coug_global_wtr, file = paste0("./Outputs/RSF_output/coug_RSF_wtr_NoHM_", Sys.Date(), ".RData"))
  save(wolf_global_smr, file = paste0("./Outputs/RSF_output/wolf_RSF_smr_NoHM_", Sys.Date(), ".RData"))
  save(wolf_global_wtr, file = paste0("./Outputs/RSF_output/wolf_RSF_wtr_NoHM_", Sys.Date(), ".RData"))
  save(bob_global_smr, file = paste0("./Outputs/RSF_output/bob_RSF_smr_NoHM_", Sys.Date(), ".RData"))
  save(bob_global_wtr, file = paste0("./Outputs/RSF_output/bob_RSF_wtr_NoHM_", Sys.Date(), ".RData"))
  save(coy_global_smr, file = paste0("./Outputs/RSF_output/coy_RSF_smr_NoHM_", Sys.Date(), ".RData"))
  save(coy_global_wtr, file = paste0("./Outputs/RSF_output/coy_RSF_wtr_NoHM_", Sys.Date(), ".RData"))
  
  save(md_SA_only_smr, file = paste0("./Outputs/RSF_output/md_RSF_smr_SAonly_", Sys.Date(), ".RData"))  
  
  
  ####  Summary tables  ####
  #'  Save model outputs in table format 
  #'  Functions extract outputs for each sub-model and appends species/season info
  
  #'  Pull out RSF results
  load("./Outputs/RSF_output/md_RSF_smr_noHM_2021-09-23.RData")  #' Make sure I read in the right dataset 2021-09-13
  load("./Outputs/RSF_output/md_RSF_wtr_noHM_2021-09-23.RData")
  load("./Outputs/RSF_output/elk_RSF_smr_noHM_2021-09-23.RData")
  load("./Outputs/RSF_output/elk_RSF_wtr_noHM_2021-09-23.RData") 
  load("./Outputs/RSF_output/wtd_RSF_smr_noHM_2021-09-23.RData")
  load("./Outputs/RSF_output/wtd_RSF_wtr_noHM_2021-09-23.RData")
  load("./Outputs/RSF_output/coug_RSF_smr_noHM_2021-09-23.RData")
  load("./Outputs/RSF_output/coug_RSF_wtr_noHM_2021-09-23.RData")
  load("./Outputs/RSF_output/wolf_RSF_smr_noHM_2021-10-30.RData") # note these are different
  load("./Outputs/RSF_output/wolf_RSF_wtr_noHM_2021-10-30.RData") # excludes dispersal events
  load("./Outputs/RSF_output/bob_RSF_smr_noHM_2021-09-23.RData")
  load("./Outputs/RSF_output/bob_RSF_wtr_noHM_2021-09-23.RData")
  load("./Outputs/RSF_output/coy_RSF_smr_noHM_2021-09-23.RData")
  load("./Outputs/RSF_output/coy_RSF_wtr_noHM_2021-09-23.RData")
  
  load("./Outputs/RSF_output/md_RSF_smr_SAonly_2021-11-10.RData")
  

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
  