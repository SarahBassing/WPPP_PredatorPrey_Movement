  #'  ===================================
  #'  Hidden Markov Movement Models 
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  February 2022
  #'  ===================================
  #'  Script to run hidden Markov movement models for deer, elk, cougars, wolves, 
  #'  coyotes, and bobcats for summer & winter, July 2018 - March 2021. 
  #'  Data were collected & generously provided by WPPP collaborators including
  #'  M. Devivo, B. Kertson, T.Ganz, T.Roussin, L.Satterfield, B.Windell, & others. 
  #'  Code adapted from momentuHMM GitHub, J.Merkel Movement Workshop, L.Satterfield, 
  #'  & R.Emmet.
  #'  
  #'  Cleaned telemetry and covariate data were prepared for HMMs with the
  #'  Collar_Movement_DataPrep.R script which took FOREVER to run so only do once.
  #'  Covariates extracted using Collar_Covariate_Extraction.R script.
  #'  ============================================
  
  #'  Clear memory
  rm(list=ls())

  #'  Load libraries
  library(momentuHMM)
  library(rgdal)
  library(ggplot2)
  library(tidyverse)

  #'  Load crwOut & covaraite data
  load("./Outputs/Telemetry_crwOut/crwOut_ALL_2022-03-14.RData")    
  load("./Outputs/Telemetry_covs/spp_telem_covs_2023-04-08.RData")  #2022-05-23 cougar/wolf RSF combined study areas
  
  #'  Merge datasets and create momentuHMMData object
  #'  Data merged and scaled by study area separately b/c different species collared
  #'  in each study area- can't test effect of study area-specific species
  #'  across both study areas.
  #'  OKANOGAN data sets
  spp_dataPrep_OK <- function(crwOut, telem_covs){
    #'  Merge crawlOut data with extracted covariate data
    crwlMerge <- crawlMerge(crwOut, telem_covs, Time.name = "time")
    #'  Make categorical variables factors
    crwlMerge$crwPredict$StudyArea <- as.factor(crwlMerge$crwPredict$StudyArea)
    crwlMerge$crwPredict$Sex <- as.factor(crwlMerge$crwPredict$Sex)
    crwlMerge$crwPredict$Season <- as.factor(crwlMerge$crwPredict$Season)
    crwlMerge$crwPredict$SnowCover <- as.factor(crwlMerge$crwPredict$SnowCover)
    crwlMerge$crwPredict$daytime <- as.factor(crwlMerge$crwPredict$daytime)
    #'  Standardize continuous variables
    crwlMerge$crwPredict$Dist2Road <- scale(crwlMerge$crwPredict$Dist2Road)
    # crwlMerge$crwPredict$NDVI <- scale(crwlMerge$crwPredict$NDVI)
    crwlMerge$crwPredict$PercOpen <- scale(crwlMerge$crwPredict$PercOpen)
    crwlMerge$crwPredict$TRI <- scale(crwlMerge$crwPredict$TRI)
    crwlMerge$crwPredict$MD_RSF <- scale(crwlMerge$crwPredict$MD_RSF)
    #' crwlMerge$crwPredict$ELK_RSF <- scale(crwlMerge$crwPredict$ELK_RSF)   # not in OK data sets
    #' crwlMerge$crwPredict$WTD_RSF <- scale(crwlMerge$crwPredict$WTD_RSF)   # not in OK data sets
    crwlMerge$crwPredict$COUG_RSF <- scale(crwlMerge$crwPredict$COUG_RSF)
    crwlMerge$crwPredict$WOLF_RSF <- scale(crwlMerge$crwPredict$WOLF_RSF)
    # crwlMerge$crwPredict$BOB_RSF <- scale(crwlMerge$crwPredict$BOB_RSF)
    # crwlMerge$crwPredict$COY_RSF <- scale(crwlMerge$crwPredict$COY_RSF)
    crwlMerge$crwPredict$hour <- as.integer(crwlMerge$crwPredict$hour)
    crwlMerge$crwPredict$hour_fix <- as.integer(crwlMerge$crwPredict$hour_fix)
    crwlMerge$crwPredict$hour3 <- as.integer(crwlMerge$crwPredict$hour3)
    #'  Prep crwlMerge data for fitHMM function
    Data <- prepData(data = crwlMerge, covNames = c("Dist2Road", "PercOpen", #"NDVI", 
                                                    "SnowCover", "TRI", "MD_RSF", 
                                                    "COUG_RSF", "WOLF_RSF", # "BOB_RSF", "COY_RSF", 
                                                    "hour", "hour_fix",
                                                    "hour3", "daytime", "Sex", 
                                                    "StudyArea", "Season"))  # "ELK_RSF", "WTD_RSF", 
    return(Data)
  }
  #'  Run season & species-specific data from the Okanogan through prep function
  #'  Warnings are due to missing data for interpolated locations. prepData 
  #'  command automatically fills in values with closest following value.
  mdData_smr <- spp_dataPrep_OK(crwOut_ALL[[1]], spp_telem_covs[[1]])      # HUGE data set- need lab computer
  mdData_wtr <- spp_dataPrep_OK(crwOut_ALL[[2]], spp_telem_covs[[2]])      # HUGE data set- need lab computer
  cougData_smr_OK <- spp_dataPrep_OK(crwOut_ALL[[7]], spp_telem_covs[[7]]) 
  cougData_wtr_OK <- spp_dataPrep_OK(crwOut_ALL[[8]], spp_telem_covs[[8]])
  wolfData_smr_OK <- spp_dataPrep_OK(crwOut_ALL[[11]], spp_telem_covs[[11]])
  wolfData_wtr_OK <- spp_dataPrep_OK(crwOut_ALL[[12]], spp_telem_covs[[12]])
  bobData_smr_OK <- spp_dataPrep_OK(crwOut_ALL[[15]], spp_telem_covs[[15]])
  bobData_wtr_OK <- spp_dataPrep_OK(crwOut_ALL[[16]], spp_telem_covs[[16]])
  coyData_smr_OK <- spp_dataPrep_OK(crwOut_ALL[[19]], spp_telem_covs[[19]])
  coyData_wtr_OK <- spp_dataPrep_OK(crwOut_ALL[[20]], spp_telem_covs[[20]])

  
  #'  NORTHEAST data sets
  spp_dataPrep_NE <- function(crwOut, telem_covs){
    #'  Merge crawlOut data with extracted covariate data
    crwlMerge <- crawlMerge(crwOut, telem_covs, Time.name = "time")
    #'  Make categorical variables factors
    crwlMerge$crwPredict$StudyArea <- as.factor(crwlMerge$crwPredict$StudyArea)
    crwlMerge$crwPredict$Sex <- as.factor(crwlMerge$crwPredict$Sex)
    crwlMerge$crwPredict$Season <- as.factor(crwlMerge$crwPredict$Season)
    crwlMerge$crwPredict$SnowCover <- as.factor(crwlMerge$crwPredict$SnowCover)
    crwlMerge$crwPredict$daytime <- as.factor(crwlMerge$crwPredict$daytime)
    #' #'  Standardize continuous variables
    crwlMerge$crwPredict$Dist2Road <- scale(crwlMerge$crwPredict$Dist2Road)
    # crwlMerge$crwPredict$NDVI <- scale(crwlMerge$crwPredict$NDVI)
    crwlMerge$crwPredict$PercOpen <- scale(crwlMerge$crwPredict$PercOpen)
    crwlMerge$crwPredict$TRI <- scale(crwlMerge$crwPredict$TRI)
    #' crwlMerge$crwPredict$MD_RSF <- scale(crwlMerge$crwPredict$MD_RSF)   # Not in NE data set  
    crwlMerge$crwPredict$ELK_RSF <- scale(crwlMerge$crwPredict$ELK_RSF)
    crwlMerge$crwPredict$WTD_RSF <- scale(crwlMerge$crwPredict$WTD_RSF)
    crwlMerge$crwPredict$COUG_RSF <- scale(crwlMerge$crwPredict$COUG_RSF)
    crwlMerge$crwPredict$WOLF_RSF <- scale(crwlMerge$crwPredict$WOLF_RSF)
    # crwlMerge$crwPredict$BOB_RSF <- scale(crwlMerge$crwPredict$BOB_RSF)
    # crwlMerge$crwPredict$COY_RSF <- scale(crwlMerge$crwPredict$COY_RSF)
    crwlMerge$crwPredict$hour <- as.integer(crwlMerge$crwPredict$hour)
    crwlMerge$crwPredict$hour_fix <- as.integer(crwlMerge$crwPredict$hour_fix)
    crwlMerge$crwPredict$hour3 <- as.integer(crwlMerge$crwPredict$hour3)
    #'  Prep crwlMerge data for fitHMM function
    Data <- prepData(data = crwlMerge, covNames = c("Dist2Road", "PercOpen", #"NDVI", 
                                                    "SnowCover", "TRI", "ELK_RSF", 
                                                    "WTD_RSF", "COUG_RSF", "WOLF_RSF", #"BOB_RSF", "COY_RSF", 
                                                    "hour", "hour_fix", "hour3", "daytime",
                                                    "Sex", "StudyArea", "Season")) # "MD_RSF",
    return(Data)
  }
  #'  Run season & species-specific data from the Northeast through prep function
  #'  Warnings are due to missing data for interpolated locations. prepData 
  #'  command automatically fills in values with closest following value.
  elkData_smr <- spp_dataPrep_NE(crwOut_ALL[[3]], spp_telem_covs[[3]])       # HUGE data set- need lab computer
  elkData_wtr <- spp_dataPrep_NE(crwOut_ALL[[4]], spp_telem_covs[[4]])
  wtdData_smr <- spp_dataPrep_NE(crwOut_ALL[[5]], spp_telem_covs[[5]])
  wtdData_wtr <- spp_dataPrep_NE(crwOut_ALL[[6]], spp_telem_covs[[6]])
  cougData_smr_NE <- spp_dataPrep_NE(crwOut_ALL[[9]], spp_telem_covs[[9]])
  cougData_wtr_NE <- spp_dataPrep_NE(crwOut_ALL[[10]], spp_telem_covs[[10]])
  wolfData_smr_NE <- spp_dataPrep_NE(crwOut_ALL[[13]], spp_telem_covs[[13]])
  wolfData_wtr_NE <- spp_dataPrep_NE(crwOut_ALL[[14]], spp_telem_covs[[14]])
  bobData_smr_NE <- spp_dataPrep_NE(crwOut_ALL[[17]], spp_telem_covs[[17]])
  bobData_wtr_NE <- spp_dataPrep_NE(crwOut_ALL[[18]], spp_telem_covs[[18]])
  coyData_smr_NE <- spp_dataPrep_NE(crwOut_ALL[[21]], spp_telem_covs[[21]])
  coyData_wtr_NE <- spp_dataPrep_NE(crwOut_ALL[[22]], spp_telem_covs[[22]])
  
  #'  Save data prepped for HMMs
  hmm_data <- list(mdData_smr, mdData_wtr, elkData_smr, elkData_wtr, wtdData_smr, 
                   wtdData_wtr, cougData_smr_OK, cougData_wtr_OK, cougData_smr_NE, 
                   cougData_wtr_NE, wolfData_smr_OK, wolfData_wtr_OK, wolfData_smr_NE, 
                   wolfData_wtr_NE, bobData_smr_OK, bobData_wtr_OK, bobData_smr_NE, 
                   bobData_wtr_NE, coyData_smr_OK, coyData_wtr_OK, coyData_smr_NE, 
                   coyData_wtr_NE)
  # save(hmm_data, file = paste0("./Outputs/Telemetry_crwOut/crwOut_ALL_wCovs_", Sys.Date(), ".RData"))
  load("./Outputs/Telemetry_crwOut/crwOut_ALL_wCovs_2023-04-08.RData") #2022-05-23
  names(hmm_data) <- c("mdData_smr", "mdData_wtr", "elkData_smr", "elkData_wtr", "wtdData_smr", "wtdData_wtr",
                  "cougData_smr_OK", "cougData_wtr_OK", "cougData_smr_NE", "cougData_wtr_NE",
                  "wolfData_smr_OK", "wolfData_wtr_OK", "wolfData_smr_NE", "wolfData_wtr_NE",
                  "bobData_smr_OK", "bobData_wtr_OK", "bobData_smr_NE", "bobData_wtr_NE",
                  "coyData_smr_OK", "coyData_wtr_OK", "coyData_smr_NE", "coyData_wtr_NE")
  
  #' #'  Remove coordinates from crwOut data for publication
  #' drop_coords <- function(dat){
  #'   skinny_df <- dat %>%
  #'     dplyr::select(-c("Long", "Lat", "nu.x", "nu.y", "se.mu.x", "se.nu.x", "se.mu.y", "se.nu.y", "x", "y", "TimeNum", "locType", "ID2"))
  #'   return(skinny_df)
  #' }
  #' hmm_data <- lapply(hmm_data, drop_coords)
  #' names(hmm_data) <- c("mdData_smr", "mdData_wtr", "elkData_smr", "elkData_wtr", "wtdData_smr", "wtdData_wtr",
  #'                 "cougData_smr_OK", "cougData_wtr_OK", "cougData_smr_NE", "cougData_wtr_NE",
  #'                 "wolfData_smr_OK", "wolfData_wtr_OK", "wolfData_smr_NE", "wolfData_wtr_NE",
  #'                 "bobData_smr_OK", "bobData_wtr_OK", "bobData_smr_NE", "bobData_wtr_NE",
  #'                 "coyData_smr_OK", "coyData_wtr_OK", "coyData_smr_NE", "coyData_wtr_NE")
  #' save(hmm_data, file = paste0("./Outputs/Telemetry_crwOut/crwOut_ALL_wCovs_for_pub_", Sys.Date(), ".RData"))
  
  
  #'  Correlation Matrix
  #'  ==================
  #'  Function to create correlation matrix for all continuous covariates at once
  cov_correlation_OK <- function(dat) {
    covs <- dat %>%
      dplyr::select(c("Dist2Road", "PercOpen", "TRI", "MD_RSF", "COUG_RSF", #"NDVI", 
                      "WOLF_RSF"))#, "BOB_RSF", "COY_RSF")) 
    cor_matrix <- cor(covs, use = "complete.obs")
    return(cor_matrix)
  }
  #'  Generate correlation matrix for each species and season
  (md_smr_corr <- cov_correlation_OK(mdData_smr)) # COY & TRI highly correlated (-0.75)
  (md_wtr_corr <- cov_correlation_OK(mdData_wtr)) # BOB & COUG correlated (0.60), COY & TRI correlated (-0.60)
  (coug_smr_OK_corr <- cov_correlation_OK(cougData_smr_OK))
  (coug_wtr_OK_corr <- cov_correlation_OK(cougData_wtr_OK)) # MD & TRI highly correlated (0.76)
  (wolf_smr_OK_corr <- cov_correlation_OK(wolfData_smr_OK))
  (wolf_wtr_OK_corr <- cov_correlation_OK(wolfData_wtr_OK)) # MD & TRI correlated (0.64)
  (bob_smr_OK_corr <- cov_correlation_OK(bobData_smr_OK))
  (bob_wtr_OK_corr <- cov_correlation_OK(bobData_wtr_OK)) # MD & TRI correlated (0.67)
  (coy_smr_OK_corr <- cov_correlation_OK(coyData_smr_OK)) 
  (coy_wtr_OK_corr <- cov_correlation_OK(coyData_wtr_OK)) # MD & TRI highly correlated (0.71)
  
  cov_correlation_NE <- function(dat) {
    covs <- dat %>%
      dplyr::select(c("Dist2Road", "PercOpen", "TRI", "ELK_RSF", "WTD_RSF", #"NDVI", 
                      "COUG_RSF", "WOLF_RSF"))#, "BOB_RSF", "COY_RSF")) 
    cor_matrix <- cor(covs, use = "complete.obs")
    return(cor_matrix)
  }
  #'  Generate correlation matrix for each species and season
  (elk_smr_corr <- cov_correlation_NE(elkData_smr)) # COY & TRI correlated (-0.67) 
  (elk_wtr_corr <- cov_correlation_NE(elkData_wtr))
  (wtd_smr_corr <- cov_correlation_NE(wtdData_smr)) # BOB & WOLF correlated (0.62)
  (wtd_wtr_corr <- cov_correlation_NE(wtdData_wtr)) 
  (coug_smr_NE_corr <- cov_correlation_NE(cougData_smr_NE)) # WTD & TRI correlated (-0.60)
  (coug_wtr_NE_corr <- cov_correlation_NE(cougData_wtr_NE))
  (wolf_smr_NE_corr <- cov_correlation_NE(wolfData_smr_NE))  
  (wolf_wtr_NE_corr <- cov_correlation_NE(wolfData_wtr_NE))
  (bob_smr_NE_corr <- cov_correlation_NE(bobData_smr_NE))
  (bob_wtr_NE_corr <- cov_correlation_NE(bobData_wtr_NE))
  (coy_smr_NE_corr <- cov_correlation_NE(coyData_smr_NE))
  (coy_wtr_NE_corr <- cov_correlation_NE(coyData_wtr_NE))
  
  
  #' #'  Visualize data to inform initial parameter specifications
  mean(mdData_smr$step, na.rm = T); sd(mdData_smr$step, na.rm = T)
  mean(mdData_wtr$step, na.rm = T); sd(mdData_wtr$step, na.rm = T)
  mean(elkData_smr$step, na.rm = T); sd(elkData_smr$step, na.rm = T)
  mean(elkData_wtr$step, na.rm = T); sd(elkData_wtr$step, na.rm = T)
  mean(wtdData_smr$step, na.rm = T); sd(wtdData_smr$step, na.rm = T)
  mean(wtdData_wtr$step, na.rm = T); sd(wtdData_wtr$step, na.rm = T)
  mean(cougData_smr_OK$step, na.rm = T); sd(cougData_smr_OK$step, na.rm = T)
  mean(cougData_wtr_OK$step, na.rm = T); sd(cougData_wtr_OK$step, na.rm = T)
  mean(cougData_smr_NE$step, na.rm = T); sd(cougData_smr_NE$step, na.rm = T)
  mean(cougData_wtr_NE$step, na.rm = T); sd(cougData_wtr_NE$step, na.rm = T)
  mean(wolfData_smr_OK$step, na.rm = T); sd(wolfData_smr_OK$step, na.rm = T)
  mean(wolfData_wtr_OK$step, na.rm = T); sd(wolfData_wtr_OK$step, na.rm = T)
  mean(wolfData_smr_NE$step, na.rm = T); sd(wolfData_smr_NE$step, na.rm = T)
  mean(wolfData_wtr_NE$step, na.rm = T); sd(wolfData_wtr_NE$step, na.rm = T)
  mean(bobData_smr_OK$step, na.rm = T); sd(bobData_smr_OK$step, na.rm = T)
  mean(bobData_wtr_OK$step, na.rm = T); sd(bobData_wtr_OK$step, na.rm = T)
  mean(bobData_smr_NE$step, na.rm = T); sd(bobData_smr_NE$step, na.rm = T)
  mean(bobData_wtr_NE$step, na.rm = T); sd(bobData_wtr_NE$step, na.rm = T)
  mean(coyData_smr_OK$step, na.rm = T); sd(coyData_smr_OK$step, na.rm = T)
  mean(coyData_wtr_OK$step, na.rm = T); sd(coyData_wtr_OK$step, na.rm = T)
  mean(coyData_smr_NE$step, na.rm = T); sd(coyData_smr_NE$step, na.rm = T)
  mean(coyData_wtr_NE$step, na.rm = T); sd(coyData_wtr_NE$step, na.rm = T)
  #' plot(mdData_smr)  #250, 500, 250, 500 -- old values
  #' plot(elkData_smr)  #500, 1000, 500, 1000 -- old values
  #' plot(wtdData_smr)  #100, 500, 100, 500 -- old values
  #' plot(cougData_smr_OK)  #500, 1500, 500, 1500 -- old values
  #' plot(cougData_smr_NE)  #500, 1500, 500, 1500 -- old values
  #' plot(wolfData_smr_OK)  #500, 3000, 500, 3000 -- old values
  #' plot(wolfData_smr_NE)  #500, 3000, 500, 3000 -- old values
  #' plot(bobData_smr_OK)  #500, 1000, 500, 1000 -- old values
  #' plot(bobData_smr_NE)  #500, 1000, 500, 1000 -- old values
  #' plot(coyData_smr_OK)  #500, 2000, 500, 2000 -- old values
  #' plot(coyData_smr_NE)  #500, 2000, 500, 2000 -- old values
  
  #'  Visualize data to identify potential temporal autocorrelation
  #'  lag.max is measured in hours
  acf(mdData_smr$step[!is.na(mdData_smr$step)],lag.max=100)
  acf(mdData_wtr$step[!is.na(mdData_wtr$step)],lag.max=100)
  acf(elkData_smr$step[!is.na(elkData_smr$step)],lag.max=100)
  acf(elkData_wtr$step[!is.na(elkData_wtr$step)],lag.max=100)
  acf(wtdData_smr$step[!is.na(wtdData_smr$step)],lag.max=100)
  acf(wtdData_wtr$step[!is.na(wtdData_wtr$step)],lag.max=100)
  acf(cougData_smr_OK$step[!is.na(cougData_smr_OK$step)],lag.max=100)
  acf(cougData_wtr_OK$step[!is.na(cougData_wtr_OK$step)],lag.max=100)
  acf(cougData_smr_NE$step[!is.na(cougData_smr_NE$step)],lag.max=100)
  acf(cougData_wtr_NE$step[!is.na(cougData_wtr_NE$step)],lag.max=100)
  acf(wolfData_smr_OK$step[!is.na(wolfData_smr_OK$step)],lag.max=100)
  acf(wolfData_wtr_OK$step[!is.na(wolfData_wtr_OK$step)],lag.max=100)
  acf(wolfData_smr_NE$step[!is.na(wolfData_smr_NE$step)],lag.max=100)
  acf(wolfData_wtr_NE$step[!is.na(wolfData_wtr_NE$step)],lag.max=100)
  acf(bobData_smr_OK$step[!is.na(bobData_smr_OK$step)],lag.max=100)
  acf(bobData_wtr_OK$step[!is.na(bobData_wtr_OK$step)],lag.max=100)
  acf(bobData_smr_NE$step[!is.na(bobData_smr_NE$step)],lag.max=100)
  acf(bobData_wtr_NE$step[!is.na(bobData_wtr_NE$step)],lag.max=100)
  acf(coyData_smr_OK$step[!is.na(coyData_smr_OK$step)],lag.max=100)
  acf(coyData_wtr_OK$step[!is.na(coyData_wtr_OK$step)],lag.max=100)
  acf(coyData_smr_NE$step[!is.na(coyData_smr_NE$step)],lag.max=100)
  acf(coyData_wtr_NE$step[!is.na(coyData_wtr_NE$step)],lag.max=100)
 
  #'  What's up with the ACF? Plot step lengths against hour to look for patterns
  mdData_smr <- hmm_data[[1]] 
  mdData_wtr <- hmm_data[[2]] 
  elkData_smr <- hmm_data[[3]]
  elkData_wtr <- hmm_data[[4]]
  wtdData_smr <- hmm_data[[5]] 
  wtdData_wtr <- hmm_data[[6]] 
  cougData_smr_OK <- hmm_data[[7]] 
  cougData_wtr_OK <- hmm_data[[8]] 
  cougData_smr_NE <- hmm_data[[9]]
  cougData_wtr_NE <- hmm_data[[10]]
  wolfData_smr_OK <- hmm_data[[11]]
  wolfData_wtr_OK <- hmm_data[[12]]
  wolfData_smr_NE <- hmm_data[[13]] 
  wolfData_wtr_NE <- hmm_data[[14]] 
  
  
  write.csv(mdData_smr, "./Outputs/Telemetry_crwOut/mdData_smr_crwOut.csv")
  write.csv(mdData_wtr, "./Outputs/Telemetry_crwOut/mdData_wtr_crwOut.csv")
  write.csv(wtdData_smr, "./Outputs/Telemetry_crwOut/wtdData_smr_crwOut.csv")
  write.csv(wtdData_wtr, "./Outputs/Telemetry_crwOut/wtdData_wtr_crwOut.csv")
  write.csv(cougData_smr_OK, "./Outputs/Telemetry_crwOut/cougData_smr_OK_crwOut.csv")
  write.csv(cougData_wtr_OK, "./Outputs/Telemetry_crwOut/cougData_wtr_OK_crwOut.csv")
  write.csv(wolfData_smr_NE, "./Outputs/Telemetry_crwOut/wolfData_smr_NE_crwOut.csv")
  write.csv(wolfData_wtr_NE, "./Outputs/Telemetry_crwOut/wolfData_wtr_NE_crwOut.csv")
  
  pdf(file = "./Outputs/Telemetry_crwOut/Step_Length_by_Hour.pdf")
  #' Boxplot of step length by hour, EXCLUDES outliers
  boxplot(step~hour_fix,data=mdData_smr, main="Mule Deer Summer Daily Step Length, no outliers",
          xlab="GPS Collar Fix Time", ylab="Step Length (m)", outline=FALSE)
  boxplot(step~hour_fix,data=mdData_wtr, main="Mule Deer Winter Daily Step Length, no outliers",
          xlab="GPS Collar Fix Time", ylab="Step Length (m)", outline=FALSE)
  boxplot(step~hour_fix,data=wtdData_smr, main="White-tailed Deer Summer Daily Step Length, no outliers",
          xlab="GPS Collar Fix Time", ylab="Step Length (m)", outline=FALSE)
  boxplot(step~hour_fix,data=wtdData_wtr, main="White-tailed Deer Winter Daily Step Length, no outliers",
          xlab="GPS Collar Fix Time", ylab="Step Length (m)", outline=FALSE)
  boxplot(step~hour_fix,data=cougData_smr_OK, main="Cougar Okanogan Summer Daily Step Length, no outliers",
          xlab="GPS Collar Fix Time", ylab="Step Length (m)", outline=FALSE)
  boxplot(step~hour_fix,data=cougData_wtr_OK, main="Cougar Okanogan Winter Daily Step Length, no outliers",
          xlab="GPS Collar Fix Time", ylab="Step Length (m)", outline=FALSE)
  boxplot(step~hour_fix,data=wolfData_smr_NE, main="Wolf Northeast Summer Daily Step Length, no outliers",
          xlab="GPS Collar Fix Time", ylab="Step Length (m)", outline=FALSE)
  boxplot(step~hour_fix,data=wolfData_wtr_NE, main="Wolf Northeast Winter Daily Step Length, no outliers",
          xlab="GPS Collar Fix Time", ylab="Step Length (m)", outline=FALSE)
  #' Boxplot of step length by hour, INCLUDES outliers
  boxplot(step~hour_fix,data=mdData_smr, main="Mule Deer Summer Daily Step Length, w/ outliers",
          xlab="GPS Collar Fix Time", ylab="Step Length (m)")
  boxplot(step~hour_fix,data=mdData_wtr, main="Mule Deer Winter Daily Step Length, w/ outliers",
          xlab="GPS Collar Fix Time", ylab="Step Length (m)")
  boxplot(step~hour_fix,data=wtdData_smr, main="White-tailed Deer Summer Daily Step Length, w/ outliers",
          xlab="GPS Collar Fix Time", ylab="Step Length (m)")
  boxplot(step~hour_fix,data=wtdData_wtr, main="White-tailed Deer Winter Daily Step Length, w/ outliers",
          xlab="GPS Collar Fix Time", ylab="Step Length (m)")
  boxplot(step~hour_fix,data=cougData_smr_OK, main="Cougar Okanogan Summer Daily Step Length, w/ outliers",
          xlab="GPS Collar Fix Time", ylab="Step Length (m)")
  boxplot(step~hour_fix,data=cougData_wtr_OK, main="Cougar Okanogan Winter Daily Step Length, w/ outliers",
          xlab="GPS Collar Fix Time", ylab="Step Length (m)")
  boxplot(step~hour_fix,data=wolfData_smr_NE, main="Wolf Northeast Summer Daily Step Length, w/ outliers",
          xlab="GPS Collar Fix Time", ylab="Step Length (m)")
  boxplot(step~hour_fix,data=wolfData_wtr_NE, main="Wolf Northeast Winter Daily Step Length, w/ outliers",
          xlab="GPS Collar Fix Time", ylab="Step Length (m)")
  dev.off()
  
  
  ####  Initial model set up  ####
  #'  ============================
  #'  Define initial parameters associated with each distribution & each state
  #'  Species-specific parameters based on viewing plotted data and mean step lengths
  #'  Providing value close to mean step length as "exploratory" mean & SD
  Par0_m1_md <- list(step = c(100, 250, 100, 250, 0.01, 0.005), angle = c(0.1, 0.5))  
  Par0_m1_elk <- list(step = c(100, 450, 100, 450, 0.01, 0.005), angle = c(0.1, 0.5))  
  Par0_m1_wtd <- list(step = c(100, 260, 100, 260, 0.01, 0.005), angle = c(0.1, 0.5))  
  Par0_m1_coug <- list(step = c(100, 650, 100, 650, 0.01, 0.005), angle = c(0.1, 0.5))  
  Par0_m1_wolf <- list(step = c(100, 1600, 100, 1600), angle = c(0.1, 0.5))  
  Par0_m1_bob <- list(step = c(100, 470, 100, 580), angle = c(0.1, 0.5))  
  Par0_m1_bob_zmass <- list(step = c(100, 470, 100, 580, 0.01, 0.005), angle = c(0.1, 0.5))
  Par0_m1_coy <- list(step = c(100, 850, 100, 850), angle = c(0.1, 0.5))  
  #'  Step arguments: report 2 means then the 2 SD for the two different states
  #'  Gamma distribution: mean & standard deviation of step lengths for each state
  #'  Michelot & Langrock 2019 recommend using same value for mean and SD per state
  #'  Wrapped Cauchy distribution: concentration of turning angles for each state
  #'  Include zero-mass parameters when there are 0s in the data w/gamma, Weibull, 
  #'  etc. distributions, e.g., zeromass0 <- c(0.1,0.05) # step zero-mass
  #'  Applies to mule deer, elk, white-tailed deer, and cougars
  
  #'  Label states
  stateNames <- c("encamped", "exploratory")
  
  ####  Models describing State-Dependent Distributions  ####
  #' Distributions for observation processes
  #' Step length: gamma or Weibull; Turning angle: von Mises or wrapped Cauchy
  #' State dwell time: geometric distribution
  #' Weibull = "weibull"; von Mises = "vm"
  dists_wc <- list(step = "gamma", angle = "wrpcauchy")  
  dists_vm <- list(step = "gamma", angle = "vm")
  
  #'  Define formula(s) to be applied to state-dependent distributions
  #'  Covariates that help describe movement patterns of a given state
  #'  Add zeromass = formula for species that need zeromass parameters above
  #'  Note that factor-level covariates must be individually specified 
  #'  (e.g., 'sexF', 'sexM') when using pseudo-design matrix (harbourSealExample)
  DM_formula_null <- ~1

  #'  Create pseudo-design matrices for state-dependent distributions
  DM_null <- list(step = list(mean = ~1, sd = ~1), angle = list(concentration = ~1))
  DM_null_ZeroMass <- list(step = list(mean = ~1, sd = ~1, zeromass = ~1), angle = list(concentration = ~1)) # includes zeromass parameters
  DM_time <- list(step = list(mean = ~daytime + cosinor(hour_fix, period = 12), sd = ~daytime + cosinor(hour_fix, period = 12)), angle = list(concentration = ~1))
  DM_Zerotime <- list(step = list(mean = ~daytime + cosinor(hour_fix, period = 12), sd = ~daytime + cosinor(hour_fix, period = 12), zeromass = ~1), angle = list(concentration = ~1)) # includes zeromass parameters
  
  # period = 4   # period = 6   # period = 12   # period = 24
  
  ####  Models describing Transition Probabilities  ####
  #'  Define formula(s) to be applied to transition probabilities
  #'  Covariates affecting probability of transitioning from one state to another
  #'  and associated with behavioral states
  #'  Univariate models
  trans_formula_null <- ~1
  trans_formula_time <- ~cosinor(hour, period = 24)
  trans_formula_TRI <- ~TRI
  trans_formula_Open <- ~PercOpen
  trans_formula_Snow <- ~SnowCover
  trans_formula_Road <- ~Dist2Road 
  trans_formula_coug <- ~COUG_RSF
  trans_formula_wolf <- ~WOLF_RSF
  trans_formula_bob <- ~BOB_RSF
  trans_formula_coy <- ~COY_RSF
  trans_formula_md <- ~MD_RSF
  trans_formula_elk <- ~ELK_RSF
  trans_formula_wtd <- ~WTD_RSF
  
  #'  Global models
  #'  For prey species
  trans_formula_smr_all <- ~TRI + PercOpen + Dist2Road + COUG_RSF + WOLF_RSF #+ BOB_RSF + COY_RSF 
  trans_formula_wtr_all <- ~TRI + PercOpen + Dist2Road + SnowCover + COUG_RSF + WOLF_RSF #+ BOB_RSF + COY_RSF 
  trans_formula_smr_all_noCoy <- ~TRI + PercOpen + Dist2Road + COUG_RSF + WOLF_RSF #+ BOB_RSF 
  trans_formula_smr_all_noTRI <- ~PercOpen + Dist2Road + COUG_RSF + WOLF_RSF #+ BOB_RSF + COY_RSF
  trans_formula_smr_all_noBob <- ~TRI + PercOpen + Dist2Road + COUG_RSF + WOLF_RSF #+ COY_RSF 
  trans_formula_wtr_all_noBob <- ~TRI + PercOpen + Dist2Road + SnowCover + COUG_RSF + WOLF_RSF #+ COY_RSF
  #'  For predator species
  trans_formula_smr_OK <- ~TRI + PercOpen + Dist2Road + MD_RSF 
  trans_formula_wtr_OK <- ~TRI + PercOpen + Dist2Road + SnowCover + MD_RSF 
  trans_formula_wtr_OK_noMD <- ~TRI + PercOpen + Dist2Road + SnowCover
  trans_formula_wtr_OK_noTRI <- ~PercOpen + Dist2Road + SnowCover + MD_RSF
  trans_formula_smr_NE <- ~TRI + PercOpen + Dist2Road + ELK_RSF + WTD_RSF 
  trans_formula_wtr_NE <- ~TRI + PercOpen + Dist2Road + SnowCover + ELK_RSF + WTD_RSF 
  
  
  
  ####  It's H[a]MM[er] Time!  ####
  #'  =============================
  #'  Keep in mind I can fit covariates on the state transition probabilities, 
  #'  meaning the variables that influence whether an animal will transition from
  #'  one state to the other, or on the state-dependent observation distributions,
  #'  meaning variables that influence step length and/or turning angle for each
  #'  of the states. 

  #'  Use retryFits argument to specify the number of attempts to minimize the 
  #'  negative log-likelihood based on random perturbations of the parameter 
  #'  estimates at the current minimum- helps ensure convergence
  
  #'  Function to run data through null and global HMM for each species
  HMM_fit <- function(Data, dists, Par0_m1, dm, tformula, fits) { 
    
    #' Fit basic model with no covariates
    m1 <- fitHMM(data = Data, nbStates = 2, dist = dists, Par0 = Par0_m1,
                 estAngleMean = list(angle = FALSE), stateNames = stateNames,
                 retryFits = fits)
    
    #'  Get new initial parameter values for global model based on nested m1 model
    Par0_m2 <- getPar0(model = m1, DM = dm, formula = tformula)   
    
    #'  Fit model with sex covariate on transition probability
    m2 <- fitHMM(data = Data, nbStates = 2, dist = dists, Par0 = Par0_m2$Par,
                 stateNames = stateNames, DM = dm, beta0 = Par0_m2$beta, formula = tformula) 
    
    #'  What proportion of the locations fall within each state?
    states <- viterbi(m2)
    print(table(states)/nrow(Data))
    
    #'  Model selection with AIC
    print(AIC(m1,m2))
    
    #'  Model summary and covariate effects
    print(m2)
    
    global_est <- CIbeta(m2, alpha = 0.95)
    print(global_est[[3]])
    
    return(m2)
    
  }
  ####  MULE DEER HMMS  ####     
  #'  Summer
  #' #'  Univariate models   
  #' md_HMM_smr_null <- HMM_fit(mdData_smr, dists_vm, Par0_m1_md, DM_Zerotime, trans_formula_null, fits = 1)
  #' md_HMM_smr_time <- HMM_fit(mdData_smr, dists_vm, Par0_m1_md, DM_Zerotime, trans_formula_time, fits = 1)
  #' md_HMM_smr_TRI <- HMM_fit(mdData_smr, dists_vm, Par0_m1_md, DM_Zerotime, trans_formula_TRI, fits = 1)
  #' md_HMM_smr_Open <- HMM_fit(mdData_smr, dists_vm, Par0_m1_md, DM_Zerotime, trans_formula_Open, fits = 1)
  #' md_HMM_smr_Road <- HMM_fit(mdData_smr, dists_vm, Par0_m1_md, DM_Zerotime, trans_formula_Road, fits = 1)
  #' md_HMM_smr_COUG <- HMM_fit(mdData_smr, dists_vm, Par0_m1_md, DM_Zerotime, trans_formula_coug, fits = 1)
  #' md_HMM_smr_WOLF <- HMM_fit(mdData_smr, dists_vm, Par0_m1_md, DM_Zerotime, trans_formula_wolf, fits = 1)
  #' md_HMM_smr_BOB <- HMM_fit(mdData_smr, dists_vm, Par0_m1_md, DM_Zerotime, trans_formula_bob, fits = 1)
  #' md_HMM_smr_COY <- HMM_fit(mdData_smr, dists_vm, Par0_m1_md, DM_Zerotime, trans_formula_coy, fits = 1)
  #' #'  Global model 
  #' #' md_HMM_smr_wc <- HMM_fit(mdData_smr, dists_wc, Par0_m1_md, DM_Zerotime, trans_formula_smr_all_noCoy, fits = 1)
  #' #'  COY & TRI highly correlated (-0.75) so running models with TRI & MD separately
  #' #'  Will use AIC to choose final model  
  #' md_HMM_smr_noCoy <- HMM_fit(mdData_smr, dists_vm, Par0_m1_md, DM_Zerotime, trans_formula_smr_all_noCoy, fits = 1)
  #' md_HMM_smr_noTRI <- HMM_fit(mdData_smr, dists_vm, Par0_m1_md, DM_Zerotime, trans_formula_smr_all_noTRI, fits = 1)
  #' AIC(md_HMM_smr_noTRI, md_HMM_smr_noCoy)
  #' #'  Final model based on AIC above (noCoy has lower AIC by 83)
  #' md_HMM_smr <- HMM_fit(mdData_smr, dists_vm, Par0_m1_md, DM_Zerotime, trans_formula_smr_all_noCoy, fits = 1)
  md_HMM_smr <- HMM_fit(mdData_smr, dists_vm, Par0_m1_md, DM_Zerotime, trans_formula_smr_all, fits = 1)
  #'  QQplot of residuals
  plotPR(md_HMM_smr, lag.max = 100, ncores = 4)
  #'  Does temporal autocorrelation look any better? NOPE looks real bad on step length
  pr_md_HMM_smr <- pseudoRes(md_HMM_smr)
  acf(pr_md_HMM_smr$stepRes[is.finite(pr_md_HMM_smr$stepRes)], lag.max = 100)
  plot(md_HMM_smr, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  
  #'  Winter
  #' #'  Univariate models
  #' md_HMM_wtr_time <- HMM_fit(mdData_wtr, dists_vm, Par0_m1_md, DM_Zerotime, trans_formula_time, fits = 1)
  #' md_HMM_wtr_TRI <- HMM_fit(mdData_wtr, dists_vm, Par0_m1_md, DM_Zerotime, trans_formula_TRI, fits = 1)
  #' md_HMM_wtr_Open <- HMM_fit(mdData_wtr, dists_vm, Par0_m1_md, DM_Zerotime, trans_formula_Open, fits = 1)
  #' md_HMM_wtr_Road <- HMM_fit(mdData_wtr, dists_vm, Par0_m1_md, DM_Zerotime, trans_formula_Road, fits = 1)
  #' md_HMM_wtr_Snow <- HMM_fit(mdData_wtr, dists_vm, Par0_m1_md, DM_Zerotime, trans_formula_Snow, fits = 1)
  #' md_HMM_wtr_COUG <- HMM_fit(mdData_wtr, dists_vm, Par0_m1_md, DM_Zerotime, trans_formula_coug, fits = 1)
  #' md_HMM_wtr_WOLF <- HMM_fit(mdData_wtr, dists_vm, Par0_m1_md, DM_Zerotime, trans_formula_wolf, fits = 1)
  #' md_HMM_wtr_BOB <- HMM_fit(mdData_wtr, dists_vm, Par0_m1_md, DM_Zerotime, trans_formula_bob, fits = 1)
  #' md_HMM_wtr_COY <- HMM_fit(mdData_wtr, dists_vm, Par0_m1_md, DM_Zerotime, trans_formula_coy, fits = 1)
  #' #'  Global model
  #' md_HMM_wtr_wc <- HMM_fit(mdData_wtr, dists_wc, Par0_m1_md, DM_Zerotime, trans_formula_wtr_all, fits = 1)
  md_HMM_wtr <- HMM_fit(mdData_wtr, dists_vm, Par0_m1_md, DM_Zerotime, trans_formula_wtr_all, fits = 1)
  #'  QQplot of residuals
  plotPR(md_HMM_wtr, lag.max = 100, ncores = 4)
  #'  Does temporal autocorrelation look any better?
  pr_md_HMM_wtr <- pseudoRes(md_HMM_wtr)
  acf(pr_md_HMM_wtr$stepRes[is.finite(pr_md_HMM_wtr$stepRes)], lag.max = 100)
  plot(md_HMM_wtr, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  
  
  ####  ELK HMMS  ####
  #'  Summer
  #' #'  Univariate models
  #' elk_HMM_smr_time <- HMM_fit(elkData_smr, dists_vm, Par0_m1_elk, DM_Zerotime, trans_formula_time, fits = 1)
  #' elk_HMM_smr_TRI <- HMM_fit(elkData_smr, dists_vm, Par0_m1_elk, DM_Zerotime, trans_formula_TRI, fits = 1)
  #' elk_HMM_smr_Open <- HMM_fit(elkData_smr, dists_vm, Par0_m1_elk, DM_Zerotime, trans_formula_Open, fits = 1)
  #' elk_HMM_smr_Road <- HMM_fit(elkData_smr, dists_vm, Par0_m1_elk, DM_Zerotime, trans_formula_Road, fits = 1)
  #' elk_HMM_smr_COUG <- HMM_fit(elkData_smr, dists_vm, Par0_m1_elk, DM_Zerotime, trans_formula_coug, fits = 1)
  #' elk_HMM_smr_WOLF <- HMM_fit(elkData_smr, dists_vm, Par0_m1_elk, DM_Zerotime, trans_formula_wolf, fits = 1)
  #' elk_HMM_smr_BOB <- HMM_fit(elkData_smr, dists_vm, Par0_m1_elk, DM_Zerotime, trans_formula_bob, fits = 1)
  #' elk_HMM_smr_COY <- HMM_fit(elkData_smr, dists_vm, Par0_m1_elk, DM_Zerotime, trans_formula_coy, fits = 1)
  #' #'  Global model
  #' elk_HMM_smr_wc <- HMM_fit(elkData_smr, dists_wc, Par0_m1_elk, DM_Zerotime, trans_formula_smr_all_noBob, fits = 1)
  elk_HMM_smr <- HMM_fit(elkData_smr, dists_vm, Par0_m1_elk, DM_Zerotime, trans_formula_smr_all, fits = 1)
  #'  QQplot of residuals
  plotPR(elk_HMM_smr, lag.max = 100, ncores = 4)
  #'  Does temporal autocorrelation look any better?
  pr_elk_HMM_smr <- pseudoRes(elk_HMM_smr)
  acf(pr_elk_HMM_smr$stepRes[!is.na(pr_elk_HMM_smr$stepRes)],lag.max = 100)
  #'  Review parameter estimates
  plot(elk_HMM_smr, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE) # decoded locations don't seem to have any encamped observations?
  
  
  #'  Winter
  #' #'  Univariate models
  #' elk_HMM_wtr_time <- HMM_fit(elkData_wtr, dists_vm, Par0_m1_elk, DM_Zerotime, trans_formula_time, fits = 1)
  #' elk_HMM_wtr_TRI <- HMM_fit(elkData_wtr, dists_vm, Par0_m1_elk, DM_Zerotime, trans_formula_TRI, fits = 1)
  #' elk_HMM_wtr_Open <- HMM_fit(elkData_wtr, dists_vm, Par0_m1_elk, DM_Zerotime, trans_formula_Open, fits = 1)
  #' elk_HMM_wtr_Road <- HMM_fit(elkData_wtr, dists_vm, Par0_m1_elk, DM_Zerotime, trans_formula_Road, fits = 1)
  #' elk_HMM_wtr_Snow <- HMM_fit(elkData_wtr, dists_vm, Par0_m1_elk, DM_Zerotime, trans_formula_Snow, fits = 1)
  #' elk_HMM_wtr_COUG <- HMM_fit(elkData_wtr, dists_vm, Par0_m1_elk, DM_Zerotime, trans_formula_coug, fits = 1)
  #' elk_HMM_wtr_WOLF <- HMM_fit(elkData_wtr, dists_vm, Par0_m1_elk, DM_Zerotime, trans_formula_wolf, fits = 1)
  #' elk_HMM_wtr_BOB <- HMM_fit(elkData_wtr, dists_vm, Par0_m1_elk, DM_Zerotime, trans_formula_bob, fits = 1)
  #' elk_HMM_wtr_COY <- HMM_fit(elkData_wtr, dists_vm, Par0_m1_elk, DM_Zerotime, trans_formula_coy, fits = 1)
  #' #'  Global model
  #' elk_HMM_wtr_wc <- HMM_fit(elkData_wtr, dists_wc, Par0_m1_elk, DM_Zerotime, trans_formula_wtr_all_noBob, fits = 1)
  elk_HMM_wtr <- HMM_fit(elkData_wtr, dists_vm, Par0_m1_elk, DM_Zerotime, trans_formula_wtr_all, fits = 1)
  #'  QQplot of residuals
  plotPR(elk_HMM_wtr, lag.max = 100, ncores = 4)
  #'  Does temporal autocorrelation look any better?
  pr_elk_HMM_wtr <- pseudoRes(elk_HMM_wtr)
  acf(pr_elk_HMM_wtr$stepRes[!is.na(pr_elk_HMM_wtr$stepRes)],lag.max = 100)
  #'  Review parameter estimates
  plot(elk_HMM_wtr, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  
  
  ####  WHITE-TAILED DEER HMMS  ####
  #'  Summer
  #' #'  Univariate models
  #' wtd_HMM_smr_time <- HMM_fit(wtdData_smr, dists_vm, Par0_m1_wtd, DM_Zerotime, trans_formula_time, fits = 1)
  #' wtd_HMM_smr_TRI <- HMM_fit(wtdData_smr, dists_vm, Par0_m1_wtd, DM_Zerotime, trans_formula_TRI, fits = 1)
  #' wtd_HMM_smr_Open <- HMM_fit(wtdData_smr, dists_vm, Par0_m1_wtd, DM_Zerotime, trans_formula_Open, fits = 1)
  #' wtd_HMM_smr_Road <- HMM_fit(wtdData_smr, dists_vm, Par0_m1_wtd, DM_Zerotime, trans_formula_Road, fits = 1)
  #' wtd_HMM_smr_COUG <- HMM_fit(wtdData_smr, dists_vm, Par0_m1_wtd, DM_Zerotime, trans_formula_coug, fits = 1)
  #' wtd_HMM_smr_WOLF <- HMM_fit(wtdData_smr, dists_vm, Par0_m1_wtd, DM_Zerotime, trans_formula_wolf, fits = 1)
  #' wtd_HMM_smr_BOB <- HMM_fit(wtdData_smr, dists_vm, Par0_m1_wtd, DM_Zerotime, trans_formula_bob, fits = 1)
  #' wtd_HMM_smr_COY <- HMM_fit(wtdData_smr, dists_vm, Par0_m1_wtd, DM_Zerotime, trans_formula_coy, fits = 1)
  #' #'  Global model
  #' wtd_HMM_smr_wc <- HMM_fit(wtdData_smr, dists_wc, Par0_m1_wtd, DM_Zerotime, trans_formula_smr_all, fits = 1)
  wtd_HMM_smr <- HMM_fit(wtdData_smr, dists_vm, Par0_m1_wtd, DM_Zerotime, trans_formula_smr_all, fits = 1)
  #'  QQplot of residuals
  plotPR(wtd_HMM_smr, lag.max = 100, ncores = 4)
  #'  Does temporal autocorrelation look any better?
  pr_wtd_HMM_smr <- pseudoRes(wtd_HMM_smr)
  acf(pr_wtd_HMM_smr$stepRes[!is.na(pr_wtd_HMM_smr$stepRes)],lag.max = 100)
  #'  Review parameter estimates
  plot(wtd_HMM_smr, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  
  #'  Winter
  #' #'  Univariate models
  #' wtd_HMM_wtr_time <- HMM_fit(wtdData_wtr, dists_vm, Par0_m1_wtd, DM_Zerotime, trans_formula_time, fits = 1)
  #' wtd_HMM_wtr_TRI <- HMM_fit(wtdData_wtr, dists_vm, Par0_m1_wtd, DM_Zerotime, trans_formula_TRI, fits = 1)
  #' wtd_HMM_wtr_Open <- HMM_fit(wtdData_wtr, dists_vm, Par0_m1_wtd, DM_Zerotime, trans_formula_Open, fits = 1)
  #' wtd_HMM_wtr_Road <- HMM_fit(wtdData_wtr, dists_vm, Par0_m1_wtd, DM_Zerotime, trans_formula_Road, fits = 1)
  #' wtd_HMM_wtr_Snow <- HMM_fit(wtdData_wtr, dists_vm, Par0_m1_wtd, DM_Zerotime, trans_formula_Snow, fits = 1)
  #' wtd_HMM_wtr_COUG <- HMM_fit(wtdData_wtr, dists_vm, Par0_m1_wtd, DM_Zerotime, trans_formula_coug, fits = 1)
  #' wtd_HMM_wtr_WOLF <- HMM_fit(wtdData_wtr, dists_vm, Par0_m1_wtd, DM_Zerotime, trans_formula_wolf, fits = 1)
  #' wtd_HMM_wtr_BOB <- HMM_fit(wtdData_wtr, dists_vm, Par0_m1_wtd, DM_Zerotime, trans_formula_bob, fits = 1)
  #' wtd_HMM_wtr_COY <- HMM_fit(wtdData_wtr, dists_vm, Par0_m1_wtd, DM_Zerotime, trans_formula_coy, fits = 1)
  #' #'  Global model
  #' wtd_HMM_wtr_wc <- HMM_fit(wtdData_wtr, dists_wc, Par0_m1_wtd, DM_Zerotime, trans_formula_wtr_all, fits = 1)
  wtd_HMM_wtr <- HMM_fit(wtdData_wtr, dists_vm, Par0_m1_wtd, DM_Zerotime, trans_formula_wtr_all, fits = 1)
  #'  QQplot of residuals
  plotPR(wtd_HMM_wtr, lag.max = 100, ncores = 4)
  #'  Does temporal autocorrelation look any better?
  pr_wtd_HMM_wtr <- pseudoRes(wtd_HMM_wtr)
  acf(pr_wtd_HMM_wtr$stepRes[!is.na(pr_wtd_HMM_wtr$stepRes)],lag.max = 100)
  #'  Review parameter estimates
  plot(wtd_HMM_wtr, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  
  
  ####  COUGAR HMMS  ####       
  #'  Okanogan Summer
  #' #'  Univariate models
  #' coug_HMM_smr_OK_time <- HMM_fit(cougData_smr_OK, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_time, fits = 1)
  #' coug_HMM_smr_OK_TRI <- HMM_fit(cougData_smr_OK, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_TRI, fits = 1)
  #' coug_HMM_smr_OK_Open <- HMM_fit(cougData_smr_OK, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_Open, fits = 1)
  #' coug_HMM_smr_OK_Road <- HMM_fit(cougData_smr_OK, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_Road, fits = 1)
  #' coug_HMM_smr_OK_MD <- HMM_fit(cougData_smr_OK, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_md, fits = 1)
  #' #'  Global model
  #' coug_HMM_smr_OK_wc <- HMM_fit(cougData_smr_OK, dists_wc, Par0_m1_coug, DM_Zerotime, trans_formula_smr_OK, fits = 1)
  coug_HMM_smr_OK <- HMM_fit(cougData_smr_OK, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_smr_OK, fits = 1)
  #'  QQplot of residuals
  plotPR(coug_HMM_smr_OK, lag.max = 100, ncores = 4)
  #'  Does temporal autocorrelation look any better?
  pr_coug_HMM_smr_OK <- pseudoRes(coug_HMM_smr_OK)
  acf(pr_coug_HMM_smr_OK$stepRes[!is.na(pr_coug_HMM_smr_OK$stepRes)],lag.max = 100)
  #'  Review parameter estimates
  plot(coug_HMM_smr_OK, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  
  #'  Okanogan Winter
  #' #'  Univariate models
  #' coug_HMM_wtr_OK_time <- HMM_fit(cougData_wtr_OK, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_time, fits = 1)
  #' coug_HMM_wtr_OK_TRI <- HMM_fit(cougData_wtr_OK, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_TRI, fits = 1)
  #' coug_HMM_wtr_OK_Open <- HMM_fit(cougData_wtr_OK, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_Open, fits = 1)
  #' coug_HMM_wtr_OK_Road <- HMM_fit(cougData_wtr_OK, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_Road, fits = 1)
  #' coug_HMM_wtr_OK_Snow <- HMM_fit(cougData_wtr_OK, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_Snow, fits = 1)
  #' coug_HMM_wtr_OK_MD <- HMM_fit(cougData_wtr_OK, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_md, fits = 1)
  #' #'  Global model
  #' coug_HMM_wtr_OK_wc <- HMM_fit(cougData_wtr_OK, dists_wc, Par0_m1_coug, DM_Zerotime, trans_formula_wtr_OK_noMD, fits = 1)
  #'  MD & TRI highly correlated (0.76) so running models with TRI & MD separately
  #'  Will use AIC to choose final model
  coug_HMM_wtr_OK_noMD <- HMM_fit(cougData_wtr_OK, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_wtr_OK_noMD, fits = 1)
  coug_HMM_wtr_OK_noTRI <- HMM_fit(cougData_wtr_OK, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_wtr_OK_noTRI, fits = 1)
  AIC(coug_HMM_wtr_OK_noMD, coug_HMM_wtr_OK_noTRI)
  #'  Final model based on AIC above (noMD has lower AIC by 12)
  coug_HMM_wtr_OK <- HMM_fit(cougData_wtr_OK, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_wtr_OK_noMD, fits = 1)
  #'  QQplot of residuals
  plotPR(coug_HMM_wtr_OK, lag.max = 100, ncores = 4)
  #'  Does temporal autocorrelation look any better?
  pr_coug_HMM_wtr_OK <- pseudoRes(coug_HMM_wtr_OK)
  acf(pr_coug_HMM_wtr_OK$stepRes[!is.na(pr_coug_HMM_wtr_OK$stepRes)],lag.max = 100)
  #'  Review parameter estimates
  plot(coug_HMM_wtr_OK, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  
  #'  Northeast Summer
  #' #'  Univariate models
  #' coug_HMM_smr_NE_time <- HMM_fit(cougData_smr_NE, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_time, fits = 1)
  #' coug_HMM_smr_NE_TRI <- HMM_fit(cougData_smr_NE, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_TRI, fits = 1)
  #' coug_HMM_smr_NE_Open <- HMM_fit(cougData_smr_NE, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_Open, fits = 1)
  #' coug_HMM_smr_NE_Road <- HMM_fit(cougData_smr_NE, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_Road, fits = 1)
  #' coug_HMM_smr_NE_ELK <- HMM_fit(cougData_smr_NE, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_elk, fits = 1)
  #' coug_HMM_smr_NE_WTD <- HMM_fit(cougData_smr_NE, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_wtd, fits = 1)
  #' #'  Global model
  #' coug_HMM_smr_NE_wc <- HMM_fit(cougData_smr_NE, dists_wc, Par0_m1_coug, DM_Zerotime, trans_formula_smr_NE, fits = 1)
  coug_HMM_smr_NE <- HMM_fit(cougData_smr_NE, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_smr_NE, fits = 1)
  #'  QQplot of residuals
  plotPR(coug_HMM_smr_NE, lag.max = 100, ncores = 4)
  #'  Does temporal autocorrelation look any better?
  pr_coug_HMM_smr_NE <- pseudoRes(coug_HMM_smr_NE)
  acf(pr_coug_HMM_smr_NE$stepRes[!is.na(pr_coug_HMM_smr_NE$stepRes)],lag.max = 100)
  #'  Review parameter estimates
  plot(coug_HMM_smr_NE, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  
  #'  Northeast Winter
  #' #'  Univariate models
  #' coug_HMM_wtr_time <- HMM_fit(cougData_wtr_NE, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_time, fits = 1)
  #' coug_HMM_wtr_TRI <- HMM_fit(cougData_wtr_NE, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_TRI, fits = 1)
  #' coug_HMM_wtr_Open <- HMM_fit(cougData_wtr_NE, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_Open, fits = 1)
  #' coug_HMM_wtr_Road <- HMM_fit(cougData_wtr_NE, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_Road, fits = 1)
  #' coug_HMM_wtr_Snow <- HMM_fit(cougData_wtr_NE, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_Snow, fits = 1)
  #' coug_HMM_wtr_ELK <- HMM_fit(cougData_wtr_NE, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_elk, fits = 1)
  #' coug_HMM_wtr_WTD <- HMM_fit(cougData_wtr_NE, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_wtd, fits = 1)
  #' #'  Global model
  #' coug_HMM_wtr_NE_wc <- HMM_fit(cougData_wtr_NE, dists_wc, Par0_m1_coug, DM_Zerotime, trans_formula_wtr_NE, fits = 1)
  coug_HMM_wtr_NE <- HMM_fit(cougData_wtr_NE, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_wtr_NE, fits = 1)
  #'  QQplot of residuals
  plotPR(coug_HMM_wtr_NE, lag.max = 100, ncores = 4)
  #'  Does temporal autocorrelation look any better?
  pr_coug_HMM_wtr_NE <- pseudoRes(coug_HMM_wtr_NE)
  acf(pr_coug_HMM_wtr_NE$stepRes[!is.na(pr_coug_HMM_wtr_NE$stepRes)],lag.max = 100)
  #'  Review parameter estimates
  plot(coug_HMM_wtr_NE, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  
  
  ####  WOLF HMMS  ####
  #'  Okanogan Summer
  #' #'  Univariate models
  #' wolf_HMM_smr_OK_time <- HMM_fit(wolfData_smr_OK, dists_vm, Par0_m1_wolf, DM_time, trans_formula_time, fits = 1)
  #' wolf_HMM_smr_OK_TRI <- HMM_fit(wolfData_smr_OK, dists_vm, Par0_m1_wolf, DM_time, trans_formula_TRI, fits = 1)
  #' wolf_HMM_smr_OK_Open <- HMM_fit(wolfData_smr_OK, dists_vm, Par0_m1_wolf, DM_time, trans_formula_Open, fits = 1)
  #' wolf_HMM_smr_OK_Road <- HMM_fit(wolfData_smr_OK, dists_vm, Par0_m1_wolf, DM_time, trans_formula_Road, fits = 1)
  #' wolf_HMM_smr_OK_MD <- HMM_fit(wolfData_smr_OK, dists_vm, Par0_m1_wolf, DM_time, trans_formula_md, fits = 1)
  #' #'  Global model
  #' wolf_HMM_smr_OK_wc <- HMM_fit(wolfData_smr_OK, dists_wc, Par0_m1_wolf, DM_time, trans_formula_smr_OK, fits = 1)
  wolf_HMM_smr_OK <- HMM_fit(wolfData_smr_OK, dists_vm, Par0_m1_wolf, DM_time, trans_formula_smr_OK, fits = 1)
  #'  QQplot of residuals
  plotPR(wolf_HMM_smr_OK, lag.max = 100, ncores = 4)
  #'  Does temporal autocorrelation look any better?
  pr_wolf_HMM_smr_OK <- pseudoRes(wolf_HMM_smr_OK)
  acf(pr_wolf_HMM_smr_OK$stepRes[!is.na(pr_wolf_HMM_smr_OK$stepRes)],lag.max = 100)
  #'  Review parameter estimates
  plot(wolf_HMM_smr_OK, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  
  #'  Okanogan Winter
  #' #'  Univariate models
  #' wolf_HMM_wtr_OK_time <- HMM_fit(wolfData_wtr_OK, dists_vm, Par0_m1_wolf, DM_time, trans_formula_time, fits = 1)
  #' wolf_HMM_wtr_OK_TRI <- HMM_fit(wolfData_wtr_OK, dists_vm, Par0_m1_wolf, DM_time, trans_formula_TRI, fits = 1)
  #' wolf_HMM_wtr_OK_Open <- HMM_fit(wolfData_wtr_OK, dists_vm, Par0_m1_wolf, DM_time, trans_formula_Open, fits = 1)
  #' wolf_HMM_wtr_OK_Road <- HMM_fit(wolfData_wtr_OK, dists_vm, Par0_m1_wolf, DM_time, trans_formula_Road, fits = 1)
  #' wolf_HMM_wtr_OK_Snow <- HMM_fit(wolfData_wtr_OK, dists_vm, Par0_m1_wolf, DM_time, trans_formula_Snow, fits = 1)
  #' wolf_HMM_wtr_OK_MD <- HMM_fit(wolfData_wtr_OK, dists_vm, Par0_m1_wolf, DM_time, trans_formula_md, fits = 1)
  #' #'  Global model
  #' wolf_HMM_wtr_OK_wc <- HMM_fit(wolfData_wtr_OK, dists_wc, Par0_m1_wolf, DM_time, trans_formula_wtr_OK, fits = 1)
  wolf_HMM_wtr_OK <- HMM_fit(wolfData_wtr_OK, dists_vm, Par0_m1_wolf, DM_time, trans_formula_wtr_OK, fits = 1)
  #'  QQplot of residuals
  plotPR(wolf_HMM_wtr_OK, lag.max = NULL, ncores = 4)
  #'  Does temporal autocorrelation look any better?
  pr_wolf_HMM_wtr_OK <- pseudoRes(wolf_HMM_wtr_OK)
  acf(pr_wolf_HMM_wtr_OK$stepRes[!is.na(pr_wolf_HMM_wtr_OK$stepRes)],lag.max = 100)
  #'  Review parameter estimates
  plot(wolf_HMM_wtr_OK, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  
  #'  Northeast Summer
  #' #'  Univariate models
  #' wolf_HMM_smr_NE_time <- HMM_fit(wolfData_smr_NE, dists_vm, Par0_m1_wolf, DM_time, trans_formula_time, fits = 1)
  #' wolf_HMM_smr_NE_TRI <- HMM_fit(wolfData_smr_NE, dists_vm, Par0_m1_wolf, DM_time, trans_formula_TRI, fits = 1)
  #' wolf_HMM_smr_NE_Open <- HMM_fit(wolfData_smr_NE, dists_vm, Par0_m1_wolf, DM_time, trans_formula_Open, fits = 1)
  #' wolf_HMM_smr_NE_Road <- HMM_fit(wolfData_smr_NE, dists_vm, Par0_m1_wolf, DM_time, trans_formula_Road, fits = 1)
  #' wolf_HMM_smr_NE_ELK <- HMM_fit(wolfData_smr_NE, dists_vm, Par0_m1_wolf, DM_time, trans_formula_elk, fits = 1)
  #' wolf_HMM_smr_NE_WTD <- HMM_fit(wolfData_smr_NE, dists_vm, Par0_m1_wolf, DM_time, trans_formula_wtd, fits = 1)
  #' #'  Global model
  #' wolf_HMM_smr_NE_wc <- HMM_fit(wolfData_smr_NE, dists_wc, Par0_m1_wolf, DM_time, trans_formula_smr_NE, fits = 1)
  wolf_HMM_smr_NE <- HMM_fit(wolfData_smr_NE, dists_vm, Par0_m1_wolf, DM_time, trans_formula_smr_NE, fits = 1)
  #'  QQplot of residuals
  plotPR(wolf_HMM_smr_NE, lag.max = NULL, ncores = 4)
  #'  Does temporal autocorrelation look any better?
  pr_wolf_HMM_smr_NE <- pseudoRes(wolf_HMM_smr_NE)
  acf(pr_wolf_HMM_smr_NE$stepRes[!is.na(pr_wolf_HMM_smr_NE$stepRes)],lag.max = 100)
  #'  Review parameter estimates
  plot(wolf_HMM_smr_NE, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  
  #'  Northeast Winter
  #' #'  Univariate models
  #' wolf_HMM_wtr_NE_time <- HMM_fit(wolfData_wtr_NE, dists_vm, Par0_m1_wolf, DM_time, trans_formula_time, fits = 1)
  #' wolf_HMM_wtr_NE_TRI <- HMM_fit(wolfData_wtr_NE, dists_vm, Par0_m1_wolf, DM_time, trans_formula_TRI, fits = 1)
  #' wolf_HMM_wtr_NE_Open <- HMM_fit(wolfData_wtr_NE, dists_vm, Par0_m1_wolf, DM_time, trans_formula_Open, fits = 1)
  #' wolf_HMM_wtr_NE_Road <- HMM_fit(wolfData_wtr_NE, dists_vm, Par0_m1_wolf, DM_time, trans_formula_Road, fits = 1)
  #' wolf_HMM_wtr_NE_Snow <- HMM_fit(wolfData_wtr_NE, dists_vm, Par0_m1_wolf, DM_time, trans_formula_Snow, fits = 1)
  #' wolf_HMM_wtr_NE_ELK <- HMM_fit(wolfData_wtr_NE, dists_vm, Par0_m1_wolf, DM_time, trans_formula_elk, fits = 1)
  #' wolf_HMM_wtr_NE_WTD <- HMM_fit(wolfData_wtr_NE, dists_vm, Par0_m1_wolf, DM_time, trans_formula_wtd, fits = 1)
  #' #'  Global model
  #' wolf_HMM_wtr_NE_wc <- HMM_fit(wolfData_wtr_NE, dists_wc, Par0_m1_wolf, DM_time, trans_formula_wtr_NE, fits = 1)
  wolf_HMM_wtr_NE <- HMM_fit(wolfData_wtr_NE, dists_vm, Par0_m1_wolf, DM_time, trans_formula_wtr_NE, fits = 1)
  #'  QQplot of residuals
  plotPR(wolf_HMM_wtr_NE, lag.max = NULL, ncores = 4)
  #'  Does temporal autocorrelation look any better?
  pr_wolf_HMM_wtr_NE <- pseudoRes(wolf_HMM_wtr_NE)
  acf(pr_wolf_HMM_wtr_NE$stepRes[!is.na(pr_wolf_HMM_wtr_NE$stepRes)],lag.max = 100)
  #'  Review parameter estimates
  plot(wolf_HMM_wtr_NE, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  
  
  ####  BOBCAT HMMS  ####
  #'  Okanogan Summer
  #' #'  Univariate models
  #' bob_HMM_smr_OK_time <- HMM_fit(bobData_smr_OK, dists_vm, Par0_m1_bob, DM_time, trans_formula_time, fits = 1)
  #' bob_HMM_smr_OK_TRI <- HMM_fit(bobData_smr_OK, dists_vm, Par0_m1_bob, DM_time, trans_formula_TRI, fits = 1)
  #' bob_HMM_smr_OK_Open <- HMM_fit(bobData_smr_OK, dists_vm, Par0_m1_bob, DM_time, trans_formula_Open, fits = 1)
  #' bob_HMM_smr_OK_Road <- HMM_fit(bobData_smr_OK, dists_vm, Par0_m1_bob, DM_time, trans_formula_Road, fits = 1)
  #' bob_HMM_smr_OK_MD <- HMM_fit(bobData_smr_OK, dists_vm, Par0_m1_bob, DM_time, trans_formula_md, fits = 1)
  #' #'  Global model
  #' bob_HMM_smr_OK_wc <- HMM_fit(bobData_smr_OK, dists_wc, Par0_m1_bob_zmass, DM_Zerotime, trans_formula_smr_OK, fits = 1)
  bob_HMM_smr_OK <- HMM_fit(bobData_smr_OK, dists_vm, Par0_m1_bob_zmass, DM_Zerotime, trans_formula_smr_OK, fits = 1)
  #'  QQplot of residuals
  plotPR(bob_HMM_smr_OK, lag.max = NULL, ncores = 4)
  #'  Does temporal autocorrelation look any better?
  pr_bob_HMM_smr_OK <- pseudoRes(bob_HMM_smr_OK)
  acf(pr_bob_HMM_smr_OK$stepRes[!is.na(pr_bob_HMM_smr_OK$stepRes)],lag.max = 100)
  #'  Review parameter estimates
  plot(bob_HMM_smr_OK, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  
  #' #'  Okanogan Winter
  #' #'  Univariate models
  #' bob_HMM_wtr_OK_time <- HMM_fit(bobData_wtr_OK, dists_vm, Par0_m1_bob, DM_time, trans_formula_time, fits = 3)
  #' bob_HMM_wtr_OK_TRI <- HMM_fit(bobData_wtr_OK, dists_vm, Par0_m1_bob, DM_time, trans_formula_TRI, fits = 3)
  #' bob_HMM_wtr_OK_Open <- HMM_fit(bobData_wtr_OK, dists_vm, Par0_m1_bob, DM_time, trans_formula_Open, fits = 3)
  #' bob_HMM_wtr_OK_Road <- HMM_fit(bobData_wtr_OK, dists_vm, Par0_m1_bob, DM_time, trans_formula_Road, fits = 3)
  #' bob_HMM_wtr_OK_Snow <- HMM_fit(bobData_wtr_OK, dists_vm, Par0_m1_bob, DM_time, trans_formula_Snow, fits = 3)
  #' bob_HMM_wtr_OK_MD <- HMM_fit(bobData_wtr_OK, dists_vm, Par0_m1_bob, DM_time, trans_formula_md, fits = 3)
  #' #'  Global model
  #' #'  Error in nlm(nLogLike, optPar, nbStates, newformula, p$bounds, p$parSize,  : 
  #' #'  non-finite value supplied by 'nlm'
  #' bob_HMM_wtr_OK_wc <- HMM_fit(bobData_wtr_OK, dists_wc, Par0_m1_bob, DM_time, trans_formula_wtr_OK, fits = 3)
  #' bob_HMM_wtr_OK <- HMM_fit(bobData_wtr_OK, dists_vm, Par0_m1_bob_zmass, DM_Zerotime, trans_formula_wtr_OK, fits = 3)
  #' #'  QQplot of residuals
  #' plotPR(bob_HMM_wtr_OK, lag.max = NULL, ncores = 4)
  #' #'  Does temporal autocorrelation look any better?
  #' pr_bob_HMM_wtr_OK <- pseudoRes(bob_HMM_wtr_OK)
  #' acf(pr_bob_HMM_wtr_OK$stepRes[!is.na(pr_bob_HMM_wtr_OK$stepRes)],lag.max = 100)
  #' #'  Review parameter estimates
  #' plot(bob_HMM_wtr_OK, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  
  #' #'  Northeast Summer
  #' #'  Univariate models
  #' bob_HMM_smr_NE_time <- HMM_fit(bobData_smr_NE, dists_vm, Par0_m1_bob, DM_time, trans_formula_time, fits = 1)
  #' bob_HMM_smr_NE_TRI <- HMM_fit(bobData_smr_NE, dists_vm, Par0_m1_bob, DM_time, trans_formula_TRI, fits = 1)
  #' bob_HMM_smr_NE_Open <- HMM_fit(bobData_smr_NE, dists_vm, Par0_m1_bob, DM_time, trans_formula_Open, fits = 1)
  #' bob_HMM_smr_NE_Road <- HMM_fit(bobData_smr_NE, dists_vm, Par0_m1_bob, DM_time, trans_formula_Road, fits = 1)
  #' bob_HMM_smr_NE_ELK <- HMM_fit(bobData_smr_NE, dists_vm, Par0_m1_bob, DM_time, trans_formula_elk, fits = 1)
  #' bob_HMM_smr_NE_WTD <- HMM_fit(bobData_smr_NE, dists_vm, Par0_m1_bob, DM_time, trans_formula_wtd, fits = 1)
  #' #'  Global model
  #' #'  THIS ONE CONVERGES BUT THE ESTIMATES SEEM WEIRD
  #' bob_HMM_smr_NE_wc <- HMM_fit(bobData_smr_NE, dists_wc, Par0_m1_bob, DM_time, trans_formula_smr_NE, fits = 3)
  #' #'  CURRENTLY THROWING AN ERROR
  #' #'  Error in nlm(nLogLike, optPar, nbStates, newformula, p$bounds, p$parSize,  : 
  #' #'  non-finite value supplied by 'nlm'
  #' bob_HMM_smr_NE <- HMM_fit(bobData_smr_NE, dists_vm, Par0_m1_bob, DM_time, trans_formula_smr_NE, fits = 3)
  #' #'  QQplot of residuals
  #' plotPR(bob_HMM_smr_NE, lag.max = NULL, ncores = 4)
  #' #'  Does temporal autocorrelation look any better?
  #' pr_bob_HMM_smr_NE <- pseudoRes(bob_HMM_smr_NE)
  #' acf(pr_bob_HMM_smr_NE$stepRes[!is.na(pr_bob_HMM_smr_NE$stepRes)],lag.max = 100)
  #' #'  Review parameter estimates
  #' plot(bob_HMM_smr_NE, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  
  #'  Northeast Winter
  #' #'  Univariate models
  #' bob_HMM_wtr_NE_time <- HMM_fit(bobData_wtr_NE, dists_vm, Par0_m1_bob, DM_time, trans_formula_time, fits = 1)
  #' bob_HMM_wtr_NE_TRI <- HMM_fit(bobData_wtr_NE, dists_vm, Par0_m1_bob, DM_time, trans_formula_TRI, fits = 1)
  #' bob_HMM_wtr_NE_Open <- HMM_fit(bobData_wtr_NE, dists_vm, Par0_m1_bob, DM_time, trans_formula_Open, fits = 1)
  #' bob_HMM_wtr_NE_Road <- HMM_fit(bobData_wtr_NE, dists_vm, Par0_m1_bob, DM_time, trans_formula_Road, fits = 1)
  #' bob_HMM_wtr_NE_Snow <- HMM_fit(bobData_wtr_NE, dists_vm, Par0_m1_bob, DM_time, trans_formula_Snow, fits = 1)
  #' bob_HMM_wtr_NE_ELK <- HMM_fit(bobData_wtr_NE, dists_vm, Par0_m1_bob, DM_time, trans_formula_elk, fits = 1)
  #' bob_HMM_wtr_NE_WTD <- HMM_fit(bobData_wtr_NE, dists_vm, Par0_m1_bob, DM_time, trans_formula_wtd, fits = 1)
  #' #'  Global model
  #' bob_HMM_wtr_NE_wc <- HMM_fit(bobData_wtr_NE, dists_wc, Par0_m1_bob, DM_time, trans_formula_wtr_NE, fits = 1)
  bob_HMM_wtr_NE <- HMM_fit(bobData_wtr_NE, dists_vm, Par0_m1_bob, DM_time, trans_formula_wtr_NE, fits = 1)
  #'  QQplot of residuals
  plotPR(bob_HMM_wtr_NE, lag.max = NULL, ncores = 4)
  #'  Does temporal autocorrelation look any better?
  pr_bob_HMM_wtr_NE <- pseudoRes(bob_HMM_wtr_NE)
  acf(pr_bob_HMM_wtr_NE$stepRes[!is.na(pr_bob_HMM_wtr_NE$stepRes)],lag.max = 100)
  #'  Review parameter estimates
  plot(bob_HMM_wtr_NE, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  
  
  ####  COYOTE HMMS  ####   
  #'  Okanogan Summer
  #' #'  Univariate models
  #' coy_HMM_smr_OK_time <- HMM_fit(coyData_smr_OK, dists_vm, Par0_m1_coy, DM_time, trans_formula_time, fits = 1)
  #' coy_HMM_smr_OK_TRI <- HMM_fit(coyData_smr_OK, dists_vm, Par0_m1_coy, DM_time, trans_formula_TRI, fits = 1)
  #' coy_HMM_smr_OK_Open <- HMM_fit(coyData_smr_OK, dists_vm, Par0_m1_coy, DM_time, trans_formula_Open, fits = 1)
  #' coy_HMM_smr_OK_Road <- HMM_fit(coyData_smr_OK, dists_vm, Par0_m1_coy, DM_time, trans_formula_Road, fits = 1)
  #' coy_HMM_smr_OK_MD <- HMM_fit(coyData_smr_OK, dists_vm, Par0_m1_coy, DM_time, trans_formula_md, fits = 1)
  #' #'  Global model
  #' coy_HMM_smr_OK_wc <- HMM_fit(coyData_smr_OK, dists_wc, Par0_m1_coy, DM_time, trans_formula_smr_OK, fits = 1)
  coy_HMM_smr_OK <- HMM_fit(coyData_smr_OK, dists_vm, Par0_m1_coy, DM_time, trans_formula_smr_OK, fits = 1)
  #'  QQplot of residuals
  plotPR(coy_HMM_smr_OK, lag.max = NULL, ncores = 4)
  #'  Does temporal autocorrelation look any better?
  pr_coy_HMM_smr_OK <- pseudoRes(coy_HMM_smr_OK)
  acf(pr_coy_HMM_smr_OK$stepRes[!is.na(pr_coy_HMM_smr_OK$stepRes)],lag.max = 100)
  #'  Review parameter estimates
  plot(coy_HMM_smr_OK, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  
  #'  Okanogan Winter
  #' #'  Univariate models
  #' coy_HMM_wtr_OK_time <- HMM_fit(coyData_wtr_OK, dists_vm, Par0_m1_coy, DM_time, trans_formula_time, fits = 1)
  #' coy_HMM_wtr_OK_TRI <- HMM_fit(coyData_wtr_OK, dists_vm, Par0_m1_coy, DM_time, trans_formula_TRI, fits = 1)
  #' coy_HMM_wtr_OK_Open <- HMM_fit(coyData_wtr_OK, dists_vm, Par0_m1_coy, DM_time, trans_formula_Open, fits = 1)
  #' coy_HMM_wtr_OK_Road <- HMM_fit(coyData_wtr_OK, dists_vm, Par0_m1_coy, DM_time, trans_formula_Road, fits = 1)
  #' coy_HMM_wtr_OK_Snow <- HMM_fit(coyData_wtr_OK, dists_vm, Par0_m1_coy, DM_time, trans_formula_Snow, fits = 1)
  #' coy_HMM_wtr_OK_MD <- HMM_fit(coyData_wtr_OK, dists_vm, Par0_m1_coy, DM_time, trans_formula_md, fits = 1)
  #' #'  Global model
  #' coy_HMM_wtr_OK_wc <- HMM_fit(coyData_wtr_OK, dists_wc, Par0_m1_coy, DM_time, trans_formula_wtr_OK_noMD, fits = 1)
  #'  MD & TRI highly correlated (0.71) so running models with TRI & MD separately
  #'  Will use AIC to choose final model 
  coy_HMM_wtr_OK_noMD <- HMM_fit(coyData_wtr_OK, dists_vm, Par0_m1_coy, DM_time, trans_formula_wtr_OK_noMD, fits = 1)
  coy_HMM_wtr_OK_noTRI <- HMM_fit(coyData_wtr_OK, dists_vm, Par0_m1_coy, DM_time, trans_formula_wtr_OK_noTRI, fits = 1)
  AIC(coy_HMM_wtr_OK_noTRI, coy_HMM_wtr_OK_noMD)
  #'  Final model based on AIC above (noMD has lower AIC by 4)
  coy_HMM_wtr_OK <- HMM_fit(coyData_wtr_OK, dists_vm, Par0_m1_coy, DM_time, trans_formula_wtr_OK_noMD, fits = 1)
  #'  QQplot of residuals
  plotPR(coy_HMM_wtr_OK, lag.max = NULL, ncores = 4)
  #'  Does temporal autocorrelation look any better?
  pr_coy_HMM_wtr_OK <- pseudoRes(coy_HMM_wtr_OK)
  acf(pr_coy_HMM_wtr_OK$stepRes[!is.na(pr_coy_HMM_wtr_OK$stepRes)],lag.max = 100)
  #'  Review parameter estimates
  plot(coy_HMM_wtr_OK, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  
  #'  Northeast Summer
  #' #'  Univariate models
  #' coy_HMM_smr_NE_time <- HMM_fit(coyData_smr_NE, dists_vm, Par0_m1_coy, DM_time, trans_formula_time, fits = 1)
  #' coy_HMM_smr_NE_TRI <- HMM_fit(coyData_smr_NE, dists_vm, Par0_m1_coy, DM_time, trans_formula_TRI, fits = 1)
  #' coy_HMM_smr_NE_Open <- HMM_fit(coyData_smr_NE, dists_vm, Par0_m1_coy, DM_time, trans_formula_Open, fits = 1)
  #' coy_HMM_smr_NE_Road <- HMM_fit(coyData_smr_NE, dists_vm, Par0_m1_coy, DM_time, trans_formula_Road, fits = 1)
  #' coy_HMM_smr_NE_ELK <- HMM_fit(coyData_smr_NE, dists_vm, Par0_m1_coy, DM_time, trans_formula_elk, fits = 1)
  #' coy_HMM_smr_NE_WTD <- HMM_fit(coyData_smr_NE, dists_vm, Par0_m1_coy, DM_time, trans_formula_wtd, fits = 1)
  #' #'  Global model
  #' coy_HMM_smr_NE_wc <- HMM_fit(coyData_smr_NE, dists_wc, Par0_m1_coy, DM_time, trans_formula_smr_NE, fits = 1)
  coy_HMM_smr_NE <- HMM_fit(coyData_smr_NE, dists_vm, Par0_m1_coy, DM_time, trans_formula_smr_NE, fits = 1)
  #'  QQplot of residuals
  plotPR(coy_HMM_smr_NE, lag.max = NULL, ncores = 4)
  #'  Does temporal autocorrelation look any better?
  pr_coy_HMM_smr_NE <- pseudoRes(coy_HMM_smr_NE)
  acf(pr_coy_HMM_smr_NE$stepRes[!is.na(pr_coy_HMM_smr_NE$stepRes)],lag.max = 100)
  #'  Review parameter estimates
  plot(coy_HMM_smr_NE, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  
  #'  Northeast Winter
  #' #'  Univariate models
  #' coy_HMM_wtr_NE_time <- HMM_fit(coyData_wtr_NE, dists_vm, Par0_m1_coy, DM_time, trans_formula_time, fits = 1)
  #' coy_HMM_wtr_NE_TRI <- HMM_fit(coyData_wtr_NE, dists_vm, Par0_m1_coy, DM_time, trans_formula_TRI, fits = 1)
  #' coy_HMM_wtr_NE_Open <- HMM_fit(coyData_wtr_NE, dists_vm, Par0_m1_coy, DM_time, trans_formula_Open, fits = 1)
  #' coy_HMM_wtr_NE_Road <- HMM_fit(coyData_wtr_NE, dists_vm, Par0_m1_coy, DM_time, trans_formula_Road, fits = 1)
  #' coy_HMM_wtr_NE_Snow <- HMM_fit(coyData_wtr_NE, dists_vm, Par0_m1_coy, DM_time, trans_formula_Snow, fits = 1)
  #' coy_HMM_wtr_NE_ELK <- HMM_fit(coyData_wtr_NE, dists_vm, Par0_m1_coy, DM_time, trans_formula_elk, fits = 1)
  #' coy_HMM_wtr_NE_WTD <- HMM_fit(coyData_wtr_NE, dists_vm, Par0_m1_coy, DM_time, trans_formula_wtd, fits = 1)
  #' #'  Global model
  #' coy_HMM_wtr_NE_wc <- HMM_fit(coyData_wtr_NE, dists_wc, Par0_m1_coy, DM_time, trans_formula_wtr_NE, fits = 1)
  coy_HMM_wtr_NE <- HMM_fit(coyData_wtr_NE, dists_vm, Par0_m1_coy, DM_time, trans_formula_wtr_NE, fits = 1)
  #'  QQplot of residuals
  plotPR(coy_HMM_wtr_NE, lag.max = NULL, ncores = 4)
  #'  Does temporal autocorrelation look any better?
  pr_coy_HMM_wtr_NE <- pseudoRes(coy_HMM_wtr_NE)
  acf(pr_coy_HMM_wtr_NE$stepRes[!is.na(pr_coy_HMM_wtr_NE$stepRes)],lag.max = 100)
  #'  Review parameter estimates
  plot(coy_HMM_wtr_NE, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)

  
  
  #'  Save model results
  ungulate_HMM_output <- list(md_HMM_smr, md_HMM_wtr, elk_HMM_smr, elk_HMM_wtr, wtd_HMM_smr, wtd_HMM_wtr)
  save(ungulate_HMM_output, file = paste0("./Outputs/HMM_output/ungulate_HMM_output_", Sys.Date(), ".RData"))
  
  #'  NOTE THE ORDER IS NOW CHANGED SINCE DROPPING SOME BOBCAT MODELS (No BOB smr NE [[17]])
  spp_HMM_output <- list(md_HMM_smr, md_HMM_wtr, elk_HMM_smr, elk_HMM_wtr, wtd_HMM_smr, 
                         wtd_HMM_wtr, coug_HMM_smr_OK, coug_HMM_wtr_OK, coug_HMM_smr_NE, 
                         coug_HMM_wtr_NE, wolf_HMM_smr_OK, wolf_HMM_wtr_OK, 
                         wolf_HMM_smr_NE, wolf_HMM_wtr_NE, bob_HMM_smr_OK, 
                         #bob_HMM_wtr_OK, #bob_HMM_smr_NE, 
                         bob_HMM_wtr_NE, coy_HMM_smr_OK, coy_HMM_wtr_OK, 
                         coy_HMM_smr_NE, coy_HMM_wtr_NE)
  save(spp_HMM_output, file = paste0("./Outputs/HMM_output/spp_HMM_output_", Sys.Date(), ".RData"))
  

  ####  Summarize Results  ####
  load("./Outputs/HMM_output/spp_HMM_output_2023-04-08.RData") #2022-06-15
  load("./Outputs/HMM_output/ungulate_HMM_output_2023-04-04.RData")
  
  #' #'  Swap original ungulate models with updated models (2023-04-04: no mesopredators)
  #' spp_HMM_output_new <- list(ungulate_HMM_output[[1]], ungulate_HMM_output[[2]], 
  #'                            ungulate_HMM_output[[3]], ungulate_HMM_output[[4]], 
  #'                            ungulate_HMM_output[[5]], ungulate_HMM_output[[6]], 
  #'                            spp_HMM_output[[7]], spp_HMM_output[[8]], 
  #'                            spp_HMM_output[[9]], spp_HMM_output[[10]], spp_HMM_output[[11]], 
  #'                            spp_HMM_output[[12]], spp_HMM_output[[13]], spp_HMM_output[[14]],
  #'                            spp_HMM_output[[15]], spp_HMM_output[[16]], spp_HMM_output[[17]], 
  #'                            spp_HMM_output[[18]], spp_HMM_output[[19]], spp_HMM_output[[20]])
  #'                            # spp_HMM_output[[21]], spp_HMM_output[[22]])
  #' #'  Rename to keep life simpler
  #' spp_HMM_output <- spp_HMM_output_new
  
  #'  Review model output
  print(spp_HMM_output[[1]]) # md_HMM_smr
  print(spp_HMM_output[[2]]) # md_HMM_wtr
  print(spp_HMM_output[[3]]) # elk_HMM_smr
  print(spp_HMM_output[[4]]) # elk_HMM_wtr
  print(spp_HMM_output[[5]]) # wtd_HMM_smr
  print(spp_HMM_output[[6]]) # wtr_HMM_wtr
  print(spp_HMM_output[[7]]) # coug_HMM_smr_OK
  print(spp_HMM_output[[8]]) # coug_HMM_wtr_OK
  print(spp_HMM_output[[9]]) # coug_HMM_smr_NE
  print(spp_HMM_output[[10]]) # coug_HMM_wtr_NE
  print(spp_HMM_output[[11]]) # wolf_HMM_smr_OK
  print(spp_HMM_output[[12]]) # wolf_HMM_wtr_OK
  print(spp_HMM_output[[13]]) # wolf_HMM_smr_NE
  print(spp_HMM_output[[14]]) # wolf_HMM_wtr_NE
  print(spp_HMM_output[[15]]) # bob_HMM_smr_OK
  #'  List order changes here if using 2022-02-22 data or later
  # print(spp_HMM_output[[16]]) # bob_HMM_wtr_OK
  # print(spp_HMM_output[[17]]) # bob_HMM_smr_NE
  print(spp_HMM_output[[16]]) # bob_HMM_wtr_NE
  print(spp_HMM_output[[17]]) # coy_HMM_smr_OK
  print(spp_HMM_output[[18]]) # coy_HMM_wtr_OK
  print(spp_HMM_output[[19]]) # coy_HMM_smr_NE
  print(spp_HMM_output[[20]]) # coy_HMM_wtr_NE
  

  ####  State-Dependent Distributions  ####
  #'  Function to report state-dependent distribution parameters, including zero-mass parameters
  step_turn_parms_zmass <- function(mod, spp, season, area){ 
    #'  Pull out turning angle parameters
    step_out <- as.data.frame(mod$mle[[1]])
    step_out$Species <- spp
    step_out$Season <- season
    step_out$StudyArea <- area    #####  UPDATE FOR DAYTIME VARIABLE  ####
    colnames(step_out) <- c("State1 Intercept_mu", "State1 Daylight_mu", "State1 Cos_mu", "State1 Sin_mu", 
                            "State2 Intercept_mu", "State2 Daylight_mu", "State2 Cos_mu", "State2 Sin_mu",
                            "State1 Intercept_sd", "State1 Daylight_sd", "State1 Cos_sd", "State1 Sin_sd", 
                            "State2 Intercept_sd", "State2 Daylight_sd", "State2 Cos_sd", "State2 Sin_sd",
                            "State1 Intercept_zmass", #"State1 Daylight_zmass", "State1 Cos_zmass", "State1 Sin_zmass", 
                            "State2 Intercept_zmass", #"State2 Daylight_zmass", "State2 Cos_zmass", "State2 Sin_zmass", 
                            "Species", "Season", "StudyArea")
    #'  Wrangle parameters into an interpret-able table
    step_table <- step_out %>%
      pivot_longer(!c(Species, Season, StudyArea), names_to = "Parameter", values_to = "Estimate") %>%
      separate(Parameter, c("State", "Parameter"), sep = " ") %>%
      pivot_wider(names_from = "State", values_from = "Estimate") %>%
      separate(Parameter, c("Coefficient", "Parameter"), sep = "_") %>%
      pivot_wider(names_from = "Parameter", values_from = c("State1", "State2"))
    #'  Create separate tables for state 1 & 2 parameters
    state1 <- step_table[,1:7]
    state1$State <- "Encamped"
    colnames(state1) <- c("Species", "Season", "Study Area", "Coefficient", "Mean", "SD", "Zeromass", "State")
    state2 <- step_table[,c(1:4,8:10)]
    state2$State <- "Exploratory"
    colnames(state2) <- c("Species", "Season", "Study Area", "Coefficient", "Mean", "SD", "Zeromass", "State")
    #'  Merge into one single table of step length parameters
    step_out_tbl <- rbind(state1, state2) %>%
      relocate(State, .before = "Coefficient")

    #'  Turning angles parameters
    turn_out <- as.data.frame(mod$mle[[2]])
    turn_out$Species <- spp
    turn_out$Season <- season
    turn_out$StudyArea <- area
    turn_out_tbl <- turn_out %>%
      relocate(Species, .before = "encamped") %>%
      relocate(Season, .after = "Species") %>%
      relocate(StudyArea, .after = "Season") %>%
      rownames_to_column(var = "Parameter") %>%
      relocate(Parameter, .after = "StudyArea") %>%
      mutate(Parameter = ifelse(Parameter == "mean", "Mean", "Concentration"))
    colnames(turn_out_tbl) <- c("Species", "Season", "Study Area", "Parameter", "Encamped", "Exploratory")
      
    #'  List parameter tables together
    params_out <- list(step_out_tbl, turn_out_tbl)
    return(params_out)

  }
  #'  Create parameter tables for species that included zero-mass parameters
  md_smr_params <- step_turn_parms_zmass(spp_HMM_output[[1]], spp = "Mule Deer", season = "Summer", area = "Okanogan")
  md_wtr_params <- step_turn_parms_zmass(spp_HMM_output[[2]], spp = "Mule Deer", season = "Winter", area = "Okanogan")
  elk_smr_params <- step_turn_parms_zmass(spp_HMM_output[[3]], spp = "Elk", season = "Summer", area = "Northeast")
  elk_wtr_params <- step_turn_parms_zmass(spp_HMM_output[[4]], spp = "Elk", season = "Winter", area = "Northeast")
  wtd_smr_params <- step_turn_parms_zmass(spp_HMM_output[[5]], spp = "White-tailed Deer", season = "Summer", area = "Northeast")
  wtd_wtr_params <- step_turn_parms_zmass(spp_HMM_output[[6]], spp = "White-tailed Deer", season = "Winter", area = "Northeast")
  coug_smr_params_OK <- step_turn_parms_zmass(spp_HMM_output[[7]], spp = "Cougar", season = "Summer", area = "Okanogan")
  coug_wtr_params_OK <- step_turn_parms_zmass(spp_HMM_output[[8]], spp = "Cougar", season = "Winter", area = "Okanogan")
  coug_smr_params_NE <- step_turn_parms_zmass(spp_HMM_output[[9]], spp = "Cougar", season = "Summer", area = "Northeast")
  coug_wtr_params_NE <- step_turn_parms_zmass(spp_HMM_output[[10]], spp = "Cougar", season = "Winter", area = "Northeast")
  #' #' And NE winter bobcat- note the list order is a little wonky starting here
  #' bob_wtr_params_NE <- step_turn_parms_zmass(spp_HMM_output[[16]], spp = "Bobcat", season = "Winter", area = "Northeast")
  

  #'  Function to report state-dependent distribution parameters, excluding zero-mass parameters
  step_turn_parms <- function(mod, spp, season, area){ 
    #'  Pull out turning angle parameters
    step_out <- as.data.frame(mod$mle[[1]])
    step_out$Species <- spp
    step_out$Season <- season
    step_out$StudyArea <- area
    colnames(step_out) <- c("State1 Intercept_mu", "State1 Daylight_mu", "State1 Cos_mu", "State1 Sin_mu", 
                            "State2 Intercept_mu", "State2 Daylight_mu", "State2 Cos_mu", "State2 Sin_mu",
                            "State1 Intercept_sd", "State1 Daylight_sd", "State1 Cos_sd", "State1 Sin_sd", 
                            "State2 Intercept_sd", "State2 Daylight_sd", "State2 Cos_sd", "State2 Sin_sd",
                            "Species", "Season", "StudyArea")
    #'  Wrangle parameters into an interpret-able table
    step_table <- step_out %>%
      pivot_longer(!c(Species, Season, StudyArea), names_to = "Parameter", values_to = "Estimate") %>%
      separate(Parameter, c("State", "Parameter"), sep = " ") %>%
      pivot_wider(names_from = "State", values_from = "Estimate") %>%
      separate(Parameter, c("Coefficient", "Parameter"), sep = "_") %>%
      pivot_wider(names_from = "Parameter", values_from = c("State1", "State2"))
    #'  Create separate tables for state 1 & 2 parameters
    state1 <- step_table[,1:6]
    state1$State <- "Encamped"
    colnames(state1) <- c("Species", "Season", "Study Area", "Coefficient", "Mean", "SD", "State")
    state2 <- step_table[,c(1:4,7:8)]
    state2$State <- "Exploratory"
    colnames(state2) <- c("Species", "Season", "Study Area", "Coefficient", "Mean", "SD", "State")
    #'  Merge into one single table of step length parameters
    step_out_tbl <- rbind(state1, state2) %>%
      relocate(State, .before = "Coefficient") %>%
      mutate(zmass = NA)
    colnames(step_out_tbl) <- c("Species", "Season", "Study Area", "State", "Coefficient", 
                                "Mean", "SD", "Zeromass")
    
    #'  Turning angles parameters
    turn_out <- as.data.frame(mod$mle[[2]])
    turn_out$Species <- spp
    turn_out$Season <- season
    turn_out$StudyArea <- area
    turn_out_tbl <- turn_out %>%
      relocate(Species, .before = "encamped") %>%
      relocate(Season, .after = "Species") %>%
      relocate(StudyArea, .after = "Season") %>%
      rownames_to_column(var = "Parameter") %>%
      relocate(Parameter, .after = "StudyArea") %>%
      mutate(Parameter = ifelse(Parameter == "mean", "Mean", "Concentration"))
    colnames(turn_out_tbl) <- c("Species", "Season", "Study Area", "Parameter", "Encamped", "Exploratory")
    
    #'  List parameter tables together
    params_out <- list(step_out_tbl, turn_out_tbl)
    return(params_out)
    
  }
  #'  Create parameter tables for species that don't include zero-mass parameters
  wolf_smr_params_OK <- step_turn_parms(spp_HMM_output[[11]], spp = "Wolf", season = "Summer", area = "Okanogan")
  wolf_wtr_params_OK <- step_turn_parms(spp_HMM_output[[12]], spp = "Wolf", season = "Winter", area = "Okanogan")
  wolf_smr_params_NE <- step_turn_parms(spp_HMM_output[[13]], spp = "Wolf", season = "Summer", area = "Northeast")
  wolf_wtr_params_NE <- step_turn_parms(spp_HMM_output[[14]], spp = "Wolf", season = "Winter", area = "Northeast")
  bob_smr_params_OK <- step_turn_parms(spp_HMM_output[[15]], spp = "Bobcat", season = "Summer", area = "Okanogan")
  #'  Note the list order gets a little wonky here b/c excluding OK winter & NE summer bob for now
  # bob_wtr_params_OK <- step_turn_parms(spp_HMM_output[[16]], spp = "Bobcat", season = "Winter", area = "Okanogan")
  # bob_smr_params_NE <- step_turn_parms(spp_HMM_output[[17]], spp = "Bobcat", season = "Summer", area = "Northeast")
  bob_wtr_params_NE <- step_turn_parms(spp_HMM_output[[16]], spp = "Bobcat", season = "Winter", area = "Northeast")
  coy_smr_params_OK <- step_turn_parms(spp_HMM_output[[17]], spp = "Coyote", season = "Summer", area = "Okanogan")
  coy_wtr_params_OK <- step_turn_parms(spp_HMM_output[[18]], spp = "Coyote", season = "Winter", area = "Okanogan")
  coy_smr_params_NE <- step_turn_parms(spp_HMM_output[[19]], spp = "Coyote", season = "Summer", area = "Northeast")
  coy_wtr_params_NE <- step_turn_parms(spp_HMM_output[[20]], spp = "Coyote", season = "Winter", area = "Northeast")
  
  #'  Make single giant table of all step length parameters
  all_steps <- bind_rows(md_smr_params[[1]], md_wtr_params[[1]], elk_smr_params[[1]], 
                         elk_wtr_params[[1]], wtd_smr_params[[1]], wtd_wtr_params[[1]],
                         coug_smr_params_OK[[1]], coug_wtr_params_OK[[1]], 
                         coug_smr_params_NE[[1]], coug_wtr_params_NE[[1]], 
                         wolf_smr_params_OK[[1]], wolf_wtr_params_OK[[1]], 
                         wolf_smr_params_NE[[1]], wolf_wtr_params_NE[[1]], 
                         bob_smr_params_OK[[1]],  
                         # bob_wtr_params_OK[[1]],bob_smr_params_NE[[1]], 
                         bob_wtr_params_NE[[1]], 
                         coy_smr_params_OK[[1]], coy_wtr_params_OK[[1]], 
                         coy_smr_params_NE[[1]], coy_wtr_params_NE[[1]]) %>%
    mutate(Mean = round(Mean, 2),
           SD = round(SD, 2),
           Zeromass = round(Zeromass, 2)) %>%
    arrange(Species)
  
  #'  Make single giant table of all turning angles parameters
  all_turns <- bind_rows(md_smr_params[[2]], md_wtr_params[[2]], elk_smr_params[[2]], 
                         elk_wtr_params[[2]], wtd_smr_params[[2]], wtd_wtr_params[[2]],
                         coug_smr_params_OK[[2]], coug_wtr_params_OK[[2]], 
                         coug_smr_params_NE[[2]], coug_wtr_params_NE[[2]], 
                         wolf_smr_params_OK[[2]], wolf_wtr_params_OK[[2]], 
                         wolf_smr_params_NE[[2]], wolf_wtr_params_NE[[2]], 
                         bob_smr_params_OK[[2]], 
                         # bob_wtr_params_OK[[2]], bob_smr_params_NE[[2]], 
                         bob_wtr_params_NE[[2]], 
                         coy_smr_params_OK[[2]], coy_wtr_params_OK[[2]], 
                         coy_smr_params_NE[[2]], coy_wtr_params_NE[[2]]) %>%
    mutate(Encamped = round(Encamped, 2),
           Exploratory = round(Exploratory, 2)) %>%
    arrange(Species)
  
  # write.csv(all_steps, paste0("./Outputs/HMM_output/HMM_Results_StepLength_", Sys.Date(), ".csv"))
  # write.csv(all_turns, paste0("./Outputs/HMM_output/HMM_Results_TurningAngle_", Sys.Date(), ".csv"))
  
  
  ####  Transition Probabilities  ####
  #'  Function to report transition probability coefficients in a table
  rounddig <- 2
  hmm_out <- function(mod, spp, season, area) {
    #'  Extract estimates, standard error, and 95% Confidence Intervals for effect
    #'  of each covariate on transition probabilities
    est_out <- CIbeta(mod, alpha = 0.95)
    beta1.2 <- formatC(round(est_out[[3]]$est[,1], rounddig), rounddig, format="f")
    beta2.1 <- formatC(round(est_out[[3]]$est[,2], rounddig), rounddig, format="f")
    se1.2 <- formatC(round(est_out[[3]]$se[,1], rounddig), rounddig, format="f")
    se2.1 <- formatC(round(est_out[[3]]$se[,2], rounddig), rounddig, format="f")
    lci1.2 <- formatC(round(est_out[[3]]$lower[,1], rounddig), rounddig, format="f")
    lci2.1 <- formatC(round(est_out[[3]]$lower[,2], rounddig), rounddig, format="f")
    uci1.2 <- formatC(round(est_out[[3]]$upper[,1], rounddig), rounddig, format="f")
    uci2.1 <- formatC(round(est_out[[3]]$upper[,2], rounddig), rounddig, format="f")
    #'  Merge into a data frame and organize
    out1.2 <- as.data.frame(cbind(beta1.2, se1.2, lci1.2, uci1.2)) %>%
      mutate(
        Parameter = row.names(est_out[[3]]$est),
        Species = rep(spp, nrow(.)),
        Season = rep(season, nrow(.)),
        StudyArea = rep(area, nrow(.)),
        Transition = rep("Trans.1->2", nrow(.))
      ) %>%
      relocate(Parameter, .before = beta1.2) %>%
      relocate(Species, .before = Parameter) %>%
      relocate(Season, .before = Parameter) %>%
      relocate(StudyArea, .before = Parameter) %>%
      relocate(Transition, .before = Parameter)
    colnames(out1.2) <- c("Species", "Season", "Study Area", "Transition", "Parameter", "Estimate", "SE", "Lower", "Upper")
    out2.1 <- as.data.frame(cbind(beta2.1, se2.1, lci2.1, uci2.1)) %>%
      mutate(
        Parameter = row.names(est_out[[3]]$est),
        Species = rep(spp, nrow(.)),
        Season = rep(season, nrow(.)),
        StudyArea = rep(area, nrow(.)),
        Transition = rep("Trans.2->1", nrow(.))
      ) %>%
      relocate(Parameter, .before = beta2.1) %>%
      relocate(Species, .before = Parameter) %>%
      relocate(Season, .before = Parameter) %>%
      relocate(StudyArea, .before = Parameter) %>%
      relocate(Transition, .before = Parameter)
    colnames(out2.1) <- c("Species", "Season", "Study Area", "Transition", "Parameter", "Estimate", "SE", "Lower", "Upper")
    out <- as.data.frame(rbind(out1.2, out2.1))
    return(out)
  }
  #'  Run each season and species-specific model through function
  md_smr_hmm <- hmm_out(spp_HMM_output[[1]], "Mule Deer", "Summer", "Okanogan") #md_HMM_smr
  md_wtr_hmm <- hmm_out(spp_HMM_output[[2]], "Mule Deer", "Winter", "Okanogan") #md_HMM_wtr
  elk_smr_hmm <- hmm_out(spp_HMM_output[[3]], "Elk", "Summer", "Northeast") #elk_HMM_smr
  elk_wtr_hmm <- hmm_out(spp_HMM_output[[4]], "Elk", "Winter", "Northeast") #elk_HMM_wtr
  wtd_smr_hmm <- hmm_out(spp_HMM_output[[5]], "White-tailed Deer", "Summer", "Northeast") #wtd_HMM_smr
  wtd_wtr_hmm <- hmm_out(spp_HMM_output[[6]], "White-tailed Deer", "Winter", "Northeast") #wtd_HMM_wtr
  coug_smr_hmm_OK <- hmm_out(spp_HMM_output[[7]], "Cougar", "Summer", "Okanogan") #coug_HMM_smr_OK
  coug_wtr_hmm_OK <- hmm_out(spp_HMM_output[[8]], "Cougar", "Winter", "Okanogan") #coug_HMM_wtr_OK
  coug_smr_hmm_NE <- hmm_out(spp_HMM_output[[9]], "Cougar", "Summer", "Northeast") #coug_HMM_smr_NE
  coug_wtr_hmm_NE <- hmm_out(spp_HMM_output[[10]], "Cougar", "Winter", "Northeast") #coug_HMM_wtr_NE
  wolf_smr_hmm_OK <- hmm_out(spp_HMM_output[[11]], "Wolf", "Summer", "Okanogan") #wolf_HMM_smr_OK
  wolf_wtr_hmm_OK <- hmm_out(spp_HMM_output[[12]], "Wolf", "Winter", "Okanogan") #wolf_HMM_wtr_OK
  wolf_smr_hmm_NE <- hmm_out(spp_HMM_output[[13]], "Wolf", "Summer", "Northeast") #wolf_HMM_smr_NE
  wolf_wtr_hmm_NE <- hmm_out(spp_HMM_output[[14]], "Wolf", "Winter", "Northeast") #wolf_HMM_wtr_NE
  bob_smr_hmm_OK <- hmm_out(spp_HMM_output[[15]], "Bobcat", "Summer", "Okanogan") #bob_HMM_smr_OK
  #'  Note the list order gets wonky here b/c excluding bob wtr OK & smr NE right now
  # bob_wtr_hmm_OK <- hmm_out(spp_HMM_output[[16]], "Bobcat", "Winter", "Okanogan") #bob_HMM_wtr_OK
  # bob_smr_hmm_NE <- hmm_out(spp_HMM_output[[17]], "Bobcat", "Summer", "Northeast") #bob_HMM_smr_NE
  bob_wtr_hmm_NE <- hmm_out(spp_HMM_output[[16]], "Bobcat", "Winter", "Northeast") #bob_HMM_wtr_NE
  coy_smr_hmm_OK <- hmm_out(spp_HMM_output[[17]], "Coyote", "Summer", "Okanogan") #coy_HMM_smr_OK
  coy_wtr_hmm_OK <- hmm_out(spp_HMM_output[[18]], "Coyote", "Winter", "Okanogan") #coy_HMM_wtr_OK
  coy_smr_hmm_NE <- hmm_out(spp_HMM_output[[19]], "Coyote", "Summer", "Northeast") #coy_HMM_smr_NE
  coy_wtr_hmm_NE <- hmm_out(spp_HMM_output[[20]], "Coyote", "Winter", "Northeast") #coy_HMM_wtr_NE
  
  #'  Gather prey and predator results to put into a single results table
  results_hmm_TransPr <- rbind(md_smr_hmm, md_wtr_hmm, elk_smr_hmm, elk_wtr_hmm, 
                               wtd_smr_hmm, wtd_wtr_hmm, coug_smr_hmm_OK, coug_wtr_hmm_OK, 
                               coug_smr_hmm_NE, coug_wtr_hmm_NE, wolf_smr_hmm_OK, 
                               wolf_wtr_hmm_OK, wolf_smr_hmm_NE, wolf_wtr_hmm_NE, 
                               bob_smr_hmm_OK, #bob_wtr_hmm_OK, #bob_smr_hmm_NE, 
                               bob_wtr_hmm_NE, coy_smr_hmm_OK, coy_wtr_hmm_OK, 
                               coy_smr_hmm_NE, coy_wtr_hmm_NE)
  
  results_hmm_TransPr_prey <- rbind(md_smr_hmm, md_wtr_hmm, elk_smr_hmm, elk_wtr_hmm, 
                                    wtd_smr_hmm, wtd_wtr_hmm) %>%
    unite(CI95, Lower, Upper, sep = ", ") %>%
    # dplyr::select(-SE) %>%
    mutate(
      Parameter = ifelse(Parameter == "(Intercept)", "Intercept", Parameter),
      Parameter = ifelse(Parameter == "TRI", "Terrain Ruggedness", Parameter),
      Parameter = ifelse(Parameter == "PercOpen", "Percent Open", Parameter),
      Parameter = ifelse(Parameter == "Dist2Road", "Nearest Road", Parameter),
      Parameter = ifelse(Parameter == "SnowCover1", "Snow Cover (Y)", Parameter),
      Parameter = ifelse(Parameter == "COUG_RSF", "Pr(Cougar)", Parameter),
      Parameter = ifelse(Parameter == "WOLF_RSF", "Pr(Wolf)", Parameter),
      Parameter = ifelse(Parameter == "BOB_RSF", "Pr(Bobcat)", Parameter),
      Parameter = ifelse(Parameter == "COY_RSF", "Pr(Coyote)", Parameter)
    ) #%>%
    # group_by(Species) %>%
    # arrange(match(`Study Area`, c("Okanogan", "Northeast")), .by_group = TRUE) %>%
    # ungroup()
  colnames(results_hmm_TransPr_prey) <- c("Species", "Season", "Study Area", 
                                          "Transition", "Parameter", "Estimate", 
                                          "SE", "CI95")
  results_hmm_TransPr_pred <- rbind(coug_smr_hmm_OK, coug_wtr_hmm_OK, coug_smr_hmm_NE, 
                                    coug_wtr_hmm_NE, wolf_smr_hmm_OK, wolf_wtr_hmm_OK, 
                                    wolf_smr_hmm_NE, wolf_wtr_hmm_NE, bob_smr_hmm_OK, 
                                    #bob_wtr_hmm_OK, #bob_smr_hmm_NE, 
                                    bob_wtr_hmm_NE, coy_smr_hmm_OK, coy_wtr_hmm_OK, 
                                    coy_smr_hmm_NE, coy_wtr_hmm_NE) %>%
    unite(CI95, Lower, Upper, sep = ", ") %>%
    # dplyr::select(-SE) %>%
    mutate(
      Parameter = ifelse(Parameter == "(Intercept)", "Intercept", Parameter),
      Parameter = ifelse(Parameter == "TRI", "Terrain Ruggedness", Parameter),
      Parameter = ifelse(Parameter == "PercOpen", "Percent Open", Parameter),
      Parameter = ifelse(Parameter == "Dist2Road", "Nearest Road", Parameter),
      Parameter = ifelse(Parameter == "SnowCover1", "Snow Cover (Y)", Parameter),
      Parameter = ifelse(Parameter == "MD_RSF", "Pr(Mule Deer)", Parameter),
      Parameter = ifelse(Parameter == "ELK_RSF", "Pr(Elk)", Parameter),
      Parameter = ifelse(Parameter == "WTD_RSF", "Pr(White-tailed Deer)", Parameter)
    ) %>%
    filter(!Species == "Bobcat")
    # group_by(Species) %>%
    # arrange(Parameter, Season, `Study Area`) %>% 
    # arrange(match(Transition, c("Trans.1->2", "Trans.2->1")), .by_group = TRUE) %>%
    # arrange(match(`Study Area`, c("Okanogan", "Northeast")), .by_group = TRUE) %>%
    # ungroup()
  colnames(results_hmm_TransPr_pred) <- c("Species", "Season", "Study Area", 
                                          "Transition", "Parameter", "Estimate",
                                          "SE", "CI95")
  
  write.csv(results_hmm_TransPr_prey, paste0("./Outputs/HMM_output/HMM_Results_TransPr_prey_long", Sys.Date(), ".csv"))
  write.csv(results_hmm_TransPr_pred, paste0("./Outputs/HMM_output/HMM_Results_TransPr_pred_long", Sys.Date(), ".csv"))
  

  #'  Spread results so the coefficient effects are easier to compare between 
  #'  transition probabilities and across species
  #'  Prey HMM results
  results_hmm_wide_TransPr_prey <- results_hmm_TransPr_prey %>%  
    # dplyr::select(-z) %>%
    mutate(
      # SE = round(SE, 2),
      SE = paste0("(", SE, ")"),
    ) %>%
    #'  Bold significant variables- doesn't work if continue manipulating data frame
    # condformat(.) %>%
    # rule_text_bold(c(Estimate, SE, Pval), expression = Pval <= 0.05) %>%
    unite(Est_SE, Estimate, SE, sep = " ") %>%
    # unite(CI95, Lower, Upper, sep = ", ") %>%
    unite(Est_SE_CI, Est_SE, CI95, sep = "_") %>%
    spread(Parameter, Est_SE_CI) %>%
    separate("Intercept", c("Intercept (SE)", "Intercept 95% CI"), sep = "_") %>%
    separate("Terrain Ruggedness", c("Terrain Ruggedness (SE)", "Terrain Ruggedness 95% CI"), sep = "_") %>%
    separate("Percent Open", c("Percent Open (SE)", "Percent Open 95% CI"), sep = "_") %>%
    separate("Nearest Road", c("Nearest Road (SE)", "Nearest Road 95% CI"), sep = "_") %>%
    separate("Snow Cover (Y)", c("Snow Cover (Y) (SE)", "Snow Cover (Y) 95% CI"), sep = "_") %>%
    separate("Pr(Cougar)", c("Pr(Cougar) (SE)", "Pr(Cougar) 95% CI"), sep = "_") %>%
    separate("Pr(Wolf)", c("Pr(Wolf) (SE)", "Pr(Wolf) 95% CI"), sep = "_") %>%
    # separate("Pr(Bobcat)", c("Pr(Bobcat) (SE)", "Pr(Bobcat) 95% CI"), sep = "_") %>%
    # separate("Pr(Coyote)", c("Pr(Coyote) (SE)", "Pr(Coyote) 95% CI"), sep = "_") %>%
    arrange(match(Species, c("Mule Deer", "Elk", "White-tailed Deer")))

  
  write.csv(results_hmm_wide_TransPr_prey, paste0("./Outputs/HMM_output/HMM_Results_TransPr_prey_wide", Sys.Date(), ".csv"))
  
  #'  Predators HMM results
  results_hmm_wide_TransPr_pred <- results_hmm_TransPr_pred %>% 
    # dplyr::select(-z) %>%
    mutate(
      # SE = round(SE, 2),
      SE = paste0("(", SE, ")"),
    ) %>%
    #'  Bold significant variables- doesn't work if continue manipulating data frame
    # condformat(.) %>%
    # rule_text_bold(c(Estimate, SE, Pval), expression = Pval <= 0.05) %>%
    unite(Est_SE, Estimate, SE, sep = " ") %>%
    # unite(CI95, Lower, Upper, sep = ", ") %>%
    unite(Est_SE_CI, Est_SE, CI95, sep = "_") %>%
    spread(Parameter, Est_SE_CI) %>%
    separate("Intercept", c("Intercept (SE)", "Intercept 95% CI"), sep = "_") %>%
    separate("Terrain Ruggedness", c("Terrain Ruggedness (SE)", "Terrain Ruggedness 95% CI"), sep = "_") %>%
    separate("Percent Open", c("Percent Open (SE)", "Percent Open 95% CI"), sep = "_") %>%
    separate("Nearest Road", c("Nearest Road (SE)", "Nearest Road 95% CI"), sep = "_") %>%
    separate("Snow Cover (Y)", c("Snow Cover (Y) (SE)", "Snow Cover (Y) 95% CI"), sep = "_") %>%
    separate("Pr(Mule Deer)", c("Pr(Mule Deer) (SE)", "Pr(Mule Deer) 95% CI"), sep = "_") %>%
    separate("Pr(Elk)", c("Pr(Elk) (SE)", "Pr(Elk) 95% CI"), sep = "_") %>%
    separate("Pr(White-tailed Deer)", c("Pr(White-tailed Deer) (SE)", "Pr(White-tailed Deer) 95% CI"), sep = "_") %>%
    group_by(Species) %>%
    arrange(match(`Study Area`, c("Okanogan", "Northeast")), .by_group = TRUE) %>%
    ungroup()
  
  write.csv(results_hmm_wide_TransPr_pred, paste0("./Outputs/HMM_output/HMM_Results_TransPr_pred_wide", Sys.Date(), ".csv"))
  
 
  ####  Back-transformed Results  ####
  #'  Back-transform HMM results to the real (natural) scale of the data
  #'  Extract parameter means, SE, and 95% CI on natural scale when all covariates
  #'  are held at their mean value (i.e., 0 since covariates are scaled)
  backtrans_params <- function(mod, spp, season, area) {
    
    #'  CIreal has 4 lists: [[1]] step length params, [[2]] turning angle concentration,
    #'  [[3]] transition probabilities, and [[4]] initial state for each track.
    #'  Step length includes 4-6 lists depending on if zeromass parameter is needed
    #'  Lists 1:3 are State1 mean, sd, zeromass, 4:6 are State2 mean, sd, zeromass
    #'  Transition probability included 4 lists: [[1]] staying in State1, [[2]] 
    #'  transition from State1 to State2, [[3]] transitioning from State2 to State1,
    #'  and [[4]] staying in State2
    ci_nat <- CIreal(mod)
    
    #'  Table of step lengths (in meters) and 95% CI
    steps_state1 <- c(ci_nat[[1]]$est[[1]], ci_nat[[1]]$lower[[1]], ci_nat[[1]]$upper[[1]])
    steps_state2 <- c(ci_nat[[1]]$est[[4]], ci_nat[[1]]$lower[[4]], ci_nat[[1]]$upper[[4]])
    steps_real <- as.data.frame(rbind(steps_state1, steps_state2))
    colnames(steps_real) <- c("Mean", "Lower", "Upper")
    steps_real <- rownames_to_column(steps_real, var = "State") %>%
      mutate(Species = spp,
             Season = season,
             StudyArea = area,
             State = ifelse(State == "steps_state1", "Encamped", "Exploratory"),
             Mean = round(Mean, 2),
             Lower = round(Lower, 2),
             Upper = round(Upper, 2)) %>%
      unite("95%CI", Lower:Upper, sep = " - ") %>%
      relocate(Species, .before = "State") %>%
      relocate(StudyArea, .after = "Species") %>%
      relocate(Season, .after = "StudyArea")
    
    #'  Table of turning angles and concentrations
    turn_matrix <- matrix(c(0, 0, ci_nat[[2]]$est[[1]], ci_nat[[2]]$est[[2]]),nrow=2,ncol=2,byrow=TRUE)
    colnames(turn_matrix) <- c("Encamped", "Exploratory")
    rownames(turn_matrix) <- c("Mean", "Concentration")
    turn_real <- rownames_to_column(as.data.frame(turn_matrix), var = "Parameter") %>%
      mutate(Species = spp,
             Season = season,
             StudyArea = area,
             Encamped = round(Encamped, 2),
             Exploratory = round(Exploratory, 2)) %>%
      relocate(Species, .before = "Parameter") %>%
      relocate(StudyArea, .after = "Species") %>%
      relocate(Season, .after = "StudyArea")
      
    #'  Table of transition probabilties and 95% CI
    trans_probs <- ci_nat[[3]]$est
    trans_real <- rownames_to_column(as.data.frame(trans_probs), var = "States") %>%
      mutate(States = ifelse(States == "encamped", "Encamped", "Exploratory"),
             Species = spp,
             Season = season,
             StudyArea = area,
             encamped = round(encamped, 2),
             exploratory = round(exploratory, 2)) %>%
      relocate(Species, .before = "States") %>%
      relocate(StudyArea, .after = "Species") %>%
      relocate(Season, .after = "StudyArea") 
    colnames(trans_real) <- c("Species", "Study Area", "Season", "Start State", "Pr(To Encamped)", "Pr(To Exploratory)")
    # trans_state11 <- c(ci_nat[[3]]$est[[1]], ci_nat[[3]]$lower[[1]], ci_nat[[3]]$upper[[1]])
    # trans_state12 <- c(ci_nat[[3]]$est[[2]], ci_nat[[3]]$lower[[2]], ci_nat[[3]]$upper[[2]])
    # trans_state22 <- c(ci_nat[[3]]$est[[4]], ci_nat[[3]]$lower[[4]], ci_nat[[3]]$upper[[4]])
    # trans_state21 <- c(ci_nat[[3]]$est[[3]], ci_nat[[3]]$lower[[3]], ci_nat[[3]]$upper[[3]])
    # trans_real <- as.data.frame(rbind(trans_state11, trans_state12, trans_state22, trans_state21))
    # colnames(trans_real) <- c("Mean", "Lower", "Upper")
    # trans_real <- rownames_to_column(trans_real, var = "States") %>%
    #   mutate(Species = spp,
    #          Season = season,
    #          StudyArea = area,
    #          States = ifelse(States == "trans_state11", "Encamped to Encamped", States),
    #          States = ifelse(States == "trans_state22", "Exploratory to Exploratory", States),
    #          States = ifelse(States == "trans_state12", "Encamped to Exploratory", States),
    #          States = ifelse(States == "trans_state21", "Exploratory to Encamped", States), 
    #          Mean = round(Mean, 2), 
    #          Lower = round(Lower, 2), 
    #          Upper = round(Upper, 2)) %>%
    #   unite("95%CI", Lower:Upper, sep = " - ") %>%
    #   relocate(Species, .before = "States") %>%
    #   relocate(Season, .after = "Species") %>%
    #   relocate(StudyArea, .after = "Season")
 
    print(round(ci_nat[[1]]$est, 2))
    print(round(ci_nat[[2]]$est, 2))
    print(round(ci_nat[[3]]$est, 2))
    
    table_list <- list(steps_real, turn_real, trans_real)
    
    return(table_list)
  }
  md_smr_backtrans <- backtrans_params(spp_HMM_output[[1]], spp = "Mule Deer", season = "Summer", area = "Okanogan")
  md_wtr_backtrans <- backtrans_params(spp_HMM_output[[2]], spp = "Mule Deer", season = "Winter", area = "Okanogan")
  elk_smr_backtrans <- backtrans_params(spp_HMM_output[[3]], spp = "Elk", season = "Summer", area = "Northeast")
  elk_wtr_backtrans <- backtrans_params(spp_HMM_output[[4]], spp = "Elk", season = "Winter", area = "Northeast")
  wtd_smr_backtrans <- backtrans_params(spp_HMM_output[[5]], spp = "White-tailed Deer", season = "Summer", area = "Northeast")
  wtd_wtr_backtrans <- backtrans_params(spp_HMM_output[[6]], spp = "White-tailed Deer", season = "Winter", area = "Northeast")
  coug_smr_backtrans_OK <- backtrans_params(spp_HMM_output[[7]], spp = "Cougar", season = "Summer", area = "Okanogan")
  coug_wtr_backtrans_OK <- backtrans_params(spp_HMM_output[[8]], spp = "Cougar", season = "Winter", area = "Okanogan")
  coug_smr_backtrans_NE <- backtrans_params(spp_HMM_output[[9]], spp = "Cougar", season = "Summer", area = "Northeast")
  coug_wtr_backtrans_NE <- backtrans_params(spp_HMM_output[[10]], spp = "Cougar", season = "Winter", area = "Northeast")
  wolf_smr_backtrans_OK <- backtrans_params(spp_HMM_output[[11]], spp = "Wolf", season = "Summer", area = "Okanogan")
  wolf_wtr_backtrans_OK <- backtrans_params(spp_HMM_output[[12]], spp = "Wolf", season = "Winter", area = "Okanogan")
  wolf_smr_backtrans_NE <- backtrans_params(spp_HMM_output[[13]], spp = "Wolf", season = "Summer", area = "Northeast")
  wolf_wtr_backtrans_NE <- backtrans_params(spp_HMM_output[[14]], spp = "Wolf", season = "Winter", area = "Northeast")
  bob_smr_backtrans_OK <- backtrans_params(spp_HMM_output[[15]], spp = "Bobcat", season = "Summer", area = "Okanogan")
  #'  List order changes when using von Mises distribution
  # bob_wtr_backtrans_OK <- backtrans_params(spp_HMM_output[[16]], spp = "Bobcat", season = "Winter", area = "Okanogan")
  # bob_smr_backtrans_NE <- backtrans_params(spp_HMM_output[[17]], spp = "Bobcat", season = "Summer", area = "Northeast")
  bob_wtr_backtrans_NE <- backtrans_params(spp_HMM_output[[16]], spp = "Bobcat", season = "Winter", area = "Northeast")
  coy_smr_backtrans_OK <- backtrans_params(spp_HMM_output[[17]], spp = "Coyote", season = "Summer", area = "Okanogan")
  coy_wtr_backtrans_OK <- backtrans_params(spp_HMM_output[[18]], spp = "Coyote", season = "Winter", area = "Okanogan")
  coy_smr_backtrans_NE <- backtrans_params(spp_HMM_output[[19]], spp = "Coyote", season = "Summer", area = "Northeast")
  coy_wtr_backtrans_NE <- backtrans_params(spp_HMM_output[[20]], spp = "Coyote", season = "Winter", area = "Northeast")
  
  
  #'  Table for back-transformed step lengths
  all_steps_backtrans <- bind_rows(md_smr_backtrans[[1]], md_wtr_backtrans[[1]], 
                                   elk_smr_backtrans[[1]], elk_wtr_backtrans[[1]], 
                                   wtd_smr_backtrans[[1]], wtd_wtr_backtrans[[1]],
                                   coug_smr_backtrans_OK[[1]], coug_wtr_backtrans_OK[[1]], 
                                   coug_smr_backtrans_NE[[1]], coug_wtr_backtrans_NE[[1]], 
                                   wolf_smr_backtrans_OK[[1]], wolf_wtr_backtrans_OK[[1]], 
                                   wolf_smr_backtrans_NE[[1]], wolf_wtr_backtrans_NE[[1]], 
                                   bob_smr_backtrans_OK[[1]], 
                                   # bob_wtr_backtrans_OK[[1]], bob_smr_backtrans_NE[[1]], 
                                   bob_wtr_backtrans_NE[[1]], 
                                   coy_smr_backtrans_OK[[1]], coy_wtr_backtrans_OK[[1]], 
                                   coy_smr_backtrans_NE[[1]], coy_wtr_backtrans_NE[[1]]) %>%
    arrange(Species) %>%
    pivot_wider(names_from = "State", values_from = c("Mean", "95%CI")) %>%
    relocate('95%CI_Encamped', .after = "Mean_Encamped") %>%
    filter(!Species == "Bobcat")
  colnames(all_steps_backtrans) <- c("Species", "Study Area", "Season", "Mean Encamped", "95% CI Encamped", "Mean Exploratory", "95% CI Exploratory")
  
  #'  Table for back-transformed turning angles
  all_turns_backtrans <- bind_rows(md_smr_backtrans[[2]], md_wtr_backtrans[[2]], 
                                   elk_smr_backtrans[[2]], elk_wtr_backtrans[[2]], 
                                   wtd_smr_backtrans[[2]], wtd_wtr_backtrans[[2]],
                                   coug_smr_backtrans_OK[[2]], coug_wtr_backtrans_OK[[2]], 
                                   coug_smr_backtrans_NE[[2]], coug_wtr_backtrans_NE[[2]], 
                                   wolf_smr_backtrans_OK[[2]], wolf_wtr_backtrans_OK[[2]], 
                                   wolf_smr_backtrans_NE[[2]], wolf_wtr_backtrans_NE[[2]], 
                                   bob_smr_backtrans_OK[[2]], 
                                   # bob_wtr_backtrans_OK[[2]], bob_smr_backtrans_NE[[2]], 
                                   bob_wtr_backtrans_NE[[2]], 
                                   coy_smr_backtrans_OK[[2]], coy_wtr_backtrans_OK[[2]], 
                                   coy_smr_backtrans_NE[[2]], coy_wtr_backtrans_NE[[2]]) %>%
    arrange(Species) %>%
    filter(!Species == "Bobcat")
  colnames(all_turns_backtrans) <- c("Species", "Study Area", "Season", "Parameter", "Encamped", "Exploratory")
  
  #'  Table for back-transformed transition probabilities
  all_TransPr_backtrans <- bind_rows(md_smr_backtrans[[3]], md_wtr_backtrans[[3]], 
                                   elk_smr_backtrans[[3]], elk_wtr_backtrans[[3]], 
                                   wtd_smr_backtrans[[3]], wtd_wtr_backtrans[[3]],
                                   coug_smr_backtrans_OK[[3]], coug_wtr_backtrans_OK[[3]], 
                                   coug_smr_backtrans_NE[[3]], coug_wtr_backtrans_NE[[3]], 
                                   wolf_smr_backtrans_OK[[3]], wolf_wtr_backtrans_OK[[3]], 
                                   wolf_smr_backtrans_NE[[3]], wolf_wtr_backtrans_NE[[3]], 
                                   bob_smr_backtrans_OK[[3]], 
                                   # bob_wtr_backtrans_OK[[3]], bob_smr_backtrans_NE[[3]], 
                                   bob_wtr_backtrans_NE[[3]], 
                                   coy_smr_backtrans_OK[[3]], coy_wtr_backtrans_OK[[3]], 
                                   coy_smr_backtrans_NE[[3]], coy_wtr_backtrans_NE[[3]]) %>%
    arrange(Species) %>%
    filter(!Species == "Bobcat")
  
  #'  Save tables
  write.csv(all_steps_backtrans, paste0("./Outputs/HMM_output/HMM_Results_StepLength_BackTrans_", Sys.Date(), ".csv"))
  write.csv(all_turns_backtrans, paste0("./Outputs/HMM_output/HMM_Results_TurningAngle_BackTrans_", Sys.Date(), ".csv"))
  write.csv(all_TransPr_backtrans, paste0("./Outputs/HMM_output/HMM_Results_TransPr_BackTrans_", Sys.Date(), ".csv"))
  
  
  
  ####  Plot Stationary-State Probabilities  ####
  #'  Functions to extract stationary state probabilities & plot predicted responses
  stay_probs_prey <- function(hmmm) {
    #'  Calculate stationary state probs. for each state based on covariate data
    #'  for each time step
    stay_pr <- stationary(hmmm)
    stay_pr <- stay_pr[[1]]
    #'  Calculate stationary state probs. for each state when covariate data are
    #'  held at their mean value (0 b/c data are centered and scaled)
    stay_mu0 <- stationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0, 
                                                   SnowCover = 0, TRI = 0, 
                                                   COUG_RSF = 0, WOLF_RSF = 0))#,
                                                   #BOB_RSF = 0, COY_RSF = 0))
    print(stay_mu0)
    #'  Plot stationary state probabilities and extract predicted estimates
    fig <- plotStationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                  SnowCover = 0, TRI = 0, 
                                                  COUG_RSF = 0, WOLF_RSF = 0),
                                                  #BOB_RSF = 0, COY_RSF = 0),
                          col = c("red", "blue"), plotCI = TRUE, alpha = 0.95, return =  TRUE) 
    stationary_probs <- list(stay_pr, fig)
    
    return(stationary_probs)
  }
  #'  Extract stationary state probabilities for deer and elk
  stay_md_smr <- stay_probs_prey(spp_HMM_output[[1]])
  stay_md_wtr <- stay_probs_prey(spp_HMM_output[[2]])
  stay_elk_smr <- stay_probs_prey(spp_HMM_output[[3]])
  stay_elk_wtr <- stay_probs_prey(spp_HMM_output[[4]])
  stay_wtd_smr <- stay_probs_prey(spp_HMM_output[[5]])
  stay_wtd_wtr <- stay_probs_prey(spp_HMM_output[[6]])
  
  #'  Stationary probabilities for predators in the Okanogan
  stay_probs_pred_OK <- function(hmmm) {
    stay_pr <- stationary(hmmm)
    stay_pr <- stay_pr[[1]]
    stay_mu0 <- stationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                   SnowCover = 0, TRI = 0, MD_RSF = 0)) 
    print(stay_mu0) 
    #'  Plot stationary state probabilities and extract predicted estimates
    fig <- plotStationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                  SnowCover = 0, TRI = 0, MD_RSF = 0),  
                          col = c("red", "blue"), plotCI = TRUE, alpha = 0.95, return =  TRUE) 
    stationary_probs <- list(stay_pr, fig)
    
    return(stationary_probs)
  }
  #'  Extract stationary state probabilities
  stay_coug_smr_OK <- stay_probs_pred_OK(spp_HMM_output[[7]])
  stay_coug_wtr_OK <- stay_probs_pred_OK(spp_HMM_output[[8]])
  stay_wolf_smr_OK <- stay_probs_pred_OK(spp_HMM_output[[11]])
  stay_wolf_wtr_OK <- stay_probs_pred_OK(spp_HMM_output[[12]])
  stay_bob_smr_OK <- stay_probs_pred_OK(spp_HMM_output[[15]])
  # stay_bob_wtr_OK <- stay_probs_pred_OK(spp_HMM_output[[16]])
  stay_coy_smr_OK <- stay_probs_pred_OK(spp_HMM_output[[17]])
  stay_coy_wtr_OK <- stay_probs_pred_OK(spp_HMM_output[[18]])
  
  #'  Stationary state probabilities for predators in the Northeast
  stay_probs_pred_NE <- function(hmmm) {
    stay_pr <- stationary(hmmm)
    stay_pr <- stay_pr[[1]]
    stay_mu0 <- stationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                   SnowCover = 0, TRI = 0, 
                                                   ELK_RSF = 0, WTD_RSF = 0))
    print(stay_mu0) 
    #'  Plot stationary state probabilities and extract predicted estimates
    fig <- plotStationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                  SnowCover = 0, TRI = 0, 
                                                  ELK_RSF = 0, WTD_RSF = 0),    
                          col = c("red", "blue"), plotCI = TRUE, alpha = 0.95, return =  TRUE)
    stationary_probs <- list(stay_pr, fig)
    
    return(stationary_probs)
  }
  #'  Extract stationary state probabilities
  stay_coug_smr_NE <- stay_probs_pred_NE(spp_HMM_output[[9]])
  stay_coug_wtr_NE <- stay_probs_pred_NE(spp_HMM_output[[10]])
  stay_wolf_smr_NE <- stay_probs_pred_NE(spp_HMM_output[[13]])
  stay_wolf_wtr_NE <- stay_probs_pred_NE(spp_HMM_output[[14]])
  # stay_bob_smr_NE <- stay_probs_pred_NE(spp_HMM_output[[17]])
  stay_bob_wtr_NE <- stay_probs_pred_NE(spp_HMM_output[[16]])
  stay_coy_smr_NE <- stay_probs_pred_NE(spp_HMM_output[[19]])
  stay_coy_wtr_NE <- stay_probs_pred_NE(spp_HMM_output[[20]])
  
  
  
  ####  Prettier Plots for Stationary State Probabilities  ####
  #'  Function to extract stationary state probabilities and plot outputs
  stay_plots <- function(stay, season, spp, area) {
    #'  Extract list of calculated stationary states for range of covariate values 
    #'  from HMM stationary output
    stay_covs <- stay[[2]]
    #'  Create empty list
    covs_out <- list()
    #'  Loop through all list elements (results for each covariate)
    for(l in 1:length(stay_covs)){
      #'  Hold list of interest
      cov <- stay_covs[[l]]
      #'  Add column indicating which behavioral state values belong to
      cov[[1]]$State <- "Encamped"
      cov[[2]]$State <- "Exploratory"
      #'  Convert to data frame instead of list
      cov <- rbind(as.data.frame(cov[[1]]), as.data.frame(cov[[2]]))
      cov$State <- as.factor(cov$State)
      #'  Append to new list of data frames
      covs_out[[l]] <- cov
    }
    #'  Rename list elements based on covariate
    names(covs_out) <- names(stay_covs)
    
    #'  Extract names of list elements and clean up for plotting
    list_names <- as.data.frame(names(covs_out))
    colnames(list_names) <- "nms"
    list_names <- list_names %>%
      mutate(nms = ifelse(nms == "PercOpen", "Habitat Openness", nms),
             nms = ifelse(nms == "Dist2Road", "Distance to Road", nms),
             nms = ifelse(nms == "SnowCover", "Snow Cover", nms),
             nms = ifelse(nms == "MD_RSF", "Mule Deer Presence", nms),
             nms = ifelse(nms == "ELK_RSF", "Elk Presence", nms),
             nms = ifelse(nms == "WTD_RSF", "White-tailed Deer Presence", nms),
             nms = ifelse(nms == "COUG_RSF", "Cougar Presence", nms),
             nms = ifelse(nms == "WOLF_RSF", "Wolf Presence", nms),
             nms = ifelse(nms == "BOB_RSF", "Bobcat Presence", nms),
             nms = ifelse(nms == "COY_RSF", "Coyote Presence", nms))
    #'  Force back to an atomic vector of characters (needed for looping below)
    list_names <- list_names$nms
    
    #'  Create empty list
    stay_figs <- list()
    #'  Loop through each data frame
    for(l in 1:length(covs_out)){
      #'  Create a figure plotting the stationary state probabilities & 95% CI
      stay_plot <- ggplot(covs_out[[l]], aes(x = cov, y = est, group = State)) + 
        geom_line(aes(color = State)) + 
        #'  Add confidence intervals
        geom_ribbon(aes(ymin = lci, ymax = uci, fill = State), alpha = 0.2) +
        #'  Get rid of lines and gray background
        theme_bw() +
        theme(panel.border = element_blank()) +
        theme(axis.line = element_line(color = 'black')) +
        #'  Force y-axis from 0 to 1
        ylim(0,1.0) +
        #'  Use list name as X-axis title
        xlab(list_names[l]) +
        ylab("Stationary State Probability") +
        labs(#title = paste(area, season, spp, "Stationary State Probabilities"), 
             fill = "Movement State", color = "Movement State") 
        #theme(legend.position="bottom")
      #'  Review figure
      plot(stay_plot)
      #'  Append figures
      stay_figs[[l]] <- stay_plot
    }
    
    return(stay_figs)
  }
  #'  Run each species through- for loops should allow the different coefficients
  #'  to still plot nicely
  md_smr_fig <- stay_plots(stay_md_smr, season = "Summer", spp = "Mule Deer", area = "Okanogan")
  md_wtr_fig <- stay_plots(stay_md_wtr, season = "Winter", spp = "Mule Deer", area = "Okanogan")
  elk_smr_fig <- stay_plots(stay_elk_smr, season = "Summer", spp = "Elk", area = "Northeast")
  elk_wtr_fig <- stay_plots(stay_elk_wtr, season = "Winter", spp = "Elk", area = "Northeast")
  wtd_smr_fig <- stay_plots(stay_wtd_smr, season = "Summer", spp = "White-tailed Deer", area = "Northeast")
  wtd_wtr_fig <- stay_plots(stay_wtd_wtr, season = "Winter", spp = "White-tailed Deer", area = "Northeast")
  coug_smr_OK_fig <- stay_plots(stay_coug_smr_OK, season = "Summer", spp = "Cougar", area = "Okanogan")
  coug_wtr_OK_fig <- stay_plots(stay_coug_wtr_OK, season = "Winter", spp = "Cougar", area = "Okanogan")
  coug_smr_NE_fig <- stay_plots(stay_coug_smr_NE, season = "Summer", spp = "Cougar", area = "Northeast")
  coug_wtr_NE_fig <- stay_plots(stay_coug_wtr_NE, season = "Winter", spp = "Cougar", area = "Northeast")
  wolf_smr_OK_fig <- stay_plots(stay_wolf_smr_OK, season = "Summer", spp = "Wolf", area = "Okanogan")
  wolf_wtr_OK_fig <- stay_plots(stay_wolf_wtr_OK, season = "Winter", spp = "Wolf", area = "Okanogan")
  wolf_smr_NE_fig <- stay_plots(stay_wolf_smr_NE, season = "Summer", spp = "Wolf", area = "Northeast")
  wolf_wtr_NE_fig <- stay_plots(stay_wolf_wtr_NE, season = "Winter", spp = "Wolf", area = "Northeast")
  coy_smr_OK_fig <- stay_plots(stay_coy_smr_OK, season = "Summer", spp = "Coyote", area = "Okanogan")
  coy_wtr_OK_fig <- stay_plots(stay_coy_wtr_OK, season = "Winter", spp = "Coyote", area = "Okanogan")
  coy_smr_NE_fig <- stay_plots(stay_coy_smr_NE, season = "Summer", spp = "Coyote", area = "Northeast")
  coy_wtr_NE_fig <- stay_plots(stay_coy_wtr_NE, season = "Winter", spp = "Coyote", area = "Northeast")
  
  
  #'  Patchwork figures together in panels
  library(patchwork)
  #'  MULE DEER panels
  length(md_smr_fig)
  (md_smr_patch <- md_smr_fig[[1]] + md_smr_fig[[2]] + theme(axis.title.y = element_blank()) + 
      md_smr_fig[[3]] + md_smr_fig[[4]] + theme(axis.title.y = element_blank()) + 
      md_smr_fig[[5]] + theme(axis.title.y = element_blank()) + # md_smr_fig +
      guide_area() + plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Summer Mule Deer Stationary State Probabilities',
                      subtitle = '     Okanogan 2018 - 2021') + plot_layout(ncol = 2))
  length(md_wtr_fig)
  (md_wtr_patch <- md_wtr_fig[[1]] + md_wtr_fig[[2]] + theme(axis.title.y = element_blank()) + 
      md_wtr_fig[[3]] + theme(axis.title.y = element_blank()) + md_wtr_fig[[4]] + 
      md_wtr_fig[[5]] + theme(axis.title.y = element_blank()) + md_wtr_fig[[6]] + 
      theme(axis.title.y = element_blank()) + #md_wtr_fig[[7]] + md_wtr_fig[[8]] + 
      theme(axis.title.y = element_blank()) + #guide_area() + 
      plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Winter Mule Deer Stationary State Probabilities',
                      subtitle = '     Okanogan 2018 - 2021') + plot_layout(ncol = 2))
  #'  ELK panels
  length(elk_smr_fig)
  (elk_smr_patch <- elk_smr_fig[[1]] + elk_smr_fig[[2]] + theme(axis.title.y = element_blank()) + 
      elk_smr_fig[[3]] + theme(axis.title.y = element_blank()) + elk_smr_fig[[4]] + 
      elk_smr_fig[[5]] + theme(axis.title.y = element_blank()) + #elk_smr_fig[[6]] + elk_smr_fig[[7]] +
      theme(axis.title.y = element_blank()) + plot_layout(guides = 'collect') + 
      guide_area() + plot_annotation(title = 'Summer Elk Stationary State Probabilities',
                                     subtitle = '     Northeast 2018 - 2021') + plot_layout(ncol = 3))
  length(elk_wtr_fig)
  (elk_wtr_patch <- elk_wtr_fig[[1]] + elk_wtr_fig[[2]] + theme(axis.title.y = element_blank()) + 
      elk_wtr_fig[[3]] + theme(axis.title.y = element_blank()) + elk_wtr_fig[[4]] + 
      elk_wtr_fig[[5]] + theme(axis.title.y = element_blank()) + elk_wtr_fig[[6]] + 
      theme(axis.title.y = element_blank()) + #elk_wtr_fig[[7]] + elk_wtr_fig[[8]] + 
      theme(axis.title.y = element_blank()) + guide_area() + plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Winter Elk Stationary State Probabilities',
                      subtitle = '    Northeast 2018 - 2021') + plot_layout(ncol = 2))
  #'  WHITE-TAILED DEER panels
  length(wtd_smr_fig)
  (wtd_smr_patch <- wtd_smr_fig[[1]] + wtd_smr_fig[[2]] + theme(axis.title.y = element_blank()) + 
      wtd_smr_fig[[3]] + theme(axis.title.y = element_blank()) + wtd_smr_fig[[4]] + 
      wtd_smr_fig[[5]] + theme(axis.title.y = element_blank()) + #wtd_smr_fig[[6]] + wtd_smr_fig[[7]] + 
      theme(axis.title.y = element_blank()) + plot_layout(guides = 'collect') + 
      guide_area() + plot_annotation(title = 'Summer White-tailed Deer Stationary State Probabilities',
                                     subtitle = '     Northeast 2018 - 2021') + plot_layout(ncol = 3))
  length(wtd_wtr_fig)
  (wtd_wtr_patch <- wtd_wtr_fig[[1]] + wtd_wtr_fig[[2]] + theme(axis.title.y = element_blank()) + 
      wtd_wtr_fig[[3]] + theme(axis.title.y = element_blank()) + wtd_wtr_fig[[4]] + 
      wtd_wtr_fig[[5]] + theme(axis.title.y = element_blank()) + wtd_wtr_fig[[6]] + 
      theme(axis.title.y = element_blank()) + #wtd_wtr_fig[[7]] + wtd_wtr_fig[[8]] + guide_area() + 
      theme(axis.title.y = element_blank()) + plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Winter White-tailed Deer Stationary State Probabilities',
                      subtitle = '     Northeast 2018 - 2021') + plot_layout(ncol = 3))
  #'  COUGAR panels
  length(coug_smr_OK_fig)
  (coug_smr_OK_patch <- coug_smr_OK_fig[[1]] + coug_smr_OK_fig[[2]] + 
      theme(axis.title.y = element_blank()) + coug_smr_OK_fig[[3]] + 
      coug_smr_OK_fig[[4]] + theme(axis.title.y = element_blank()) + 
      plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Summer Cougar Stationary State Probabilities',
                      subtitle = '     Okanogan 2018 - 2021') + plot_layout(ncol = 2))
  length(coug_wtr_OK_fig)
  (coug_wtr_OK_patch <- coug_wtr_OK_fig[[1]] + coug_wtr_OK_fig[[2]] + 
      theme(axis.title.y = element_blank()) + coug_wtr_OK_fig[[3]] + 
      coug_wtr_OK_fig[[4]] + theme(axis.title.y = element_blank()) + 
      plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Winter Cougar Stationary State Probabilities',
                      subtitle = '     Okanogan 2018 - 2021') + plot_layout(ncol = 2))
  length(coug_smr_NE_fig)
  (coug_smr_NE_patch <- coug_smr_NE_fig[[1]] + coug_smr_NE_fig[[2]] + 
      theme(axis.title.y = element_blank()) + coug_smr_NE_fig[[3]] + 
      coug_smr_NE_fig[[4]] + theme(axis.title.y = element_blank()) + 
      coug_smr_NE_fig[[5]] + guide_area() + plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Summer Cougar Stationary State Probabilities',
                      subtitle = '     Northeast 2018 - 2021') + plot_layout(ncol = 2))
  length(coug_wtr_NE_fig)
  (coug_wtr_NE_patch <- coug_wtr_NE_fig[[1]] + coug_wtr_NE_fig[[2]] + 
      theme(axis.title.y = element_blank()) + coug_wtr_NE_fig[[3]] + 
      coug_wtr_NE_fig[[4]] + theme(axis.title.y = element_blank()) + 
      coug_wtr_NE_fig[[5]] + coug_wtr_NE_fig[[6]] + 
      theme(axis.title.y = element_blank()) + plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Winter Cougar Stationary State Probabilities',
                      subtitle = '     Northeast 2018 - 2021') + plot_layout(ncol = 2))
  #'  WOLF panels
  length(wolf_smr_OK_fig)
  (wolf_smr_OK_patch <- wolf_smr_OK_fig[[1]] + wolf_smr_OK_fig[[2]] + 
      theme(axis.title.y = element_blank()) + wolf_smr_OK_fig[[3]] + 
      wolf_smr_OK_fig[[4]] + theme(axis.title.y = element_blank()) + 
      plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Summer Wolf Stationary State Probabilities',
                    subtitle = '     Okanogan 2018 - 2021') + plot_layout(ncol = 2))
  length(wolf_wtr_OK_fig)
  (wolf_wtr_OK_patch <- wolf_wtr_OK_fig[[1]] + wolf_wtr_OK_fig[[2]] + 
      theme(axis.title.y = element_blank()) + wolf_wtr_OK_fig[[3]] + 
      wolf_wtr_OK_fig[[4]] + theme(axis.title.y = element_blank()) + 
      wolf_wtr_OK_fig[[5]] + guide_area() + plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Winter Wolf Stationary State Probabilities',
                    subtitle = '     Okanogan 2018 - 2021') + plot_layout(ncol = 2))
  length(wolf_smr_NE_fig)
  (wolf_smr_NE_patch <- wolf_smr_NE_fig[[1]] + wolf_smr_NE_fig[[2]] + 
      theme(axis.title.y = element_blank()) + wolf_smr_NE_fig[[3]] + 
      wolf_smr_NE_fig[[4]] + theme(axis.title.y = element_blank()) + 
      wolf_smr_NE_fig[[5]] + guide_area() + plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Summer Wolf Stationary State Probabilities',
                    subtitle = '     Northeast 2018 - 2021') + plot_layout(ncol = 2))
  length(wolf_wtr_NE_fig)
  (wolf_wtr_NE_patch <- wolf_wtr_NE_fig[[1]] + wolf_wtr_NE_fig[[2]] + 
      theme(axis.title.y = element_blank()) + wolf_wtr_NE_fig[[3]] + 
      wolf_wtr_NE_fig[[4]] + theme(axis.title.y = element_blank()) + 
      wolf_wtr_NE_fig[[5]] + wolf_wtr_NE_fig[[6]] + theme(axis.title.y = element_blank()) + 
      plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Winter Wolf Stationary State Probabilities',
                    subtitle = '     Northeast 2018 - 2021') + plot_layout(ncol = 2))
  #'  COYOTE panels
  length(coy_smr_OK_fig)
  (coy_smr_OK_patch <- coy_smr_OK_fig[[1]] + coy_smr_OK_fig[[2]] + 
      theme(axis.title.y = element_blank()) + coy_smr_OK_fig[[3]] + 
      coy_smr_OK_fig[[4]] + theme(axis.title.y = element_blank()) + 
      plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Summer Coyote Stationary State Probabilities',
                      subtitle = '     Okanogan 2018 - 2021') + plot_layout(ncol = 2))
  length(coy_wtr_OK_fig)
  (coy_wtr_OK_patch <- coy_wtr_OK_fig[[1]] + coy_wtr_OK_fig[[2]] + 
      theme(axis.title.y = element_blank()) + coy_wtr_OK_fig[[3]] + 
      coy_wtr_OK_fig[[4]] + theme(axis.title.y = element_blank()) + 
      plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Winter Coyote Stationary State Probabilities',
                      subtitle = '     Okanogan 2018 - 2021') + plot_layout(ncol = 2))
  length(coy_smr_NE_fig)
  (coy_smr_NE_patch <- coy_smr_NE_fig[[1]] + coy_smr_NE_fig[[2]] + 
      theme(axis.title.y = element_blank()) + coy_smr_NE_fig[[3]] + 
      coy_smr_NE_fig[[4]] + theme(axis.title.y = element_blank()) + 
      coy_smr_NE_fig[[5]] + guide_area() + plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Summer Coyote Stationary State Probabilities',
                      subtitle = '     Northeast 2018 - 2021') + plot_layout(ncol = 2))
  length(coy_wtr_NE_fig)
  (coy_wtr_NE_patch <- coy_wtr_NE_fig[[1]] + coy_wtr_NE_fig[[2]] + 
      theme(axis.title.y = element_blank()) + coy_wtr_NE_fig[[3]] + 
      coy_wtr_NE_fig[[4]] + theme(axis.title.y = element_blank()) + 
      coy_wtr_NE_fig[[5]] + coy_wtr_NE_fig[[6]] + plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Winter Coyote Stationary State Probabilities',
                      subtitle = '     Northeast 2018 - 2021') + plot_layout(ncol = 2))
  
  
  pdf(file = "./Outputs/HMM_output/Stationary_State_Prob_Plots_04.08.23.pdf")
  plot(md_smr_patch, main = "Stationary State Probabilties for Summer Mule Deer")
  plot(md_wtr_patch, main = "Stationary State Probabilties for Winter Mule Deer")
  plot(elk_smr_patch, main = "Stationary State Probabilties for Summer Elk")
  plot(elk_wtr_patch, main = "Stationary State Probabilties for Winter Elk")
  plot(wtd_smr_patch, main = "Stationary State Probabilties for Summer White-tailed Deer")
  plot(wtd_wtr_patch, main = "Stationary State Probabilties for Winter White-tailed Deer")
  plot(coug_smr_OK_patch, main = "Stationary State Probabilties for Summer Cougar, OK")
  plot(coug_wtr_OK_patch, main = "Stationary State Probabilties for Winter Cougar, OK")
  plot(coug_smr_NE_patch, main = "Stationary State Probabilties for Summer Cougar, NE")
  plot(coug_wtr_NE_patch, main = "Stationary State Probabilties for Winter Cougar, NE")
  plot(wolf_smr_OK_patch, main = "Stationary State Probabilties for Summer Wolf, OK")
  plot(wolf_wtr_OK_patch, main = "Stationary State Probabilties for Winter Wolf, OK")
  plot(wolf_smr_NE_patch, main = "Stationary State Probabilties for Summer Wolf, NE")
  plot(wolf_wtr_NE_patch, main = "Stationary State Probabilties for Winter Wolf, NE")
  plot(coy_smr_OK_patch, main = "Stationary State Probabilties for Summer Coyote, OK")
  plot(coy_wtr_OK_patch, main = "Stationary State Probabilties for Winter Coyote, OK")
  plot(coy_smr_NE_patch, main = "Stationary State Probabilties for Summer Coyote, NE")
  plot(coy_wtr_NE_patch, main = "Stationary State Probabilties for Winter Coyote, NE")
  dev.off()
  
  #'  Save individual plots
  png(file="./Outputs/Figures for ms/MD_smr_TRI.png", width = 700, height = 500)
  (md_smr_tri <- md_smr_fig[[1]] + plot_annotation(title = 'Summer Mule Deer Stationary State Probabilities', subtitle = '     Okanogan 2018 - 2021'))
  dev.off()
  png(file="./Outputs/Figures for ms/MD_smr_WOLF.png", width = 700, height = 500)
  (md_smr_wolf <- md_smr_fig[[5]] + plot_annotation(title = 'Summer Mule Deer Stationary State Probabilities', subtitle = '     Okanogan 2018 - 2021'))
  dev.off()
  png(file="./Outputs/Figures for ms/MD_smr_BOB.png", width = 700, height = 500)
  (md_smr_bob <- md_smr_fig[[6]] + plot_annotation(title = 'Summer Mule Deer Stationary State Probabilities', subtitle = '     Okanogan 2018 - 2021'))
  dev.off()
  png(file="./Outputs/Figures for ms/MD_wtr_COUG.png", width = 700, height = 500)
  (md_wtr_coug <- md_wtr_fig[[5]] + plot_annotation(title = 'Winter Mule Deer Stationary State Probabilities', subtitle = '     Okanogan 2018 - 2021'))
  dev.off()
  png(file="./Outputs/Figures for ms/MD_wtr_COY.png", width = 700, height = 500)
  (md_wtr_coy <- md_wtr_fig[[8]] + plot_annotation(title = 'Winter Mule Deer Stationary State Probabilities', subtitle = '     Okanogan 2018 - 2021'))
  dev.off()
  
  
  png(file="./Outputs/Figures for ms/WTD_smr_COY.png", width = 700, height = 500)
  (wtd_smr_coy <- wtd_smr_fig[[7]] + plot_annotation(title = 'Summer White-tailed Deer Stationary State Probabilities', subtitle = '     Okanogan 2018 - 2021'))
  dev.off()
  png(file="./Outputs/Figures for ms/ELK_wtr_WOLF.png", width = 700, height = 500)
  (elk_wtr_wolf <- elk_wtr_fig[[6]] + plot_annotation(title = 'Winter Elk Stationary State Probabilities', subtitle = '     Northeast 2018 - 2021'))
  dev.off()

  

  #' ####  Tables for manuscript  ####
  #' #'  Reformat transition probability tables for manuscript
  #' #'  Prey HMM results
  #' results_hmm_TransPr_prey_ms <- results_hmm_TransPr_prey %>%  
  #'   unite(CI95, Lower, Upper, sep = ", ") %>%
  #'   unite(Est_CI, Estimate, CI95, sep = "_") %>%
  #'   dplyr::select(-SE) %>%
  #'   spread(Parameter, Est_CI) %>%
  #'   separate("(Intercept)", c("Intercept", "Intercept 95% CI"), sep = "_") %>%
  #'   separate("TRI", c("TRI", "TRI 95% CI"), sep = "_") %>%
  #'   separate("PercOpen", c("Percent Open", "Percent Open 95% CI"), sep = "_") %>%
  #'   separate("Dist2Road", c("Nearest Road", "Nearest Road 95% CI"), sep = "_") %>%
  #'   separate("SnowCover1", c("Snow Cover (Y)", "Snow Cover (Y) 95% CI"), sep = "_") %>%
  #'   separate("COUG_RSF", c("Pr(Cougar)", "Pr(Cougar) 95% CI"), sep = "_") %>%
  #'   separate("WOLF_RSF", c("Pr(Wolf)", "Pr(Wolf) 95% CI"), sep = "_") %>%
  #'   # separate("BOB_RSF", c("Pr(Bobcat)", "Pr(Bobcat) 95% CI"), sep = "_") %>%
  #'   # separate("COY_RSF", c("Pr(Coyote)", "Pr(Coyote) 95% CI"), sep = "_") %>%
  #'   arrange(match(Species, c("Mule Deer", "Elk", "White-tailed Deer")))
  #' 
  #' write.csv(results_hmm_TransPr_prey_ms, paste0("./Outputs/HMM_output/HMM_Results_TransPr_prey_forMS_", Sys.Date(), ".csv"))
  #' 
  #' #'  Predators HMM results
  #' results_hmm_TransPr_pred_ms <- results_hmm_TransPr_pred %>% 
  #'   unite(CI95, Lower, Upper, sep = ", ") %>%
  #'   unite(Est_CI, Estimate, CI95, sep = "_") %>%
  #'   dplyr::select(-SE) %>%
  #'   spread(Parameter, Est_CI) %>%
  #'   separate("(Intercept)", c("Intercept (SE)", "Intercept 95% CI"), sep = "_") %>%
  #'   separate("TRI", c("TRI (SE)", "TRI 95% CI"), sep = "_") %>%
  #'   separate("PercOpen", c("Percent Open (SE)", "Percent Open 95% CI"), sep = "_") %>%
  #'   separate("Dist2Road", c("Nearest Road (SE)", "Nearest Road 95% CI"), sep = "_") %>%
  #'   separate("SnowCover1", c("Snow Cover (Y) (SE)", "Snow Cover (Y) 95% CI"), sep = "_") %>%
  #'   separate("MD_RSF", c("Pr(Mule Deer) (SE)", "Pr(Mule Deer) 95% CI"), sep = "_") %>%
  #'   separate("ELK_RSF", c("Pr(Elk) (SE)", "Pr(Elk) 95% CI"), sep = "_") %>%
  #'   separate("WTD_RSF", c("Pr(White-tailed Deer) (SE)", "Pr(White-tailed Deer) 95% CI"), sep = "_") %>%
  #'   filter(!Species == "Bobcat") %>%
  #'   group_by(Species) %>%
  #'   arrange(match(`Study Area`, c("Okanogan", "Northeast")), .by_group = TRUE) %>%
  #'   ungroup()
  #' 
  #' write.csv(results_hmm_TransPr_pred_ms, paste0("./Outputs/HMM_output/HMM_Results_TransPr_pred_forMS_", Sys.Date(), ".csv"))
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ####  Viterbi Algorithm  ####
  #'  Function to extract most likely state sequence for all locations based on
  #'  the Viterbi algorithm and the fitted HMM
  loc_states <- function(mod, locs) {
    #'  Decode most likely state for each observation
    states <- viterbi(mod)
    #'  Append state classification to location data
    dat <- cbind(locs, states)

    #'  Print derived percentage of time spent in each state
    print(table(states)/nrow(locs))

    return(dat)
  }
  #'  Run fitted HMMs and location data through function
  md_state_smr <- loc_states(spp_HMM_output[[1]], hmm_data[[1]])
  #' md_state_wtr <- loc_states(md_HMM_wtr, mdData_wtr)
  #' elk_state_smr <- loc_states(elk_HMM_smr, elkData_smr)
  #' elk_state_wtr <- loc_states(elk_HMM_wtr, elkData_wtr)
  #' wtd_state_smr <- loc_states(wtd_HMM_smr, wtdData_smr)
  #' wtd_state_wtr <- loc_states(wtd_HMM_wtr, wtdData_wtr)
  #' coug_state_smr <- loc_states(coug_HMM_smr, cougData_smr)
  #' coug_state_wtr <- loc_states(coug_HMM_wtr, cougData_wtr)
  #' wolf_state_smr <- loc_states(wolf_HMM_smr, wolfData_smr)
  #' wolf_state_wtr <- loc_states(wolf_HMM_wtr, wolfData_wtr)
  #' bob_state_smr <- loc_states(bob_HMM_smr, bobData_smr)
  #' bob_state_wtr <- loc_states(bob_HMM_wtr, bobData_wtr)
  #' coy_state_smr <- loc_states(coy_HMM_smr, coyData_smr)
  #' coy_state_wtr <- loc_states(coy_HMM_wtr, coyData_wtr)
  #' 
  #' #'  Save state sequences for external analyses
  #' spp_state_output <- list(md_state_smr, md_state_wtr, elk_state_smr, elk_state_wtr, 
  #'                          wtd_state_smr, wtd_state_wtr, coug_state_smr, coug_state_wtr, 
  #'                          wolf_state_smr, wolf_state_wtr, bob_state_smr, 
  #'                          bob_state_wtr, coy_state_smr, coy_state_wtr)
  #' #'  Make sure to note whether covariates were included on transition probabilities
  #' # save(spp_state_output, file = paste0("./Outputs/spp_state_output_", Sys.Date(), ".RData"))
  #' save(spp_state_output, file = paste0("./Outputs/spp_state_output_NULLtrans_", Sys.Date(), ".RData"))
  
 
  
  #'  Make panel of figures
  #'  https://www.benjaminbell.co.uk/2018/02/creating-multi-panel-plots-and-figures.html
  layout(matrix(1:8, ncol=2, byrow=TRUE))
  #'  Adjusts the margins
  par(oma=c(4, 4, 4, 4), mar=c(4, 4, 4, 4))
  stay_md_smr <- stay_probs(spp_HMM_output[[1]])
  
  layout(matrix(1:8, ncol=2, byrow=TRUE))
  par(oma=c(4, 4, 4, 4), mar=c(4, 4, 4, 4))
  stay_md_wtr <- stay_probs(spp_HMM_output[[2]])
  
  layout(matrix(1:8, ncol=2, byrow=TRUE))
  par(oma=c(4, 4, 4, 4), mar=c(4, 4, 4, 4))
  stay_coug_smr <- stay_probs(spp_HMM_output[[7]])
  
  layout(matrix(1:8, ncol=2, byrow=TRUE))
  par(oma=c(4, 4, 4, 4), mar=c(4, 4, 4, 4))
  stay_coug_wtr <- stay_probs(spp_HMM_output[[8]])
  
  

  


  
  
  