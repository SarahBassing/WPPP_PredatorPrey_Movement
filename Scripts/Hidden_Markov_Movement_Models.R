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
  library(tidyverse)

  #'  Load crwOut & covaraite data
  load("./Outputs/Telemetry_crwOut/crwOut_ALL_2022-02-03.RData") 
  load("./Outputs/Telemetry_covs/spp_telem_covs_noNDVI_2022-02-07.RData") 
  
  
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
    crwlMerge$crwPredict$BOB_RSF <- scale(crwlMerge$crwPredict$BOB_RSF)
    crwlMerge$crwPredict$COY_RSF <- scale(crwlMerge$crwPredict$COY_RSF)
    crwlMerge$crwPredict$hour <- as.integer(strftime(crwlMerge$crwPredict$time, format = "%H", tz="Etc/GMT+8"))
    #'  Prep crwlMerge data for fitHMM function
    Data <- prepData(data = crwlMerge, covNames = c("Dist2Road", "PercOpen", #"NDVI", 
                                                    "SnowCover", "TRI", "MD_RSF", 
                                                    "COUG_RSF", "WOLF_RSF", "BOB_RSF", 
                                                    "COY_RSF", "hour", "Sex", 
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
    crwlMerge$crwPredict$BOB_RSF <- scale(crwlMerge$crwPredict$BOB_RSF)
    crwlMerge$crwPredict$COY_RSF <- scale(crwlMerge$crwPredict$COY_RSF)
    crwlMerge$crwPredict$hour <- as.integer(strftime(crwlMerge$crwPredict$time, format = "%H", tz="Etc/GMT+8"))
    #'  Prep crwlMerge data for fitHMM function
    Data <- prepData(data = crwlMerge, covNames = c("Dist2Road", "PercOpen", #"NDVI", 
                                                    "SnowCover", "TRI", "ELK_RSF", 
                                                    "WTD_RSF", "COUG_RSF", "WOLF_RSF",
                                                    "BOB_RSF", "COY_RSF", "hour",
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
  save(hmm_data, file = paste0("./Outputs/Telemetry_crwOut/crwOut_ALL_wCovs_", Sys.Date(), ".RData"))
  
  
  
  #'  Correlation Matrix
  #'  ==================
  #'  Function to create correlation matrix for all continuous covariates at once
  cov_correlation_OK <- function(dat) {
    covs <- dat %>%
      dplyr::select(c("Dist2Road", "PercOpen", "TRI", "MD_RSF", "COUG_RSF", 
                      "WOLF_RSF", "BOB_RSF", "COY_RSF")) #"NDVI",
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
      dplyr::select(c("Dist2Road", "PercOpen", "TRI", "ELK_RSF", "WTD_RSF", 
                      "COUG_RSF", "WOLF_RSF", "BOB_RSF", "COY_RSF")) #"NDVI",
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
  #' plot(mdData_smr)  #250, 500, 250, 500
  #' plot(elkData_smr)  #500, 1000, 500, 1000
  #' plot(wtdData_smr)  #100, 500, 100, 500
  #' plot(cougData_smr_OK)  #500, 1500, 500, 1500
  #' plot(cougData_smr_NE)  #500, 1500, 500, 1500
  #' plot(wolfData_smr_OK)  #500, 3000, 500, 3000
  #' plot(wolfData_smr_NE)  #500, 3000, 500, 3000
  #' plot(bobData_smr_OK)  #500, 1000, 500, 1000
  #' plot(bobData_smr_NE)  #500, 1000, 500, 1000
  #' plot(coyData_smr_OK)  #500, 2000, 500, 2000
  #' plot(coyData_smr_NE)  #500, 2000, 500, 2000
  
  
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
  
  
  
  ####  Initial model set up  ####
  #'  ============================
  #'  Define initial parameters associated with each distribution & each state
  #'  Species-specific parameters based on viewing plotted data
  Par0_m1_md <- list(step = c(250, 500, 250, 500, 0.01, 0.005), angle = c(0.3, 0.7))  #zero-mass params needed
  Par0_m1_elk <- list(step = c(500, 1000, 500, 1000, 0.01, 0.005), angle = c(0.3, 0.7))  #zero-mass params needed
  Par0_m1_wtd <- list(step = c(100, 500, 100, 500, 0.01, 0.005), angle = c(0.3, 0.7))  #zero-mass params needed
  Par0_m1_coug <- list(step = c(500, 1500, 500, 1500, 0.01, 0.005), angle = c(0.3, 0.7))  #zero-mass params needed
  Par0_m1_wolf <- list(step = c(500, 3000, 500, 3000), angle = c(0.3, 0.7))  
  Par0_m1_bob <- list(step = c(500, 1000, 500, 1000), angle = c(0.3, 0.7))  
  Par0_m1_coy <- list(step = c(500, 2000, 500, 2000), angle = c(0.3, 0.7))  
  #'  Step arguments: report 2 means then the 2 SD for the two different states
  #'  Gamma distribution: mean & standard deviation of step lengths for each state
  #'  Michelot & Langrock 2019 recommend using same value for mean and SD per state
  #'  Wrapped Cauchy distribution: concentration of turning angles for each state
  #'  Include zero-mass parameters when there are 0s in the data w/gamma, Weibull, etc. distributions
  #'  e.g., zeromass0 <- c(0.1,0.05) # step zero-mass
  
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
  DM_time <- list(step = list(mean = ~cosinor(hour, period = 24), sd = ~cosinor(hour, period = 24)), angle = list(concentration = ~1))
  DM_null_Zerotime <- list(step = list(mean = ~cosinor(hour, period = 24), sd = ~cosinor(hour, period = 24), zeromass = ~cosinor(hour, period = 24)), angle = list(concentration = ~1)) # includes zeromass parameters
  
    
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
  trans_formula_smr_hab_pred <- ~TRI + PercOpen + Dist2Road + COUG_RSF + WOLF_RSF + BOB_RSF + COY_RSF + cosinor(hour, period = 24)
  trans_formula_wtr_hab_pred <- ~TRI + PercOpen + Dist2Road + SnowCover + COUG_RSF + WOLF_RSF + BOB_RSF + COY_RSF + cosinor(hour, period = 24)
  trans_formula_smr_hab_pred_noCoy <- ~TRI + PercOpen + Dist2Road + COUG_RSF + WOLF_RSF + BOB_RSF + cosinor(hour, period = 24)
  trans_formula_wtr_hab_pred_noTime <- ~TRI + PercOpen + Dist2Road + SnowCover + COUG_RSF + WOLF_RSF + BOB_RSF + COY_RSF
  #'  For predator species
  trans_formula_smr_OK <- ~TRI + PercOpen + Dist2Road + MD_RSF + cosinor(hour, period = 24)
  trans_formula_wtr_OK <- ~TRI + PercOpen + Dist2Road + SnowCover + MD_RSF + cosinor(hour, period = 24)
  trans_formula_smr_NE <- ~TRI + PercOpen + Dist2Road + ELK_RSF + WTD_RSF + cosinor(hour, period = 24)
  trans_formula_wtr_NE <- ~TRI + PercOpen + Dist2Road + SnowCover + ELK_RSF + WTD_RSF + cosinor(hour, period = 24)
  
  
  
  
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
  HMM_fit <- function(Data, dists, Par0_m1, dm, tformula) { 
    
    #' Fit basic model with no covariates
    m1 <- fitHMM(data = Data, nbStates = 2, dist = dists, Par0 = Par0_m1,
                 estAngleMean = list(angle = FALSE), stateNames = stateNames)

    #' #'  Compute the most likely state sequence
    #' states <- viterbi(m1)
    #' #'  Derive percentage of time spent in each state
    #' table(states)/nrow(Data)
    
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
  #'  Univariate models   
  md_HMM_smr_time <- HMM_fit(mdData_smr, dists_wc, Par0_m1_md, DM_null_ZeroMass, trans_formula_time)
  md_HMM_smr_TRI <- HMM_fit(mdData_smr, dists_wc, Par0_m1_md, DM_null_ZeroMass, trans_formula_TRI)
  md_HMM_smr_Open <- HMM_fit(mdData_smr, dists_wc, Par0_m1_md, DM_null_ZeroMass, trans_formula_Open)
  md_HMM_smr_Road <- HMM_fit(mdData_smr, dists_wc, Par0_m1_md, DM_null_ZeroMass, trans_formula_Road)
  md_HMM_smr_COUG <- HMM_fit(mdData_smr, dists_wc, Par0_m1_md, DM_null_ZeroMass, trans_formula_coug)
  md_HMM_smr_WOLF <- HMM_fit(mdData_smr, dists_wc, Par0_m1_md, DM_null_ZeroMass, trans_formula_wolf)
  md_HMM_smr_BOB <- HMM_fit(mdData_smr, dists_wc, Par0_m1_md, DM_null_ZeroMass, trans_formula_bob)
  md_HMM_smr_COY <- HMM_fit(mdData_smr, dists_wc, Par0_m1_md, DM_null_ZeroMass, trans_formula_coy)
  #'  COY & TRI highly correlated (-0.75) so identify which is more supported
  AIC(md_HMM_smr_TRI, md_HMM_smr_COY)
  #'  Global model 
  #'  Excluding COY_RSF from md_HMM_smr based on univariate TRI model having lower AIC than COY model
  #'  Including cosinor parameters on DM step length to help with autocorrelation  
  md_HMM_smr <- HMM_fit(mdData_smr, dists_wc, Par0_m1_md, DM_null_Zerotime, trans_formula_smr_hab_pred_noCoy)  
  #'  QQplot of residuals
  plotPR(md_HMM_smr, lag.max = NULL, ncores = 4)
  #'  Does temporal autocorrelation look any better? NOPE looks real bad on step length
  pr_md_HMM_smr <- pseudoRes(md_HMM_smr)
  acf(pr_md_HMM_smr$stepRes[is.finite(pr_md_HMM_smr$stepRes)], lag.max = 300)
  plot(md_HMM_smr, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  #'  Winter
  #'  Univariate models
  #'  SE, 95% CI on Sin effect on Pr(2 -> 1) in time model is NA
  md_HMM_wtr_time <- HMM_fit(mdData_wtr, dists_wc, Par0_m1_md, DM_null_ZeroMass, trans_formula_time)
  md_HMM_wtr_TRI <- HMM_fit(mdData_wtr, dists_wc, Par0_m1_md, DM_null_ZeroMass, trans_formula_TRI)
  md_HMM_wtr_Open <- HMM_fit(mdData_wtr, dists_wc, Par0_m1_md, DM_null_ZeroMass, trans_formula_Open)
  md_HMM_wtr_Road <- HMM_fit(mdData_wtr, dists_wc, Par0_m1_md, DM_null_ZeroMass, trans_formula_Road)
  md_HMM_wtr_Snow <- HMM_fit(mdData_wtr, dists_wc, Par0_m1_md, DM_null_ZeroMass, trans_formula_Snow)
  md_HMM_wtr_COUG <- HMM_fit(mdData_wtr, dists_wc, Par0_m1_md, DM_null_ZeroMass, trans_formula_coug)
  md_HMM_wtr_WOLF <- HMM_fit(mdData_wtr, dists_wc, Par0_m1_md, DM_null_ZeroMass, trans_formula_wolf)
  md_HMM_wtr_BOB <- HMM_fit(mdData_wtr, dists_wc, Par0_m1_md, DM_null_ZeroMass, trans_formula_bob)
  md_HMM_wtr_COY <- HMM_fit(mdData_wtr, dists_wc, Par0_m1_md, DM_null_ZeroMass, trans_formula_coy)
  #'  Global model
  #'  Including cosinor parameters on DM step length to help with autocorrelation
  #'  Removed cosinor parameters on trans prob. due to poor convergence of SE and 
  #'  95% CI on Sin effect on Pr(2 -> 1) in univariate time model
  md_HMM_wtr <- HMM_fit(mdData_wtr, dists_wc, Par0_m1_md, DM_null_Zerotime, trans_formula_wtr_hab_pred)
  #'  QQplot of residuals
  plotPR(md_HMM_wtr, lag.max = NULL, ncores = 4)
  #'  Does temporal autocorrelation look any better?
  pr_md_HMM_wtr <- pseudoRes(md_HMM_wtr)
  acf(pr_md_HMM_wtr$stepRes[is.finite(pr_md_HMM_wtr$stepRes)], lag.max = 300)
  plot(md_HMM_wtr, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  
  ####  ELK HMMS  ####
  #'  Summer
  #'  Univariate models
  elk_HMM_smr_time <- HMM_fit(elkData_smr, dists_wc, Par0_m1_elk, DM_null_ZeroMass, trans_formula_time)
  elk_HMM_smr_TRI <- HMM_fit(elkData_smr, dists_wc, Par0_m1_elk, DM_null_ZeroMass, trans_formula_TRI)
  elk_HMM_smr_Open <- HMM_fit(elkData_smr, dists_wc, Par0_m1_elk, DM_null_ZeroMass, trans_formula_Open)
  elk_HMM_smr_Road <- HMM_fit(elkData_smr, dists_wc, Par0_m1_elk, DM_null_ZeroMass, trans_formula_Road)
  elk_HMM_smr_COUG <- HMM_fit(elkData_smr, dists_wc, Par0_m1_elk, DM_null_ZeroMass, trans_formula_coug)
  elk_HMM_smr_WOLF <- HMM_fit(elkData_smr, dists_wc, Par0_m1_elk, DM_null_ZeroMass, trans_formula_wolf)
  elk_HMM_smr_BOB <- HMM_fit(elkData_smr, dists_wc, Par0_m1_elk, DM_null_ZeroMass, trans_formula_bob)
  elk_HMM_smr_COY <- HMM_fit(elkData_smr, dists_wc, Par0_m1_elk, DM_null_ZeroMass, trans_formula_coy)
  #'  Global model
  elk_HMM_smr <- HMM_fit(elkData_smr, dists_wc, Par0_m1_elk, DM_null_ZeroMass, trans_formula_smr_hab_pred)
  pr <- pseudoRes(elk_HMM_smr)
  plotPR(elk_HMM_smr, lag.max = NULL, ncores = 4)
  acf(pr$stepRes[!is.na(pr$stepRes)],lag.max = 100)
  
  #'  Winter
  #'  Univariate models
  elk_HMM_wtr_time <- HMM_fit(elkData_wtr, dists_wc, Par0_m1_elk, DM_null_ZeroMass, trans_formula_time)
  elk_HMM_wtr_TRI <- HMM_fit(elkData_wtr, dists_wc, Par0_m1_elk, DM_null_ZeroMass, trans_formula_TRI)
  elk_HMM_wtr_Open <- HMM_fit(elkData_wtr, dists_wc, Par0_m1_elk, DM_null_ZeroMass, trans_formula_Open)
  elk_HMM_wtr_Road <- HMM_fit(elkData_wtr, dists_wc, Par0_m1_elk, DM_null_ZeroMass, trans_formula_Road)
  elk_HMM_wtr_Snow <- HMM_fit(elkData_wtr, dists_wc, Par0_m1_elk, DM_null_ZeroMass, trans_formula_Snow)
  elk_HMM_wtr_COUG <- HMM_fit(elkData_wtr, dists_wc, Par0_m1_elk, DM_null_ZeroMass, trans_formula_coug)
  elk_HMM_wtr_WOLF <- HMM_fit(elkData_wtr, dists_wc, Par0_m1_elk, DM_null_ZeroMass, trans_formula_wolf)
  elk_HMM_wtr_BOB <- HMM_fit(elkData_wtr, dists_wc, Par0_m1_elk, DM_null_ZeroMass, trans_formula_bob)
  elk_HMM_wtr_COY <- HMM_fit(elkData_wtr, dists_wc, Par0_m1_elk, DM_null_ZeroMass, trans_formula_coy)
  #'  Global model
  elk_HMM_wtr <- HMM_fit(elkData_wtr, dists_wc, Par0_m1_elk, DM_null_ZeroMass, trans_formula_wtr_hab_pred)
  
  
  ####  WHITE-TAILED DEER HMMS  ####
  #'  Summer
  #'  Univariate models
  wtd_HMM_smr_time <- HMM_fit(wtdData_smr, dists_wc, Par0_m1_wtd, DM_null_ZeroMass, trans_formula_time)
  wtd_HMM_smr_TRI <- HMM_fit(wtdData_smr, dists_wc, Par0_m1_wtd, DM_null_ZeroMass, trans_formula_TRI)
  wtd_HMM_smr_Open <- HMM_fit(wtdData_smr, dists_wc, Par0_m1_wtd, DM_null_ZeroMass, trans_formula_Open)
  wtd_HMM_smr_Road <- HMM_fit(wtdData_smr, dists_wc, Par0_m1_wtd, DM_null_ZeroMass, trans_formula_Road)
  wtd_HMM_smr_COUG <- HMM_fit(wtdData_smr, dists_wc, Par0_m1_wtd, DM_null_ZeroMass, trans_formula_coug)
  wtd_HMM_smr_WOLF <- HMM_fit(wtdData_smr, dists_wc, Par0_m1_wtd, DM_null_ZeroMass, trans_formula_wolf)
  wtd_HMM_smr_BOB <- HMM_fit(wtdData_smr, dists_wc, Par0_m1_wtd, DM_null_ZeroMass, trans_formula_bob)
  wtd_HMM_smr_COY <- HMM_fit(wtdData_smr, dists_wc, Par0_m1_wtd, DM_null_ZeroMass, trans_formula_coy)
  #'  Global model
  wtd_HMM_smr <- HMM_fit(wtdData_smr, dists_wc, Par0_m1_wtd, DM_null_ZeroMass, trans_formula_smr_hab_pred)
  
  #'  Winter
  #'  Univariate models
  wtd_HMM_wtr_time <- HMM_fit(wtdData_wtr, dists_wc, Par0_m1_wtd, DM_null_ZeroMass, trans_formula_time)
  wtd_HMM_wtr_TRI <- HMM_fit(wtdData_wtr, dists_wc, Par0_m1_wtd, DM_null_ZeroMass, trans_formula_TRI)
  wtd_HMM_wtr_Open <- HMM_fit(wtdData_wtr, dists_wc, Par0_m1_wtd, DM_null_ZeroMass, trans_formula_Open)
  wtd_HMM_wtr_Road <- HMM_fit(wtdData_wtr, dists_wc, Par0_m1_wtd, DM_null_ZeroMass, trans_formula_Road)
  wtd_HMM_wtr_Snow <- HMM_fit(wtdData_wtr, dists_wc, Par0_m1_wtd, DM_null_ZeroMass, trans_formula_Snow)
  wtd_HMM_wtr_COUG <- HMM_fit(wtdData_wtr, dists_wc, Par0_m1_wtd, DM_null_ZeroMass, trans_formula_coug)
  wtd_HMM_wtr_WOLF <- HMM_fit(wtdData_wtr, dists_wc, Par0_m1_wtd, DM_null_ZeroMass, trans_formula_wolf)
  wtd_HMM_wtr_BOB <- HMM_fit(wtdData_wtr, dists_wc, Par0_m1_wtd, DM_null_ZeroMass, trans_formula_bob)
  wtd_HMM_wtr_COY <- HMM_fit(wtdData_wtr, dists_wc, Par0_m1_wtd, DM_null_ZeroMass, trans_formula_coy)
  #'  Global model
  wtd_HMM_wtr <- HMM_fit(wtdData_wtr, dists_wc, Par0_m1_wtd, DM_null_ZeroMass, trans_formula_wtr_hab_pred)
  
  
  ####  COUGAR HMMS  ####       
  #'  Okanogan Summer
  #'  Univariate models
  coug_HMM_smr_OK_time <- HMM_fit(cougData_smr_OK, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_time)
  coug_HMM_smr_OK_TRI <- HMM_fit(cougData_smr_OK, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_TRI)
  coug_HMM_smr_OK_Open <- HMM_fit(cougData_smr_OK, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_Open)
  coug_HMM_smr_OK_Road <- HMM_fit(cougData_smr_OK, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_Road)
  coug_HMM_smr_OK_MD <- HMM_fit(cougData_smr_OK, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_md)
  #'  Global model
  coug_HMM_smr_OK <- HMM_fit(cougData_smr_OK, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_smr_OK)
  
  #'  Okanogan Winter
  #'  Univariate models
  coug_HMM_wtr_OK_time <- HMM_fit(cougData_wtr_OK, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_time)
  coug_HMM_wtr_OK_TRI <- HMM_fit(cougData_wtr_OK, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_TRI)
  coug_HMM_wtr_OK_Open <- HMM_fit(cougData_wtr_OK, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_Open)
  coug_HMM_wtr_OK_Road <- HMM_fit(cougData_wtr_OK, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_Road)
  coug_HMM_wtr_OK_Snow <- HMM_fit(cougData_wtr_OK, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_Snow)
  coug_HMM_wtr_OK_MD <- HMM_fit(cougData_wtr_OK, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_md)
  #'  MD & TRI highly correlated (0.76) so identify which is more supported
  AIC(coug_HMM_wtr_OK_TRI, coug_HMM_wtr_OK_Open, coug_HMM_wtr_OK_MD)
  #'  Global model   
  #'  Exluding XXX from coug_HMM_wtr_OK model based on univariate XXX model having lower AIC than XXX model
  coug_HMM_wtr_OK <- HMM_fit(cougData_wtr_OK, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_wtr_OK)
  
  #'  Northeast Summer
  #'  Univariate models
  coug_HMM_smr_NE_time <- HMM_fit(cougData_smr_NE, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_time)
  coug_HMM_smr_NE_TRI <- HMM_fit(cougData_smr_NE, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_TRI)
  coug_HMM_smr_NE_Open <- HMM_fit(cougData_smr_NE, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_Open)
  coug_HMM_smr_NE_Road <- HMM_fit(cougData_smr_NE, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_Road)
  coug_HMM_smr_NE_ELK <- HMM_fit(cougData_smr_NE, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_elk)
  coug_HMM_smr_NE_WTD <- HMM_fit(cougData_smr_NE, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_wtd)
  #'  Global model
  coug_HMM_smr_NE <- HMM_fit(cougData_smr_NE, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_smr_NE)
  #
  #'  Northeast Winter
  #'  Univariate models
  coug_HMM_wtr_time <- HMM_fit(cougData_wtr_NE, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_time)
  coug_HMM_wtr_TRI <- HMM_fit(cougData_wtr_NE, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_TRI)
  coug_HMM_wtr_Open <- HMM_fit(cougData_wtr_NE, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_Open)
  coug_HMM_wtr_Road <- HMM_fit(cougData_wtr_NE, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_Road)
  coug_HMM_wtr_Snow <- HMM_fit(cougData_wtr_NE, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_Snow)
  coug_HMM_wtr_ELK <- HMM_fit(cougData_wtr_NE, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_elk)
  coug_HMM_wtr_WTD <- HMM_fit(cougData_wtr_NE, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_wtd)
  #'  Global model
  coug_HMM_wtr_NE <- HMM_fit(cougData_wtr_NE, dists_wc, Par0_m1_coug, DM_null_ZeroMass, trans_formula_wtr_NE)
  
  
  ####  WOLF HMMS  ####
  #'  Okanogan Summer
  #'  Univariate models
  wolf_HMM_smr_OK_time <- HMM_fit(wolfData_smr_OK, dists_wc, Par0_m1_wolf, DM_null, trans_formula_time)
  wolf_HMM_smr_OK_TRI <- HMM_fit(wolfData_smr_OK, dists_wc, Par0_m1_wolf, DM_null, trans_formula_TRI)
  wolf_HMM_smr_OK_Open <- HMM_fit(wolfData_smr_OK, dists_wc, Par0_m1_wolf, DM_null, trans_formula_Open)
  wolf_HMM_smr_OK_Road <- HMM_fit(wolfData_smr_OK, dists_wc, Par0_m1_wolf, DM_null, trans_formula_Road)
  wolf_HMM_smr_OK_MD <- HMM_fit(wolfData_smr_OK, dists_wc, Par0_m1_wolf, DM_null, trans_formula_md)
  #'  Global model
  wolf_HMM_smr_OK <- HMM_fit(wolfData_smr_OK, dists_wc, Par0_m1_wolf, DM_null, trans_formula_smr_OK)
  
  #'  Okanogan Winter
  #'  Univariate models
  wolf_HMM_wtr_OK_time <- HMM_fit(wolfData_wtr_OK, dists_wc, Par0_m1_wolf, DM_null, trans_formula_time)
  wolf_HMM_wtr_OK_TRI <- HMM_fit(wolfData_wtr_OK, dists_wc, Par0_m1_wolf, DM_null, trans_formula_TRI)
  wolf_HMM_wtr_OK_Open <- HMM_fit(wolfData_wtr_OK, dists_wc, Par0_m1_wolf, DM_null, trans_formula_Open)
  wolf_HMM_wtr_OK_Road <- HMM_fit(wolfData_wtr_OK, dists_wc, Par0_m1_wolf, DM_null, trans_formula_Road)
  wolf_HMM_wtr_OK_Snow <- HMM_fit(wolfData_wtr_OK, dists_wc, Par0_m1_wolf, DM_null, trans_formula_Snow)
  wolf_HMM_wtr_OK_MD <- HMM_fit(wolfData_wtr_OK, dists_wc, Par0_m1_wolf, DM_null, trans_formula_md)
  #'  Global model
  wolf_HMM_wtr_OK <- HMM_fit(wolfData_wtr_OK, dists_wc, Par0_m1_wolf, DM_null, trans_formula_wtr_OK)
  
  #'  Northeast Summer
  #'  Univariate models
  wolf_HMM_smr_NE_time <- HMM_fit(wolfData_smr_NE, dists_wc, Par0_m1_wolf, DM_null, trans_formula_time)
  wolf_HMM_smr_NE_TRI <- HMM_fit(wolfData_smr_NE, dists_wc, Par0_m1_wolf, DM_null, trans_formula_TRI)
  wolf_HMM_smr_NE_Open <- HMM_fit(wolfData_smr_NE, dists_wc, Par0_m1_wolf, DM_null, trans_formula_Open)
  wolf_HMM_smr_NE_Road <- HMM_fit(wolfData_smr_NE, dists_wc, Par0_m1_wolf, DM_null, trans_formula_Road)
  wolf_HMM_smr_NE_ELK <- HMM_fit(wolfData_smr_NE, dists_wc, Par0_m1_wolf, DM_null, trans_formula_elk)
  wolf_HMM_smr_NE_WTD <- HMM_fit(wolfData_smr_NE, dists_wc, Par0_m1_wolf, DM_null, trans_formula_wtd)
  #'  Global model
  wolf_HMM_smr_NE <- HMM_fit(wolfData_smr_NE, dists_wc, Par0_m1_wolf, DM_null, trans_formula_smr_NE)
  
  #'  Northeast Winter
  #'  Univariate models
  wolf_HMM_wtr_NE_time <- HMM_fit(wolfData_wtr_NE, dists_wc, Par0_m1_wolf, DM_null, trans_formula_time)
  wolf_HMM_wtr_NE_TRI <- HMM_fit(wolfData_wtr_NE, dists_wc, Par0_m1_wolf, DM_null, trans_formula_TRI)
  wolf_HMM_wtr_NE_Open <- HMM_fit(wolfData_wtr_NE, dists_wc, Par0_m1_wolf, DM_null, trans_formula_Open)
  wolf_HMM_wtr_NE_Road <- HMM_fit(wolfData_wtr_NE, dists_wc, Par0_m1_wolf, DM_null, trans_formula_Road)
  wolf_HMM_wtr_NE_Snow <- HMM_fit(wolfData_wtr_NE, dists_wc, Par0_m1_wolf, DM_null, trans_formula_Snow)
  wolf_HMM_wtr_NE_ELK <- HMM_fit(wolfData_wtr_NE, dists_wc, Par0_m1_wolf, DM_null, trans_formula_elk)
  wolf_HMM_wtr_NE_WTD <- HMM_fit(wolfData_wtr_NE, dists_wc, Par0_m1_wolf, DM_null, trans_formula_wtd)
  #'  Global model
  wolf_HMM_wtr_NE <- HMM_fit(wolfData_wtr_NE, dists_wc, Par0_m1_wolf, DM_null, trans_formula_wtr_NE)
  
  
  ####  BOBCAT HMMS  ####
  #'  Okanogan Summer
  #'  Univariate models
  bob_HMM_smr_OK_time <- HMM_fit(bobData_smr_OK, dists_wc, Par0_m1_bob, DM_null, trans_formula_time)
  bob_HMM_smr_OK_TRI <- HMM_fit(bobData_smr_OK, dists_wc, Par0_m1_bob, DM_null, trans_formula_TRI)
  bob_HMM_smr_OK_Open <- HMM_fit(bobData_smr_OK, dists_wc, Par0_m1_bob, DM_null, trans_formula_Open)
  bob_HMM_smr_OK_Road <- HMM_fit(bobData_smr_OK, dists_wc, Par0_m1_bob, DM_null, trans_formula_Road)
  bob_HMM_smr_OK_MD <- HMM_fit(bobData_smr_OK, dists_wc, Par0_m1_bob, DM_null, trans_formula_md)
  #'  Global model
  bob_HMM_smr_OK <- HMM_fit(bobData_smr_OK, dists_wc, Par0_m1_bob, DM_null, trans_formula_smr_OK)
  
  #'  Okanogan Winter
  #'  Univariate models
  bob_HMM_wtr_OK_time <- HMM_fit(bobData_wtr_OK, dists_wc, Par0_m1_bob, DM_null, trans_formula_time)
  bob_HMM_wtr_OK_TRI <- HMM_fit(bobData_wtr_OK, dists_wc, Par0_m1_bob, DM_null, trans_formula_TRI)
  bob_HMM_wtr_OK_Open <- HMM_fit(bobData_wtr_OK, dists_wc, Par0_m1_bob, DM_null, trans_formula_Open)
  bob_HMM_wtr_OK_Road <- HMM_fit(bobData_wtr_OK, dists_wc, Par0_m1_bob, DM_null, trans_formula_Road)
  bob_HMM_wtr_OK_Snow <- HMM_fit(bobData_wtr_OK, dists_wc, Par0_m1_bob, DM_null, trans_formula_Snow)
  bob_HMM_wtr_OK_MD <- HMM_fit(bobData_wtr_OK, dists_wc, Par0_m1_bob, DM_null, trans_formula_md)
  #'  Global model
  bob_HMM_wtr_OK <- HMM_fit(bobData_wtr_OK, dists_wc, Par0_m1_bob, DM_null, trans_formula_wtr_OK)
  
  #'  Northeast Summer
  #'  Univariate models
  bob_HMM_smr_NE_time <- HMM_fit(bobData_smr_NE, dists_wc, Par0_m1_bob, DM_null, trans_formula_time)
  bob_HMM_smr_NE_TRI <- HMM_fit(bobData_smr_NE, dists_wc, Par0_m1_bob, DM_null, trans_formula_TRI)
  bob_HMM_smr_NE_Open <- HMM_fit(bobData_smr_NE, dists_wc, Par0_m1_bob, DM_null, trans_formula_Open)
  bob_HMM_smr_NE_Road <- HMM_fit(bobData_smr_NE, dists_wc, Par0_m1_bob, DM_null, trans_formula_Road)
  bob_HMM_smr_NE_ELK <- HMM_fit(bobData_smr_NE, dists_wc, Par0_m1_bob, DM_null, trans_formula_elk)
  bob_HMM_smr_NE_WTD <- HMM_fit(bobData_smr_NE, dists_wc, Par0_m1_bob, DM_null, trans_formula_wtd)
  #'  Global model
  bob_HMM_smr_NE <- HMM_fit(bobData_smr_NE, dists_wc, Par0_m1_bob, DM_null, trans_formula_smr_NE)
  
  #'  Northeast Winter
  #'  Univariate models
  bob_HMM_wtr_NE_time <- HMM_fit(bobData_wtr_NE, dists_wc, Par0_m1_bob, DM_null, trans_formula_time)
  bob_HMM_wtr_NE_TRI <- HMM_fit(bobData_wtr_NE, dists_wc, Par0_m1_bob, DM_null, trans_formula_TRI)
  bob_HMM_wtr_NE_Open <- HMM_fit(bobData_wtr_NE, dists_wc, Par0_m1_bob, DM_null, trans_formula_Open)
  bob_HMM_wtr_NE_Road <- HMM_fit(bobData_wtr_NE, dists_wc, Par0_m1_bob, DM_null, trans_formula_Road)
  bob_HMM_wtr_NE_Snow <- HMM_fit(bobData_wtr_NE, dists_wc, Par0_m1_bob, DM_null, trans_formula_Snow)
  bob_HMM_wtr_NE_ELK <- HMM_fit(bobData_wtr_NE, dists_wc, Par0_m1_bob, DM_null, trans_formula_elk)
  bob_HMM_wtr_NE_WTD <- HMM_fit(bobData_wtr_NE, dists_wc, Par0_m1_bob, DM_null, trans_formula_wtd)
  #'  Global model
  bob_HMM_wtr_NE <- HMM_fit(bobData_wtr_NE, dists_wc, Par0_m1_bob, DM_null, trans_formula_wtr_NE)
  
  
  ####  COYOTE HMMS  ####   
  #'  Okanogan Summer
  #'  Univariate models
  coy_HMM_smr_OK_time <- HMM_fit(coyData_smr_OK, dists_wc, Par0_m1_coy, DM_null, trans_formula_time)
  coy_HMM_smr_OK_TRI <- HMM_fit(coyData_smr_OK, dists_wc, Par0_m1_coy, DM_null, trans_formula_TRI)
  coy_HMM_smr_OK_Open <- HMM_fit(coyData_smr_OK, dists_wc, Par0_m1_coy, DM_null, trans_formula_Open)
  coy_HMM_smr_OK_Road <- HMM_fit(coyData_smr_OK, dists_wc, Par0_m1_coy, DM_null, trans_formula_Road)
  coy_HMM_smr_OK_MD <- HMM_fit(coyData_smr_OK, dists_wc, Par0_m1_coy, DM_null, trans_formula_md)
  #'  Global model
  coy_HMM_smr_OK <- HMM_fit(coyData_smr_OK, dists_wc, Par0_m1_coy, DM_null, trans_formula_smr_OK)
  
  #'  Okanogan Winter
  #'  Univariate models
  coy_HMM_wtr_OK_time <- HMM_fit(coyData_wtr_OK, dists_wc, Par0_m1_coy, DM_null, trans_formula_time)
  coy_HMM_wtr_OK_TRI <- HMM_fit(coyData_wtr_OK, dists_wc, Par0_m1_coy, DM_null, trans_formula_TRI)
  coy_HMM_wtr_OK_Open <- HMM_fit(coyData_wtr_OK, dists_wc, Par0_m1_coy, DM_null, trans_formula_Open)
  coy_HMM_wtr_OK_Road <- HMM_fit(coyData_wtr_OK, dists_wc, Par0_m1_coy, DM_null, trans_formula_Road)
  coy_HMM_wtr_OK_Snow <- HMM_fit(coyData_wtr_OK, dists_wc, Par0_m1_coy, DM_null, trans_formula_Snow)
  coy_HMM_wtr_OK_MD <- HMM_fit(coyData_wtr_OK, dists_wc, Par0_m1_coy, DM_null, trans_formula_md)
  #'  MD & TRI highly correlated (0.71) so identify which is more supported
  AIC(coy_HMM_wtr_OK_TRI, coy_HMM_wtr_OK_Open, coy_HMM_wtr_OK_MD)
  #'  Global model 
  #'  Exluding XXX from coy_HMM_wtr_OK model based on univariate XXX model having lower AIC than XXX model  
  coy_HMM_wtr_OK <- HMM_fit(coyData_wtr_OK, dists_wc, Par0_m1_coy, DM_null, trans_formula_wtr_OK)
  
  #'  Northeast Summer
  #'  Univariate models
  coy_HMM_smr_NE_time <- HMM_fit(coyData_smr_NE, dists_wc, Par0_m1_coy, DM_null, trans_formula_time)
  coy_HMM_smr_NE_TRI <- HMM_fit(coyData_smr_NE, dists_wc, Par0_m1_coy, DM_null, trans_formula_TRI)
  coy_HMM_smr_NE_Open <- HMM_fit(coyData_smr_NE, dists_wc, Par0_m1_coy, DM_null, trans_formula_Open)
  coy_HMM_smr_NE_Road <- HMM_fit(coyData_smr_NE, dists_wc, Par0_m1_coy, DM_null, trans_formula_Road)
  coy_HMM_smr_NE_ELK <- HMM_fit(coyData_smr_NE, dists_wc, Par0_m1_coy, DM_null, trans_formula_elk)
  coy_HMM_smr_NE_WTD <- HMM_fit(coyData_smr_NE, dists_wc, Par0_m1_coy, DM_null, trans_formula_wtd)
  #'  Global model
  coy_HMM_smr_NE <- HMM_fit(coyData_smr_NE, dists_wc, Par0_m1_coy, DM_null, trans_formula_smr_NE)
  
  #'  Northeast Winter
  #'  Univariate models
  coy_HMM_wtr_NE_time <- HMM_fit(coyData_wtr_NE, dists_wc, Par0_m1_coy, DM_null, trans_formula_time)
  coy_HMM_wtr_NE_TRI <- HMM_fit(coyData_wtr_NE, dists_wc, Par0_m1_coy, DM_null, trans_formula_TRI)
  coy_HMM_wtr_NE_Open <- HMM_fit(coyData_wtr_NE, dists_wc, Par0_m1_coy, DM_null, trans_formula_Open)
  coy_HMM_wtr_NE_Road <- HMM_fit(coyData_wtr_NE, dists_wc, Par0_m1_coy, DM_null, trans_formula_Road)
  coy_HMM_wtr_NE_Snow <- HMM_fit(coyData_wtr_NE, dists_wc, Par0_m1_coy, DM_null, trans_formula_Snow)
  coy_HMM_wtr_NE_ELK <- HMM_fit(coyData_wtr_NE, dists_wc, Par0_m1_coy, DM_null, trans_formula_elk)
  coy_HMM_wtr_NE_WTD <- HMM_fit(coyData_wtr_NE, dists_wc, Par0_m1_coy, DM_null, trans_formula_wtd)
  #'  Global model
  coy_HMM_wtr_NE <- HMM_fit(coyData_wtr_NE, dists_wc, Par0_m1_coy, DM_null, trans_formula_wtr_NE)
  
  
  
  
  
  
  
  
  
  
  
  
  #'  Save model results
  spp_HMM_output <- list(md_HMM_smr, md_HMM_wtr, elk_HMM_smr, elk_HMM_wtr, wtd_HMM_smr, 
                         wtd_HMM_wtr, coug_HMM_smr_OK, coug_HMM_wtr_OK, coug_HMM_smr_NE, 
                         coug_HMM_wtr_NE, wolf_HMM_smr_OK, wolf_HMM_wtr_OK, 
                         wolf_HMM_smr_NE, wolf_HMM_wtr_NE, bob_HMM_smr_OK, 
                         bob_HMM_wtr_OK, bob_HMM_smr_NE, bob_HMM_wtr_NE, 
                         coy_HMM_smr_OK, coy_HMM_wtr_OK, coy_HMM_smr_NE, coy_HMM_wtr_NE)
  #'  Make sure to note whether covariates were included on transition probabilities
  save(spp_HMM_output, file = paste0("./Outputs/spp_HMM_output_", Sys.Date(), ".RData"))
  
  load("./Outputs/spp_HMM_output_2021-05-18.RData")

  
 
  
  #' #'  Function to extract most likely state sequence for all locations based on
  #' #'  the Viterbi algorithm and the fitted HMM
  #' #'  MAKE SURE YOU KNOW WHICH TRANSITION FORMULA INFORMED THESE CLASSIFICATIONS
  #' loc_states <- function(mod, locs) {
  #'   #'  Decode most likely state for each observation
  #'   states <- viterbi(mod)
  #'   #'  Append state classification to location data
  #'   dat <- cbind(locs, states)
  #'   
  #'   #'  Print derived percentage of time spent in each state
  #'   print(table(states)/nrow(locs))
  #'   
  #'   return(dat)
  #' }
  #' #'  Run fitted HMMs and location data through function
  #' md_state_smr <- loc_states(md_HMM_smr, mdData_smr)
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
  


  #'  Function to extract stationary state probabilities & plot predicted responses
  stay_probs <- function(hmmm) {
    stay_pr <- stationary(hmmm)
    stay_pr <- stay_pr[[1]]
    stay_mu0 <- stationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0, 
                                                   SnowCover = 0, TRI = 0, MD_RSF = 0, 
                                                   ELK_RSF = 0, WTD_RSF = 0, 
                                                   COUG_RSF = 0, WOLF_RSF = 0,
                                                   BOB_RSF = 0, COY_RSF = 0))
    print(stay_mu0)
    
    #'  Plot stationary state probabilities and extract predicted estimates
    fig <- plotStationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0, 
                                                  SnowCover = 0, TRI = 0, MD_RSF = 0, 
                                                  ELK_RSF = 0, WTD_RSF = 0, 
                                                  COUG_RSF = 0, WOLF_RSF = 0,
                                                  BOB_RSF = 0, COY_RSF = 0),
                          col = c("red", "blue"), plotCI = TRUE, alpha = 0.95, return =  TRUE)
    stationary_probs <- list(stay_pr, fig)
    
    return(stationary_probs)
  }
  #'  Extract stationary state probabilities
  stay_md_smr <- stay_probs(md_HMM_smr)
  stay_md_wtr <- stay_probs(md_HMM_wtr)
  stay_elk_smr <- stay_probs(elk_HMM_smr)
  stay_elk_wtr <- stay_probs(elk_HMM_wtr)
  stay_wtd_smr <- stay_probs(wtd_HMM_smr)
  stay_wtd_wtr <- stay_probs(wtd_HMM_wtr)
  stay_coug_smr_OK <- stay_probs(coug_HMM_smr_OK)
  stay_coug_wtr_OK <- stay_probs(coug_HMM_wtr_OK)
  stay_coug_smr_NE <- stay_probs(coug_HMM_smr_NE)
  stay_coug_wtr_NE <- stay_probs(coug_HMM_wtr_NE)
  stay_wolf_smr_OK <- stay_probs(wolf_HMM_smr_OK)
  stay_wolf_wtr_OK <- stay_probs(wolf_HMM_wtr_OK)
  stay_wolf_smr_NE <- stay_probs(wolf_HMM_smr_NE)
  stay_wolf_wtr_NE <- stay_probs(wolf_HMM_wtr_NE)
  stay_bob_smr_OK <- stay_probs(bob_HMM_smr_OK)
  stay_bob_wtr_OK <- stay_probs(bob_HMM_wtr_OK)
  stay_bob_smr_NE <- stay_probs(bob_HMM_smr_NE)
  stay_bob_wtr_NE <- stay_probs(bob_HMM_wtr_NE)
  stay_coy_smr_OK <- stay_probs(coy_HMM_smr_OK)
  stay_coy_wtr_OK <- stay_probs(coy_HMM_wtr_OK)
  stay_coy_smr_NE <- stay_probs(coy_HMM_smr_NE)
  stay_coy_wtr_NE <- stay_probs(coy_HMM_wtr_NE)
  
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
  
  #'  Function to extract covariate effects on transition-probabilities
  rounddig <- 2
  hmm_out <- function(mod, spp, season) {
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
        Transition = rep("Trans.1->2", nrow(.))
      ) %>%
      relocate(Parameter, .before = beta1.2) %>%
      relocate(Species, .before = Parameter) %>%
      relocate(Season, .before = Parameter) %>%
      relocate(Transition, .before = Parameter)
    colnames(out1.2) <- c("Species", "Season", "Transition", "Parameter", "Estimate", "SE", "Lower", "Upper")
    out2.1 <- as.data.frame(cbind(beta2.1, se2.1, lci2.1, uci2.1)) %>%
      mutate(
        Parameter = row.names(est_out[[3]]$est),
        Species = rep(spp, nrow(.)),
        Season = rep(season, nrow(.)),
        Transition = rep("Trans.2->1", nrow(.))
      ) %>%
      relocate(Parameter, .before = beta2.1) %>%
      relocate(Species, .before = Parameter) %>%
      relocate(Season, .before = Parameter) %>%
      relocate(Transition, .before = Parameter)
    colnames(out2.1) <- c("Species", "Season", "Transition", "Parameter", "Estimate", "SE", "Lower", "Upper")
    out <- as.data.frame(rbind(out1.2, out2.1))
    return(out)
  }
  #'  Run each season and species-specific model through function
  md_s1819_hmm <- hmm_out(md_HMM_smr, "Mule Deer", "Summer")
  md_w1820_hmm <- hmm_out(md_HMM_wtr, "Mule Deer", "Winter")
  elk_s1819_hmm <- hmm_out(elk_HMM_smr, "Elk", "Summer")
  elk_w1820_hmm <- hmm_out(elk_HMM_wtr, "Elk", "Winter")
  wtd_s1819_hmm <- hmm_out(wtd_HMM_smr, "White-tailed Deer", "Summer")
  wtd_w1820_hmm <- hmm_out(wtd_HMM_wtr, "White-tailed Deer", "Winter")
  coug_s1819_hmm <- hmm_out(coug_HMM_smr, "Cougar", "Summer")
  coug_w1820_hmm <- hmm_out(coug_HMM_wtr, "Cougar", "Winter")
  wolf_s1819_hmm <- hmm_out(wolf_HMM_smr, "Wolf", "Summer")
  wolf_w1820_hmm <- hmm_out(wolf_HMM_wtr, "Wolf", "Winter")
  bob_s1819_hmm <- hmm_out(bob_HMM_smr, "Bobcat", "Summer")
  bob_w1820_hmm <- hmm_out(bob_HMM_wtr, "Bobcat", "Winter")
  coy_s1819_hmm <- hmm_out(coy_HMM_smr, "Coyote", "Summer")
  coy_w1820_hmm <- hmm_out(coy_HMM_wtr, "Coyote", "Winter")
  
  #'  Gather prey and predator results to put into a single results table
  results_hmm_trpr <- rbind(bob_s1819_hmm, bob_w1820_hmm, coug_s1819_hmm, coug_w1820_hmm, 
                              coy_s1819_hmm, coy_w1820_hmm, md_s1819_hmm, md_w1820_hmm, 
                              elk_s1819_hmm, elk_w1820_hmm, wtd_s1819_hmm, wtd_w1820_hmm,
                              wolf_s1819_hmm, wolf_w1820_hmm)
  results_hmm_trpr_prey <- rbind(md_s1819_hmm, md_w1820_hmm, elk_s1819_hmm, elk_w1820_hmm, 
                            wtd_s1819_hmm, wtd_w1820_hmm)
  results_hmm_trpr_pred <- rbind(bob_s1819_hmm, bob_w1820_hmm, coug_s1819_hmm, coug_w1820_hmm, 
                            coy_s1819_hmm, coy_w1820_hmm, wolf_s1819_hmm, wolf_w1820_hmm)
  
  #'  Spread results so the coefficient effects are easier to compare between 
  #'  transition probabilities and across species
  #'  Ungulates (no study area covariate included)
  results_hmm_wide_trpr <- results_hmm_trpr %>%  #results_hmm_prey
    # dplyr::select(-z) %>%
    mutate(
      # SE = round(SE, 2),
      SE = paste0("(", SE, ")"),
    ) %>%
    #'  Bold significant variables- doesn't work if continue manipulating data frame
    # condformat(.) %>%
    # rule_text_bold(c(Estimate, SE, Pval), expression = Pval <= 0.05) %>%
    unite(Est_SE, Estimate, SE, sep = " ") %>%
    unite(CI95, Lower, Upper, sep = ", ") %>%
    unite(Est_SE_CI, Est_SE, CI95, sep = "_") %>%
    spread(Parameter, Est_SE_CI) %>%
    separate("(Intercept)", c("Intercept (SE)", "Intercept 95% CI"), sep = "_") %>%
    separate("Elev", c("Elev (SE)", "Elev 95% CI"), sep = "_") %>%
    separate("Slope", c("Slope (SE)", "Slope 95% CI"), sep = "_") %>%
    separate("PercForMix", c("PercForMix (SE)", "PercForMix 95% CI"), sep = "_") %>%
    separate("PercXGrass", c("PercXGrass (SE)", "PercXGrass 95% CI"), sep = "_") %>%
    separate("PercXShrub", c("PercXShrub (SE)", "PercXShrub 95% CI"), sep = "_") %>%
    separate("NearestRd", c("NearestRd (SE)", "NearestRd 95% CI"), sep = "_") %>%
    separate("HumanMod", c("HumanMod (SE)", "HumanMod 95% CI"), sep = "_") %>%
    # mutate(
    #   AreaOK = rep("NA", nrow(.)),
    #   AreaCI = rep("NA", nrow(.))
    # ) %>%
    # relocate(AreaOK, .before = "Elev (SE)") %>%
    # relocate(AreaCI, .before = "Elev (SE)") %>%
    arrange(match(Species, c("Mule Deer", "Elk", "White-tailed Deer")))
  # names(results_hmm_wide)[names(results_hmm_wide) == "AreaOK"] <- "AreaOK (SE)"
  # names(results_hmm_wide)[names(results_hmm_wide) == "AreaCI"] <- "AreaOK 95% CI"

  #'  Predators (study area covariate included)
  results_hmm_wide_pred <- results_hmm_pred %>% 
    # dplyr::select(-z) %>%
    mutate(
      # SE = round(SE, 2),
      SE = paste0("(", SE, ")"),
    ) %>%
    #'  Bold significant variables- doesn't work if continue manipulating data frame
    # condformat(.) %>%
    # rule_text_bold(c(Estimate, SE, Pval), expression = Pval <= 0.05) %>%
    unite(Est_SE, Estimate, SE, sep = " ") %>%
    unite(CI95, Lower, Upper, sep = ", ") %>%
    unite(Est_SE_CI, Est_SE, CI95, sep = "_") %>%
    spread(Parameter, Est_SE_CI) %>%
    separate("(Intercept)", c("Intercept (SE)", "Intercept 95% CI"), sep = "_") %>%
    separate("AreaOK", c("AreaOK (SE)", "AreaOK 95% CI"), sep = "_") %>%
    separate("Elev", c("Elev (SE)", "Elev 95% CI"), sep = "_") %>%
    separate("Slope", c("Slope (SE)", "Slope 95% CI"), sep = "_") %>%
    separate("PercForMix", c("PercForMix (SE)", "PercForMix 95% CI"), sep = "_") %>%
    separate("PercXGrass", c("PercXGrass (SE)", "PercXGrass 95% CI"), sep = "_") %>%
    separate("PercXShrub", c("PercXShrub (SE)", "PercXShrub 95% CI"), sep = "_") %>%
    separate("NearestRd", c("NearestRd (SE)", "NearestRd 95% CI"), sep = "_") %>%
    separate("HumanMod", c("HumanMod (SE)", "HumanMod 95% CI"), sep = "_") %>%
    arrange(match(Species, c("Bobcat", "Cougar", "Coyote", "Wolf"))) 
  # arrange(match(Season, c("Summer", "Winter")))
  
  results_hmm_wide <- rbind(results_hmm_wide_pred, results_hmm_wide_prey)
  
  write.csv(results_hmm_wide_trpr, paste0("./Outputs/HMM_Results_TrnsPrb_wide", Sys.Date(), ".csv"))
  
  #'  THIS ISN'T WORKING- CAN'T FIGURE OUT HOW TO EXTRACT BETA EFFECTS OF STUDY AREA & SEX ON STEP LENGTH
  #'  Function to extract 95%CI for covariate effects on state-dependent distributions
  #'  For species with ZeroMass functions (ungulates and cougars)
  #' state_dist_zm <- function(hmmm) {
  #'   #'  Calculate 95%CI on beta coefficients
  #'   global_est <- CIbeta(hmmm, alpha = 0.95)
  #'   #'  Grab names of each parameter
  #'   rnames <- colnames(global_est[[1]][[1]])
  #'   #'  Organize into a dataframe
  #'   step_out <- as.data.frame(matrix(unlist(global_est[[1]]), ncol = 4))
  #'   colnames(step_out) <- c("est", "se", "lower", "upper")
  #'   row.names(step_out) <- rnames
  #'   vrbls <- c("Intercept", "AreaOK", "SexM")  
  #'   params <- rep(vrbls, 4)
  #'   distp <- c("Mean1", "Mean2", "SD1", "SD2", "ZeroMass1", "ZeroMass2")
  #'   dist_params <- rep(distp, each = 3) 
  #'   step_out <- cbind(step_out, params)
  #'   step_out <- cbind(step_out, dist_params)
  #'   step_out <- step_out %>%
  #'     mutate(
  #'       params = as.factor(as.character(params)),
  #'       dist_params = as.factor(as.character(dist_params))) %>%
  #'     arrange(match(params, vrbls))
  #'   
  #'   return(step_out)
  #' }
  #' #'  Run each season and species-specific model through function
  #' md_s1819_hmm <- state_dist_zm(md_HMM_smr)

  #' #' Plot estimates and CIs for Pr(exploratory) at each time step
  #' plot(trProbs$est[1,2,], type="l", ylim=c(0,1), ylab="Pr(exploratory)", xlab="t", col=c("#E69F00", "#56B4E9")[coy_HMM_smr$miSum$Par$states])
  #' arrows(1:dim(trProbs$est)[3],
  #'        trProbs$lower[1,2,],
  #'        1:dim(trProbs$est)[3],
  #'        trProbs$upper[1,2,],
  #'        length=0.025, angle=90, code=3, col=c("#E69F00", "#56B4E9")[coy_HMM_smr$miSum$Par$states], lwd=1.3)
  #' abline(h=0.5,lty=2)
  #' 
  #' # proportion of entire time series spent in each state
  #' coy_HMM_smr$miSum$Par$timeInStates
  #' 
  #' # histograms of distance to water by state
  #' par(mfrow=c(2,1))
  #' hist(coy_HMM_smr$miSum$data$dist2sabie[which(coy_HMM_smr$miSum$Par$states==1)],main=stateNames[1],xlab="distance to water (m)")
  #' hist(coy_HMM_smr$miSum$data$dist2sabie[which(coy_HMM_smr$miSum$Par$states==2)],main=stateNames[2],xlab="distance to water (m)")
  
  
  
  
  #' #'  Testing with one species
  #' #'  Create momentuHMMData object from crwData object and covariates
  #' #'  Missing values due to sex missing from interpolated locations
  #' coyData_wtr <- prepData(data = coyMerge_wtr,
  #'                         covNames = c("Elev", "Slope", "HumanMod", "NearstRd",
  #'                                      "PercForMix", "PercXGrass", "PercXShrub",
  #'                                      "Year", "Sex", "Area"))
  #' # mdData <- prepData(data = mdMerge, covNames = c("Elev", "Slope", "HumanMod", "Sex", "Season"))
  #' # wolfData <- prepData(data = wolfMerge, covNames = c("Elev", "Slope", "HumanMod", "dist2road")) #"Sex"
  #' # wolfData <- prepData(data = crwOut_WOLF, covNames = "Sex") #covNames = c("Elev", "Slope", "HumanMod")
  #' 
  #' #' Fit basic model with no covariates
  #' m1 <- fitHMM(data = coyData_wtr, nbStates = 2, dist = dist, Par0 = Par0_m1_coy_wtr,
  #'              estAngleMean = list(angle = FALSE), stateNames = stateNames)
  #' #'  Compute the most likely state sequence
  #' states <- viterbi(m1)
  #' #'  Derive percentage of time spent in each state
  #' table(states)/nrow(coyData_wtr)
  #' 
  #' #'  Adding complexity
  #' formula <- ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + NearstRd + HumanMod + Area + Sex  # Remove Area + Sex if using ungulate data!
  #' #'  Consider putting sex (predators only) or season on the state-dependent distributions
  #' #'  I could see season influencing step length at the very least
  #' DM <- list(step = list(mean = ~Season, sd = ~Season), angle = list(concentration = ~1)) # zeromass = ~Season
  #' 
  #' #'  Get new initial parameter values based on nested m1 model
  #' Par0_m2_coy_wtr <- getPar0(model = m1, formula = formula)  #DM = DM
  #' Par0_m2_coy_wtr$beta  # should the covariate betas be 0.00000???
  #' 
  #' #'  Fit model with all covariates on transition probability
  #' m2 <- fitHMM(data = coyData_wtr, nbStates = 2, dist = dist, Par0 = Par0_m2_coy_wtr$Par,
  #'              beta0 = Par0_m2_coy_wtr$beta, stateNames = stateNames, formula = formula, DM = DM) #
  #' states <- viterbi(m2)
  #' table(states)/nrow(coyData_wtr)
  #' 
  #' #'  Model selection with AIC
  #' AIC(m1,m2)
  #' 
  #' Par0_m3_coy_wtr <- getPar0(model = m2, formula = formula)
  #' Par0_m3_coy_wtr$beta
  #' 
  #' #'  Plot model- this will plot every individual track....
  #' # plot(m1, plotCI = TRUE)
  


  
  
  