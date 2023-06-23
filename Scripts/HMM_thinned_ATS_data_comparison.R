  #'  ===================================
  #'  Hidden Markov Movement Models - thinned accelerometer data
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  April 2023
  #'  ===================================
  #'  Script to run hidden Markov movement models on thinned accelerometer data 
  #'  as a post-hoc analysis to assess how fix interval influences results. Original
  #'  data supplied by L. Satterfield and included relocations of four cougars in
  #'  Okanogan study area from Dec. 2019 - Feb. 2020. Collars recorded a location
  #'  every 10 minutes. These data were subsampled to a relocation every 30-min,
  #'  1 hour, 2 hours, and 4 hours.
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
  load("./Outputs/Telemetry_crwOut/crwOut_ATS_wtr.RData")
  load("./Outputs/Telemetry_covs/ats_telem_covs_2023-04-03.RData")
  
  
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
    crwlMerge$crwPredict$BOB_RSF <- scale(crwlMerge$crwPredict$BOB_RSF)
    crwlMerge$crwPredict$COY_RSF <- scale(crwlMerge$crwPredict$COY_RSF)
    crwlMerge$crwPredict$hour <- as.integer(crwlMerge$crwPredict$hour)
    crwlMerge$crwPredict$hour_fix <- as.integer(crwlMerge$crwPredict$hour_fix)
    crwlMerge$crwPredict$hour3 <- as.integer(crwlMerge$crwPredict$hour3)
    #'  Prep crwlMerge data for fitHMM function
    Data <- prepData(data = crwlMerge, covNames = c("Dist2Road", "PercOpen", #"NDVI", 
                                                    "SnowCover", "TRI", "MD_RSF", 
                                                    "COUG_RSF", "WOLF_RSF", "BOB_RSF", 
                                                    "COY_RSF", "hour", "hour_fix",
                                                    "hour3", "daytime", "Sex", 
                                                    "StudyArea", "Season"))  # "ELK_RSF", "WTD_RSF", 
    return(Data)
  }
  #'  Run data from the Okanogan through prep function
  #'  Warnings are due to missing data for interpolated locations. prepData 
  #'  command automatically fills in values with closest following value.
  ats_full <- spp_dataPrep_OK(crwOut_ATS_wtr[[1]], ats_covs[[1]])
  ats_30m <- spp_dataPrep_OK(crwOut_ATS_wtr[[2]], ats_covs[[2]])
  ats_1hr <- spp_dataPrep_OK(crwOut_ATS_wtr[[3]], ats_covs[[3]])
  ats_2hr <- spp_dataPrep_OK(crwOut_ATS_wtr[[4]], ats_covs[[4]])
  ats_4hr <- spp_dataPrep_OK(crwOut_ATS_wtr[[5]], ats_covs[[5]])
  
  ats_hmm_data <- list(ats_full, ats_30m, ats_1hr, ats_2hr, ats_4hr)   
  names(ats_hmm_data) <- c("atsData_full", "atsData_30m", "atsData_1hr", "atsData_2hr", "atsData_4hr")
  # save(ats_hmm_data, file = paste0("./Outputs/Telemetry_crwOut/crwOut_ATS_wCovs_", Sys.Date(), ".RData"))
  
  load("./Outputs/Telemetry_crwOut/crwOut_ATS_wCovs_2023-04-03.RData")
  names(ats_hmm_data) <- c("atsData_full", "atsData_30m", "atsData_1hr", "atsData_2hr", "atsData_4hr")
  
  #'  Visualize data to inform initial parameter specifications
  mean(ats_hmm_data[[1]]$step, na.rm = T); sd(ats_hmm_data[[1]]$step, na.rm = T) # 36, 77
  mean(ats_hmm_data[[2]]$step, na.rm = T); sd(ats_hmm_data[[2]]$step, na.rm = T) # 83, 181
  mean(ats_hmm_data[[3]]$step, na.rm = T); sd(ats_hmm_data[[3]]$step, na.rm = T) # 146, 302
  mean(ats_hmm_data[[4]]$step, na.rm = T); sd(ats_hmm_data[[4]]$step, na.rm = T) # 260, 478
  mean(ats_hmm_data[[5]]$step, na.rm = T); sd(ats_hmm_data[[5]]$step, na.rm = T) # 454, 714
  
  #'  Visualize data to identify potential temporal autocorrelation
  #'  lag.max is measured in hours
  acf(ats_hmm_data[[1]]$step[!is.na(ats_hmm_data[[1]]$step)],lag.max=100)
  acf(ats_hmm_data[[2]]$step[!is.na(ats_hmm_data[[2]]$step)],lag.max=100)
  acf(ats_hmm_data[[3]]$step[!is.na(ats_hmm_data[[3]]$step)],lag.max=100)
  acf(ats_hmm_data[[4]]$step[!is.na(ats_hmm_data[[4]]$step)],lag.max=100)
  acf(ats_hmm_data[[5]]$step[!is.na(ats_hmm_data[[5]]$step)],lag.max=100)
  
  ####  Initial model set up  ####
  #'  ============================
  #'  Define initial parameters associated with each distribution & each state
  #'  Species-specific parameters based on viewing plotted data and mean step lengths
  #'  Providing value close to mean step length as "exploratory" mean & SD
  Par0_m1_ats_full <- list(step = c(30, 75, 30, 75, 0.01, 0.005), angle = c(0.1, 0.5))
  Par0_m1_ats_30m <- list(step = c(80, 200, 80, 200, 0.01, 0.005), angle = c(0.1, 0.5))
  Par0_m1_ats_1hr <- list(step = c(100, 300, 100, 300, 0.01, 0.005), angle = c(0.1, 0.5))
  Par0_m1_ats_2hr <- list(step = c(100, 500, 100, 500, 0.01, 0.005), angle = c(0.1, 0.5))
  Par0_m1_ats_4hr <- list(step = c(100, 650, 100, 650, 0.01, 0.005), angle = c(0.1, 0.5))
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
  
  ####  Models describing Transition Probabilities  ####
  #'  Define formula(s) to be applied to transition probabilities
  #'  Covariates affecting probability of transitioning from one state to another
  #'  and associated with behavioral states
  #'  Univariate models
  trans_formula_null <- ~1
  trans_formula_time <- ~cosinor(hour, period = 24)
  
  #'  Covariate model reflects cougar wtr OK noMD model in original analysis but
  #'  excludes snow cover b/c almost no snow = 0 observations in this data set
  trans_formula_wtr_OK_noMDnoSnow <- ~TRI + PercOpen + Dist2Road
  
  
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
  
  ####  ATS COUGAR WINTER HMMS  ####
  #'  10-min fix intervals - note no snow covariate owing to no snow = 0 observations in data set
  ats_HMM_full <- HMM_fit(ats_hmm_data[[1]], dists_vm, Par0_m1_ats_full, DM_Zerotime, trans_formula_wtr_OK_noMDnoSnow, fits = 1)
  #'  QQplot of residuals
  plotPR(ats_HMM_full, lag.max = 100, ncores = 4)
  #'  Does temporal autocorrelation look any better?
  pr_ats_HMM_full <- pseudoRes(ats_HMM_full)
  acf(pr_ats_HMM_full$stepRes[!is.na(pr_ats_HMM_full$stepRes)],lag.max = 100)
  #'  Review parameter estimates
  plot(ats_HMM_full, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  #'  30-min fix interval - note no snow covariate owing to no snow = 0 observations in data set
  ats_HMM_30m <- HMM_fit(ats_hmm_data[[2]], dists_vm, Par0_m1_ats_30m, DM_Zerotime, trans_formula_wtr_OK_noMDnoSnow, fits = 1)
  #'  QQplot of residuals
  plotPR(ats_HMM_30m, lag.max = 100, ncores = 4)
  #'  Does temporal autocorrelation look any better?
  pr_ats_HMM_30m <- pseudoRes(ats_HMM_30m)
  acf(pr_ats_HMM_30m$stepRes[!is.na(pr_ats_HMM_30m$stepRes)],lag.max = 100)
  #'  Review parameter estimates
  plot(ats_HMM_30m, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  #'  1-hr fix interval - note no snow covariate owing to no snow = 0 observations in data set
  ats_HMM_1hr <- HMM_fit(ats_hmm_data[[3]], dists_vm, Par0_m1_ats_1hr, DM_Zerotime, trans_formula_wtr_OK_noMDnoSnow, fits = 1)
  #'  QQplot of residuals
  plotPR(ats_HMM_1hr, lag.max = 100, ncores = 4)
  #'  Does temporal autocorrelation look any better?
  pr_ats_HMM_1hr <- pseudoRes(ats_HMM_1hr)
  acf(pr_ats_HMM_1hr$stepRes[!is.na(pr_ats_HMM_1hr$stepRes)],lag.max = 100)
  #'  Review parameter estimates
  plot(ats_HMM_1hr, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  #'  2-hr fix interval - note no snow covariate owing to no snow = 0 observations in data set
  ats_HMM_2hr <- HMM_fit(ats_hmm_data[[4]], dists_vm, Par0_m1_ats_2hr, DM_Zerotime, trans_formula_wtr_OK_noMDnoSnow, fits = 1)
  #'  QQplot of residuals
  plotPR(ats_HMM_2hr, lag.max = 100, ncores = 4)
  #'  Does temporal autocorrelation look any better?
  pr_ats_HMM_2hr <- pseudoRes(ats_HMM_2hr)
  acf(pr_ats_HMM_2hr$stepRes[!is.na(pr_ats_HMM_2hr$stepRes)],lag.max = 100)
  #'  Review parameter estimates
  plot(ats_HMM_2hr, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  #'  4-hr fix interval - note no snow covariate owing to no snow = 0 observations in data set
  ats_HMM_4hr <- HMM_fit(ats_hmm_data[[5]], dists_vm, Par0_m1_ats_4hr, DM_Zerotime, trans_formula_wtr_OK_noMDnoSnow, fits = 1)
  #'  QQplot of residuals
  plotPR(ats_HMM_4hr, lag.max = 100, ncores = 4)
  #'  Does temporal autocorrelation look any better?
  pr_ats_HMM_4hr <- pseudoRes(ats_HMM_4hr)
  acf(pr_ats_HMM_4hr$stepRes[!is.na(pr_ats_HMM_4hr$stepRes)],lag.max = 100)
  #'  Review parameter estimates
  plot(ats_HMM_4hr, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  #'  Save!
  ats_HMM_output <- list(ats_HMM_full, ats_HMM_30m, ats_HMM_1hr, ats_HMM_2hr, ats_HMM_4hr)
  save(ats_HMM_output, file = paste0("./Outputs/HMM_output/ats_HMM_output_", Sys.Date(), ".RData"))
  
  
  ####  Summarize Results  ####
  load("./Outputs/HMM_output/ats_HMM_output_2023-04-03.RData")
  
  #'  Review model output
  print(ats_HMM_output[[1]]) # 10-min interval
  print(ats_HMM_output[[2]]) # 30-min interval
  print(ats_HMM_output[[3]]) # 1-hr interval
  print(ats_HMM_output[[4]]) # 2-hr interval
  print(ats_HMM_output[[5]]) # 4-hr interval
  
  
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
  ats_full_params <- step_turn_parms_zmass(ats_HMM_output[[1]], spp = "Cougar ATS 10min", season = "Winter", area = "Okanogan")
  ats_30m_params <- step_turn_parms_zmass(ats_HMM_output[[2]], spp = "Cougar ATS 30min", season = "Winter", area = "Okanogan")
  ats_1hr_params <- step_turn_parms_zmass(ats_HMM_output[[3]], spp = "Cougar ATS 1hr", season = "Winter", area = "Okanogan")
  ats_2hr_params <- step_turn_parms_zmass(ats_HMM_output[[4]], spp = "Cougar ATS 2hr", season = "Winter", area = "Okanogan")
  ats_4hr_params <- step_turn_parms_zmass(ats_HMM_output[[5]], spp = "Cougar ATS 4hr", season = "Winter", area = "Okanogan")
  

  #'  Make single giant table of all step length parameters
  ats_steps <- bind_rows(ats_full_params[[1]], ats_30m_params[[1]], ats_1hr_params[[1]], 
                         ats_2hr_params[[1]], ats_4hr_params[[1]]) %>%
    mutate(Mean = round(Mean, 2),
           SD = round(SD, 2),
           Zeromass = round(Zeromass, 2)) %>%
    arrange(Species)
  
  #'  Make single giant table of all turning angles parameters
  ats_turns <- bind_rows(ats_full_params[[2]], ats_30m_params[[2]], ats_1hr_params[[2]], 
                         ats_2hr_params[[2]], ats_4hr_params[[2]]) %>%
    mutate(Encamped = round(Encamped, 2),
           Exploratory = round(Exploratory, 2)) %>%
    arrange(Species)
  
  # write.csv(ats_steps, paste0("./Outputs/HMM_output/HMM_Results_StepLength_ATS_", Sys.Date(), ".csv"))
  # write.csv(ats_turns, paste0("./Outputs/HMM_output/HMM_Results_TurningAngle_ATS_", Sys.Date(), ".csv"))
  
  
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
  ats_full_hmm <- hmm_out(ats_HMM_output[[1]], "Cougar ATS 10min", "Winter", "Okanogan")
  ats_30m_hmm <- hmm_out(ats_HMM_output[[2]], "Cougar ATS 30min", "Winter", "Okanogan")
  ats_1hr_hmm <- hmm_out(ats_HMM_output[[3]], "Cougar ATS 1hr", "Winter", "Okanogan")
  ats_2hr_hmm <- hmm_out(ats_HMM_output[[4]], "Cougar ATS 2hr", "Winter", "Okanogan")
  ats_4hr_hmm <- hmm_out(ats_HMM_output[[5]], "Cougar ATS 4hr", "Winter", "Okanogan")
  
  #'  Gather results to put into a single results table
  results_hmm_TransPr_ats <- rbind(ats_full_hmm, ats_30m_hmm, ats_1hr_hmm, 
                                   ats_2hr_hmm, ats_4hr_hmm) %>%
    unite(CI95, Lower, Upper, sep = ", ") %>%
    mutate(
      Parameter = ifelse(Parameter == "(Intercept)", "Intercept", Parameter),
      Parameter = ifelse(Parameter == "TRI", "Terrain Ruggedness", Parameter),
      Parameter = ifelse(Parameter == "PercOpen", "Percent Open", Parameter),
      Parameter = ifelse(Parameter == "Dist2Road", "Nearest Road", Parameter),
      Parameter = ifelse(Parameter == "SnowCover1", "Snow Cover (Y)", Parameter),
      Parameter = ifelse(Parameter == "MD_RSF", "Pr(Mule Deer)", Parameter),
      Parameter = ifelse(Parameter == "ELK_RSF", "Pr(Elk)", Parameter),
      Parameter = ifelse(Parameter == "WTD_RSF", "Pr(White-tailed Deer)", Parameter)
    )
  colnames(results_hmm_TransPr_ats) <- c("Species", "Season", "Study Area", 
                                         "Transition", "Parameter", "Estimate",
                                         "SE", "CI95")
  
  # write.csv(results_hmm_TransPr_ats, paste0("./Outputs/HMM_output/HMM_Results_TransPr_ATS_long", Sys.Date(), ".csv"))
  
  
  
  #'  Spread results so the coefficient effects are easier to compare between 
  #'  transition probabilities and across species
  results_hmm_wide_TransPr_ats <- results_hmm_TransPr_ats %>% 
    mutate(
      SE = paste0("(", SE, ")"),
    ) %>%
    unite(Est_SE, Estimate, SE, sep = " ") %>%
    unite(Est_SE_CI, Est_SE, CI95, sep = "_") %>%
    spread(Parameter, Est_SE_CI) %>%
    separate("Intercept", c("Intercept (SE)", "Intercept 95% CI"), sep = "_") %>%
    separate("Terrain Ruggedness", c("Terrain Ruggedness (SE)", "Terrain Ruggedness 95% CI"), sep = "_") %>%
    separate("Percent Open", c("Percent Open (SE)", "Percent Open 95% CI"), sep = "_") %>%
    separate("Nearest Road", c("Nearest Road (SE)", "Nearest Road 95% CI"), sep = "_") %>%
    group_by(Species) %>%
    arrange(match(`Study Area`, c("Okanogan", "Northeast")), .by_group = TRUE) %>%
    ungroup()
  
  # write.csv(results_hmm_wide_TransPr_ats, paste0("./Outputs/HMM_output/HMM_Results_TransPr_ATS_wide", Sys.Date(), ".csv"))
  
  
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
    
    
    print(round(ci_nat[[1]]$est, 2))
    print(round(ci_nat[[2]]$est, 2))
    print(round(ci_nat[[3]]$est, 2))
    
    table_list <- list(steps_real, turn_real, trans_real)
    
    return(table_list)
  }
  ats_full_backtrans <- backtrans_params(ats_HMM_output[[1]], spp = "Cougar ATS 10min", season = "Winter", area = "Okanogan")
  ats_30m_backtrans <- backtrans_params(ats_HMM_output[[2]], spp = "Cougar ATS 30min", season = "Winter", area = "Okanogan")
  ats_1hr_backtrans <- backtrans_params(ats_HMM_output[[3]], spp = "Cougar ATS 1hr", season = "Winter", area = "Okanogan")
  ats_2hr_backtrans <- backtrans_params(ats_HMM_output[[4]], spp = "Cougar ATS 2hr", season = "Winter", area = "Okanogan")
  ats_3hr_backtrans <- backtrans_params(ats_HMM_output[[5]], spp = "Cougar ATS 4hr", season = "Winter", area = "Okanogan")
  
  #'  Table for back-transformed step lengths
  ats_steps_backtrans <- bind_rows(ats_full_backtrans[[1]], ats_30m_backtrans[[1]], 
                                   ats_1hr_backtrans[[1]], ats_2hr_backtrans[[1]], 
                                   ats_3hr_backtrans[[1]]) %>%
    pivot_wider(names_from = "State", values_from = c("Mean", "95%CI")) %>%
    relocate('95%CI_Encamped', .after = "Mean_Encamped")
  colnames(ats_steps_backtrans) <- c("Species", "Study Area", "Season", "Mean Encamped", "95% CI Encamped", "Mean Exploratory", "95% CI Exploratory")
  
  #'  Table for back-transformed turning angles
  ats_turns_backtrans <- bind_rows(ats_full_backtrans[[2]], ats_30m_backtrans[[2]], 
                                   ats_1hr_backtrans[[2]], ats_2hr_backtrans[[2]], 
                                   ats_3hr_backtrans[[2]]) 
  colnames(ats_turns_backtrans) <- c("Species", "Study Area", "Season", "Parameter", "Encamped", "Exploratory")  
  
  #'  Table for back-transformed transition probabilities
  ats_TransPr_backtrans <- bind_rows(ats_full_backtrans[[3]], ats_30m_backtrans[[3]], 
                                     ats_1hr_backtrans[[3]], ats_2hr_backtrans[[3]], 
                                     ats_3hr_backtrans[[3]]) 
  
  # write.csv(ats_steps_backtrans, paste0("./Outputs/HMM_output/HMM_Results_StepLength_BackTrans_ATS_", Sys.Date(), ".csv"))
  # write.csv(ats_turns_backtrans, paste0("./Outputs/HMM_output/HMM_Results_TurningAngle_BackTrans_ATS_", Sys.Date(), ".csv"))
  # write.csv(ats_TransPr_backtrans, paste0("./Outputs/HMM_output/HMM_Results_TransPr_BackTrans_ATS_", Sys.Date(), ".csv"))
  
  
  ####  Plot Stationary-State Probabilities  ####
  #'  Functions to extract stationary state probabilities & plot predicted responses
  stay_probs_pred_ats <- function(hmmm) {
    stay_pr <- stationary(hmmm)
    stay_pr <- stay_pr[[1]]
    stay_mu0 <- stationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0, TRI = 0)) 
    print(stay_mu0) 
    #'  Plot stationary state probabilities and extract predicted estimates
    fig <- plotStationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0, TRI = 0),  
                          col = c("red", "blue"), plotCI = TRUE, alpha = 0.95, return =  TRUE) 
    stationary_probs <- list(stay_pr, fig)
    
    return(stationary_probs)
  }
  stay_ats_full <- stay_probs_pred_ats(ats_HMM_output[[1]])
  stay_ats_30m <- stay_probs_pred_ats(ats_HMM_output[[2]])
  stay_ats_1hr <- stay_probs_pred_ats(ats_HMM_output[[3]])
  stay_ats_2hr <- stay_probs_pred_ats(ats_HMM_output[[4]])
  stay_ats_4hr <- stay_probs_pred_ats(ats_HMM_output[[5]])
  
  
  ####  Prettier Plots for Stationary State Probabilities  ####
  #'  Function to extract stationary state probabilities and plot outputs
  stay_plots <- function(stay, fixinterval, season, spp, area) {
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
      cov[[1]]$State <- "Slower" #"Encamped"
      cov[[2]]$State <- "Faster" #"Exploratory"
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
      mutate(nms = ifelse(nms == "TRI", "Terrain ruggedness", nms),
             nms = ifelse(nms == "PercOpen", "Open habitat", nms),
             nms = ifelse(nms == "Dist2Road", "Distance to road", nms),
             nms = ifelse(nms == "SnowCover", "Snow cover", nms),
             nms = ifelse(nms == "MD_RSF", "Mule deer RSF", nms),
             nms = ifelse(nms == "ELK_RSF", "Elk RSF", nms),
             nms = ifelse(nms == "WTD_RSF", "White-tailed deer RSF", nms),
             nms = ifelse(nms == "COUG_RSF", "Cougar RSF", nms),
             nms = ifelse(nms == "WOLF_RSF", "Wolf RSF", nms),
             nms = ifelse(nms == "BOB_RSF", "Bobcat RSF", nms),
             nms = ifelse(nms == "COY_RSF", "Coyote RSF", nms))
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
        ylab(paste(fixinterval, "fix \nStationary state probability")) +
        labs(#title = paste("Fix interval:", fixinterval), #title = paste(area, season, spp, "Stationary State Probabilities"), 
          fill = "Movement State", color = "Movement State") 
      # theme(legend.position="bottom")
      #'  Review figure
      plot(stay_plot)
      #'  Append figures
      stay_figs[[l]] <- stay_plot
    }
    
    return(stay_figs)
  }
  #'  Run each species through- for loops should allow the different coefficients
  #'  to still plot nicely
  ats_full_fig <- stay_plots(stay_ats_full, fixinterval = "10-minute", season = "Winter", spp = "Cougar ATS 10-min", area = "Okanogan")
  ats_30m_fig <- stay_plots(stay_ats_30m, fixinterval = "30-minute", season = "Winter", spp = "Cougar ATS 30-min", area = "Okanogan")
  ats_1hr_fig <- stay_plots(stay_ats_1hr, fixinterval = "1-hour", season = "Winter", spp = "Cougar ATS 1-hr", area = "Okanogan")
  ats_2hr_fig <- stay_plots(stay_ats_2hr, fixinterval = "2-hour", season = "Winter", spp = "Cougar ATS 2-hr", area = "Okanogan")
  ats_4hr_fig <- stay_plots(stay_ats_4hr, fixinterval = "4-hour", season = "Winter", spp = "Cougar ATS 4-hr", area = "Okanogan")
  
  
  #'  Patchwork figures together in panels
  library(patchwork)
  length(ats_full_fig)
  (ats_full_patch <- ats_full_fig[[1]] + theme(axis.title.x = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))) + 
      ats_full_fig[[2]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank()) + 
      ats_full_fig[[3]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank()) + 
      # plot_layout(guides = 'collect') + #& theme(legend.position = 'bottom') + 
      plot_annotation(#title = 'Winter Cougar Stationary State Probabilities',
                      subtitle = 'Fix interval: 10-min',
                      theme = theme(plot.subtitle = element_text(size = 18))) +
      plot_layout(guides = 'collect') & theme(legend.position = "none") & theme(text = element_text(size = 16)))
  length(ats_30m_fig)
  (ats_30m_patch <- ats_30m_fig[[1]] + theme(axis.title.x = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))) + 
      ats_30m_fig[[2]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank()) + 
      ats_30m_fig[[3]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank()) + 
      # plot_layout(guides = 'collect') + 
      plot_annotation(#title = 'Winter Cougar Stationary State Probabilities',
                      subtitle = 'Fix interval: 30-min',
                      theme = theme(plot.subtitle = element_text(size = 18))) +
      plot_layout(guides = 'collect') & theme(legend.position = "none") & theme(text = element_text(size = 16)))
  length(ats_1hr_fig)
  (ats_1hr_patch <- ats_1hr_fig[[1]] + theme(axis.title.x = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))) + 
      ats_1hr_fig[[2]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank()) + 
      ats_1hr_fig[[3]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank()) + 
      # plot_layout(guides = 'collect') + 
      plot_annotation(#title = 'Winter Cougar Stationary State Probabilities',
                      subtitle = 'Fix interval: 1-hour',
                      theme = theme(plot.subtitle = element_text(size = 18))) +
      plot_layout(guides = 'collect') & theme(legend.position = "none") & theme(text = element_text(size = 16)))
  length(ats_2hr_fig)
  (ats_2hr_patch <- ats_2hr_fig[[1]] + theme(axis.title.x = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))) + 
      ats_2hr_fig[[2]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank()) + 
      ats_2hr_fig[[3]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank()) + 
      # plot_layout(guides = 'collect') + 
      plot_annotation(#title = 'Winter Cougar Stationary State Probabilities',
                      subtitle = 'Fix interval: 2-hour',
                      theme = theme(plot.subtitle = element_text(size = 18))) +
      plot_layout(guides = 'collect') & theme(legend.position = "none") & theme(text = element_text(size = 16)))
  length(ats_4hr_fig)
  (ats_4hr_patch <- ats_4hr_fig[[1]] + theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))) +
      ats_4hr_fig[[2]] + theme(axis.title.y = element_blank()) + 
      ats_4hr_fig[[3]] + theme(axis.title.y = element_blank()) + 
      # plot_layout(guides = 'collect') + 
      plot_annotation(#title = 'Winter Cougar Stationary State Probabilities',
                      subtitle = 'Fix interval: 4-hour',
                      theme = theme(plot.subtitle = element_text(size = 18))) +
      plot_layout(guides = 'collect') & theme(legend.position = 'bottom') & theme(text = element_text(size = 16)))
  
  (ATS_comparision <- ats_full_patch / ats_30m_patch / ats_1hr_patch / ats_2hr_patch / ats_4hr_patch
    + plot_layout(guides = 'collect') & theme(text = element_text(size = 16))) 
  
  png(file="./Outputs/HMM_Output/ATS_fixrate_comparison.png", width = 750, height = 1200)
  (ATS_comparision <- ats_full_patch / ats_30m_patch / ats_1hr_patch / ats_2hr_patch / ats_4hr_patch
    + plot_annotation(title = 'Effect of fix interval on predicted stationary state probabilities', 
                      # subtitle = '     Relcations thinned from 10-min to 30-min, 1-hr, 2-hr, & 4-hr intervals',
                      theme = theme(plot.title = element_text(size = 18), plot.subtitle = element_text(size = 18)))
    + plot_layout(guides = 'collect') & theme(text = element_text(size = 16)))
  dev.off()

  
  
