  #'  ============================================
  #'  Plotting Plot Stationary-State Probabilities
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  June 2022
  #'  ============================================
  
  #'  Load libraries
  library(momentuHMM)
  library(ggplot2)
  library(tidyverse)  
  library(purrr)
  #devtools::install_github("thomasp85/scico")
  library(scico)
  
  #'  Load raw data with standardized covariates
  load("./Outputs/Telemetry_crwOut/crwOut_ALL_wCovs_2022-05-23.RData")

  #'  Load HMM results
  load("./Outputs/HMM_output/spp_HMM_output_2022-05-27.RData")
  
  #'  Identify range of valued for each covariated included in species-specific HMMs
  cov_range <- function(hmm_dat) {
    #'  Pull out min and max values of each standardized covariate
    dist2rd <- range(hmm_dat$Dist2Road)
    percopen <- range(hmm_dat$PercOpen)
    #snow <- unique(hmm_dat$SnowCover)
    tri <- range(hmm_dat$TRI)
    md_rsf <- range(hmm_dat$MD_RSF)
    elk_rsf <- range(hmm_dat$ELK_RSF)
    wtd_rsf <- range(hmm_dat$WTD_RSF)
    coug_rsf <- range(hmm_dat$COUG_RSF)
    wolf_rsf <- range(hmm_dat$WOLF_RSF)
    bob_rsf <- range(hmm_dat$BOB_RSF)
    coy_rsf <- range(hmm_dat$COY_RSF)
    
    #'  Merge into a single data frame
    cov_range <- as.data.frame(rbind(dist2rd, percopen, tri, md_rsf, elk_rsf, #snow, 
                                     wtd_rsf, coug_rsf, wolf_rsf, bob_rsf, coy_rsf))
    cov_range$cov <- rownames(cov_range)
    cov_range <- cov_range %>%
      relocate(cov, .before = V1) %>%
      mutate(V1 = ifelse(cov == "snow", 0, V1),
             V2 = ifelse(cov == "snow", 1, V2),
             V1 = round(V1, 2),
             V2 = round(V2, 2))
    
    #'  Remove rows where covariate not included in model and +/-Inf introduced
    cov_range <- filter(cov_range, is.finite(V1)) 
    colnames(cov_range) <- c("Covariate", "Min_val", "Max_val")
    rownames(cov_range) <- NULL
  
    #'  Make sequence of covariate values spanning range of each covariate
    #'  Note: number of values will differ with each covariate depending on the
    #'  covariate's range 
    cov_seq <- list()
    for(i in 1:nrow(cov_range)) {
      cov_seq[[i]] <- seq(from = cov_range[i,2], to = cov_range[i,3], by = 0.1)
    }
    #'  Name list based on included covariates
    list_names <- cov_range[,1]
    names(cov_seq) <- list_names
    
    return(cov_seq)
  }
  #'  Generate sequence of covariate values for each species and season
  md_smr_range <- cov_range(hmm_data[[1]])
  md_wtr_range <- cov_range(hmm_data[[2]])
  elk_smr_range <- cov_range(hmm_data[[3]])
  elk_wtr_range <- cov_range(hmm_data[[4]])
  wtd_smr_range <- cov_range(hmm_data[[5]])
  wtd_wtr_range <- cov_range(hmm_data[[6]])
  coug_smr_OK_range <- cov_range(hmm_data[[7]])
  coug_wtr_OK_range <- cov_range(hmm_data[[8]])
  coug_smr_NE_range <- cov_range(hmm_data[[9]])
  coug_wtr_NE_range <- cov_range(hmm_data[[10]])
  wolf_smr_OK_range <- cov_range(hmm_data[[11]])
  wolf_wtr_OK_range <- cov_range(hmm_data[[12]])
  wolf_smr_NE_range <- cov_range(hmm_data[[13]])
  wolf_wtr_NE_range <- cov_range(hmm_data[[14]])
  bob_smr_OK_range <- cov_range(hmm_data[[15]])
  bob_wtr_OK_range <- cov_range(hmm_data[[16]])
  bob_smr_NE_range <- cov_range(hmm_data[[17]])
  bob_wtr_NE_range <- cov_range(hmm_data[[18]])
  coy_smr_OK_range <- cov_range(hmm_data[[19]])
  coy_wtr_OK_range <- cov_range(hmm_data[[20]])
  coy_smr_NE_range <- cov_range(hmm_data[[21]])
  coy_wtr_NE_range <- cov_range(hmm_data[[22]])
  
  stay_probs_prey <- function(hmmm, ranges) {
    #'  Calculate stationary state probs. for each state based on covariate data
    #'  for each time step
    stay_pr <- stationary(hmmm)
    stay_pr <- stay_pr[[1]]
    
    #'  Covariate effects on stationary state as distance to road changes
    stay_Dist2Rd_seq <- fig_Dist2Rd_seq <- c()
    for(i in 1:length(ranges)) {
      stay_Dist2Rd_seq[i] <- stationary(hmmm, covs = data.frame(Dist2Road = ranges[[1]][i], 
                                                                PercOpen = 0, SnowCover = 0, 
                                                                TRI = 0, COUG_RSF = 0, 
                                                                WOLF_RSF = 0, BOB_RSF = 0, COY_RSF = 0))
      fig_Dist2Rd_seq[i] <- plotStationary(hmmm, covs = data.frame(Dist2Road = ranges[[1]][i], 
                                                                   PercOpen = 0, SnowCover = 0, 
                                                                   TRI = 0, COUG_RSF = 0, 
                                                                   WOLF_RSF = 0, BOB_RSF = 0, COY_RSF = 0),
                                           col = c("red", "blue"), plotCI = TRUE, alpha = 0.95, return =  TRUE)
    }
    stationary_probs_Dist2Rd <- list(stay_Dist2Rd_seq, fig_Dist2Rd_seq)
    
    #'  Covariate effects on stationary state as TRI changes
    stay_TRI_seq <- fig_TRI_seq <- c()
    for(i in 1:length(ranges)) {
      stay_TRI_seq[i] <- stationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                            SnowCover = 0, TRI = ranges[[3]][i],
                                                            COUG_RSF = 0, WOLF_RSF = 0,
                                                            BOB_RSF = 0, COY_RSF = 0))
      fig_TRI_seq[i] <- plotStationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                               SnowCover = 0, TRI = ranges[[3]][i],
                                                               COUG_RSF = 0, WOLF_RSF = 0,
                                                               BOB_RSF = 0, COY_RSF = 0),
                                       col = c("red", "blue"), plotCI = TRUE, alpha = 0.95, return =  TRUE)
    }
    stationary_probs_TRI <- list(stay_TRI_seq, fig_TRI_seq)
    
    #'  Covariate effects on stationary state when snow cover is present (winter models only)
    stay_snow_seq <- stationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                        SnowCover = 1, TRI = 0,
                                                        COUG_RSF = 0, WOLF_RSF = 0,
                                                        BOB_RSF = 0, COY_RSF = 0))
    fig_snow_seq <- plotStationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                           SnowCover = 1, TRI = 0,
                                                           COUG_RSF = 0, WOLF_RSF = 0,
                                                           BOB_RSF = 0, COY_RSF = 0),
                                   col = c("red", "blue"), plotCI = TRUE, alpha = 0.95, return =  TRUE)
    stationary_probs_SnowCover <- list(stay_snow_seq, fig_snow_seq)
    
    stay_probs <- list(stationary_probs_Dist2Rd, stationary_probs_TRI, stationary_probs_SnowCover)
    
    return(stay_probs)
  }
  stay_md_smr <- stay_probs_prey(spp_HMM_output[[1]], md_smr_range)
  stay_md_wtr <- stay_probs_prey(spp_HMM_output[[2]])
  
  
  
  #'  Calculate stationary state probs. for each state based on covariate data
  #'  for each time step
  stay_pr <- stationary(hmmm)
  stay_pr <- stay_pr[[1]]
  
  #'  Covariate effects on stationary state as distance to road changes
  stay_Dist2Rd_seq <- fig_Dist2Rd_seq <- c()
  for(i in 1:length(md_smr_range[[1]])) {
    stay_Dist2Rd_seq[i] <- stationary(hmmm, covs = data.frame(Dist2Road = md_smr_range[[1]][i], 
                                                              PercOpen = 0, SnowCover = 0, 
                                                              TRI = 0, COUG_RSF = 0, 
                                                              WOLF_RSF = 0, BOB_RSF = 0, COY_RSF = 0))
    fig_Dist2Rd_seq[i] <- plotStationary(hmmm, covs = data.frame(Dist2Road = md_smr_range[[1]][i], 
                                                                 PercOpen = 0, SnowCover = 0, 
                                                                 TRI = 0, COUG_RSF = 0, 
                                                                 WOLF_RSF = 0, BOB_RSF = 0, COY_RSF = 0),
                          col = c("red", "blue"), plotCI = TRUE, alpha = 0.95, return =  TRUE)
  }
  stationary_probs_Dist2Rd <- list(stay_Dist2Rd_seq, fig_Dist2Rd_seq)
  
  #'  Covariate effects on stationary state as TRI changes
  stay_TRI_seq <- fig_TRI_seq <- c()
  for(i in 1:length(md_smr_range[[3]])) {
    stay_TRI_seq[i] <- stationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                         SnowCover = 0, TRI = md_smr_range[[3]][i],
                                                         COUG_RSF = 0, WOLF_RSF = 0,
                                                         BOB_RSF = 0, COY_RSF = 0))
    fig_TRI_seq[i] <- plotStationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                            SnowCover = 0, TRI = md_smr_range[[3]][i],
                                                            COUG_RSF = 0, WOLF_RSF = 0,
                                                            BOB_RSF = 0, COY_RSF = 0),
                                    col = c("red", "blue"), plotCI = TRUE, alpha = 0.95, return =  TRUE)
  }
  stationary_probs_TRI <- list(stay_TRI_seq, fig_TRI_seq)
  
  #'  Covariate effects on stationary state when snow cover is present (winter models only)
  stay_snow_seq <- stationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                      SnowCover = 1, TRI = 0,
                                                      COUG_RSF = 0, WOLF_RSF = 0,
                                                      BOB_RSF = 0, COY_RSF = 0))
  fig_snow_seq <- plotStationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                         SnowCover = 1, TRI = 0,
                                                         COUG_RSF = 0, WOLF_RSF = 0,
                                                         BOB_RSF = 0, COY_RSF = 0),
                                 col = c("red", "blue"), plotCI = TRUE, alpha = 0.95, return =  TRUE)
  stationary_probs_SnowCover <- list(stay_snow_seq, fig_snow_seq)
  
  
  
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
                                                   COUG_RSF = 0, WOLF_RSF = 0,
                                                   BOB_RSF = 0, COY_RSF = 0))
    print(stay_mu0)
    #'  Plot stationary state probabilities and extract predicted estimates
    fig <- plotStationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                  SnowCover = 0, TRI = 0,
                                                  COUG_RSF = 0, WOLF_RSF = 0,
                                                  BOB_RSF = 0, COY_RSF = 0),
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
  
  stay_covs <- function(stay, season, spp, area) {
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
      cov$Species <- spp
      cov$Season <- season
      cov$StudyArea <- area
      #'  Append to new list of data frames
      covs_out[[l]] <- cov
    }
    #'  Rename list elements based on covariate
    names(covs_out) <- names(stay_covs)
    return(covs_out)
  }
  md_smr_PrStay <- stay_covs(stay_md_smr, season = "Summer", spp = "Mule Deer", area = "Okanogan")
  md_wtr_PrStay <- stay_covs(stay_md_wtr, season = "Winter", spp = "Mule Deer", area = "Okanogan")
  elk_smr_PrStay <- stay_covs(stay_elk_smr, season = "Summer", spp = "Elk", area = "Northeast")
  elk_wtr_PrStay <- stay_covs(stay_elk_wtr, season = "Winter", spp = "Elk", area = "Northeast")
  wtd_smr_PrStay <- stay_covs(stay_wtd_smr, season = "Summer", spp = "White-tailed Deer", area = "Northeast")
  wtd_wtr_PrStay <- stay_covs(stay_wtd_wtr, season = "Winter", spp = "White-tailed Deer", area = "Northeast")
  coug_smr_OK_PrStay <- stay_covs(stay_coug_smr_OK, season = "Summer", spp = "Cougar", area = "Okanogan")
  coug_wtr_OK_PrStay <- stay_covs(stay_coug_wtr_OK, season = "Winter", spp = "Cougar", area = "Okanogan")
  coug_smr_NE_PrStay <- stay_covs(stay_coug_smr_NE, season = "Summer", spp = "Cougar", area = "Northeast")
  coug_wtr_NE_PrStay <- stay_covs(stay_coug_wtr_NE, season = "Winter", spp = "Cougar", area = "Northeast")
  wolf_smr_OK_PrStay <- stay_covs(stay_wolf_smr_OK, season = "Summer", spp = "Wolf", area = "Okanogan")
  wolf_wtr_OK_PrStay <- stay_covs(stay_wolf_wtr_OK, season = "Winter", spp = "Wolf", area = "Okanogan")
  wolf_smr_NE_PrStay <- stay_covs(stay_wolf_smr_NE, season = "Summer", spp = "Wolf", area = "Northeast")
  wolf_wtr_NE_PrStay <- stay_covs(stay_wolf_wtr_NE, season = "Winter", spp = "Wolf", area = "Northeast")
  coy_smr_OK_PrStay <- stay_covs(stay_coy_smr_OK, season = "Summer", spp = "Coyote", area = "Okanogan")
  coy_wtr_OK_PrStay <- stay_covs(stay_coy_wtr_OK, season = "Winter", spp = "Coyote", area = "Okanogan")
  coy_smr_NE_PrStay <- stay_covs(stay_coy_smr_NE, season = "Summer", spp = "Coyote", area = "Northeast")
  coy_wtr_NE_PrStay <- stay_covs(stay_coy_wtr_NE, season = "Winter", spp = "Coyote", area = "Northeast")
  
  ####  Predator-Prey Stationary State Plots  ####
  #'  ----------------------------------------
  #'  Ungulate stationary state ~ Cougar RSF
  coug_effects <- rbind(md_smr_PrStay$COUG_RSF, md_wtr_PrStay$COUG_RSF,
                       elk_smr_PrStay$COUG_RSF, elk_wtr_PrStay$COUG_RSF,
                       wtd_smr_PrStay$COUG_RSF, wtd_wtr_PrStay$COUG_RSF) %>%
    filter(!State == "Encamped") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  prey_coug_plot <- ggplot(coug_effects, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "longdash")) +
    scale_color_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) + 
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(legend.position="bottom") +
    xlim(-2, 2.5) + ylim(0, 1.0) +
    xlab("Scaled cougar RSF value") +
    ylab("Probability of exploratory state") +
    labs(title = "Ungulate movement in response to relative cougar site use")
  
  #'  Ungulate stationary state ~ Wolf RSF
  wolf_effects <- rbind(md_smr_PrStay$WOLF_RSF, md_wtr_PrStay$WOLF_RSF,
                       elk_smr_PrStay$WOLF_RSF, elk_wtr_PrStay$WOLF_RSF,
                       wtd_smr_PrStay$WOLF_RSF, wtd_wtr_PrStay$WOLF_RSF) %>%
    filter(!State == "Encamped") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  prey_wolf_plot <- ggplot(wolf_effects, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "longdash")) +
    scale_color_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) + 
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(legend.position="bottom") +
    xlim(-1.5, 2) + ylim(0, 1.0) +
    xlab("Scaled wolf RSF value") +
    ylab("Probability of exploratory state") +
    labs(title = "Ungulate movement in response to relative wolf site use")
  
  #'  Ungulate stationary state ~ Bobcat RSF
  bob_effects <- rbind(md_smr_PrStay$BOB_RSF, md_wtr_PrStay$BOB_RSF,
                       #elk_smr_PrStay$BOB_RSF, elk_wtr_PrStay$BOB_RSF,
                       wtd_smr_PrStay$BOB_RSF, wtd_wtr_PrStay$BOB_RSF) %>%
    filter(!State == "Encamped") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  prey_bob_plot <- ggplot(bob_effects, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "longdash")) +
    scale_color_manual(values=c("#E66100", "#5D3A9B")) + 
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#E66100", "#5D3A9B")) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(legend.position="bottom") +
    xlim(-2, 2.5) + ylim(0, 1.0) +
    xlab("Scaled bobcat RSF value") +
    ylab("Probability of exploratory state") +
    labs(title = "Ungulate movement in response to relative bobcat site use")
  
  #'  Ungulate stationary state ~ Coyote RSF
  coy_effects <- rbind(md_smr_PrStay$COY_RSF, md_wtr_PrStay$COY_RSF,
                       elk_smr_PrStay$COY_RSF, elk_wtr_PrStay$COY_RSF,
                       wtd_smr_PrStay$COY_RSF, wtd_wtr_PrStay$COY_RSF) %>%
    filter(!State == "Encamped") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  prey_coy_plot <- ggplot(coy_effects, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "longdash")) +
    scale_color_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) +  #"#332288", "#44AA99", "#CC6677"
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) + #"#40B0A6", "#E66100", "#5D3A9B"
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(legend.position="bottom") +
    xlim(-2, 2.5) + ylim(0, 1.0) +
    xlab("Scaled coyote RSF value") +
    ylab("Probability of exploratory state") +
    labs(title = "Ungulate movement in response to relative coyote site use")
  
  #'  Predator stationary state ~ Mule Deer RSF
  md_effects <- rbind(coug_smr_OK_PrStay$MD_RSF, coug_wtr_OK_PrStay$MD_RSF,
                       wolf_smr_OK_PrStay$MD_RSF, wolf_wtr_OK_PrStay$MD_RSF,
                       coy_smr_OK_PrStay$MD_RSF, coy_wtr_OK_PrStay$MD_RSF) %>%
    filter(!State == "Encamped") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_md_plot <- ggplot(md_effects, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "longdash")) +
    scale_color_manual(values=c("#D41159", "#FFC20A", "#0C7BDC")) +  
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#FFC20A", "#0C7BDC")) + 
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(legend.position="bottom") +
    xlim(-2.5, 2) + ylim(0, 1.0) +
    xlab("Scaled mule deer RSF value") +
    ylab("Probability of exploratory state") +
    labs(title = "Predator movement in response to relative mule deer site use")
  
  #'  Predator stationary state ~ Elk RSF
  elk_effects <- rbind(coug_smr_NE_PrStay$ELK_RSF, coug_wtr_NE_PrStay$ELK_RSF,
                       wolf_smr_NE_PrStay$ELK_RSF, wolf_wtr_NE_PrStay$ELK_RSF,
                       coy_smr_NE_PrStay$ELK_RSF, coy_wtr_NE_PrStay$ELK_RSF) %>%
    filter(!State == "Encamped") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_elk_plot <- ggplot(elk_effects, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "longdash")) +
    scale_color_manual(values=c("#D41159", "#FFC20A", "#0C7BDC")) +  
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#FFC20A", "#0C7BDC")) +     
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(legend.position="bottom") +
    xlim(-2, 2.5) + ylim(0, 1.0) +
    xlab("Scaled elk RSF value") +
    ylab("Probability of exploratory state") +
    labs(title = "Predator movement in response to relative elk site use")
  
  #'  Predator stationary state ~ White-tailed Deer RSF
  wtd_effects <- rbind(coug_smr_NE_PrStay$WTD_RSF, coug_wtr_NE_PrStay$WTD_RSF,
                       wolf_smr_NE_PrStay$WTD_RSF, wolf_wtr_NE_PrStay$WTD_RSF,
                       coy_smr_NE_PrStay$WTD_RSF, coy_wtr_NE_PrStay$WTD_RSF) %>%
    filter(!State == "Encamped") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_wtd_plot <- ggplot(wtd_effects, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "longdash")) +
    scale_color_manual(values=c("#D41159", "#FFC20A", "#0C7BDC")) +  
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#FFC20A", "#0C7BDC")) +     
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(legend.position="bottom") +
    xlim(-2.5, 2) + ylim(0, 1.0) +
    xlab("Scaled white-tailed deer RSF value") +
    ylab("Probability of exploratory state") +
    labs(title = "Predator movement in response to relative white-tailed deer site use")
  
  
  ####  Landscape Effect Stationary State Plots  ####
  #'  ----------------------------------------
  #'  Ungulate stationary state ~ TRI
  tri_effects_prey <- rbind(md_smr_PrStay$TRI, md_wtr_PrStay$TRI,
                       elk_smr_PrStay$TRI, elk_wtr_PrStay$TRI,
                       wtd_smr_PrStay$TRI, wtd_wtr_PrStay$TRI) %>%
    filter(!State == "Encamped") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  prey_tri_plot <- ggplot(tri_effects_prey, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "longdash")) +
    scale_color_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) + 
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(legend.position="bottom") +
    xlim(-1, 4) + ylim(0, 1.0) +
    xlab("Scaled terrain ruggedness index (TRI)") +
    ylab("Probability of exploratory state") +
    labs(title = "Ungulate movement in response to habitat complexity")
  
  #'  Ungulate stationary state ~ Distance to Road
  dist2rd_effects_prey <- rbind(md_smr_PrStay$Dist2Road, md_wtr_PrStay$Dist2Road,
                           elk_smr_PrStay$Dist2Road, elk_wtr_PrStay$Dist2Road,
                           wtd_smr_PrStay$Dist2Road, wtd_wtr_PrStay$Dist2Road) %>%
    filter(!State == "Encamped") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  prey_dist2rd_plot <- ggplot(dist2rd_effects_prey, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "longdash")) +
    scale_color_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) + 
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(legend.position="bottom") +
    xlim(-1, 3) + ylim(0, 1.0) +
    xlab("Scaled distance to nearest road") +
    ylab("Probability of exploratory state") +
    labs(title = "Ungulate movement in response to distance to nearest road")
  
  #'  Ungulate stationary state ~ Habitat Openness
  percopen_effects_prey <- rbind(md_smr_PrStay$PercOpen, md_wtr_PrStay$PercOpen,
                                elk_smr_PrStay$PercOpen, elk_wtr_PrStay$PercOpen,
                                wtd_smr_PrStay$PercOpen, wtd_wtr_PrStay$PercOpen) %>%
    filter(!State == "Encamped") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  prey_open_plot <- ggplot(percopen_effects_prey, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "longdash")) +
    scale_color_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) + 
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(legend.position="bottom") +
    xlim(-2, 2) + ylim(0, 1.0) +
    xlab("Scaled percent open habitat") +
    ylab("Probability of exploratory state") +
    labs(title = "Ungulate movement in response to percentage of open habitat")
  
  #'  Predator stationary state ~ TRI
  tri_effects_pred_OK <- rbind(coug_smr_OK_PrStay$TRI, coug_wtr_OK_PrStay$TRI,
                               wolf_smr_OK_PrStay$TRI, wolf_wtr_OK_PrStay$TRI,
                               coy_smr_OK_PrStay$TRI, coy_wtr_OK_PrStay$TRI) %>%
    filter(!State == "Encamped") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_OK_tri_plot <- ggplot(tri_effects_pred_OK, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "longdash")) +
    scale_color_manual(values=c("#D41159", "#FFC20A", "#0C7BDC")) +  
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#FFC20A", "#0C7BDC")) +     
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(legend.position="bottom") +
    xlim(-2, 5) + ylim(0, 1.0) +
    xlab("Scaled terrain ruggedness index (TRI)") +
    ylab("Probability of exploratory state") +
    labs(title = "Predator movement in response to habitat complexity, Okanogan")
  
  tri_effects_pred_NE <- rbind(coug_smr_NE_PrStay$TRI, coug_wtr_NE_PrStay$TRI,
                               wolf_smr_NE_PrStay$TRI, wolf_wtr_NE_PrStay$TRI,
                               coy_smr_NE_PrStay$TRI, coy_wtr_NE_PrStay$TRI) %>%
    filter(!State == "Encamped") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_NE_tri_plot <- ggplot(tri_effects_pred_NE, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "longdash")) +
    scale_color_manual(values=c("#D41159", "#FFC20A", "#0C7BDC")) +  
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#FFC20A", "#0C7BDC")) +     
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(legend.position="bottom") +
    xlim(-2, 4) + ylim(0, 1.0) +
    xlab("Scaled terrain ruggedness index (TRI)") +
    ylab("Probability of exploratory state") +
    labs(title = "Predator movement in response to habitat complexity, Northeast")
  
  #'  Predator stationary state ~ Distance to Road
  dist2rd_effects_pred_OK <- rbind(coug_smr_OK_PrStay$Dist2Road, coug_wtr_OK_PrStay$Dist2Road,
                           wolf_smr_OK_PrStay$Dist2Road, wolf_wtr_OK_PrStay$Dist2Road,
                           coy_smr_OK_PrStay$Dist2Road, coy_wtr_OK_PrStay$Dist2Road) %>%
    filter(!State == "Encamped") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_OK_dist2rd_plot <- ggplot(dist2rd_effects_pred_OK, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "longdash")) +
    scale_color_manual(values=c("#D41159", "#FFC20A", "#0C7BDC")) +  
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#FFC20A", "#0C7BDC")) +     
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(legend.position="bottom") +
    xlim(-0.5, 4.5) + ylim(0, 1.0) +
    xlab("Scaled distance to nearest road") +
    ylab("Probability of exploratory state") +
    labs(title = "Predator movement in response to distance to nearest road, Okanogan")
  
  dist2rd_effects_pred_NE <- rbind(coug_smr_NE_PrStay$Dist2Road, coug_wtr_NE_PrStay$Dist2Road,
                                   wolf_smr_NE_PrStay$Dist2Road, wolf_wtr_NE_PrStay$Dist2Road,
                                   coy_smr_NE_PrStay$Dist2Road, coy_wtr_NE_PrStay$Dist2Road) %>%
    filter(!State == "Encamped") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_NE_dist2rd_plot <- ggplot(dist2rd_effects_pred_NE, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "longdash")) +
    scale_color_manual(values=c("#D41159", "#FFC20A", "#0C7BDC")) +  
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#FFC20A", "#0C7BDC")) +     
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(legend.position="bottom") +
    xlim(-1, 3.5) + ylim(0, 1.0) +
    xlab("Scaled distance to nearest road") +
    ylab("Probability of exploratory state") +
    labs(title = "Predator movement in response to distance to nearest road, Northeast")
  
  #'  Predator stationary state ~ Open Habitat
  percopen_effects_pred_OK <- rbind(coug_smr_OK_PrStay$PercOpen, coug_wtr_OK_PrStay$PercOpen,
                                    wolf_smr_OK_PrStay$PercOpen, wolf_wtr_OK_PrStay$PercOpen,
                                    coy_smr_OK_PrStay$PercOpen, coy_wtr_OK_PrStay$PercOpen) %>%
    filter(!State == "Encamped") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_OK_open_plot <- ggplot(percopen_effects_pred_OK, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "longdash")) +
    scale_color_manual(values=c("#D41159", "#FFC20A", "#0C7BDC")) +  
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#FFC20A", "#0C7BDC")) +     
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(legend.position="bottom") +
    xlim(-1.5, 1.5) + ylim(0, 1.0) +
    xlab("Scaled percent open habitat") +
    ylab("Probability of exploratory state") +
    labs(title = "Predator movement in response to percentage of open habitat, Okanogan")
  
  percopen_effects_pred_NE <- rbind(coug_smr_NE_PrStay$PercOpen, coug_wtr_NE_PrStay$PercOpen,
                                    wolf_smr_NE_PrStay$PercOpen, wolf_wtr_NE_PrStay$PercOpen,
                                    coy_smr_NE_PrStay$PercOpen, coy_wtr_NE_PrStay$PercOpen) %>%
    filter(!State == "Encamped") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_NE_open_plot <- ggplot(percopen_effects_pred_NE, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "longdash")) +
    scale_color_manual(values=c("#D41159", "#FFC20A", "#0C7BDC")) +  
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#FFC20A", "#0C7BDC")) +     
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(legend.position="bottom") +
    xlim(-1, 3) + ylim(0, 1.0) +
    xlab("Scaled percent open habitat") +
    ylab("Probability of exploratory state") +
    labs(title = "Predator movement in response to percentage of open habitat, Northeast")
  
  #'  Save
  ggsave("./Outputs/Figures for ms/HMM Stationary States/prey_coug_plot.tiff", prey_coug_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/prey_wolf_plot.tiff", prey_wolf_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/prey_bob_plot.tiff", prey_bob_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/prey_coy_plot.tiff", prey_coy_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_md_plot.tiff", pred_md_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_elk_plot.tiff", pred_elk_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_wtd_plot.tiff", pred_wtd_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/prey_tri_plot.tiff", prey_tri_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/prey_dsit2rd_plot.tiff", prey_dist2rd_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/prey_percopen_plot.tiff", prey_open_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_OK_tri_plot.tiff", pred_OK_tri_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_OK_dsit2rd_plot.tiff", pred_OK_dist2rd_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_OK_percopen_plot.tiff", pred_OK_open_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_NE_tri_plot.tiff", pred_NE_tri_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_NE_dsit2rd_plot.tiff", pred_NE_dist2rd_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_NE_percopen_plot.tiff", pred_NE_open_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  
  
  #'  Patchwork figures together in panels
  library(patchwork)
  
  (prey_pred_fig <- prey_coug_plot + prey_wolf_plot +
    plot_layout(guides = 'collect') + 
    plot_layout(ncol = 2) +
    plot_annotation(title = 'Ungulate stationary state probabilities',
                    subtitle = 'Effect of relative probability of predator use'))
  
  
  
  
  #'  Patchwork figures together in panels
  library(patchwork)
  #'  MULE DEER panels
  length(md_smr_fig)
  (md_smr_patch <- md_smr_fig[[1]] + md_smr_fig[[2]] + theme(axis.title.y = element_blank()) + 
      md_smr_fig[[3]] + md_smr_fig[[4]] + theme(axis.title.y = element_blank()) + 
      md_smr_fig[[5]] + md_smr_fig[[6]] + theme(axis.title.y = element_blank()) + 
      plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Summer Mule Deer Stationary State Probabilities',
                      subtitle = '     Okanogan 2018 - 2021') + plot_layout(ncol = 2))
  length(md_wtr_fig)
  (md_wtr_patch <- md_wtr_fig[[1]] + md_wtr_fig[[2]] + theme(axis.title.y = element_blank()) + 
      md_wtr_fig[[3]] + theme(axis.title.y = element_blank()) + md_wtr_fig[[4]] + 
      md_wtr_fig[[5]] + theme(axis.title.y = element_blank()) + md_wtr_fig[[6]] + 
      theme(axis.title.y = element_blank()) + md_wtr_fig[[7]] + 
      md_wtr_fig[[8]] + theme(axis.title.y = element_blank()) + guide_area() + 
      plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Winter Mule Deer Stationary State Probabilities',
                      subtitle = '     Okanogan 2018 - 2021') + plot_layout(ncol = 3))
  #'  ELK panels
  length(elk_smr_fig)
  (elk_smr_patch <- elk_smr_fig[[1]] + elk_smr_fig[[2]] + theme(axis.title.y = element_blank()) + 
      elk_smr_fig[[3]] + theme(axis.title.y = element_blank()) + elk_smr_fig[[4]] + 
      elk_smr_fig[[5]] + theme(axis.title.y = element_blank()) + elk_smr_fig[[6]] + 
      theme(axis.title.y = element_blank()) + elk_smr_fig[[7]] + plot_layout(guides = 'collect') + 
      guide_area() + plot_annotation(title = 'Summer Elk Stationary State Probabilities',
                                     subtitle = '     Northeast 2018 - 2021') + plot_layout(ncol = 3))
  length(elk_wtr_fig)
  (elk_wtr_patch <- elk_wtr_fig[[1]] + elk_wtr_fig[[2]] + theme(axis.title.y = element_blank()) + 
      elk_wtr_fig[[3]] + theme(axis.title.y = element_blank()) + elk_wtr_fig[[4]] + 
      elk_wtr_fig[[5]] + theme(axis.title.y = element_blank()) + elk_wtr_fig[[6]] + 
      theme(axis.title.y = element_blank()) + elk_wtr_fig[[7]] + elk_wtr_fig[[8]] + 
      theme(axis.title.y = element_blank()) + guide_area() + plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Winter Elk Stationary State Probabilities',
                      subtitle = '    Northeast 2018 - 2021') + plot_layout(ncol = 3))
  #'  WHITE-TAILED DEER panels
  length(wtd_smr_fig)
  (wtd_smr_patch <- wtd_smr_fig[[1]] + wtd_smr_fig[[2]] + theme(axis.title.y = element_blank()) + 
      wtd_smr_fig[[3]] + theme(axis.title.y = element_blank()) + wtd_smr_fig[[4]] + 
      wtd_smr_fig[[5]] + theme(axis.title.y = element_blank()) + wtd_smr_fig[[6]] + 
      theme(axis.title.y = element_blank()) + wtd_smr_fig[[7]] + plot_layout(guides = 'collect') + 
      guide_area() + plot_annotation(title = 'Summer White-tailed Deer Stationary State Probabilities',
                                     subtitle = '     Northeast 2018 - 2021') + plot_layout(ncol = 3))
  length(wtd_wtr_fig)
  (wtd_wtr_patch <- wtd_wtr_fig[[1]] + wtd_wtr_fig[[2]] + theme(axis.title.y = element_blank()) + 
      wtd_wtr_fig[[3]] + theme(axis.title.y = element_blank()) + wtd_wtr_fig[[4]] + 
      wtd_wtr_fig[[5]] + theme(axis.title.y = element_blank()) + wtd_wtr_fig[[6]] + 
      theme(axis.title.y = element_blank()) + wtd_wtr_fig[[7]] + wtd_wtr_fig[[8]] + 
      theme(axis.title.y = element_blank()) + guide_area() + plot_layout(guides = 'collect') + 
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
  
  
  pdf(file = "./Outputs/HMM_output/Stationary_State_Prob_Plots_03.15.22.pdf")
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
  