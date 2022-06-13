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
  library(patchwork)
  
  #devtools::install_github("thomasp85/scico")
  # library(scico)
  
  #'  Load raw data with standardized covariates
  load("./Outputs/Telemetry_crwOut/crwOut_ALL_wCovs_2022-05-23.RData")

  #'  Load HMM results
  load("./Outputs/HMM_output/spp_HMM_output_2022-05-27.RData")
  
  
  ####  Stationary State Probs with MEAN Covariate Values  ####
  #'  -----------------------------------------------------
  #'  Functions to extract stationary state probabilities & plot predicted responses
  stay_probs_prey <- function(hmmm, snow) {
    #'  Calculate stationary state probs. for each state based on covariate data
    #'  for each time step
    stay_pr <- stationary(hmmm)
    stay_pr <- stay_pr[[1]]
    #'  Calculate stationary state probs. for each state when covariate data are
    #'  held at their mean value (0 b/c data are centered and scaled)
    stay_mu0 <- stationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                   SnowCover = snow, TRI = 0,
                                                   COUG_RSF = 0, WOLF_RSF = 0,
                                                   BOB_RSF = 0, COY_RSF = 0))
    print(stay_mu0)
    #'  Plot stationary state probabilities and extract predicted estimates
    fig <- plotStationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                  SnowCover = snow, TRI = 0,
                                                  COUG_RSF = 0, WOLF_RSF = 0,
                                                  BOB_RSF = 0, COY_RSF = 0),
                          col = c("red", "blue"), plotCI = TRUE, alpha = 0.95, return =  TRUE)
    stationary_probs <- list(stay_pr, fig)
    
    return(stationary_probs)
  }
  #'  Extract stationary state probabilities for deer and elk
  stay_md_smr <- stay_probs_prey(spp_HMM_output[[1]], snow = 0)
  stay_md_wtr <- stay_probs_prey(spp_HMM_output[[2]], snow = 1)
  stay_elk_smr <- stay_probs_prey(spp_HMM_output[[3]], snow = 0)
  stay_elk_wtr <- stay_probs_prey(spp_HMM_output[[4]], snow = 1)
  stay_wtd_smr <- stay_probs_prey(spp_HMM_output[[5]], snow = 0)
  stay_wtd_wtr <- stay_probs_prey(spp_HMM_output[[6]], snow = 1)
  
  #'  Stationary probabilities for predators in the Okanogan
  stay_probs_pred_OK <- function(hmmm, snow) {
    stay_pr <- stationary(hmmm)
    stay_pr <- stay_pr[[1]]
    stay_mu0 <- stationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                   SnowCover = snow, TRI = 0, MD_RSF = 0)) 
    print(stay_mu0) 
    #'  Plot stationary state probabilities and extract predicted estimates
    fig <- plotStationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                  SnowCover = snow, TRI = 0, MD_RSF = 0),  
                          col = c("red", "blue"), plotCI = TRUE, alpha = 0.95, return =  TRUE) 
    stationary_probs <- list(stay_pr, fig)
    
    return(stationary_probs)
  }
  #'  Extract stationary state probabilities
  stay_coug_smr_OK <- stay_probs_pred_OK(spp_HMM_output[[7]], snow = 0)
  stay_coug_wtr_OK <- stay_probs_pred_OK(spp_HMM_output[[8]], snow = 1)
  stay_wolf_smr_OK <- stay_probs_pred_OK(spp_HMM_output[[11]], snow = 0)
  stay_wolf_wtr_OK <- stay_probs_pred_OK(spp_HMM_output[[12]], snow = 1)
  stay_bob_smr_OK <- stay_probs_pred_OK(spp_HMM_output[[15]], snow = 0)
  # stay_bob_wtr_OK <- stay_probs_pred_OK(spp_HMM_output[[16]], snow = 1)
  stay_coy_smr_OK <- stay_probs_pred_OK(spp_HMM_output[[17]], snow = 0)
  stay_coy_wtr_OK <- stay_probs_pred_OK(spp_HMM_output[[18]], snow = 1)
  
  #'  Stationary state probabilities for predators in the Northeast
  stay_probs_pred_NE <- function(hmmm, snow) {
    stay_pr <- stationary(hmmm)
    stay_pr <- stay_pr[[1]]
    stay_mu0 <- stationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                   SnowCover = snow, TRI = 0, 
                                                   ELK_RSF = 0, WTD_RSF = 0))
    print(stay_mu0) 
    #'  Plot stationary state probabilities and extract predicted estimates
    fig <- plotStationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                  SnowCover = snow, TRI = 0, 
                                                  ELK_RSF = 0, WTD_RSF = 0),    
                          col = c("red", "blue"), plotCI = TRUE, alpha = 0.95, return =  TRUE)
    stationary_probs <- list(stay_pr, fig)
    
    return(stationary_probs)
  }
  #'  Extract stationary state probabilities
  stay_coug_smr_NE <- stay_probs_pred_NE(spp_HMM_output[[9]], snow = 0)
  stay_coug_wtr_NE <- stay_probs_pred_NE(spp_HMM_output[[10]], snow = 1)
  stay_wolf_smr_NE <- stay_probs_pred_NE(spp_HMM_output[[13]], snow = 0)
  stay_wolf_wtr_NE <- stay_probs_pred_NE(spp_HMM_output[[14]], snow = 1)
  # stay_bob_smr_NE <- stay_probs_pred_NE(spp_HMM_output[[17]], snow = 0)
  stay_bob_wtr_NE <- stay_probs_pred_NE(spp_HMM_output[[16]], snow = 1)
  stay_coy_smr_NE <- stay_probs_pred_NE(spp_HMM_output[[19]], snow = 0)
  stay_coy_wtr_NE <- stay_probs_pred_NE(spp_HMM_output[[20]], snow = 1)
  
  
  #' #'  Functions to extract stationary state probabilities & plot predicted responses
  #' stay_probs_prey <- function(hmmm) {
  #'   #'  Calculate stationary state probs. for each state based on covariate data
  #'   #'  for each time step
  #'   stay_pr <- stationary(hmmm)
  #'   stay_pr <- stay_pr[[1]]
  #'   #'  Calculate stationary state probs. for each state when covariate data are
  #'   #'  held at their mean value (0 b/c data are centered and scaled)
  #'   stay_mu0 <- stationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
  #'                                                  SnowCover = 0, TRI = 0,
  #'                                                  COUG_RSF = 0, WOLF_RSF = 0,
  #'                                                  BOB_RSF = 0, COY_RSF = 0))
  #'   print(stay_mu0)
  #'   #'  Plot stationary state probabilities and extract predicted estimates
  #'   fig <- plotStationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
  #'                                                 SnowCover = 0, TRI = 0,
  #'                                                 COUG_RSF = 0, WOLF_RSF = 0,
  #'                                                 BOB_RSF = 0, COY_RSF = 0),
  #'                         col = c("red", "blue"), plotCI = TRUE, alpha = 0.95, return =  TRUE)
  #'   stationary_probs <- list(stay_pr, fig)
  #'   
  #'   return(stationary_probs)
  #' }
  #' #'  Extract stationary state probabilities for deer and elk
  #' stay_md_smr <- stay_probs_prey(spp_HMM_output[[1]])
  #' stay_md_wtr <- stay_probs_prey(spp_HMM_output[[2]])
  #' stay_elk_smr <- stay_probs_prey(spp_HMM_output[[3]])
  #' stay_elk_wtr <- stay_probs_prey(spp_HMM_output[[4]])
  #' stay_wtd_smr <- stay_probs_prey(spp_HMM_output[[5]])
  #' stay_wtd_wtr <- stay_probs_prey(spp_HMM_output[[6]])
  #' 
  #' #'  Stationary probabilities for predators in the Okanogan
  #' stay_probs_pred_OK <- function(hmmm) {
  #'   stay_pr <- stationary(hmmm)
  #'   stay_pr <- stay_pr[[1]]
  #'   stay_mu0 <- stationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
  #'                                                  SnowCover = 0, TRI = 0, MD_RSF = 0)) 
  #'   print(stay_mu0) 
  #'   #'  Plot stationary state probabilities and extract predicted estimates
  #'   fig <- plotStationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
  #'                                                 SnowCover = 0, TRI = 0, MD_RSF = 0),  
  #'                         col = c("red", "blue"), plotCI = TRUE, alpha = 0.95, return =  TRUE) 
  #'   stationary_probs <- list(stay_pr, fig)
  #'   
  #'   return(stationary_probs)
  #' }
  #' #'  Extract stationary state probabilities
  #' stay_coug_smr_OK <- stay_probs_pred_OK(spp_HMM_output[[7]])
  #' stay_coug_wtr_OK <- stay_probs_pred_OK(spp_HMM_output[[8]])
  #' stay_wolf_smr_OK <- stay_probs_pred_OK(spp_HMM_output[[11]])
  #' stay_wolf_wtr_OK <- stay_probs_pred_OK(spp_HMM_output[[12]])
  #' stay_bob_smr_OK <- stay_probs_pred_OK(spp_HMM_output[[15]])
  #' # stay_bob_wtr_OK <- stay_probs_pred_OK(spp_HMM_output[[16]])
  #' stay_coy_smr_OK <- stay_probs_pred_OK(spp_HMM_output[[17]])
  #' stay_coy_wtr_OK <- stay_probs_pred_OK(spp_HMM_output[[18]])
  #' 
  #' #'  Stationary state probabilities for predators in the Northeast
  #' stay_probs_pred_NE <- function(hmmm) {
  #'   stay_pr <- stationary(hmmm)
  #'   stay_pr <- stay_pr[[1]]
  #'   stay_mu0 <- stationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
  #'                                                  SnowCover = 0, TRI = 0, 
  #'                                                  ELK_RSF = 0, WTD_RSF = 0))
  #'   print(stay_mu0) 
  #'   #'  Plot stationary state probabilities and extract predicted estimates
  #'   fig <- plotStationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
  #'                                                 SnowCover = 0, TRI = 0, 
  #'                                                 ELK_RSF = 0, WTD_RSF = 0),    
  #'                         col = c("red", "blue"), plotCI = TRUE, alpha = 0.95, return =  TRUE)
  #'   stationary_probs <- list(stay_pr, fig)
  #'   
  #'   return(stationary_probs)
  #' }
  #' #'  Extract stationary state probabilities
  #' stay_coug_smr_NE <- stay_probs_pred_NE(spp_HMM_output[[9]])
  #' stay_coug_wtr_NE <- stay_probs_pred_NE(spp_HMM_output[[10]])
  #' stay_wolf_smr_NE <- stay_probs_pred_NE(spp_HMM_output[[13]])
  #' stay_wolf_wtr_NE <- stay_probs_pred_NE(spp_HMM_output[[14]])
  #' # stay_bob_smr_NE <- stay_probs_pred_NE(spp_HMM_output[[17]])
  #' stay_bob_wtr_NE <- stay_probs_pred_NE(spp_HMM_output[[16]])
  #' stay_coy_smr_NE <- stay_probs_pred_NE(spp_HMM_output[[19]])
  #' stay_coy_wtr_NE <- stay_probs_pred_NE(spp_HMM_output[[20]]) 
  
  #'  Function to extract stationary state probabilities for prettier plotting
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
    #theme(legend.position="bottom") +
    xlim(-2, 2.5) + ylim(0, 1.0) +
    xlab("Scaled cougar RSF value") +
    ylab("Probability of exploratory state") #+
    #labs(title = "Ungulate movement in response to relative probability of use by cougars")
  
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
    #theme(legend.position="bottom") +
    xlim(-1.5, 2) + ylim(0, 1.0) +
    xlab("Scaled wolf RSF value") +
    ylab("Probability of exploratory state") #+
    #labs(title = "Ungulate movement in response to relative wolf site use")
  
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
    #theme(legend.position="bottom") +
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
    #theme(legend.position="bottom") +
    xlim(-2, 2.5) + ylim(0, 1.0) +
    xlab("Scaled coyote RSF value") +
    ylab("Probability of exploratory state") #+
    #labs(title = "Ungulate movement in response to relative coyote site use")
  
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
    #theme(legend.position="bottom") +
    xlim(-2.5, 2) + ylim(0, 1.0) +
    xlab("Scaled mule deer RSF value") +
    ylab("Probability of exploratory state") #+
    #labs(title = "Predator movement in response to relative mule deer site use")
  
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
    #theme(legend.position="bottom") +
    xlim(-2, 2.5) + ylim(0, 1.0) +
    xlab("Scaled elk RSF value") +
    ylab("Probability of exploratory state") #+
    #labs(title = "Predator movement in response to relative elk site use")
  
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
    #theme(legend.position="bottom") +
    xlim(-2.5, 2) + ylim(0, 1.0) +
    xlab("Scaled white-tailed deer RSF value") +
    ylab("Probability of exploratory state") #+
    #labs(title = "Predator movement in response to relative white-tailed deer site use")
  
  #'  Save
  ggsave("./Outputs/Figures for ms/HMM Stationary States/prey_coug_plot.tiff", prey_coug_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/prey_wolf_plot.tiff", prey_wolf_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/prey_bob_plot.tiff", prey_bob_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/prey_coy_plot.tiff", prey_coy_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_md_plot.tiff", pred_md_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_elk_plot.tiff", pred_elk_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_wtd_plot.tiff", pred_wtd_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  
  
  #'  Patchwork figures together in panels
  prey_pred_fig <- prey_coug_plot + prey_wolf_plot + prey_coy_plot +
      plot_layout(guides = 'collect') + 
      plot_layout(ncol = 3) +
      plot_annotation(title = 'Effect of relative probability of use by predators on ungulate stationary state probabilities')
  pred_prey_fig <- pred_elk_plot + pred_md_plot + pred_wtd_plot +
    plot_layout(guides = 'collect') + 
    plot_layout(ncol = 3) +
    plot_annotation(title = 'Effect of relative probability of use by prey on predator stationary state probabilities')
  
  ggsave("./Outputs/Figures for ms/HMM Stationary States/PreyEffect_StationaryProb_plot.tiff", prey_pred_fig, width = 11, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/PredEffect_StationaryProb_plot.tiff", pred_prey_fig, width = 11, height = 7, dpi = 800, units = "in", device = 'tiff')
  
  
  
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
    #theme(legend.position="bottom") +
    xlim(-1, 4) + ylim(0, 1.0) +
    xlab("Scaled terrain ruggedness index (TRI)") +
    ylab("Probability of exploratory state") #+
    #labs(title = "Ungulate movement in response to habitat complexity")
  
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
    #theme(legend.position="bottom") +
    xlim(-1, 3) + ylim(0, 1.0) +
    xlab("Scaled distance to nearest road") +
    ylab("Probability of exploratory state") #+
    #labs(title = "Ungulate movement in response to distance to nearest road")
  
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
    #theme(legend.position="bottom") +
    xlim(-1.5, 2.5) + ylim(0, 1.0) +
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
    #theme(legend.position="bottom") +
    xlim(-2, 5) + ylim(0, 1.0) +
    xlab("Scaled terrain ruggedness index (TRI), \nOkanogan") +
    ylab("Probability of exploratory state") #+
    #labs(title = "Predator movement in response to habitat complexity, Okanogan")
  
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
    #theme(legend.position="bottom") +
    xlim(-2, 4) + ylim(0, 1.0) +
    xlab("Scaled terrain ruggedness index (TRI), \nNortheast") +
    ylab("Probability of exploratory state") #+
    #labs(title = "Predator movement in response to habitat complexity, Northeast")
  
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
    #theme(legend.position="bottom") +
    xlim(-0.5, 4.5) + ylim(0, 1.0) +
    xlab("Scaled distance to nearest road, \nOkanogan") +
    ylab("Probability of exploratory state") #+
    #labs(title = "Predator movement in response to distance to nearest road, Okanogan")
  
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
    #theme(legend.position="bottom") +
    xlim(-1, 3.5) + ylim(0, 1.0) +
    xlab("Scaled distance to nearest road, \nNortheast") +
    ylab("Probability of exploratory state") #+
    #labs(title = "Predator movement in response to distance to nearest road, Northeast")
  
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
    #theme(legend.position="bottom") +
    xlim(-1.5, 1.5) + ylim(0, 1.0) +
    xlab("Scaled percent open habitat") +
    ylab("Probability of exploratory state, \nOkanogan") +
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
    #theme(legend.position="bottom") +
    xlim(-1, 3) + ylim(0, 1.0) +
    xlab("Scaled percent open habitat, \nNortheast") +
    ylab("Probability of exploratory state") +
    labs(title = "Predator movement in response to percentage of open habitat, Northeast")
  
  #'  Save
  ggsave("./Outputs/Figures for ms/HMM Stationary States/prey_tri_plot.tiff", prey_tri_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/prey_dsit2rd_plot.tiff", prey_dist2rd_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/prey_percopen_plot.tiff", prey_open_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_OK_tri_plot.tiff", pred_OK_tri_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_OK_dsit2rd_plot.tiff", pred_OK_dist2rd_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_OK_percopen_plot.tiff", pred_OK_open_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_NE_tri_plot.tiff", pred_NE_tri_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_NE_dsit2rd_plot.tiff", pred_NE_dist2rd_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_NE_percopen_plot.tiff", pred_NE_open_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  
  
  ##'  Patchwork figures together in panels
  tri_fig <- prey_tri_plot + pred_OK_tri_plot + pred_NE_tri_plot +
      plot_layout(guides = 'collect') + 
      plot_layout(ncol = 3) +
      plot_annotation(title = 'Effect of Terrain ruggedness on stationary state probabilities')
  road_fig <- prey_dist2rd_plot + pred_OK_dist2rd_plot + pred_NE_dist2rd_plot +
      plot_layout(guides = 'collect') + 
      plot_layout(ncol = 3) +
      plot_annotation(title = 'Effect of distance to nearest road on stationary state probabilities')
  
  
  ggsave("./Outputs/Figures for ms/HMM Stationary States/TRI_StationaryProb_plot.tiff", tri_fig, width = 11, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/RoadDist_StationaryProb_plot.tiff", road_fig, width = 11, height = 7, dpi = 800, units = "in", device = 'tiff')
  
  
  
  ####  Stationary State Probs over a RANGE of Covatiate Values  ####
  #'  -----------------------------------------------------------
  #'  Identify range of valued for each covariate included in species-specific HMMs
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
  stay_md_wtr <- stay_probs_prey(spp_HMM_output[[2]], md_wtr_range)
  
  
  
  
  
  
  
  
  
 
  