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
  load("./Outputs/Telemetry_crwOut/crwOut_ALL_wCovs_2023-04-08.RData") #2022-05-23

  #'  Load HMM results
  load("./Outputs/HMM_output/spp_HMM_output_2023-04-08.RData") #2022-06-15
  
  
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
                                                   COUG_RSF = 0, WOLF_RSF = 0))
    print(stay_mu0)
    #'  Plot stationary state probabilities and extract predicted estimates
    fig <- plotStationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                  SnowCover = snow, TRI = 0,
                                                  COUG_RSF = 0, WOLF_RSF = 0),
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
      cov[[1]]$State <- "State 1"
      cov[[2]]$State <- "State 2"
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
  #'  Color schemes
  #'  Mule deer, elk, white-tailed deer ["#40B0A6", "#E66100", "#5D3A9B"]
  #'  Cougar, coyote, wolf ["#D41159", "#FFC20A", "#0C7BDC"]
  
  #'  Ungulate stationary state ~ Cougar RSF
  coug_effects <- rbind(md_smr_PrStay$COUG_RSF, md_wtr_PrStay$COUG_RSF,
                        elk_smr_PrStay$COUG_RSF, elk_wtr_PrStay$COUG_RSF,
                        wtd_smr_PrStay$COUG_RSF, wtd_wtr_PrStay$COUG_RSF) %>%
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  prey_coug_plot <- ggplot(coug_effects, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) + 
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(text = element_text(size = 14)) +
    #theme(legend.position="bottom") +
    xlim(-2, 2.5) + ylim(0, 1.0) +
    xlab("Scaled cougar RSF value") +
    ylab("Probability of state 2") #+
    #labs(title = "Ungulate movement in response to relative probability of use by cougars")
  
  #'  Ungulate stationary state ~ Wolf RSF
  wolf_effects <- rbind(md_smr_PrStay$WOLF_RSF, #md_wtr_PrStay$WOLF_RSF,
                        elk_smr_PrStay$WOLF_RSF, elk_wtr_PrStay$WOLF_RSF,
                        wtd_smr_PrStay$WOLF_RSF, wtd_wtr_PrStay$WOLF_RSF) %>%
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  prey_wolf_plot <- ggplot(wolf_effects, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) + 
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(axis.title.y = element_blank()) +
    theme(text = element_text(size = 14)) +
    #theme(legend.position="bottom") +
    xlim(-1.5, 2) + ylim(0, 1.0) +
    xlab("Scaled wolf RSF value") #+
    #ylab("Probability of exploratory state") #+
    #labs(title = "Ungulate movement in response to relative wolf site use")
  
  #' #'  Ungulate stationary state ~ Bobcat RSF
  #' bob_effects <- rbind(md_smr_PrStay$BOB_RSF, md_wtr_PrStay$BOB_RSF,
  #'                      wtd_smr_PrStay$BOB_RSF, wtd_wtr_PrStay$BOB_RSF) %>%
  #'   filter(!State == "Encamped") %>%
  #'   mutate(Species = as.factor(Species),
  #'          Season = as.factor(Season),
  #'          StudyArea = as.factor(StudyArea)) %>% 
  #'   dplyr::select(-c(StudyArea, State))
  #' prey_bob_plot <- ggplot(bob_effects, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
  #'   geom_line(size = 0.75) + 
  #'   scale_linetype_manual(values=c("solid", "dashed")) +
  #'   scale_color_manual(values=c("#E66100", "#5D3A9B")) + 
  #'   #'  Add confidence intervals
  #'   geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
  #'   scale_fill_manual(values=c("#E66100", "#5D3A9B")) +
  #'   #'  Get rid of lines and gray background
  #'   theme_bw() +
  #'   theme(panel.border = element_blank()) +
  #'   theme(axis.line = element_line(color = 'black')) +
  #'   theme(text = element_text(size = 14)) +
  #'   #theme(legend.position="bottom") +
  #'   xlim(-2, 2.5) + ylim(0, 1.0) +
  #'   xlab("Scaled bobcat RSF value") +
  #'   ylab("Probability of exploratory state") +
  #'   labs(title = "Ungulate stationary state probabilities in response to relative probability of use by bobcats")
  
  #' #'  Ungulate stationary state ~ Coyote RSF
  #' coy_effects <- rbind(md_wtr_PrStay$COY_RSF, #md_smr_PrStay$COY_RSF, 
  #'                      elk_wtr_PrStay$COY_RSF, #elk_smr_PrStay$COY_RSF, 
  #'                      wtd_smr_PrStay$COY_RSF, wtd_wtr_PrStay$COY_RSF) %>%
  #'   filter(!State == "Encamped") %>%
  #'   mutate(Species = as.factor(Species),
  #'          Season = as.factor(Season),
  #'          StudyArea = as.factor(StudyArea)) %>% 
  #'   dplyr::select(-c(StudyArea, State))
  #' prey_coy_plot <- ggplot(coy_effects, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
  #'   geom_line(size = 0.75) + 
  #'   scale_linetype_manual(values=c("solid", "dashed")) +
  #'   scale_color_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) +  #"#332288", "#44AA99", "#CC6677"
  #'   #'  Add confidence intervals
  #'   geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
  #'   scale_fill_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) + #"#40B0A6", "#E66100", "#5D3A9B"
  #'   #'  Get rid of lines and gray background
  #'   theme_bw() +
  #'   theme(panel.border = element_blank()) +
  #'   theme(axis.line = element_line(color = 'black')) +
  #'   theme(axis.title.y = element_blank()) +
  #'   theme(text = element_text(size = 14)) +
  #'   #theme(legend.position="bottom") +
  #'   xlim(-2, 2.5) + ylim(0, 1.0) +
  #'   xlab("Scaled coyote RSF value") #+
  #'   #ylab("Probability of exploratory state") #+
  #'   #labs(title = "Ungulate movement in response to relative coyote site use")
  
  #'  Predator stationary state ~ Mule Deer RSF
  md_effects <- rbind(coug_smr_OK_PrStay$MD_RSF, #coug_wtr_OK_PrStay$MD_RSF,
                      wolf_smr_OK_PrStay$MD_RSF) %>%  #wolf_wtr_OK_PrStay$MD_RSF,
                      #coy_smr_OK_PrStay$MD_RSF, coy_wtr_OK_PrStay$MD_RSF
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_md_plot <- ggplot(md_effects, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#D41159", "#0C7BDC")) +  #, "#FFC20A"
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#0C7BDC")) +  #, "#FFC20A"
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(axis.title.y = element_blank()) +
    theme(text = element_text(size = 14)) +
    theme(legend.position = "none") +
    xlim(-2.5, 2) + ylim(0, 1.0) +
    xlab("Scaled mule deer RSF value") +
    ylab("Probability of state 2") #+
    #labs(title = "Predator movement in response to relative mule deer site use")
  
  #'  Predator stationary state ~ Elk RSF
  elk_effects <- rbind(coug_smr_NE_PrStay$ELK_RSF, coug_wtr_NE_PrStay$ELK_RSF,
                       wolf_smr_NE_PrStay$ELK_RSF) %>% #, wolf_wtr_NE_PrStay$ELK_RSF,
                       #coy_smr_NE_PrStay$ELK_RSF, coy_wtr_NE_PrStay$ELK_RSF) %>%
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_elk_plot <- ggplot(elk_effects, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#D41159", "#0C7BDC")) +  #"#FFC20A",  
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#0C7BDC")) +    #"#FFC20A",  
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(text = element_text(size = 14)) +
    theme(legend.position = "none") +
    xlim(-2, 2.5) + ylim(0, 1.0) +
    xlab("Scaled elk RSF value") +
    ylab("Probability of state 2") #+
    #labs(title = "Predator movement in response to relative elk site use")
  
  #'  Predator stationary state ~ White-tailed Deer RSF
  wtd_effects <- rbind(coug_smr_NE_PrStay$WTD_RSF, #coug_wtr_NE_PrStay$WTD_RSF,
                       wolf_wtr_NE_PrStay$WTD_RSF) %>% #wolf_smr_NE_PrStay$WTD_RSF, 
                        #coy_wtr_NE_PrStay$WTD_RSF, coy_smr_NE_PrStay$WTD_RSF, 
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_wtd_plot <- ggplot(wtd_effects, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#D41159", "#0C7BDC")) +  #, "#FFC20A"
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#0C7BDC")) +    #, "#FFC20A"
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(text = element_text(size = 14)) +
    theme(axis.title.y = element_blank()) +
    #theme(legend.position="bottom") +
    xlim(-2.5, 2) + ylim(0, 1.0) +
    xlab("Scaled white-tailed deer RSF value") #+
    #ylab("Probability of state 2") #+
    #labs(title = "Predator movement in response to relative white-tailed deer site use")
  
  #'  Save
  ggsave("./Outputs/Figures for ms/HMM Stationary States/prey_coug_plot.tiff", prey_coug_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/prey_wolf_plot.tiff", prey_wolf_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  # ggsave("./Outputs/Figures for ms/HMM Stationary States/prey_bob_plot.tiff", prey_bob_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  # ggsave("./Outputs/Figures for ms/HMM Stationary States/prey_coy_plot.tiff", prey_coy_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_md_plot.tiff", pred_md_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_elk_plot.tiff", pred_elk_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_wtd_plot.tiff", pred_wtd_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  
  
  #'  Patchwork figures together in panels
  PredEffect_onPrey_fig <- prey_coug_plot + prey_wolf_plot + #prey_coy_plot +
      plot_layout(guides = 'collect') + 
      plot_layout(ncol = 2) +
      plot_annotation(tag_levels = 'a', 
                      title = 'Ungulate stationary state probabilities') & 
    theme(plot.tag = element_text(size = 12)) 
  PreyEffect_onPred_fig <- pred_elk_plot + pred_md_plot + pred_wtd_plot +
    plot_layout(guides = 'collect') + 
    plot_layout(ncol = 3) +
    plot_annotation(tag_levels = 'a', 
                    title = 'Predator stationary state probabilities') & 
    theme(plot.tag = element_text(size = 12)) 
  
  ggsave("./Outputs/Figures for ms/HMM Stationary States/PredEffect_onPrey_StationaryProb_plot_041223.tiff", PredEffect_onPrey_fig, width = 11, height = 7, dpi = 600, units = "in", device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/PreyEffect_onPred_StationaryProb_plot_041223.tiff", PreyEffect_onPred_fig, width = 11, height = 7, dpi = 600, units = "in", device = 'tiff', compression = 'lzw')
  
  
  PredPrey_patchwork <- PreyEffect_onPred_fig / PredEffect_onPrey_fig +
    plot_annotation(tag_levels = 'a') &
    theme(axis.title = element_text(size = 16)) &
    theme(axis.text = element_text(size = 16)) &
    theme(legend.text = element_text(size = 16)) &
    theme(legend.title = element_text(size = 16)) &
    theme(plot.tag = element_text(size = 16)) 
  ggsave("./Outputs/Figures for ms/HMM Stationary States/PredPrey_StationaryProb_plot_041223.tiff", PredPrey_patchwork, width = 13, height = 14, dpi = 600, units = "in", device = 'tiff', compression = 'lzw')
  
  
  
  ####  Landscape Effect Stationary State Plots  ####
  #'  ----------------------------------------
  #'  Color schemes
  #'  Mule deer, elk, white-tailed deer ["#40B0A6", "#E66100", "#5D3A9B"]
  #'  Cougar, coyote, wolf ["#D41159", "#FFC20A", "#0C7BDC"]
  
  #'  Ungulate stationary state ~ TRI
  tri_effects_prey <- rbind(md_smr_PrStay$TRI, md_wtr_PrStay$TRI,
                            elk_smr_PrStay$TRI, elk_wtr_PrStay$TRI,
                            wtd_smr_PrStay$TRI) %>% # wtd_wtr_PrStay$TRI
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  prey_tri_plot <- ggplot(tri_effects_prey, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) + 
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    #theme(legend.position="bottom") +
    xlim(-2, 4) + ylim(0, 1.0) +
    xlab("Scaled TRI") +
    ylab("Probability of state 2") #+
    #labs(title = "Ungulate movement in response to habitat complexity")
  
  #'  Ungulate stationary state ~ Distance to Road
  dist2rd_effects_prey <- rbind(md_smr_PrStay$Dist2Road, md_wtr_PrStay$Dist2Road,
                                elk_smr_PrStay$Dist2Road, #elk_wtr_PrStay$Dist2Road,
                                wtd_smr_PrStay$Dist2Road) %>% #, wtd_wtr_PrStay$Dist2Road
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  prey_dist2rd_plot <- ggplot(dist2rd_effects_prey, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) + 
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    #theme(legend.position="bottom") +
    xlim(-1, 4.5) + ylim(0, 1.0) +
    guides(fill = guide_legend(title = "Prey Species"), color = guide_legend(title = "Prey Species")) +
    xlab("Scaled distance to road") +
    ylab("Probability of state 2") #+
    #labs(title = "Ungulate stationary state probabilities in response to distance to nearest road")
  
  #'  Ungulate stationary state ~ Habitat Openness
  percopen_effects_prey <- rbind(md_smr_PrStay$PercOpen, md_wtr_PrStay$PercOpen,  
                                 elk_smr_PrStay$PercOpen) %>% #, elk_wtr_PrStay$PercOpen,
                                 #wtd_smr_PrStay$PercOpen, wtd_wtr_PrStay$PercOpen) %>%
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  prey_open_plot <- ggplot(percopen_effects_prey, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#40B0A6", "#E66100")) +  #, "#5D3A9B"
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#40B0A6", "#E66100")) + #, "#5D3A9B"
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(legend.position="none") +
    #theme(legend.position="bottom") +
    xlim(-1.5, 2.5) + ylim(0, 1.0) +
    guides(fill = guide_legend(title = "Prey Species"), color = guide_legend(title = "Prey Species")) +
    xlab("Scaled percent open habitat") +
    ylab("Probability of state 2") #+
    #labs(title = "Ungulate stationary state probabilities in response to percentage of open habitat")
  
  #'  Predator stationary state ~ TRI
  tri_effects_pred_OK <- rbind(coug_smr_OK_PrStay$TRI, coug_wtr_OK_PrStay$TRI,
                               wolf_smr_OK_PrStay$TRI) %>% #, wolf_wtr_OK_PrStay$TRI,
                               # coy_smr_OK_PrStay$TRI, coy_wtr_OK_PrStay$TRI) %>%
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_OK_tri_plot <- ggplot(tri_effects_pred_OK, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#D41159", "#0C7BDC")) +  #, "#FFC20A"
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#0C7BDC")) +   #, "#FFC20A"  
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(axis.title.y = element_blank()) +
    #theme(legend.position="bottom") +
    xlim(-2, 4) + ylim(0, 1.0) +
    xlab("Scaled TRI, Okanogan") +
    ylab("Probability of state 2") #+
    #labs(title = "Predator stationary state probabilities in response to habitat complexity, Okanogan")
  
  tri_effects_pred_NE <- rbind(coug_smr_NE_PrStay$TRI, coug_wtr_NE_PrStay$TRI,
                               wolf_smr_NE_PrStay$TRI, wolf_wtr_NE_PrStay$TRI) %>%
                               #coy_wtr_NE_PrStay$TRI) %>% #coy_smr_NE_PrStay$TRI, 
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_NE_tri_plot <- ggplot(tri_effects_pred_NE, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#D41159", "#0C7BDC")) +  #, "#FFC20A"
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#0C7BDC")) +   #, "#FFC20A"  
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(axis.title.y = element_blank()) +
    #theme(legend.position="bottom") +
    xlim(-2, 4) + ylim(0, 1.0) +
    xlab("Scaled TRI, Northeast") +
    ylab("Probability of state 2") #+
    #labs(title = "Predator stationary state probabilities in response to habitat complexity, Northeast")
  
  #'  Predator stationary state ~ Distance to Road
  dist2rd_effects_pred_OK <- rbind(coug_smr_OK_PrStay$Dist2Road, coug_wtr_OK_PrStay$Dist2Road,
                                   wolf_wtr_OK_PrStay$Dist2Road) %>% #wolf_smr_OK_PrStay$Dist2Road, 
                                   #coy_smr_OK_PrStay$Dist2Road, coy_wtr_OK_PrStay$Dist2Road) %>%
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_OK_dist2rd_plot <- ggplot(dist2rd_effects_pred_OK, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#D41159", "#0C7BDC")) +  #, "#FFC20A"
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#0C7BDC")) +  #, "#FFC20A"
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(axis.title.y = element_blank()) +
    theme(legend.position="none") +
    xlim(-1, 4.5) + ylim(0, 1.0) +
    guides(fill = guide_legend(title = "Predator Species"), color = guide_legend(title = "Predator Species")) +
    xlab("Scaled distance to road, OK") +
    ylab("Probability of state 2") #+
    #labs(title = "Predator stationary state probabilities in response to distance to nearest road, Okanogan")
  
  dist2rd_effects_pred_NE <- rbind(coug_smr_NE_PrStay$Dist2Road, coug_wtr_NE_PrStay$Dist2Road,
                                   wolf_smr_NE_PrStay$Dist2Road) %>% #, wolf_wtr_NE_PrStay$Dist2Road,
                                   #coy_smr_NE_PrStay$Dist2Road) %>% #coy_wtr_NE_PrStay$Dist2Road) %>%
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_NE_dist2rd_plot <- ggplot(dist2rd_effects_pred_NE, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#D41159", "#0C7BDC")) +  #, "#FFC20A" 
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#0C7BDC")) +  #, "#FFC20A"     
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(axis.title.y = element_blank()) +
    #theme(legend.position="bottom") +
    xlim(-1, 4.5) + ylim(0, 1.0) +
    guides(fill = guide_legend(title = "Predator Species"), color = guide_legend(title = "Predator Species")) +
    xlab("Scaled distance to road, NE") +
    ylab("Probability of state 2") #+
    #labs(title = "Predator stationary state probabilities in response to distance to nearest road, Northeast")
  
  #'  Predator stationary state ~ Open Habitat
  percopen_effects_pred_OK <- rbind(coug_wtr_OK_PrStay$PercOpen, #coug_smr_OK_PrStay$PercOpen, 
                                    wolf_smr_OK_PrStay$PercOpen, wolf_wtr_OK_PrStay$PercOpen) %>% #,
                                    #coy_wtr_OK_PrStay$PercOpen) %>% #coy_smr_OK_PrStay$PercOpen, 
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_OK_open_plot <- ggplot(percopen_effects_pred_OK, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#D41159", "#0C7BDC")) +  #, "#FFC20A"
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#0C7BDC")) +  #, "#FFC20A"
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(axis.title.y = element_blank()) +
    theme(legend.position="none") +
    #theme(legend.position="bottom") +
    xlim(-1.5, 1.5) + ylim(0, 1.0) +
    guides(fill = guide_legend(title = "Predator Species"), color = guide_legend(title = "Predator Species")) +
    xlab("Scaled percent open, OK") +
    ylab("Probability of state 2") #+
    #labs(title = "Predator stationary state probabilities in response to percentage of open habitat, Okanogan")
  
  percopen_effects_pred_NE <- rbind(coug_smr_NE_PrStay$PercOpen, coug_wtr_NE_PrStay$PercOpen) %>% #,
                                    #wolf_smr_NE_PrStay$PercOpen, wolf_wtr_NE_PrStay$PercOpen,
                                    #coy_wtr_NE_PrStay$PercOpen) %>% #coy_smr_NE_PrStay$PercOpen, 
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_NE_open_plot <- ggplot(percopen_effects_pred_NE, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#D41159", "#0C7BDC")) +  #, "#FFC20A"
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#0C7BDC")) +     #, "#FFC20A"
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(axis.title.y = element_blank()) +
    theme(legend.position="none") +
    xlim(-1, 3) + ylim(0, 1.0) +
    guides(fill = guide_legend(title = "Predator Species"), color = guide_legend(title = "Predator Species")) +
    xlab("Scaled percent open, NE") +
    ylab("Probability of state 2") #+
    #labs(title = "Predator stationary state probabilities in response to percentage of open habitat, Northeast")
  
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
      plot_annotation(tag_levels = 'a') &
    theme(axis.title = element_text(size = 16)) &
    theme(axis.text = element_text(size = 16)) &
    theme(legend.text = element_text(size = 16)) &
    theme(legend.title = element_text(size = 16))
  # road_fig <- prey_dist2rd_plot + pred_OK_dist2rd_plot + pred_NE_dist2rd_plot +
  #     plot_layout(guides = 'collect') + 
  #     plot_layout(ncol = 3) +
  #     plot_annotation(tag_levels = 'a',
  #                     title = 'Effect of distance to nearst road on stationary state probabilities') & 
  #   theme(axis.title = element_text(size = 16)) &
  #   theme(axis.text = element_text(size = 16)) &
  #   theme(legend.text = element_text(size = 16)) &
  #   theme(legend.title = element_text(size = 16)) &
  #   theme(plot.tag = element_text(size = 16)) 
  # open_fig <- prey_open_plot + pred_OK_open_plot + pred_NE_open_plot +
  #   plot_layout(guides = 'collect') + 
  #   plot_layout(ncol = 3) +
  #   plot_annotation(tag_levels = 'a', 
  #                   title = 'Effect of open habitat on stationary state probabilities') & 
  #   theme(axis.title = element_text(size = 16)) &
  #   theme(axis.text = element_text(size = 16)) &
  #   theme(legend.text = element_text(size = 16)) &
  #   theme(legend.title = element_text(size = 16)) &
  #   theme(plot.tag = element_text(size = 16))
  
  open_road_fig <- prey_dist2rd_plot + pred_OK_dist2rd_plot + pred_NE_dist2rd_plot +
    prey_open_plot + pred_OK_open_plot + pred_NE_open_plot +
    plot_layout(guides = 'collect') + 
    plot_layout(ncol = 3) +
    plot_annotation(tag_levels = 'a',
                    title = 'Effect of distance to nearst road and open habitat on stationary state probabilities') & 
    theme(title = element_text(size = 18)) &
    theme(axis.title = element_text(size = 16)) &
    theme(axis.text = element_text(size = 16)) &
    theme(legend.text = element_text(size = 16)) &
    theme(legend.title = element_text(size = 16)) &
    theme(plot.tag = element_text(size = 18)) 
  
  ggsave("./Outputs/Figures for ms/HMM Stationary States/TRI_StationaryProb_plot_041223.tiff", tri_fig, width = 13, height = 7, dpi = 800, units = "in", device = 'tiff', compression = 'lzw')
  # ggsave("./Outputs/Figures for ms/HMM Stationary States/RoadDist_StationaryProb_plot_041223.tiff", road_fig, width = 13, height = 7, dpi = 800, units = "in", device = 'tiff', compression = 'lzw')
  # ggsave("./Outputs/Figures for ms/HMM Stationary States/OpenHabitat_StationaryProb_plot_041223.tiff", open_fig, width = 13, height = 7, dpi = 800, units = "in", device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/OpenHabitat-RoadDist_StationaryProb_plot_041223.tiff", open_road_fig, width = 13, height = 13, dpi = 800, units = "in", device = 'tiff', compression = 'lzw')
  
  
  
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
  
  hmmm <- spp_HMM_output[[2]]
  ranges <- md_wtr_range
  snow <- 1
  stay_pr <- stationary(hmmm)
  
  
  # THIS WORKS FOR ONE VALUE OF DIST2RD AND SUMMER DATA (SNOW = 0)
  stay_tst <- plotStationary(hmmm, covs = data.frame(Dist2Road = ranges[["dist2rd"]][1], PercOpen = 0, #ranges[["dist2rd"]][1]
                                                     SnowCover = snow, TRI = 0,
                                                     COUG_RSF = 0, WOLF_RSF = 0,
                                                     BOB_RSF = 0, COY_RSF = 0),
                             col = c("red", "blue"), plotCI = TRUE, alpha = 0.95, return =  TRUE)
  for(i in 1:length(stay_tst)){
    stay_tst[[i]]$encamped$State <- "Encamped"
    stay_tst[[i]]$encamped$Dist2Rd <- ranges[["dist2rd"]][1]
    stay_tst[[i]]$exploratory$State <- "Exploratory"
    stay_tst[[i]]$exploratory$Dist2Rd <- ranges[["dist2rd"]][1]
    stay_tst[[i]] <- rbind(stay_tst[[i]]$encamped, stay_tst[[i]]$exploratory)
  }
  
  # THIS WORKED WITH THE INTERNAL FOR LOOP
  blah <- function(t, hmmm, snow) {
    stay_tst2 <- plotStationary(hmmm, covs = data.frame(Dist2Road = t, PercOpen = 0,
                                                        SnowCover = snow, TRI = 0,
                                                        COUG_RSF = 0, WOLF_RSF = 0,
                                                        BOB_RSF = 0, COY_RSF = 0),
                                col = c("red", "blue"),  plotCI = TRUE, alpha = 0.95, return =  TRUE)
    for(i in 1:length(stay_tst2)){
      stay_tst2[[i]]$encamped$State <- "Encamped"
      stay_tst2[[i]]$encamped$Dist2Rd <- t
      stay_tst2[[i]]$exploratory$State <- "Exploratory"
      stay_tst2[[i]]$exploratory$Dist2Rd <- t
      stay_tst2[[i]] <- rbind(stay_tst2[[i]]$encamped, stay_tst2[[i]]$exploratory)
    }
    stay_tst2_df <- bind_rows(stay_tst2)
    return(stay_tst2_df)
  }
  tst <- as.list(ranges[["dist2rd"]])
  tst2 <- lapply(tst, blah, hmmm = hmmm, snow = 0)
 
  
  
  
  # THIS ONLY WORKS WHEN SNOW = 0.... WHY NOT SNOW = 1
  
  
  #'  Function to calculate stationary state probabilities across range of 
  #'  Distance to nearest road covariate and each of the other covariates in the
  #'  model. Requires feeding range of individual Dist2Road values through function
  #'  with lapply().
  SSProb_dist2rd_range <- function(cov_val, hmmm, snow) {
    #'  Calculate SS Prob across range of covariate values and a specific 
    #'  Dist2Rd value while holding other variables at their mean (which is 0 
    #'  when standardize). plotStationary() does this for all covs in the model.
    stay_probs <- plotStationary(hmmm, covs = data.frame(Dist2Road = cov_val, PercOpen = 0,
                                                         SnowCover = snow, TRI = 0,
                                                         COUG_RSF = 0, WOLF_RSF = 0,
                                                         BOB_RSF = 0, COY_RSF = 0),
                                 col = c("red", "blue"),  plotCI = TRUE, alpha = 0.95, return =  TRUE)
    #'  Loop through list of SS Probs based on each covariate and add info about
    #'  state and specific covariate value used to calculate SS Probs
    for(i in 1:length(stay_probs)){
      stay_probs[[i]]$encamped$State <- "Encamped"
      stay_probs[[i]]$encamped$Dist2Rd <- cov_val
      stay_probs[[i]]$exploratory$State <- "Exploratory"
      stay_probs[[i]]$exploratory$Dist2Rd <- cov_val
      #'  Merge encamped and exploratory state estimates for each covariate into
      #'  one data frame (should produce one data frame per covariate in the model)
      stay_probs[[i]] <- rbind(stay_probs[[i]]$encamped, stay_probs[[i]]$exploratory)
    }
    #'  Merge SS Probs calculated across range of covariates into a single data 
    #'  frame per Dist2Rd value of interest
    stay_probs_df <- bind_rows(stay_probs)
    return(stay_probs)
  }
  #### GETS STUCK WHEN SNOW = 1 FOR SOME REASON. DIG INTO THIS.
  
  stay_md_smr_dist2rd <- lapply(as.list(md_smr_range[["dist2rd"]]), SSProb_dist2rd_range, hmmm = spp_HMM_output[[1]], snow = 0) #%>% bind_rows()
  stay_md_wtr_dist2rd <- lapply(as.list(md_wtr_range[["dist2rd"]]), SSProb_dist2rd_range, hmmm = spp_HMM_output[[2]], snow = 1) %>% bind_rows()
  stay_elk_smr_dist2rd <- lapply(as.list(elk_smr_range[["dist2rd"]]), SSProb_dist2rd_range, hmmm = spp_HMM_output[[3]], snow = 0) %>% bind_rows()
  stay_elk_wtr_dist2rd <- lapply(as.list(elk_wtr_range[["dist2rd"]]), SSProb_dist2rd_range, hmmm = spp_HMM_output[[4]], snow = 1) %>% bind_rows()
  # stay_wtd_smr_dist2rd <- lapply(as.list(wtd_smr_range[["dist2rd"]]), SSProb_dist2rd_range, hmmm = spp_HMM_output[[5]], snow = 0) %>% bind_rows()
  # stay_wtd_wtr_dist2rd <- lapply(as.list(wtd_wtr_range[["dist2rd"]]), SSProb_dist2rd_range, hmmm = spp_HMM_output[[6]], snow = 1) %>% bind_rows()
  # stay_coug_OK_smr_dist2rd <- lapply(as.list(coug_smr_OK_range[["dist2rd"]]), SSProb_dist2rd_range, hmmm = spp_HMM_output[[7]], snow = 0) %>% bind_rows()
  # stay_coug_OK_wtr_dist2rd <- lapply(as.list(coug_wtr_OK_range[["dist2rd"]]), SSProb_dist2rd_range, hmmm = spp_HMM_output[[8]], snow = 1) %>% bind_rows()
  # stay_coug_NE_smr_dist2rd <- lapply(as.list(coug_smr_NE_range[["dist2rd"]]), SSProb_dist2rd_range, hmmm = spp_HMM_output[[9]], snow = 0) %>% bind_rows()
  # stay_coug_NE_wtr_dist2rd <- lapply(as.list(coug_wtr_NE_range[["dist2rd"]]), SSProb_dist2rd_range, hmmm = spp_HMM_output[[10]], snow = 1) %>% bind_rows()
  # stay_wolf_OK_smr_dist2rd <- lapply(as.list(wolf_smr_OK_range[["dist2rd"]]), SSProb_dist2rd_range, hmmm = spp_HMM_output[[11]], snow = 0) %>% bind_rows()
  # stay_wolf_OK_wtr_dist2rd <- lapply(as.list(wolf_wtr_OK_range[["dist2rd"]]), SSProb_dist2rd_range, hmmm = spp_HMM_output[[12]], snow = 1) %>% bind_rows()
  stay_wolf_NE_smr_dist2rd <- lapply(as.list(wolf_smr_NE_range[["dist2rd"]]), SSProb_dist2rd_range, hmmm = spp_HMM_output[[13]], snow = 0) %>% bind_rows()
  stay_wolf_NE_wtr_dist2rd <- lapply(as.list(wolf_wtr_NE_range[["dist2rd"]]), SSProb_dist2rd_range, hmmm = spp_HMM_output[[14]], snow = 1) %>% bind_rows()
  # stay_bob_OK_smr_dist2rd <- lapply(as.list(bob_smr_OK_range[["dist2rd"]]), SSProb_dist2rd_range, hmmm = spp_HMM_output[[15]], snow = 0) %>% bind_rows()
  # # stay_bob_OK_wtr_dist2rd <- lapply(as.list(bob_wtr_OK_range[["dist2rd"]]), SSProb_dist2rd_range, hmmm = spp_HMM_output[[16]], snow = 1) %>% bind_rows()
  # # stay_bob_NE_smr_dist2rd <- lapply(as.list(bob_smr_NE_range[["dist2rd"]]), SSProb_dist2rd_range, hmmm = spp_HMM_output[[17]], snow = 0) %>% bind_rows()
  # stay_bob_NE_wtr_dist2rd <- lapply(as.list(bob_wtr_NE_range[["dist2rd"]]), SSProb_dist2rd_range, hmmm = spp_HMM_output[[16]], snow = 1) %>% bind_rows()
  # stay_coy_OK_smr_dist2rd <- lapply(as.list(coy_smr_OK_range[["dist2rd"]]), SSProb_dist2rd_range, hmmm = spp_HMM_output[[17]], snow = 0) %>% bind_rows()
  # stay_coy_OK_wtr_dist2rd <- lapply(as.list(coy_wtr_OK_range[["dist2rd"]]), SSProb_dist2rd_range, hmmm = spp_HMM_output[[18]], snow = 1) %>% bind_rows()
  # stay_coy_NE_smr_dist2rd <- lapply(as.list(coy_smr_NE_range[["dist2rd"]]), SSProb_dist2rd_range, hmmm = spp_HMM_output[[19]], snow = 0) %>% bind_rows()
  # stay_coy_NE_wtr_dist2rd <- lapply(as.list(coy_wtr_NE_range[["dist2rd"]]), SSProb_dist2rd_range, hmmm = spp_HMM_output[[20]], snow = 1) %>% bind_rows()
  
  
  
  # DOESN'T SEEM TO WORK THE SAME WAY AS ABOVE- DATA COMES OUT FORMATTED DIFFERENTLY WHYYYYYYYYYY
  
  #'  Function to calculate stationary state probabilities across range of 
  #'  Distance to nearest road covariate and each of the other covariates in the
  #'  model. Requires feeding range of individual Dist2Road values through function
  #'  with lapply().
  SSProb_tri_range <- function(cov_val, hmmm, snow) {
    #'  Calculate SS Prob across range of covariate values and a specific TRI
    #'  value while holding other variables at their mean (which is 0 when
    #'  standardize). plotStationary() does this for all covariates in the model.
    stay_probs <- plotStationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                         SnowCover = snow, TRI = cov_val,
                                                         COUG_RSF = 0, WOLF_RSF = 0,
                                                         BOB_RSF = 0, COY_RSF = 0),
                                 col = c("red", "blue"),  plotCI = TRUE, alpha = 0.95, return =  TRUE)
    #'  Loop through list of SS Probs based on each covariate and add info about
    #'  state and specific covariate value used to calculate SS Probs
    for(i in 1:length(stay_probs)){
      stay_probs[[i]]$encamped$State <- "Encamped"
      stay_probs[[i]]$encamped$TRI <- cov_val
      stay_probs[[i]]$exploratory$State <- "Exploratory"
      stay_probs[[i]]$exploratory$TRI <- cov_val
      #'  Merge encamped and exploratory state estimates for each covariate into
      #'  one data frame (should produce one data frame per covariate in the model)
      stay_probs[[i]] <- rbind(stay_probs[[i]]$encamped, stay_probs[[i]]$exploratory)
    }
    #'  Merge SS Probs calculated across range of covariates into a single data 
    #'  frame per Dist2Rd value of interest
    stay_probs_df <- bind_rows(stay_probs)
    return(stay_probs_df)
  }
  # stay_md_smr_tri <- lapply(as.list(md_smr_range[["tri"]]), SSProb_tri_range, hmmm = spp_HMM_output[[1]], snow = 0) %>% bind_rows()
  # stay_md_wtr_tri <- lapply(as.list(md_wtr_range[["tri"]]), SSProb_tri_range, hmmm = spp_HMM_output[[2]], snow = 1) %>% bind_rows()
  stay_elk_smr_tri <- lapply(as.list(elk_smr_range[["tri"]]), SSProb_tri_range, hmmm = spp_HMM_output[[3]], snow = 0) #%>% bind_rows()
  stay_elk_wtr_tri <- lapply(as.list(elk_wtr_range[["tri"]]), SSProb_tri_range, hmmm = spp_HMM_output[[4]], snow = 1) %>% bind_rows()
  # stay_wtd_smr_tri <- lapply(as.list(wtd_smr_range[["tri"]]), SSProb_tri_range, hmmm = spp_HMM_output[[5]], snow = 0) %>% bind_rows()
  # stay_wtd_wtr_tri <- lapply(as.list(wtd_wtr_range[["tri"]]), SSProb_tri_range, hmmm = spp_HMM_output[[6]], snow = 1) %>% bind_rows()
  stay_coug_OK_smr_tri <- lapply(as.list(coug_smr_OK_range[["tri"]]), SSProb_tri_range, hmmm = spp_HMM_output[[7]], snow = 0) %>% bind_rows()
  # stay_coug_OK_wtr_tri <- lapply(as.list(coug_wtr_OK_range[["tri"]]), SSProb_tri_range, hmmm = spp_HMM_output[[8]], snow = 1) %>% bind_rows()
  # stay_coug_NE_smr_tri <- lapply(as.list(coug_smr_NE_range[["tri"]]), SSProb_tri_range, hmmm = spp_HMM_output[[9]], snow = 0) %>% bind_rows()
  stay_coug_NE_wtr_tri <- lapply(as.list(coug_wtr_NE_range[["tri"]]), SSProb_tri_range, hmmm = spp_HMM_output[[10]], snow = 1) %>% bind_rows()
  # stay_wolf_OK_smr_tri <- lapply(as.list(wolf_smr_OK_range[["tri"]]), SSProb_tri_range, hmmm = spp_HMM_output[[11]], snow = 0) %>% bind_rows()
  # stay_wolf_OK_wtr_tri <- lapply(as.list(wolf_wtr_OK_range[["tri"]]), SSProb_tri_range, hmmm = spp_HMM_output[[12]], snow = 1) %>% bind_rows()
  # stay_wolf_NE_smr_tri <- lapply(as.list(wolf_smr_NE_range[["tri"]]), SSProb_tri_range, hmmm = spp_HMM_output[[13]], snow = 0) %>% bind_rows()
  # stay_wolf_NE_wtr_tri <- lapply(as.list(wolf_wtr_NE_range[["tri"]]), SSProb_tri_range, hmmm = spp_HMM_output[[14]], snow = 1) %>% bind_rows()
  # stay_bob_OK_smr_tri <- lapply(as.list(bob_smr_OK_range[["tri"]]), SSProb_tri_range, hmmm = spp_HMM_output[[15]], snow = 0) %>% bind_rows()
  # # stay_bob_OK_wtr_tri <- lapply(as.list(bob_wtr_OK_range[["tri"]]), SSProb_tri_range, hmmm = spp_HMM_output[[16]], snow = 1) %>% bind_rows()
  # # stay_bob_NE_smr_tri <- lapply(as.list(bob_smr_NE_range[["tri"]]), SSProb_tri_range, hmmm = spp_HMM_output[[17]], snow = 0) %>% bind_rows()
  # stay_bob_NE_wtr_tri <- lapply(as.list(bob_wtr_NE_range[["tri"]]), SSProb_tri_range, hmmm = spp_HMM_output[[16]], snow = 1) %>% bind_rows()
  # stay_coy_OK_smr_tri <- lapply(as.list(coy_smr_OK_range[["tri"]]), SSProb_tri_range, hmmm = spp_HMM_output[[17]], snow = 0) %>% bind_rows()
  # stay_coy_OK_wtr_tri <- lapply(as.list(coy_wtr_OK_range[["tri"]]), SSProb_tri_range, hmmm = spp_HMM_output[[18]], snow = 1) %>% bind_rows()
  # stay_coy_NE_smr_tri <- lapply(as.list(coy_smr_NE_range[["tri"]]), SSProb_tri_range, hmmm = spp_HMM_output[[19]], snow = 0) %>% bind_rows()
  # stay_coy_NE_wtr_tri <- lapply(as.list(coy_wtr_NE_range[["tri"]]), SSProb_tri_range, hmmm = spp_HMM_output[[20]], snow = 1) %>% bind_rows()
  
  
  save(stay_elk_smr_tri, file = "./Outputs/Figures for ms/HMM Stationary States/SSPlot_Elk_smr_TRIrange.RData")
  save(stay_coug_OK_smr_tri, file = "./Outputs/Figures for ms/HMM Stationary States/SSPlot_Cougar_OK_smr_TRIrange.RData")
  
  
  
  # CAN'T FIGURE OUT HOW TO EXTRACT DATA I NEED TO MAKE FUCKING PLOT
  # 
  # 
  
  ####  Stationary State Plots Over Range of TRI Values  ####
  #'  ---------------------------------------------------
  #'  Mule Deer stationary state across varying Distance to Road and Cougar RSF values
  coug_Dist2Rd_effects_md <- cbind(stay_md_smr_dist2rd$COUG_RSF$est, stay_md_smr_dist2rd$COUG_RSF$cov, stay_md_smr_dist2rd$COUG_RSF$Dist2Rd, stay_md_smr_dist2rd$COUG_RSF$State)
    as.data.frame(tst) %>% #stay_md_smr_dist2rd
    dplyr::select(c(tst$COUG_RSF$est, tst$COUG_RSF$cov, tst$COUG_RSF$Dist2Rd, tst$COUG_RSF$State)) %>%
    filter(!tst$COUG_RSF$State == "Encamped") 
  coug_tri_elk_plot <- ggplot(coug_tri_effects_elk, aes(x = COUG.cov, y = COUG.TRI, fill = COUG.est)) + 
    geom_tile() + 
    # theme(panel.border = element_blank()) +
    # theme(axis.line = element_line(color = 'black')) +
    # xlim(-1, 4) + ylim(0, 1.0) +
    xlab("Scaled cougar RSF value") +
    ylab("Scaled terrain ruggedness index (TRI)") +
    labs(title = "Probability of exploratory state for elk \n with varying terrain ruggedness and cougar RSF values")
  
  
  
  
  
  
 
  