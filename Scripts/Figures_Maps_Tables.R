  #'  =================================================
  #'  Maps and figures for predator-prey movement
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  March 2022
  #'  =================================================
  #'  1. Plot study area map with inset map showing study area locations in  
  #'     relation to Washington State
  #'  2. Plot stationary-state probabilities derived from HMM results
  #'  =================================================
  
  #'  Clear memory
  rm(list=ls())
  
  #'  Load libraries
  library(momentuHMM)
  library(ggplot2)
  library(ggspatial)
  library(cowplot)
  library(patchwork)
  library(grid)
  library(png)
  library(RCurl)
  library(RColorBrewer)
  library(rphylopic)
  library(sf)
  library(raster)
  library(tidyverse)
  
  #'  Load HMM results
  load("./Outputs/HMM_output/spp_HMM_output_2022-03-15.RData")
  load("./Outputs/Telemetry_crwOut/crwOut_ALL_wCovs_2022-03-16.RData")
  
  #'  Get some basics pulled together to be used across most figures
  #'  -----------------------------
  ####  Spatial Data for Mapping  ####
  #'  -----------------------------
  #'  Define projections
  wgs84 <- projection("+proj=longlat +datum=WGS84 +no_defs")
  sa_proj <- projection("+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs ")
  
  #'  Read in study area data and reproject
  WA <- st_read("./Shapefiles/Washington_State_Boundary/WA_State_Geospatial_Open_Data", layer = "WA_State_Boundary") %>%
    st_transform(crs = sa_proj)
  OK_SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "METHOW_SA") %>%
    st_transform(crs = sa_proj)
  OK_SA$NAME <- "Okanogan"
  NE_SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "NE_SA") %>%
    st_transform(crs = sa_proj)
  NE_SA$NAME <- "Northeast"
  
  projection(WA)
  projection(OK_SA)
  extent(OK_SA)
  extent(NE_SA)
  
  #'  Centroid of each polygon
  st_centroid(OK_SA)
  st_centroid(NE_SA)
  
  #'  Reduce DEM raster resolution and prep new raster for ggplot
  # dem <- raster("./Shapefiles/WA DEM rasters/WPPP_DEM_30m_reproj.tif")
  # projection(dem)
  # dem_low <- aggregate(dem, fact = 10)
  # writeRaster(dem_low, file = "./Shapefiles/WA DEM rasters/dem_reproj_low", format = "GTiff")
  dem_low <- raster("./Shapefiles/WA DEM rasters/dem_reproj_low.tif")
  dem_p_low <- rasterToPoints(dem_low)
  dem_p_df <- as.data.frame(dem_p_low)
  colnames(dem_p_df) <- c("x", "y", "value")
  
  
  ####  1. Map study area and WA State inset  ####
  #'  ============================================
  #'  Plot state of WA with study areas
  #'  https://r-spatial.org/r/2018/10/25/ggplot2-sf.html
  WA_SA_map <- ggplot() + 
    geom_sf(data = WA, fill = "gray95") +
    #' #'  Label map of WA with "Washington State"
    #' geom_sf_text(data = WA, aes(label = JURISDIC_3, hjust = 0.8, vjust = 3), size = 12) +
    geom_sf(data = OK_SA, fill = "grey25", color = "grey20") +
    geom_sf(data = NE_SA, fill = "grey25", color = "grey20") +
    #'  Get rid of lines and axis labels
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),
          #'  No margins around figure
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
  plot(WA_SA_map)
  
  #'  Plot study areas against DEM
  SA_map <- ggplot() +
    geom_raster(data = dem_p_df, aes(x = x, y = y, fill = value, alpha = value), show.legend = FALSE) + 
    #'  alpha adjusts transparency of the raster (can also just set it range = 0.7)
    scale_alpha(range = c(0.3, 0.8)) +
    #'  Change colors of the raster
    scale_fill_gradient2(low = "grey95", high = "tan4") + #gray20
    #'  Add study area outlines and label with their names
    geom_sf(data = OK_SA, fill = NA, color="black", size = 0.80) +
    #'  Note the hjust & vjust need to change based on font size and coord_sf
    #'  DEPENDS ON HOW BIG YOUR PLOT WINDOW IS TOO!!!!
    geom_sf_text(data = OK_SA, aes(label = NAME, hjust = 1.3, vjust = 7), size = 7) + #vjust = -6.90
    geom_sf(data = NE_SA, fill = NA, color="black", size = 0.80) +
    geom_sf_text(data = NE_SA, aes(label = NAME, hjust = 1.3, vjust = -6.5), size = 7) +
    #' #'  Add camera locations and vary color by deployment year
    #' geom_sf(data = cams_reproj, aes(color = Year), shape = 16, size = 3) +
    #' #'  Change camera data aesthetics (make sure it's colorblind friendly)
    #' scale_discrete_manual(aesthetics = "color", values = c("#a6611a", "#018571")) + #c("#dfc27d", "#80cdc1") #c("#601A4A", "#63ACBE")
    #' labs(colour = "Camera\nlocations") +
    #'  Constrain plot to two study areas plus some room on the side & bottom
    coord_sf(xlim = c(480000.0, 810000.0), ylim = c(39000.0, 218000.0), expand = FALSE) +
    #'  Constrain map to just the two study areas only
    # coord_sf(xlim = c(504659.0, 781979.9), ylim = c(102808.3, 211000.4)) +
    #'  Get rid of lines and axis names
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16)) +
    #'  Add north arrow
    annotation_north_arrow(location = "bl", which_north = "true", 
                           pad_x = unit(0.25, "in"), pad_y = unit(0.3, "in"),
                           style = north_arrow_fancy_orienteering) +
    #'  Add scale bar (be sure to double check the scale)
    annotation_scale(location = "bl", width_hint = 0.5)
  plot(SA_map)
  
  #'  Build plot with map of study areas and inset map of WA
  #'  https://upgo.lab.mcgill.ca/2019/12/13/making-beautiful-maps/
  #'  Requires "cowplot" package
  #'  Don't use png or other calls to save while plotting- formatting gets messed up
  #'  Use export option in Plot window and formatting holds
  png(file = "./Outputs/Figures for ms/StudyAreas.png",
      width = 1000, height = 691)
  StudyArea_Map <- ggdraw(SA_map) +
    draw_plot(
      {
        WA_SA_map +
          #'  Label map of WA with "Washington State"
          #'  hjust & vjust will depend on inset map's width/height specified below
          geom_sf_text(data = WA, aes(label = JURISDIC_3, hjust = 0.5, vjust = 2)) 
      },
      #'  Distance along a (0,1) x-axis to draw the left edge of the plot
      x = 0.60,
      #'  Distance along a (0,1) y-axis to draw the bottom edge of the plot
      y = 0.20,
      #'  Width & height of the plot expressed as proportion of the entire ggdraw object
      #'  THIS DEPENDS ON HOW BIG YOUR PLOT WINDOW IS TOO!!!!
      width = 0.22,
      height = 0.22)
  plot(StudyArea_Map)
  dev.off()

  
  ####  2. Plot Stationary-State Probabilities  ####
  #'  ==============================================
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
  
  
  #'  Patchwork figures together in panels
  #'  Only plotting significant relationships based on which covariates were
  #'  significant on the transition probabilities
  library(patchwork)
  #'  MULE DEER panels
  (md_smr_patch <- md_smr_fig[[1]] + 
      md_smr_fig[[3]] + theme(axis.title.y = element_blank()) + 
      md_smr_fig[[5]] + md_smr_fig[[6]] + theme(axis.title.y = element_blank()) + 
      plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Summer Mule Deer Stationary State Probabilities',
                      subtitle = '     Okanogan 2018 - 2021') + plot_layout(ncol = 2))
  (md_wtr_patch <- md_wtr_fig[[1]] + md_wtr_fig[[2]] + theme(axis.title.y = element_blank()) + 
      md_wtr_fig[[3]] + theme(axis.title.y = element_blank()) + #md_wtr_fig[[4]] + 
      md_wtr_fig[[5]] + theme(axis.title.y = element_blank()) + md_wtr_fig[[6]] + 
      theme(axis.title.y = element_blank()) + md_wtr_fig[[7]] + theme(axis.title.y = element_blank()) + 
      md_wtr_fig[[8]] + guide_area() + 
      plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Winter Mule Deer Stationary State Probabilities',
                      subtitle = '     Okanogan 2018 - 2021') + plot_layout(ncol = 3))
  
  #'  ELK panels
  (elk_smr_patch <- elk_smr_fig[[1]] +
      elk_smr_fig[[3]] + theme(axis.title.y = element_blank()) + elk_smr_fig[[4]] + 
      elk_smr_fig[[6]] + theme(axis.title.y = element_blank()) + plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Summer Elk Stationary State Probabilities',
                                     subtitle = '     Northeast 2018 - 2021') + plot_layout(ncol = 2))
  (elk_wtr_patch <- elk_wtr_fig[[1]] + theme(axis.title.y = element_blank()) + 
      elk_wtr_fig[[3]] + theme(axis.title.y = element_blank()) + #elk_wtr_fig[[4]] + 
      elk_wtr_fig[[5]] + theme(axis.title.y = element_blank()) + elk_wtr_fig[[6]] + 
      theme(axis.title.y = element_blank()) + elk_wtr_fig[[8]] + 
      theme(axis.title.y = element_blank()) + guide_area() + plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Winter Elk Stationary State Probabilities',
                      subtitle = '    Northeast 2018 - 2021') + plot_layout(ncol = 3))
  
  #'  WHITE-TAILED DEER panels
  (wtd_smr_patch <- wtd_smr_fig[[1]] + 
      wtd_smr_fig[[4]] + theme(axis.title.y = element_blank()) + wtd_smr_fig[[5]] + 
      theme(axis.title.y = element_blank()) + wtd_smr_fig[[6]] + 
      wtd_smr_fig[[7]] + theme(axis.title.y = element_blank()) + plot_layout(guides = 'collect') + 
      guide_area() + plot_annotation(title = 'Summer White-tailed Deer Stationary State Probabilities',
                                     subtitle = '     Northeast 2018 - 2021') + plot_layout(ncol = 3))
  (wtd_wtr_patch <- wtd_wtr_fig[[1]] + wtd_wtr_fig[[5]] + 
      theme(axis.title.y = element_blank()) + wtd_wtr_fig[[6]] + 
      theme(axis.title.y = element_blank()) + wtd_wtr_fig[[7]] + wtd_wtr_fig[[8]] + 
      theme(axis.title.y = element_blank()) + guide_area() + plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Winter White-tailed Deer Stationary State Probabilities',
                      subtitle = '     Northeast 2018 - 2021') + plot_layout(ncol = 3))
  
  #'  COUGAR panels
  (coug_smr_OK_patch <- coug_smr_OK_fig[[1]] +  
      coug_smr_OK_fig[[3]] + theme(axis.title.y = element_blank()) + 
      coug_smr_OK_fig[[4]] + guide_area() + plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Summer Cougar Stationary State Probabilities',
                      subtitle = '     Okanogan 2018 - 2021') + plot_layout(ncol = 2))
  (coug_wtr_OK_patch <- coug_wtr_OK_fig[[1]] + coug_wtr_OK_fig[[2]] + 
      theme(axis.title.y = element_blank()) + #coug_wtr_OK_fig[[4]] + 
      coug_wtr_OK_fig[[3]] + guide_area() + plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Winter Cougar Stationary State Probabilities',
                      subtitle = '     Okanogan 2018 - 2021') + plot_layout(ncol = 2))
  (coug_smr_NE_patch <- coug_smr_NE_fig[[1]] + coug_smr_NE_fig[[2]] + 
      theme(axis.title.y = element_blank()) + coug_smr_NE_fig[[3]] + 
      theme(axis.title.y = element_blank()) + coug_smr_NE_fig[[4]] + 
      theme(axis.title.y = element_blank()) + coug_smr_NE_fig[[5]] + 
      guide_area() + plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Summer Cougar Stationary State Probabilities',
                      subtitle = '     Northeast 2018 - 2021') + plot_layout(ncol = 2))
  (coug_wtr_NE_patch <- coug_wtr_NE_fig[[1]] + coug_wtr_NE_fig[[2]] + 
      theme(axis.title.y = element_blank()) + coug_wtr_NE_fig[[3]] + 
      coug_wtr_NE_fig[[5]] + #coug_wtr_NE_fig[[4]] + 
      theme(axis.title.y = element_blank()) + plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Winter Cougar Stationary State Probabilities',
                      subtitle = '     Northeast 2018 - 2021') + plot_layout(ncol = 2))
  #'  WOLF panels
  (wolf_smr_OK_patch <- wolf_smr_OK_fig[[1]] + wolf_smr_OK_fig[[2]] + 
      theme(axis.title.y = element_blank()) + 
      wolf_smr_OK_fig[[4]] + guide_area() + plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Summer Wolf Stationary State Probabilities',
                      subtitle = '     Okanogan 2018 - 2021') + plot_layout(ncol = 2))
  (wolf_wtr_OK_patch <- wolf_wtr_OK_fig[[2]] + wolf_wtr_OK_fig[[3]] + 
      plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Winter Wolf Stationary State Probabilities',
                      subtitle = '     Okanogan 2018 - 2021') + plot_layout(ncol = 2))
  (wolf_smr_NE_patch <- wolf_smr_NE_fig[[1]] +  
      wolf_smr_NE_fig[[3]] + theme(axis.title.y = element_blank()) + 
      wolf_smr_NE_fig[[4]] + guide_area() + plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Summer Wolf Stationary State Probabilities',
                      subtitle = '     Northeast 2018 - 2021') + plot_layout(ncol = 2))
  (wolf_wtr_NE_patch <-  wolf_wtr_NE_fig[[6]] +
      plot_annotation(title = 'Winter Wolf Stationary State Probabilities',
                      subtitle = '     Northeast 2018 - 2021')) #+ plot_layout(ncol = 2)
  #'  COYOTE panels
  (coy_smr_OK_patch <- coy_smr_OK_fig[[1]] +  coy_smr_OK_fig[[4]] + 
      theme(axis.title.y = element_blank()) + plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Summer Coyote Stationary State Probabilities',
                      subtitle = '     Okanogan 2018 - 2021') + plot_layout(ncol = 2))
  (coy_wtr_OK_patch <- coy_wtr_OK_fig[[1]] + coy_wtr_OK_fig[[2]] + 
      theme(axis.title.y = element_blank()) + 
      plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Winter Coyote Stationary State Probabilities',
                      subtitle = '     Okanogan 2018 - 2021') + plot_layout(ncol = 2))
  (coy_smr_NE_patch <- coy_smr_NE_fig[[3]] + 
      plot_annotation(title = 'Summer Coyote Stationary State Probabilities',
                      subtitle = '     Northeast 2018 - 2021')) #+ plot_layout(ncol = 2)
  (coy_wtr_NE_patch <- coy_wtr_NE_fig[[1]] + coy_wtr_NE_fig[[2]] + 
      theme(axis.title.y = element_blank()) +  
      # coy_wtr_NE_fig[[4]] + theme(axis.title.y = element_blank()) + 
      coy_wtr_NE_fig[[6]] + guide_area() + plot_layout(guides = 'collect') + 
      plot_annotation(title = 'Winter Coyote Stationary State Probabilities',
                      subtitle = '     Northeast 2018 - 2021') + plot_layout(ncol = 2))
  

  
  ####  Save patches for manuscript  ####
  png(file="./Outputs/Figures for ms/HMM Stationary States/MD_smr_OK.png", width = 700, height = 500)
  plot(md_smr_patch)
  dev.off()
  png(file="./Outputs/Figures for ms/HMM Stationary States/MD_wtr_OK.png", width = 700, height = 500)
  plot(md_wtr_patch)
  dev.off()
  png(file="./Outputs/Figures for ms/HMM Stationary States/ELK_smr_NE.png", width = 700, height = 500)
  plot(elk_smr_patch)
  dev.off()
  png(file="./Outputs/Figures for ms/HMM Stationary States/ELK_wtr_NE.png", width = 700, height = 500)
  plot(elk_wtr_patch)
  dev.off()
  png(file="./Outputs/Figures for ms/HMM Stationary States/WTD_smr_NE.png", width = 700, height = 500)
  plot(wtd_smr_patch)
  dev.off()
  png(file="./Outputs/Figures for ms/HMM Stationary States/WTD_wtr_NE.png", width = 700, height = 500)
  plot(wtd_wtr_patch)
  dev.off()
  png(file="./Outputs/Figures for ms/HMM Stationary States/COUG_smr_OK.png", width = 700, height = 500)
  plot(coug_smr_OK_patch)
  dev.off()
  png(file="./Outputs/Figures for ms/HMM Stationary States/COUG_wtr_OK.png", width = 700, height = 500)
  plot(coug_wtr_OK_patch)
  dev.off()
  png(file="./Outputs/Figures for ms/HMM Stationary States/COUG_smr_NE.png", width = 700, height = 500)
  plot(coug_smr_NE_patch)
  dev.off()
  png(file="./Outputs/Figures for ms/HMM Stationary States/COUG_wtr_NE.png", width = 700, height = 500)
  plot(coug_wtr_NE_patch)
  dev.off()
  png(file="./Outputs/Figures for ms/HMM Stationary States/WOLF_smr_OK.png", width = 700, height = 500)
  plot(wolf_smr_OK_patch)
  dev.off()
  png(file="./Outputs/Figures for ms/HMM Stationary States/WOLF_wtr_OK.png", width = 700, height = 500)
  plot(wolf_wtr_OK_patch)
  dev.off()
  png(file="./Outputs/Figures for ms/HMM Stationary States/WOLF_smr_NE.png", width = 700, height = 500)
  plot(wolf_smr_NE_patch)
  dev.off()
  png(file="./Outputs/Figures for ms/HMM Stationary States/WOLF_wtr_NE.png", width = 700, height = 500)
  plot(wolf_wtr_NE_patch)
  dev.off()
  png(file="./Outputs/Figures for ms/HMM Stationary States/COY_smr_OK.png", width = 700, height = 500)
  plot(coy_smr_OK_patch)
  dev.off()
  png(file="./Outputs/Figures for ms/HMM Stationary States/COY_wtr_OK.png", width = 700, height = 500)
  plot(coy_wtr_OK_patch)
  dev.off()
  png(file="./Outputs/Figures for ms/HMM Stationary States/COY_smr_NE.png", width = 700, height = 500)
  plot(coy_smr_NE_patch)
  dev.off()
  png(file="./Outputs/Figures for ms/HMM Stationary States/COY_wtr_NE.png", width = 700, height = 500)
  plot(coy_wtr_NE_patch)
  dev.off()
  
  
  
  
  ####  Summary table for GPS collar data  ####
  #'  =========================================
  #'  Load track and crwOut data
  load("./Outputs/Telemetry_tracks/spp_all_tracks_noDis_noMig_SAspecific_2022-03-14.RData")
  load("./Outputs/Telemetry_crwOut/crwOut_ALL_wCovs_2022-03-16.RData")
  
  #'  Function to generate summary info about collar data 
  #'  Unique collars included in study
  collar_deets <- function(smr_tracks, wtr_tracks, smr_crw, wtr_crw, spp) { #, area
    #'  Merge summer and winter track data
    locs <- rbind(smr_tracks, wtr_tracks) 
    # locs <- rbind(spp_all_tracks[[1]], spp_all_tracks[[2]])
    
    #' Identify unique number of collars per species
    smr_collars <- length(unique(locs$AnimalID[locs$Season == "Summer18" | locs$Season == "Summer19" | locs$Season == "Summer20"]))
    wtr_collars <- length(unique(locs$AnimalID[locs$Season == "Winter1819" | locs$Season == "Winter1920" | locs$Season == "Winter2021"]))
    print("number of collars in summer")
    print(smr_collars)
    print("number of collars in winter")
    print(wtr_collars)
    
    #' Identify unique number of collars per season
    ncollars <- locs %>%
      group_by(Season, StudyArea) %>%
      summarise(collars = n_distinct(AnimalID)) %>%
      ungroup
    
    #'  Count the number of GPS locations generated per season
    nlocs <- locs %>%
      group_by(Season, StudyArea) %>%
      summarise(n = n()) %>%
      ungroup()
    colnames(nlocs) <- c("Season", "StudyArea", "nlocs")
    
    #'  Count the number of animal tracks per season
    ntracks <- locs %>%
      group_by(Season, StudyArea) %>%
      summarise(tracks = n_distinct(ID)) %>%
      ungroup()
    colnames(ntracks) <- c("Season", "StudyArea", "ntracks")
    
    #'  Merge summer and winter crwOut data
    spp_crw <- rbind(smr_crw, wtr_crw) 
    # spp_crw <- rbind(hmm_data[[1]], hmm_data[[2]])
    # spp <- "Mule Deer"
    #'  Count the number of GPS locations + interpolated locations
    ninterp <- spp_crw %>%
      group_by(Season, StudyArea) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      full_join(nlocs, by = c("Season", "StudyArea")) %>%
      #'  Calculate the proportion of interpolated locations per season
      mutate(interpolated = n - nlocs,
             prop = interpolated/n) %>%
      dplyr::select(c("Season", "StudyArea", "nlocs", "interpolated", "prop")) %>%
      #'  Add additional location and track counts to make a single data frame
      full_join(ncollars, by = c("Season", "StudyArea")) %>%
      full_join(ntracks, by = c("Season", "StudyArea")) %>%
      mutate(Species = spp, 
             Season2 = ifelse(Season == "Summer18" | Season == "Summer19" | Season == "Summer20", "Summer", "Winter")) %>%
      relocate(Species, .before = Season) %>%
      relocate(StudyArea, .after = Species) %>%
      relocate(collars, .after = Season) %>%
      relocate(ntracks, .after = collars)
    
    return(ninterp)
    
  }
  #'  Compile collar data for each species
  md_collars <- collar_deets(spp_all_tracks[[1]], spp_all_tracks[[2]], hmm_data[[1]], hmm_data[[2]], spp = "Mule Deer")
  elk_collars <- collar_deets(spp_all_tracks[[3]], spp_all_tracks[[4]], hmm_data[[3]], hmm_data[[4]], spp = "Elk")
  wtd_collars <- collar_deets(spp_all_tracks[[5]], spp_all_tracks[[6]], hmm_data[[5]], hmm_data[[6]], spp = "White-tailed Deer")
  coug_OK_collars <- collar_deets(spp_all_tracks[[7]], spp_all_tracks[[8]], hmm_data[[7]], hmm_data[[8]], spp = "Cougar")
  coug_NE_collars <- collar_deets(spp_all_tracks[[9]], spp_all_tracks[[10]], hmm_data[[9]], hmm_data[[10]], spp = "Cougar")
  wolf_OK_collars <- collar_deets(spp_all_tracks[[11]], spp_all_tracks[[12]], hmm_data[[11]], hmm_data[[12]], spp = "Wolf")
  wolf_NE_collars <- collar_deets(spp_all_tracks[[13]], spp_all_tracks[[14]], hmm_data[[13]], hmm_data[[14]], spp = "Wolf")
  bob_OK_collars <- collar_deets(spp_all_tracks[[15]], spp_all_tracks[[16]], hmm_data[[15]], hmm_data[[16]], spp = "Bobcat")
  bob_NE_collars <- collar_deets(spp_all_tracks[[17]], spp_all_tracks[[18]], hmm_data[[17]], hmm_data[[18]], spp = "Bobcat")
  coy_OK_collars <- collar_deets(spp_all_tracks[[19]], spp_all_tracks[[20]], hmm_data[[19]], hmm_data[[20]], spp = "Coyote")
  coy_NE_collars <- collar_deets(spp_all_tracks[[21]], spp_all_tracks[[22]], hmm_data[[21]], hmm_data[[22]], spp = "Coyote")
  
  #'  Merge into one large table
  collar_table <- rbind(md_collars, elk_collars, wtd_collars, coug_OK_collars, 
                        coug_NE_collars, wolf_OK_collars, wolf_NE_collars, 
                        bob_OK_collars, bob_NE_collars, coy_OK_collars, coy_NE_collars)
  
  #'  Summarize collar table for manuscript
  #'  Average across years
  mu_collars <- collar_table %>%
    group_by(Species, StudyArea, Season2) %>%
    summarise(across(c("collars", "ntracks", "nlocs", "interpolated", "prop"), ~ mean(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(
      collars = round(collars, 0),
      ntracks = round(ntracks, 0),
      nlocs = round(nlocs, 0),
      interpolated = round(interpolated, 0),
      prop = round(prop, 2)
    )
  colnames(mu_collars) <- c("Species", "Study Area", "Season", "Mean_Collars", "Mean_Tracks", "Mean_Locations", "Mean_Interpolated", "Mean_Proportion_Interpolated")
  
  #'  Calculate standard error for averages
  se_collars <- collar_table %>%
    group_by(Species, StudyArea, Season2) %>%
    summarise(across(c("collars", "ntracks", "nlocs", "interpolated", "prop"), ~sd(.x)/sqrt(length(.x)))) %>% 
    ungroup() %>%
    transmute(
      Species= Species,
      StudyArea = StudyArea,
      Season = Season2,
      #'  Round values and put inside parentheses
      se_collars = paste0("(",round(collars, 0), ")"),
      se_ntracks = paste0("(",round(ntracks, 0), ")"),
      se_nlocs = paste0("(",round(nlocs, 0), ")"),
      se_interpolated = paste0("(",round(interpolated, 0), ")"),
      se_prop = paste0("(",round(prop, 3), ")")
    )
  colnames(se_collars) <- c("Species", "Study Area", "Season", "SE_Collars", "SE_Tracks", "SE_Locations", "SE_Interpolated", "SE_Proportion_Interpolated")
  
  #'  Create table of means and SE
  dat_table <- full_join(mu_collars, se_collars, by = c("Species", "Study Area", "Season")) %>%
    unite(Collars, Mean_Collars, SE_Collars, sep = " ") %>%
    unite(Tracks, Mean_Tracks, SE_Tracks, sep = " ") %>%
    unite(Locations, Mean_Locations, SE_Locations, sep = " ") %>%
    unite(Interpolated, Mean_Interpolated, SE_Interpolated, sep = " ") %>%
    unite(Proportion, Mean_Proportion_Interpolated, SE_Proportion_Interpolated, sep = " ")
  colnames(dat_table) <- c("Species", "Study Area", "Season", "Mean Collars (SE)", "Mean Tracks (SE)", "Mean Locations (SE)", "Mean Interpolated Locations (SE)", "Mean Proportion Interpolated Locations (SE)")
  
  summary_table <- as.data.frame(dat_table)
  
  #'  Save
  write.csv(summary_table, file = paste0("./Outputs/Figures for ms/collar_summary_table_", Sys.Date(), ".csv"))
  
  #'  Quick summary stats for in-text reporting
  season_means <- collar_table %>%
    group_by(Season2) %>%
    summarise(across(c("collars", "ntracks", "nlocs", "interpolated", "prop"), ~ mean(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(
      collars = round(collars, 0),
      ntracks = round(ntracks, 0),
      nlocs = round(nlocs, 0),
      interpolated = round(interpolated, 0),
      prop = round(prop, 2)
    )
  colnames(season_means) <- c("Season", "Mean_Collars", "Mean_Tracks", "Mean_Locations", "Mean_Interpolated", "Mean_Proportion_Interpolated")
  
  season_se <- collar_table %>%
    group_by(Season2) %>%
    summarise(across(c("collars", "ntracks", "nlocs", "interpolated", "prop"), ~sd(.x)/sqrt(length(.x)))) %>% 
    ungroup() %>%
    transmute(
      Season = Season2,
      se_collars = round(collars, 0),
      se_ntracks = round(ntracks, 0),
      se_nlocs = round(nlocs, 0),
      se_interpolated = round(interpolated, 0),
      se_prop = round(prop, 3)
    )
  colnames(se_collars) <- c("Species", "Study Area", "Season", "SE_Collars", "SE_Tracks", "SE_Locations", "SE_Interpolated", "SE_Proportion_Interpolated")
  
  season_summary <- full_join(season_means, season_se, by = "Season")
  
  #'  Save
  write.csv(season_summary, file = paste0("./Outputs/Figures for ms/seasonal_summary_table_", Sys.Date(), ".csv"))
  
