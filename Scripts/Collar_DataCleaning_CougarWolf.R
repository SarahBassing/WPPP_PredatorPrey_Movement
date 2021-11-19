  ##  Wolf & Cougar Telemetry Data Cleaning
  ##  Washington Predator-Prey Project
  ##  July 2021
  ##  Sarah Bassing
  ##  =========================================================
  ##  These data were previously cleaned, thinned, floored, tz adjusted, and
  ##  AnimalID attached by L.Satterfield. Now looking to see if any dispersals
  ##  need to be excluded.
  ##  Note: L.Satterfield already created a cougar dataset with dispersals
  ##  excluded (Cougar_Vectronic_ATS_Spring2021_4hrs.csv) but I want to make 
  ##  sure these exclusions work for my purposes as well.
  ##
  ##  Update: L.Satterfield's dispersal calls work for me so will continue to use
  ##  the Cougar_Vectronic_ATS_Spring2021_4hrs.csv dataset that excludes major
  ##  dispersals where the cougar leaves the study area and never returns (and
  ##  generally does not set up a new territory elsewhere... probably died).
  ##  Only wolf dispersal I'm excluding is a series of extra-territorial forays
  ##  that end in a final dispersal (and likely death) to Idaho for W91F.
  ##  =========================================================
  
  #  Load libraries
  library(sf)
  library(ggplot2)
  library(tidyverse)

  #  Read in data
  wolf_skinny <- read.csv("./Data/Wolf_Vectronic_Spring2021_4hrs.csv") %>%
    mutate(ID = as.factor(as.character(ID)),
           CollarID = Collar,
           Sex = str_sub(ID, -1),
           Latitude = Lat,
           Longitude = Long,
           StudyArea = ifelse(grepl("W61M", ID), "OK", "NE"),  
           StudyArea = ifelse(grepl("W71F", ID), "OK", StudyArea),
           StudyArea = ifelse(grepl("W88M", ID), "OK", StudyArea),
           StudyArea = ifelse(grepl("W93M", ID), "OK", StudyArea),
           StudyArea = ifelse(grepl("W94M", ID), "OK", StudyArea),
           daytime = as.POSIXct(LMT_DateTime, format = "%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+8"),
           Finaldt = daytime,
           Floordt = daytime) %>%
    dplyr::select("No", "ID", "CollarID", "Sex", "Latitude", "Longitude", "LMT_DateTime", 
                  "StudyArea", "daytime", "Finaldt", "Floordt")
  # W91 has several extraterritorial forays starting 2019-10-31 & never returns
  wolf_skinny <- wolf_skinny[!(wolf_skinny$ID == "W91F" & wolf_skinny$Floordt > "2019-10-31 02:00:00"),]
  # W88M goes on a W Cascades walkabout (5/13/19 - 5/31/19)- Nix anything west of 120.75W Long
  wolf_skinny <- wolf_skinny[!(wolf_skinny$ID == "W88M" & wolf_skinny$Longitude <= -120.75),]
  # W88M goes on a secondary walkabout north- Nix anything north of 49.3 Lat
  wolf_skinny <- wolf_skinny[!(wolf_skinny$ID == "W88M" & wolf_skinny$Latitude >= 49.3),]
  # W70F dispersed & established Onion Crk pack - Nix anything north of 48.7 Lat
  wolf_skinny <- wolf_skinny[!(wolf_skinny$ID == "W70F" & wolf_skinny$Latitude >= 48.7),]
  
  
  coug_skinny <- read.csv("./Data/Cougar_Vectronic_ATS_Spring2021_4hrs.csv") %>% #Cougar_Vectronic_ATS_Spring2021_4hrs_wDispersal.csv
    mutate(ID = as.factor(as.character(ID)),
           CollarID = Collar,
           Sex = str_sub(ID, -1),
           Latitude = Lat,
           Longitude = Long,
           StudyArea = ifelse(grepl("NE", ID), "NE", "OK"),
           daytime = as.POSIXct(LMT_DateTime, format = "%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+8"),
           Finaldt = daytime,
           Floordt = daytime) %>%
    dplyr::select("No", "ID", "CollarID", "Sex", "Latitude", "Longitude", "LMT_DateTime", 
                  "StudyArea", "daytime", "Finaldt", "Floordt")
   

  ####  ============================================
  ####  Review individual collars for oddities  ####
  
  #  Plot all locations for a given species and look at their distribution
  #  Plot all locations for a given individual and look at their distribution
  #  Take special note of individuals with odd locations or distributions
  
  #  Read in study area shapefiles and reproject to lat/long
  OK_SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "METHOW_SA")
  NE_SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "NE_SA")
  wgs84 <- st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  OK_wgs84 <- st_transform(OK_SA, wgs84)
  NE_wgs84 <- st_transform(NE_SA, wgs84)
  
  #  Make collar location data spatial
  coug_spdf <- st_as_sf(coug_skinny, coords = c("Longitude", "Latitude"), crs = wgs84) 
  wolf_spdf <- st_as_sf(wolf_skinny, coords = c("Longitude", "Latitude"), crs = wgs84)
  
  
  #  Plot all locations (takes a hot minute)
  ggplot() +
    geom_sf(data = coug_spdf, aes(colour = ID)) 
  ggplot() +
    geom_sf(data = wolf_spdf, aes(colour = ID)) 

  #  Plot an individual collar
  ind_coug <- group_split(coug_spdf, coug_spdf$ID)
  ind_wolf <- group_split(wolf_spdf, wolf_spdf$ID)

  #  Plot by study area with study area boundary for context
  #  Function to plot locations from individual animals in the NORTHEAST study area
  plot_telem_NE <- function(spdf){
    #  Split out spatial points df by individual animal ID
    ind_animal <- group_split(spdf, spdf$ID)
    #  Place holder for unique animal ID
    names <- c()
    #  Empty list to hold individual maps
    plot_list <- list()
    #  Loop through all animals one at a time to create maps of their locations
    for(i in 1:length(unique(ind_animal))) {
      names <- c(names, unique(as.character(ind_animal[[i]]$ID)))
      plot <- ggplot() +
        geom_sf(data = NE_SA, fill = NA) +
        geom_sf(data = ind_animal[[i]], aes(color = ID)) +
        labs(title = paste(names[i], "Locations", sep = " "), x = "Longitude",
             y = "Latitude") +
        theme(legend.position = "none")
      plot_list[[i]] <- plot
    }
    return(plot_list)
  }
  
  #  Function to plot locations from individual animals in the OKANOGAN study area
  plot_telem_OK <- function(spdf){
    #  Split out spatial points df by individual animal ID
    ind_animal <- group_split(spdf, spdf$ID)
    #  Place holder for unique animal ID
    names <- c()
    #  Empty list to hold individual maps
    plot_list <- list()
    #  Loop through all animals one at a time to create maps of their locations
    for(i in 1:length(unique(ind_animal))) {
      names <- c(names, unique(as.character(ind_animal[[i]]$ID)))
      plot <- ggplot() +
        geom_sf(data = OK_SA, fill = NA) +
        geom_sf(data = ind_animal[[i]], aes(color = ID)) +
        labs(title = paste(names[i], "Locations", sep = " "), x = "Longitude",
             y = "Latitude") +
        theme(legend.position = "none")
      plot_list[[i]] <- plot
    }
    return(plot_list)
  }
  
  #  Feed wolf/cougar data through function to map individual telemetry data in NE
  coug_NE_maps <- plot_telem_NE(coug_spdf[coug_spdf$StudyArea == "NE",])
  wolf_NE_maps <- plot_telem_NE(wolf_spdf[wolf_spdf$StudyArea == "NE",])
  
  #  Feed wolf/cougar through function to map individual telemetry data in OK
  coug_OK_maps <- plot_telem_OK(coug_spdf[coug_spdf$StudyArea == "OK",])
  wolf_OK_maps <- plot_telem_OK(wolf_spdf[wolf_spdf$StudyArea == "OK",])
  
  #  Save individual plots in a single pdf for each species
  #  With NE or OK study area boundary for context
  #  Cougar
  pdf("./Outputs/coug_NE_maps2.pdf")
  for (i in 1:length(unique(coug_NE_maps))) {
    print(coug_NE_maps[[i]])
  }
  dev.off()
  pdf("./Outputs/coug_OK_maps2.pdf")
  for (i in 1:length(unique(coug_OK_maps))) {
    print(coug_OK_maps[[i]])
  }
  dev.off()
  #  WOlf
  pdf("./Outputs/wolf_NE_maps2.pdf")
  for (i in 1:length(unique(wolf_NE_maps))) {
    print(wolf_NE_maps[[i]])
  }
  dev.off()
  pdf("./Outputs/wolf_OK_maps2.pdf")
  for (i in 1:length(unique(wolf_OK_maps))) {
    print(wolf_OK_maps[[i]])
  }
  dev.off()
  
  #'  Consider trimming W91F, W88M, and W70F- looks like there were some big dispersal events
  #'  Nix W88M west of 120.75W Long- extra-territorial adventure
  #'  Nix W91F south of 48.2N Lat & east of 117.75W Long - dispersed & never returned
  #'  Nix W70F north of 48.7 Lat - she dispersed starting May 20, 2017 & est. Onion Crk pack
  
  #'  Next step is Collar_Truncating&Filtering.R script
  
  
