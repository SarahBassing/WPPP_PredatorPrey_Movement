  ##  Wolf & Cougar Telemetry Data Cleaning
  ##  Washington Predator-Prey Project
  ##  July 2021
  ##  Sarah Bassing
  ##  =========================================================
  ##  Most data were previously cleaned, thinned, floored, tz adjusted, and
  ##  AnimalID attached by L.Satterfield. Now looking to see if any dispersals
  ##  need to be excluded.
  ##  Note: L.Satterfield already created a cougar dataset with dispersals
  ##  excluded (Cougar_Vectronic_ATS_Spring2021_4hrs.csv) but I want to make 
  ##  sure these exclusions work for my purposes as well.
  ##  Last chunk of data (Nov. 2020 - Dec. 2021) must be cleaned and integrated
  ##  with previously cleaned data here.
  ##
  ##  Update: L.Satterfield's dispersal calls work for me so will continue to use
  ##  the Cougar_Vectronic_ATS_Spring2021_4hrs.csv dataset that excludes major
  ##  dispersals where the cougar leaves the study area and never returns (and
  ##  generally does not set up a new territory elsewhere... probably died).
  ##  Wolf dispersals and extra-territorial forays are excluded below.
  ##  =========================================================
  
  #  Load libraries
  library(sf)
  library(ggplot2)
  library(tidyverse)
  library(lubridate)
  library(zoo)
  library(purrr)
  
  #  Turn off scientific notation
  options(scipen = 999) 
  #  Set digits to 15 to ensure GPS coordinates aren't truncated
  options(digits=15) 

  #  Read in data
  #  L. Satterfield data
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
  
  #  Additional location data in need of thorough cleaning
  #  Note that seconds are included in the ObservationDateTimePST for these data
  #  so mdy_hm is slightly different from the mdy_hms used above
  wolf_2021 <- read.csv("./Data/dev_telem_wolf_11.2020-12.2021.csv") %>%
    mutate(daytime = mdy_hm(ObservationDateTimePST, tz = "America/Los_Angeles"),
           UTCdt = with_tz(daytime, "UTC"),
           Finaldt = with_tz(UTCdt, tzone = "Etc/GMT+8"),
           Floordt = floor_date(Finaldt, unit = "hour"))
  #  W107F dispersed out of the study area Nov. 2020 and established new pack- Nix all locations
  wolf_2021 <- wolf_2021[!(wolf_2021$IndividualName == "W107F"),]
  
  coug_2021 <- read.csv("./Data/dev_telem_coug_11.2020-12.2021.csv") %>%
    mutate(daytime = mdy_hm(ObservationDateTimePST, tz = "America/Los_Angeles"),
           UTCdt = with_tz(daytime, "UTC"),
           Finaldt = with_tz(UTCdt, tzone = "Etc/GMT+8"),
           Floordt = floor_date(Finaldt, unit = "hour"))
  #  NEC117F points look like the collar sat in Chewelah & Colville for awhile- Nix anything north of 48.25 Lat
  coug_2021 <- coug_2021[!(coug_2021$IndividualName == "NEC117F" & coug_2021$Latitude >= 48.25),]
  #  NEC144M points look like the collar sat in Chewelah & Colville for awhile- Nix anything north of 48.4 Lat
  #  And anything east of 117.8W Long
  coug_2021 <- coug_2021[!(coug_2021$IndividualName == "NEC144M" & coug_2021$Latitude >= 48.4),]
  coug_2021 <- coug_2021[!(coug_2021$IndividualName == "NEC144M" & coug_2021$Longitude >= -117.8),]
  #  NEC160F has all kinds of weird locations- some in Colville but others way N, S, and W
  coug_2021 <- coug_2021[!(coug_2021$IndividualName == "NEC160F" & coug_2021$Latitude >= 48.5),]
  coug_2021 <- coug_2021[!(coug_2021$IndividualName == "NEC160F" & coug_2021$Latitude <= 48.2),]
  coug_2021 <- coug_2021[!(coug_2021$IndividualName == "NEC160F" & coug_2021$Longitude <= -119.0),]
  #  MVC227F appears to have been dispersing in Nov 2020- Nix all locations west of 120.18W Long
  coug_2021 <- coug_2021[!(coug_2021$IndividualName == "MVC227F" & coug_2021$Longitude <= -120.18),]
  
  #'  Clean up columns to match previous datasets
  clean_2021 <- function(locs) {
    clean <- locs %>%
      arrange(IndividualName, Finaldt) %>%
      transmute(
        OBJECTID = OBJECTID,
        PositionID = PositionID,
        CollarID = CollarID,
        Latitude = Latitude,
        Longitude = Longitude,
        ObservationDateTimePST = ObservationDateTimePST,
        TransmissionDateTimePST = TransmissionDateTimePST,
        DbLoadedDateTimePST = DbLoadedDateTimePST,
        ValidLocation = as.numeric(ValidLocation),
        ValidDate = as.numeric(ValidDate),
        VEC_MortalityStatus = as.factor(as.character(VEC_MortalityStatus)),
        VEC_FixType = as.factor(as.character(VEC_FixType)),
        VEC_Origin = as.factor(as.character(VEC_Origin)),
        VEC_DOP = as.numeric(VEC_DOP),
        VEC_Height = as.numeric(VEC_Height),
        Project = Project,
        CaptureID = CaptureID,
        DeploymentID = DeploymentID,
        TransmitterID = TransmitterID,
        Species = Species,
        Sex = as.factor(as.character(Sex)),
        IndividualName = IndividualName,
        Fate = Fate,
        SerialNumber = SerialNumber,
        CaptureDate = CaptureDate,
        FateDate = FateDate,
        DateMortality = DateMortality,
        DaysDelta = DaysDelta,
        daytime = daytime,
        UTCdt = UTCdt,
        Finaldt = Finaldt,
        Floordt = Floordt,
        ID = IndividualName
      ) %>%
      filter(VEC_FixType != "No fix") %>%
      filter(VEC_FixType != "GPS-2D") %>%
      filter(VEC_DOP <= 8)
    return(clean)
  }
  coug_clean <- clean_2021(coug_2021)
  wolf_clean <- clean_2021(wolf_2021)
  
  #  VEH_Height also associated with accuracy so need to see if there are any
  #  odd outliers that should be excluded due to low accuracy fixes
  hist(coug_clean$VEC_Height, breaks = c(100), main = "Cougar locations (ALL VEC_Height values)")
  summary(coug_clean$VEC_Height)
  oddball_coug <- coug_clean[coug_clean$VEC_Height >= 2000 | coug_clean$VEC_Height < 0,]
  
  hist(wolf_clean$VEC_Height, breaks = c(100), main = "Wolf locations (ALL VEC_Height values)")
  summary(wolf_clean$VEC_Height)
  oddball_wolf <- wolf_clean[wolf_clean$VEC_Height >= 2000 | wolf_clean$VEC_Height < 0,]
  
  
  #  Remove oddball locations and drop unnecessary columns so 2021 data match the
  #  rest of the cougar and wolf data cleaned by L. Satterfield
  coug21_clean <- coug_clean %>%
    anti_join(oddball_coug) %>%
    transmute(No = PositionID,
              ID = ID,
              CollarID = CollarID,
              Sex = str_sub(ID, -1),
              Latitude = Latitude,
              Longitude = Longitude,
              LMT_DateTime = as.POSIXct(ObservationDateTimePST, format = "%m/%d/%Y %H:%M", tz = "America/Los_Angeles"),
              StudyArea = ifelse(grepl("NE", ID), "NE", "OK"),
              daytime = daytime,
              Finaldt = Finaldt,
              Floordt = Floordt)
  
  wolf21_clean <- wolf_clean %>%
    anti_join(oddball_wolf) %>%
    transmute(No = PositionID,
              ID = ID,
              CollarID = CollarID,
              Sex = str_sub(ID, -1),
              Latitude = Latitude,
              Longitude = Longitude,
              LMT_DateTime = as.POSIXct(ObservationDateTimePST, format = "%m/%d/%Y %H:%M", tz = "America/Los_Angeles"),
              StudyArea = ifelse(grepl("W107F", ID), "NE", "OK"),  
              StudyArea = ifelse(grepl("W108F", ID), "NE", StudyArea),
              daytime = daytime,
              Finaldt = Finaldt,
              Floordt = Floordt)
  
  
  #  Merge into a single dataset for each species
  coug_cleaned <- rbind(coug_skinny, coug21_clean)
  wolf_cleaned <- rbind(wolf_skinny, wolf21_clean)
  
  write.csv(coug_cleaned, paste0("coug_clean ", Sys.Date(), ".csv"))
  write.csv(wolf_cleaned, paste0("wolf_clean ", Sys.Date(), ".csv"))
   

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
  
  coug_spdf <- st_as_sf(coug21_clean, coords = c("Longitude", "Latitude"), crs = wgs84) 
  wolf_spdf <- st_as_sf(wolf21_clean, coords = c("Longitude", "Latitude"), crs = wgs84)
  
  
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
  pdf("./Outputs/coug21_NE_maps2.pdf")
  for (i in 1:length(unique(coug_NE_maps))) {
    print(coug_NE_maps[[i]])
  }
  dev.off()
  pdf("./Outputs/coug21_OK_maps2.pdf")
  for (i in 1:length(unique(coug_OK_maps))) {
    print(coug_OK_maps[[i]])
  }
  dev.off()
  #  WOlf
  pdf("./Outputs/wolf21_NE_maps2.pdf")
  for (i in 1:length(unique(wolf_NE_maps))) {
    print(wolf_NE_maps[[i]])
  }
  dev.off()
  pdf("./Outputs/wolf21_OK_maps2.pdf")
  for (i in 1:length(unique(wolf_OK_maps))) {
    print(wolf_OK_maps[[i]])
  }
  dev.off()
  
  
  #'  Next step is Collar_Truncating&Filtering.R script
  
  
