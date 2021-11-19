  ##  Mesopredator Telemetry Data Cleaning
  ##  Washington Predator-Prey Project
  ##  April 2021
  ##  Sarah Bassing
  ##  =========================================================
  ##  Script to create MASTER & CLEANED data sets for WPPP mesopredator data.
  ##    1. Combines GPS collar location data with unique animal IDs, & other
  ##       capture/mortality information provided by B. Windell
  ##    2. Creates a MASTER data set that includes all fixes (missed & successful) 
  ##       and makes sure date/times are in correct time zone (Pacific Standard 
  ##       Time all year round). 
  ##    3. Creates CLEAN data set that excludes missed locations and low accuracy 
  ##       fixes. Removes extra locations generated when fix rate increased beyond 
  ##       4-hr fix schedule and on wrong 4-hr schedule. Truncate by capture date  
  ##       and mortality date (if applicable). 
  ##    4. Maps telemetry locations by individual animal for more thorough review
  ##       of location data.
  ##    Bonus: Calculates summary stats about problem locations by species and 
  ##       individual animal.
  ##  Script was originally written by Taylor Ganz in the Prugh Lab to prepare
  ##  mule deer location data for the WPPP. I have expanded it.
  ##  =========================================================
  
  #  Clean work space and load libraries
  rm(list = ls())
  
  library(lubridate)
  library(zoo)
  library(purrr)
  library(tidyverse)
  
  #  Turn off scientific notation
  options(scipen = 999) 
  #  Set digits to 15 to ensure GPS coordinates aren't truncated
  options(digits=15) 
  
  
  ####  ====================================================
  ####  Read in Capture, Mortality, and Telemetry data  ####
  
  #  B.Windell combined capture & mortality data (provided Nov. 16 2020)
  #  L.Prugh sent udpated data (Nov. 10, 2021)
  meso_capmort <- read.csv("./Data/MesocarnivoreOverview_11112021.csv", stringsAsFactors = FALSE) #04192021
  
  #  Telemetry data (MesocarnivoreData_AllLocations_041521.csv provided Nov. 16 2020)
  #  Updated data (WPPP_meso_GSPlocs_11.10.21.csv provided Nov. 10, 2021) but missing
  #  some of the original collars from AllLocations
  #  NOTE: Acquisition time is UTC time without daylight savings!!!
  #  Add column with date/time in usable format & set timezone that data downloaded in
  #  UTC column kept for consistency across data sets
  #  Add column that converts times to Pacific Standard Time (add 8 hours), ignoring daylight savings
  #  Note: Etc/GMT+8 is UTC -8 and outputs in Pacific standard time only (so no daylight savings time)
  #  Add column that floors the time to the nearest hour not ahead
  meso_tel <- read.csv("./Data/WPPP_meso_GPSlocs_11.10.21.csv") %>%    #MesocarnivoreData_AllLocations_041521.csv
    mutate(daytime = as.POSIXct(Acquisition.Time, format = "%m/%d/%Y %H:%M", tz = "UTC"), #"%Y.%m.%d %H:%M:%S"
           UTCdt = daytime,
           Finaldt = with_tz(UTCdt, tzone = "Etc/GMT+8"),
           Floordt = floor_date(Finaldt, unit = "hour")) #%>%
    # dplyr::select(-X)
  meso_tel_missing <- read.csv("./Data/WPPP_meso_MISSING_GPSlocs_11.10.21.csv") %>%    
    mutate(daytime = as.POSIXct(Acquisition.Time, format = "%Y.%m.%d %H:%M:%S", tz = "UTC"), 
           UTCdt = daytime,
           Finaldt = with_tz(UTCdt, tzone = "Etc/GMT+8"),
           Floordt = floor_date(Finaldt, unit = "hour"))
  #  Combine datasets (keep in mind Acquisition.Time are in two different date 
  #  formats so this column will be useless from here-on out)
  meso_tel <- rbind(meso_tel, meso_tel_missing)
  
  #  Check out the time data
  #  What's the timezone?
  tz(meso_tel$daytime); tz(meso_tel$UTCdt); tz(meso_tel$Finaldt)
  #  Does it account for daylight savings time?
  head(dst(meso_tel$daytime))
  head(dst(meso_tel$Finaldt))  # better say FALSE
  
  #str(elk_tel)#; head(elk_tel)
  
  ####  ======================================================
  ####  Combine Capture & Mortality data for each animal  ####

  #  Combine capture and mortality data by unique animal ID               
  meso_info <- meso_capmort %>%
    transmute(
      IndividualIdentifier = as.factor(AnimalID),
      IndividualSpecies = Species,
      CaptureDate = mdy(CaptureDate),
      GPSCollarSerialNumber = CollarID,
      EndDate = mdy(EndDate),
      EndCause = EndCause,
      MortalityType = MortalityType,
      MortalitySubType = MortalitySubType,
      LastTransmission = mdy(LastTransmission),
      Notes = Notes
    ) %>%
    #  Fill in EndDate with last transmission date if animal did not die/collar did not fail
    mutate(
      EndDate = ifelse(is.na(EndDate), as.Date(LastTransmission), as.Date(EndDate)),
      EndDate = as.Date(EndDate)
    )
  #  Remove individuals that don't have a corresponding GPSCollarSerialNumber
  which(is.na(meso_info$GPSCollarSerialNumber))
  meso_info <- droplevels(meso_info[!is.na(meso_info$GPSCollarSerialNumber),])
  #  Remove individuals where collar failed on capture date (MVCOY98M, MVCOY53M)
  deadcap <- as.character(meso_info$IndividualIdentifier[which(meso_info$CaptureDate == meso_info$EndDate)])
  if(is_empty(deadcap) != TRUE) {
    meso_info <- meso_info[!(meso_info$IndividualIdentifier %in% deadcap),]
  } else {
    meso_info <- meso_info
  }
  #  Identify any GPS collars on captured mesopredators that never generated telemetry data
  notel <- meso_info$GPSCollarSerialNumber[!(meso_info$GPSCollarSerialNumber %in% meso_tel$CollarID)]
  if(is_empty(notel) != TRUE) {
    meso_info <- meso_info[!(meso_info$GPSCollarSerialNumber %in% notel),]
  } else {
    meso_info <- meso_info
  }
  #  Drops NEBOB9F (707810A) & MVCOY72F (699719B) due to mortalities within a few days of capture

  
  #  Save for later use
  # write.csv(meso_info, paste0('./Data/meso_info ', Sys.Date(), '.csv'))
  
  ####  Running list of collars that got nixed in this stage  ####
  #  MVCOY53M, MCVOY98M, MVCOY72F, NEBOB9F
  
  ####  =====================================================
  ####  Create MASTER data set for each ungulate species ####
  
  #  Connect location data to individual animal IDs:
  #  Attach animal IDs to their respective telemetry locations.
  #  This truncates the data so that day of capture and day of mortality 
  #  (or today's locations) are excluded from the data set. Further truncating 
  #  can happen for individual analyses.
  
  # #  Make sure the same collars are in the meso_info & meso_tel data frames
  # ninfo <- meso_info %>%
  #   arrange(IndividualIdentifier)
  # ninfo <- unique(as.character(ninfo$IndividualIdentifier))
  # ntel <- meso_tel %>%
  #   arrange(AnimalID)
  # ntel <- unique(as.character(ntel$AnimalID))
  # count <- 1:length(ninfo)
  # diff <- as.data.frame(cbind(count, ninfo, ntel))
  
  #  Create empty data frame to fill iteratively
  clean <- data.frame()
  #  How many individuals are looped over?
  nrow(meso_info)
  #  Loop over every unique individual animal and...
  for(i in 1:nrow(meso_info)){
    #  Take the individual animal ID
    ID <- droplevels(meso_info$IndividualIdentifier[i])
    #  Take the animal's GPS collar serial number 
    SN <- meso_info$GPSCollarSerialNumber[i]
    #  Buffer capture date to remove locations affected by capture event
    #  Suggested to only use data from 2 weeks after the capture data
    #  Currently not buffering but can change the +/- values in the future
    start <- meso_info$CaptureDate[i] + 1
    #  Exclude locations 1 day before estimated mortality date
    end <- meso_info$EndDate[i] - 1
    
    #  Subset telemetry data to the specific individual
    collar <- subset(meso_tel, CollarID == SN)
    #  Add a new column to the telemetry data with the animal's individual ID
    collar$ID <- ID
    #  truncate telemetry data by new start and end dates for that individual
    collarlive <- subset(collar, Finaldt >= start & Finaldt <= end)
    
    #  Append each unique animal's locations to a clean data frame
    clean <- rbind(clean, collarlive)
  }
    
  # Organize by individual ID and chronological order of locations
  # Format data fields and retain data on fix quality to help with data cleaning
  meso_master <- clean %>%
    arrange(ID, Finaldt) %>%
    transmute(
      CollarID = CollarID,
      Acquisition.Time.UTC = Acquisition.Time,
      Latitude = GPS.Latitude,
      Longitude = GPS.Longitude,
      # Iridium.Latitude = Iridium.Latitude,
      # Iridium.Longitude = Iridium.Longitude,
      UTM.Zone = GPS.UTM.Zone,
      UTM.Northing = GPS.UTM.Northing,
      UTM.Easting = GPS.UTM.Easting,
      Species = Species,
      Sex = str_sub(ID, -1),
      Fix.Quality = GPS.Fix.Attempt,
      # GPS.Altitude = GPS.Altitude,
      GPS.Horizontal.Error = GPS.Horizontal.Error,
      GPS.Horizontal.Dilution = GPS.Horizontal.Dilution,
      daytime = daytime,
      UTCdt = UTCdt,
      Finaldt = Finaldt,
      Floordt = Floordt,
      ID = ID
    )
  #  Pacific Standard Time for location data!
  
  #  Remove obviously wrong locations for individual collars
  meso_master <- meso_master %>%
    filter(ID != "NEBOB10F" | Longitude < -117.2) %>%
    filter(ID != "NEBOB37M" | Latitude < 48.2) %>%
    filter(ID != "NECOY44M" | Longitude < -117.95) %>%
    filter(ID != "NECOY44M" | Latitude < 48.6) %>%
    filter(ID != "NECOY43F" | Longitude < -117.2) %>%
    filter(ID != "NEBOB23M" | Latitude > 48.2) %>%
    filter(ID != "MVBOB90M"  | Latitude < 49.0) %>%
    filter(ID != "MVBOB69F" | Longitude < -119.8) %>%
    filter(ID != "MVBOB69F" | Longitude > -120.3) %>%
    filter(ID != "MVBOB69F" | Latitude < 48.7) %>%
    filter(ID != "MVCOY59F" | Latitude < 48.6) %>%
    filter(ID != "MVCOY18F" | Latitude > 47.8) %>%
    filter(ID != "MVBOB51M" | Latitude > 48.4) %>%
    filter(ID != "MVBOB55M" | Latitude < 49.0) %>%
    filter(ID != "MVBOB55M" | Latitude > 48.0) %>%
    filter(ID != "MVCOY79F" | Longitude < -120.0) %>%
    filter(ID != "MVBOB52M" | Longitude < -120.0) %>%
    filter(ID != "MVBOB52M" | Longitude > -121.0)
  
  
  #  Save MASTER data files
  # write.csv(meso_master, paste0('meso_master ', Sys.Date(), '.csv'))

  
  ####  ====================================================
  ####  Create CLEAN data set for each meso species ####
  
  #  Further cleaning of telemetry data ready for analyses-
  #  Drop missing fixes and locations with poor accuracy
  #  Drop extra locations when fix schedule increases
  
  #  Drop missing fixes 
  #  Drop locations with poor accuracy for a cleaner data set
  #  Fix Quality: Resolved QFP (Quick Fix Points) accurate as 3D GPS fix;
  #               Resolved QFP (Uncertain) are least accurate (btwn 100 to 1000s m off);
  #               Succeeded is a regular GPS with associated Horizontal Error- unusable
  #               ---> Retaining ONLY Resolved QFP locations
  #  GPS Altitude: geometric altitude above mean sea level (meters); generally accurate to 10-20m
  #               ---> Excluding locations with GPS Altitude < 0 or > 3000
  #  GPS Horizontal Error: error metric for Succeeded points (meters)
  #               ---> Excluding locations with any GPS Horizontal error
  #  GPS Horizontal Dilution: error metric for QFP (meters), similar to Vectronic DOP???
  #  DOP >10 is not good (~20-30 meter accuracy)--- using DOP = 8 as cutoff
  meso_clean <- meso_master %>%
    # filter(is.na(Iridium.Latitude)) %>%
    filter(!is.na(Longitude)) %>%
    filter(Fix.Quality == "Resolved QFP") %>%
    filter(GPS.Horizontal.Dilution < 8) #%>%
    # filter(GPS.Altitude > 0) %>%
    # filter(GPS.Altitude < 3000)
  
  # #  Altitude & Horizontal Dilution are associated with accuracy so need to see 
  # #  if there are any outliers that should be excluded due to low accuracy fixes
  # hist(meso_clean$GPS.Altitude, breaks = c(100), main = "Meso locations (ALL GPS.Altitude values)")
  # plot(meso_clean$GPS.Altitude)  # takes awhile
  # oddball <- meso_clean[meso_clean$GPS.Altitude >= 3000,] %>%
  #   filter(!is.na(CollarID))
  # hist(meso_clean$GPS.Horizontal.Dilution, breaks = c(100), main = "Meso locations (ALL GPS.Horizontal.Dilution values)")
  # plot(meso_clean$GPS.Horizontal.Dilution)  # takes awhile
  # oddball <- meso_clean[meso_clean$GPS.Horizontal.Dilution >= 9,] %>%
  #   filter(!is.na(CollarID))

  #  Save locations with high accuracy fixes only
  # write.csv(meso_clean, paste0('meso_clean ', Sys.Date(), '.csv'))

  
  
  ####  ==========================
  ####  Rarefy location data  ####
  #  KEEP locations with 4-hr fix schedule with 2, 6, 10, 14, 18, 22 hr fixes
  #  DISCARD: 
  #    -Extra locations from collars in mortality-mode where collars increase fix
  #     schedule to 20-min. interval.
  #    -Initial deployment data on the wrong fix schedule (0, 4, 8 hr)- NOT DOING currently.
  
  thin_locs <- function(clean) {
    #  Keep only locations on correct 4-hr fix schedule
    skinny <- with(clean, clean[hour(Floordt) == 2 | hour(Floordt) == 6 | 
                                  hour(Floordt) == 10 | hour(Floordt) == 14 | 
                                  hour(Floordt) == 18 | hour(Floordt) == 22,])
    #  Also keep locations on wrong 4-hr fix schedule 
    #  This effectively does NOT thin extra 2-hr elk locations
    skinny2 <- with(clean, clean[hour(Floordt) == 0 | hour(Floordt) == 4 |
                                   hour(Floordt) == 8 | hour(Floordt) == 12 |
                                   hour(Floordt) == 16 | hour(Floordt) == 20,])
    skinny <- as.data.frame(rbind(skinny, skinny2)) %>%
      arrange(ID, Floordt) %>%
      #  Make sure lat/long are in a numeric format
      mutate(
        Latitude = as.numeric(Latitude),
        Longitude = as.numeric(Longitude),
        StudyArea = ifelse(grepl("NE", ID), "NE", "OK")
      ) %>%
      #  Remove empty columns
      dplyr::select(-GPS.Horizontal.Error) %>% #c(Iridium.Latitude, Iridium.Longitude, GPS.Horizontal.Error)
      #  Thin data to retain only the 1st location on the hour
      #  Important when collar goes into mortality mode but animal is still alive- 
      #  Fix rate increases but flooring process puts all those times on the hour
      group_by(ID) %>%
      distinct(Floordt, .keep_all = TRUE) %>%   # .keep_all = TRUE saves all columns
      ungroup()
    return(skinny)
  }
  
  meso_skinny <- thin_locs(meso_clean)
  
  #  Gut check- did I drop too many locations?
  nrow(meso_clean) - nrow(meso_skinny)
  
  # #  Identify potentially problematic locations
  # meso_skinny <- meso_skinny %>%
  #   mutate(
  #     badAlt = ifelse(GPS.Altitude < 0 | GPS.Altitude > 3000, 1, 0),
  #     StudyArea = ifelse(grepl("NE", ID), "NE", "OK")
  #   )
  
  #  Drop oddball location that's clearly outside NEBOB33M home range
  meso_skinny <- meso_skinny[!(meso_skinny$ID == "NEBOB33M" & meso_skinny$Longitude > -117.4),]
  
  #  Exclude locations associated with translocation, dispersal, extra-territorial forays
  #  NEBOB13F: Translocation (3/28/19 - 8/8/19)
  meso_skinny <- meso_skinny[!(meso_skinny$ID == "NEBOB13F" & meso_skinny$Floordt < "2019-10-12 00:00:00"),]
  # #  MVBOB71M: Dispersal (9/22/19 - 10/11/19) and remains far outside study area
  # meso_skinny <- meso_skinny[!(meso_skinny$ID == "MVBOB71M" & meso_skinny$Floordt > "2019-09-22 00:00:00" & meso_skinny$Floordt < "2019-10-12 00:00:00"),]
  # # Cut all locations after dispersal since MVBOB71M moves so far outside study area
  meso_skinny <- meso_skinny[!(meso_skinny$ID == "MVBOB71M" & meso_skinny$Floordt > "2019-09-22 00:00:00"),]
  # #  MVBOB66M: Extra-territorial foray (11/28/19 - 3/6/20)
  # #  Appears to be transient in general (used different home ranges each season & year)
  # meso_skinny <- meso_skinny[!(meso_skinny$ID == "MVBOB66M" & meso_skinny$Floordt > "2019-11-28 00:00:00" & meso_skinny$Floordt < "2020-03-21 00:00:00"),]
  # #  NEBOB6F: Extra-territorial foray (3/8/19 - 3/31/19 & 1/15/20 - 3/5/20)
  # #  She does it both winters so I'm not sure it should be excluded- maybe normal extension of winter home range
  # meso_skinny <- meso_skinny[!(meso_skinny$ID == "NEBOB6F" & meso_skinny$Floordt > "2019-03-08 00:00:00" & meso_skinny$Floordt < "2019-04-01 00:00:00"),]
  # meso_skinny <- meso_skinny[!(meso_skinny$ID == "NEBOB6F" & meso_skinny$Floordt > "2020-01-15 00:00:00" & meso_skinny$Floordt < "2020-03-06 00:00:00"),]
  # #  Others to consider: MVBOB69F, NECOY20F, NECOY12F, MVBOB80M, MVBOB54F
  
  # # Rename data to acknowledge dispersals were excluded from this dataset
  # meso_skinny_noDispersal <- meso_skinny
  
  #  Save locations thinned to correct fix schedule
  # write.csv(meso_skinny, paste0('meso_skinny ', Sys.Date(), '.csv'))
  # write.csv(meso_skinny_noDispersal, paste0('meso_skinny_noDispersal', Sys.Date(), '.csv'))
  
  #  FYI: meso_skinny 2021-07-22.csv and 2021-11-12.csv excludes translocation and MVBOB71M's
  #  dispersal/new territory on the Colville but does not exclude other apparent
  #  extra-territorial movements or general transient behavior.


  
  
  ####  ============================================
  ####  Review individual collars for oddities  ####
  
  #  Plot all locations for a given species and look at their distribution
  #  Plot all locations for a given individual and look at their distribution
  #  Take special note of individuals with odd locations or distributions
  
  #  Load required packages for spatial work
  require(sf)
  require(ggplot2)
  
  #  Read in study area shapefiles and reproject to lat/long
  OK_SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "METHOW_SA")
  NE_SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "NE_SA")
  wgs84 <- st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  OK_wgs84 <- st_transform(OK_SA, wgs84)
  NE_wgs84 <- st_transform(NE_SA, wgs84)
  
  #  Make collar location data spatial
  meso_spdf <- st_as_sf(meso_skinny, coords = c("Longitude", "Latitude"), crs = wgs84) 


  #  Plot all locations (takes forever!)
  ggplot() +
    geom_sf(data = meso_spdf, aes(colour = ID)) #(colour = badAlt)

  
  #  Plot an individual collar
  ind_animal <- group_split(meso_spdf, meso_spdf$ID)
  
  #  Plot one animal
  ggplot() +
    geom_sf(data = NE_SA, fill = NA) +
    geom_sf(data = OK_SA, fill = NA) +
    geom_sf(data = ind_animal[[24]], aes(color = ID))
  
  #  Plotting by study area with study area boundary for context
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
  
  #  Feed meso data through function to map individual telemetry data in NE
  #  Plot an example map
  meso_NE_maps <- plot_telem_NE(meso_spdf[meso_spdf$StudyArea == "NE",])
  plot(meso_NE_maps[[24]])
  
  #  Feed meso through function to map individual telemetry data in OK
  #  Plot an example map
  meso_OK_maps <- plot_telem_OK(meso_spdf[meso_spdf$StudyArea == "OK",])
  plot(meso_OK_maps[[2]])
  
  
  #  Plotting zoomed in location data without study area boundary for context
  plot_telem <- function(spdf){
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
        geom_sf(data = ind_animal[[i]], aes(color = Floordt)) +  #badAlt
        labs(title = paste(names[i], "Locations", sep = " "), x = "Longitude",
             y = "Latitude") +
        #  Keep the legend but drop the title
        theme(legend.title = element_blank())
      plot_list[[i]] <- plot
      #print(plot)
    }
    return(plot_list)
  }
  
  #  Feed all species through function to map zoomed in telemetry data
  meso_maps <- plot_telem(meso_spdf)
  
  print(meso_maps[[24]])
  
  #  Save individual plots in a single pdf for each species
  #  With NE or OK study area boundary for context
  pdf("./Outputs/meso_NE_maps2.pdf")
  for (i in 1:length(unique(meso_NE_maps))) {
    print(meso_NE_maps[[i]])
  }
  dev.off()
  pdf("./Outputs/meso_OK_maps2.pdf")
  for (i in 1:length(unique(meso_OK_maps))) {
    print(meso_OK_maps[[i]])
  }
  dev.off()

  #  Without study area boundary for context
  pdf("./Outputs/meso_maps2.pdf")
  for (i in 1:length(unique(meso_maps))) {
    print(meso_maps[[i]])
  }
  dev.off()
 

  #  Fin
  #  Next step is Collar_Truncating&Filtering.R script
  
