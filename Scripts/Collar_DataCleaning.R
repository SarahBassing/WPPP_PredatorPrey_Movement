  ##  Ungulate Telemetry Data Cleaning
  ##  Washington Predator-Prey Project
  ##  Oct. 14, 2020
  ##  Sarah Bassing
  ##  =========================================================
  ##  Script to create MASTER & CLEANED data sets for WPPP ungulate data.
  ##    1. Combines WPPP GPS collar location data with unique animal IDs, 
  ##       and other capture/mortality information. 
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
  
  #  Capture data (latest download: 11.02.21; 11.16.20)
  #  Note the "ï.."
  md_cap <- read.csv("./Data/Capture (MD)_11.02.21.csv", stringsAsFactors = FALSE) %>% #_11.16.20
    mutate(IndividualIdentifier = ï..IndividualIdentifier) %>%
    select(-ï..IndividualIdentifier)
  elk_cap <- read.csv("./Data/Capture (Elk)_11.02.21.csv", stringsAsFactors = FALSE) %>%
    mutate(IndividualIdentifier = ï..IndividualIdentifier) %>%
    select(-ï..IndividualIdentifier)
  wtd_cap <- read.csv("./Data/Capture (WTD)_11.02.21.csv", stringsAsFactors = FALSE) %>%
    mutate(IndividualIdentifier = ï..IndividualIdentifier) %>%
    select(-ï..IndividualIdentifier) 
  
  #str(md_cap)#; head(md_cap)
  
  #  Mortality data (latest download: 11.02.21; 11.16.20)
  md_mort <- read.csv("./Data/Mortality (MD)_11.02.21.csv", stringsAsFactors = FALSE)  #_11.16.20
  elk_mort <- read.csv("./Data/Mortality (Elk)_11.02.21.csv", stringsAsFactors = FALSE) 
  wtd_mort <- read.csv("./Data/Mortality (WTD)_11.02.21.csv", stringsAsFactors = FALSE) 
  
  #str(md_mort)#; head(md_mort)
  
  #  Taylor combined capture & mortality data (last downloaded 11.02.21)
  md_capmort <- read.csv("./Data/WPPP_md_2021-07-12.csv", stringsAsFactors = FALSE) #2020-09-28
  elk_capmort <- read.csv("./Data/WPPP_elk_2021-11-2.csv", stringsAsFactors = FALSE) #2020-10-5
  wtd_capmort <- read.csv("./Data/WPPP_wtd_2021-08-10.csv", stringsAsFactors = FALSE) #2020-09-30
  
  #  Telemetry data (latest download: 11.02.21; 11.16.20)
  #  Add column with date/time in usable format & set timezone that data downloaded in (includes daylight savings)
  #  Add column that converts times to UTC
  #  Add column that converts times to Pacific Standard Time (add 8 hours), ignoring daylight savings
  #  Note: Etc/GMT+8 is UTC -8 and outputs in Pacific standard time only (so no daylight savings time)
  #  Add column that floors the time to the nearest hour not ahead
  #  Note the "ï.." with the latest download (11.02.21)
  md_tel <- read.csv("./Data/dev_telem_md_11.02.21.csv") %>%    # 11.16.20
    mutate(OBJECTID = ï..OBJECTID, # required for 11.02.21 data
           daytime = mdy_hms(ObservationDateTimePST, tz = "America/Los_Angeles"),
           UTCdt = with_tz(daytime, "UTC"),
           Finaldt = with_tz(UTCdt, tzone = "Etc/GMT+8"),
           Floordt = floor_date(Finaldt, unit = "hour")) %>%
    select(-c(ï..OBJECTID, InWashingtonJurisdiction)) # remove i..OBJECTID for 11.02.21 data
  elk_tel <- read.csv("./Data/dev_telem_elk_11.02.21.csv") %>%  
    mutate(OBJECTID = ï..OBJECTID,
           daytime = mdy_hms(ObservationDateTimePST, tz = "America/Los_Angeles"),
           UTCdt = with_tz(daytime, "UTC"),
           Finaldt = with_tz(UTCdt, tzone = "Etc/GMT+8"),
           Floordt = floor_date(Finaldt, unit = "hour")) %>%
    select(-c(ï..OBJECTID, InWashingtonJurisdiction))
  wtd_tel <- read.csv("./Data/dev_telem_wtd_11.02.21.csv") %>%
    mutate(OBJECTID = ï..OBJECTID,
           daytime = mdy_hms(ObservationDateTimePST, tz = "America/Los_Angeles"),
           UTCdt = with_tz(daytime, "UTC"),
           Finaldt = with_tz(UTCdt, tzone = "Etc/GMT+8"),  
           Floordt = floor_date(Finaldt, unit = "hour")) %>%
    select(-c(ï..OBJECTID, InWashingtonJurisdiction))
  no_fix <- read.csv("./Data/dev_telem_vec_nofix_11.02.21.csv") %>%
    #  Add extra columns to match the telemetry data with successful fixes
    mutate(Latitude = "NA",
           Longitude = "NA",
           ValidLocation = "NA",
           VEC_DOP = "NA",
           VEC_Height = "NA",
           OBJECTID = ï..OBJECTID, # required for 11.02.21 data
           daytime = mdy_hms(ObservationDateTimePST, tz = "America/Los_Angeles"),
           UTCdt = with_tz(daytime, "UTC"),
           Finaldt = with_tz(UTCdt, tzone = "Etc/GMT+8"),
           Floordt = floor_date(Finaldt, unit = "hour")) %>%
    select(-ï..OBJECTID) %>% # remove i..OBJECTID for 11.02.21 data
    #select(-OBJECTID) %>%
    #  Reorganize columns to match column order in telemetry data
    relocate(c("Latitude", "Longitude"), .after = CollarID) %>%
    relocate("ValidLocation", .after = DbLoadedDateTimePST) %>%
    relocate(c("VEC_DOP", "VEC_Height"), .after = VEC_Origin)
  #colnames(no_fix)[1] <- "OBJECTID"
  
  #  Extract species-specific missing fixes
  md_nofix <- no_fix[no_fix$Species == "Mule Deer",]
  elk_nofix <- no_fix[no_fix$Species == "Elk",]
  wtd_nofix <- no_fix[no_fix$Species == "White-tailed Deer",]
  
  #  Merge location and missed fixes data
  md_tel <- rbind(md_tel, md_nofix) %>%
    relocate(OBJECTID, .before = PositionID)
  elk_tel <- rbind(elk_tel, elk_nofix) %>%
    relocate(OBJECTID, .before = PositionID)
  wtd_tel <- rbind(wtd_tel, wtd_nofix) %>%
    relocate(OBJECTID, .before = PositionID)
  
  # #  Check out available timezones
  # OlsonNames()
  
  #  Check out the time data
  #  What's the timezone?
  tz(elk_tel$daytime); tz(elk_tel$UTCdt); tz(elk_tel$Finaldt)
  #  Does it account for daylight savings time?
  head(dst(elk_tel$daytime))
  head(dst(elk_tel$Finaldt))  # better say FALSE
  
  #str(elk_tel)#; head(elk_tel)
  
  ####  ======================================================
  ####  Combine Capture & Mortality data for each animal  ####
  
  spp_info <- function(cap, mort, capmort, tel) {
  #  Choose an end date if there is no mortality (chooses today's date)
    lastend <- ymd(Sys.Date())
  #  Combine capture and mortality data by unique animal ID               
    info <- full_join(cap, mort, by = "IndividualIdentifier") %>%
      #  Add in additional info about mortality & collar status
      right_join(capmort, by = c("IndividualIdentifier")) %>%
      transmute(
        IndividualIdentifier = as.factor(IndividualIdentifier),
        IndividualSpecies = IndividualSpecies.x,
        IndividualSex = IndividualSex.x,
        IndividualID = IndividualID,
        CaptureID =  CaptureID,
        LifeStage = LifeStage.x,
        CaptureDate = mdy(CaptureDate.x),
        GPSCollarSerialNumber = GPSCollarSerialNumber,
        MortalityID = ï..MortalityID,                        # note the "ï.."
        MortalityLifeStage = LifeStage.y,
        # EndDate = mdy(EstimatedMortalityDate),
        EndDate = mdy(EndDate),
        EndCause = EndCause,
        MortalityType = MortalityType.x,
        MortalitySubType = MortalitySubType.x,
        LastTransmission = mdy(LastTransmission),
        Notes = Notes
      )
  #  Fill in EndDate with chosen date if animal did not die/collar did not fail
    info$EndDate[is.na(info$EndDate)] <- lastend
  #  Remove individuals that don't have a corresponding GPSCollarSerialNumber
    which(is.na(info$GPSCollarSerialNumber))
    info <- droplevels(info[!is.na(info$GPSCollarSerialNumber),])
  #  Remove individuals that died on capture date
    deadcap <- as.character(info$IndividualIdentifier[which(info$CaptureDate == info$EndDate)])
    if(is_empty(deadcap) != TRUE) {
      info <- info[!(info$IndividualIdentifier %in% deadcap),]
    } else {
      info <- info
    }
  #  Identify any GPS collars on captured deer that never generated telemetry data
    notel <- info$GPSCollarSerialNumber[!(info$GPSCollarSerialNumber %in% tel$CollarID)]
    if(is_empty(notel) != TRUE) {
      info <- droplevels(info[info$GPSCollarSerialNumber != notel,])
    } else {
      info <- info
    }
    
  }
  
  #  Run function to merge capture and mortality data for individual animals
  md_info <- spp_info(md_cap, md_mort, md_capmort, md_tel)
  elk_info <- spp_info(elk_cap, elk_mort, elk_capmort, elk_tel)  
  wtd_info <- spp_info(wtd_cap, wtd_mort, wtd_capmort, wtd_tel)  
  
  #  Save for later use
  write.csv(md_info, paste0('./Data/md_info ', Sys.Date(), '.csv'))
  write.csv(elk_info, paste0('./Data/elk_info ', Sys.Date(), '.csv'))
  write.csv(wtd_info, paste0('./Data/wtd_info ', Sys.Date(), '.csv'))
  
  ####  Running list of collars that got nixed in this stage  ####
  #  Mule deer
  #  3969MD17 (male), 3958MD17 (capture-related mortality)
  
  #  Elk
  #  4836ELK20 (male)
  
  #  White-tailed deer
  #  019WTD17 (capture-related mortality), 023WTD17 (collar never transmitted), 
  #  70WTD18 (capture-related mortality), 85WTD19 (collar never transmitted), 
  #  90WTD19 (collar never transmitted)
  
  ####  =====================================================
  ####  Create MASTER data set for each ungulate species ####
  
  #  Connect location data to individual animal IDs:
  #  IDtelem function attaches animal IDs to their respective telemetry locations.
  #  This function truncates the data so that day of capture and day of mortality 
  #  (or today's locations) are excluded from the dataset. Further truncating can 
  #  happen for individual analyses.
    
  IDtelem <- function(info, telem) {
    #  Create empty dataframe to fill iteratively
    clean <- data.frame()
    #  How many individuals are looped over?
    nrow(info)
    #  Loop over every unique individual animal and...
    for(i in 1:nrow(info)){
      #  Take the individual animal ID
      ID <- droplevels(info$IndividualIdentifier[i])
      #  Take the animal's GPS collar serial number 
      SN <- info$GPSCollarSerialNumber[i]
      #  Buffer capture date to remove locations affected by capture event
      #  Suggested to only use data from 2 weeks after the capture data
      #  Currently not buffering but can change the +/- values in the future
      start <- info$CaptureDate[i] + 1
      #  Exclude locations 1 day before estimated mortality date
      end <- info$EndDate[i] - 1
      
      #  Subset telemetry data to the specific individual
      collar <- subset(telem, CollarID == SN)
      #  Add a new column to the telemetry data with the animal's individual ID
      collar$ID <- ID
      #  truncate telemetry data by new start and end dates for that individual
      collarlive <- subset(collar, Finaldt >= start & Finaldt <= end)
      
      #  Append each unique animal's locations to a clean dataframe
      clean <- rbind(clean, collarlive)
    }
    
    #  Organize by individual ID and chronological order of locations
    #  Format data fields
    clean <- clean %>%
      arrange(ID, Finaldt) %>%
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
        ID = ID
      )
    
    return(clean)
  }
  
  #  Run species-specific ID and telemetry data through the function
  md_master <- IDtelem(md_info, md_tel)
  elk_master <- IDtelem(elk_info, elk_tel)
  wtd_master <- IDtelem(wtd_info, wtd_tel)  
  #  IGNORE the warnings! Just saying that there are NAs in these columns 
  #  due to missed fixes- these NAs are expected
  
  #'  Remove spurious locations that snuck through the cleaning process
  md_master <- md_master %>%
    filter(IndividualName != "187MD20" | Latitude < 48.47) # mortality

  #  Pacific Standard Time for location data!

  #  Save MASTER data files
  write.csv(md_master, paste0('md_master ', Sys.Date(), '.csv'))
  write.csv(elk_master, paste0('elk_master ', Sys.Date(), '.csv'))
  write.csv(wtd_master, paste0('wtd_master ', Sys.Date(), '.csv'))
  
  
  ####  NEED TO ADD SOMETHING THAT PLOTS DISTANCE BTWN TIMES TO MAKE SURE I'M NOT MISSING LOCATIONS/ATTEMPTED FIXES
  
  
  ####  ====================================================
  ####  Create CLEAN data set for each ungulate species ####
  
  #  Further cleaning of telemetry data ready for analyses-
  #  Drop missing fixes and locations with poor accuracy
  #  Drop extra locations when fix schedule increases
  
  #  Drop missing fixes 
  #  Drop locations with poor accuracy for a cleaner data set
  #  No 2D locations and no locations with VEC_Height <0 or >2000 (arbitrary cutoffs)
  #  DOP >10 is not good (~20-30 meter accuracy)--- using DOP = 8 as cutoff
  #  From Vectronic: Dilution of precision (DOP); "High values indicate that the 
  #  location is likely to be inaccurate while low values indicate a better precision."
  md_clean <- md_master %>%
    filter(VEC_FixType != "No fix") %>%
    filter(VEC_FixType != "GPS-2D") %>%
    #filter(VEC_Height < 2000 & VEC_Height > 0) %>%
    filter(VEC_DOP <= 8)
  elk_clean <- elk_master %>%
    filter(VEC_FixType != "No fix") %>%
    filter(VEC_FixType != "GPS-2D") %>%
    filter(VEC_Height < 2000 & VEC_Height > 0) %>%
    filter(VEC_DOP <= 8)
  wtd_clean <- wtd_master %>%
    filter(VEC_FixType != "No fix") %>%
    filter(VEC_FixType != "GPS-2D") %>%
    filter(VEC_Height < 3000 & VEC_Height > 0) %>%
    filter(VEC_DOP <= 8) 
  
  #  VEH_Height also associated with accuracy so need to see if there are any
  #  odd outliers that should be excluded due to low accuracy fixes
  hist(md_clean$VEC_Height, breaks = c(100), main = "Mule Deer locations (ALL VEC_Height values)")
  # plot(md_clean$VEC_Height)  # takes awhile
  summary(md_clean$VEC_Height)
  oddball_md <- md_clean[md_clean$VEC_Height >= 2000 | md_clean$VEC_Height < 0,]
  hist(elk_clean$VEC_Height, breaks = c(100), main = "Elk locations (ALL VEC_Height values)")
  # plot(elk_clean$VEC_Height)  # takes awhile
  summary(elk_clean$VEC_Height)
  oddball_elk <- elk_clean[elk_clean$VEC_Height >= 2000 | elk_clean$VEC_Height < 0,]
  hist(wtd_clean$VEC_Height, breaks = c(100), main = "WTD locations (ALL VEC_Height values)")
  # plot(wtd_clean$VEC_Height)  # takes awhile
  summary(wtd_clean$VEC_Height)
  oddball_wtd <- wtd_clean[wtd_clean$VEC_Height >= 2000 | wtd_clean$VEC_Height < 0,]
  
  #  Save locations with high accuracy fixes only
  write.csv(md_clean, paste0('md_clean ', Sys.Date(), '.csv'))
  write.csv(elk_clean, paste0('elk_clean ', Sys.Date(), '.csv'))
  write.csv(wtd_clean, paste0('wtd_clean ', Sys.Date(), '.csv'))
  
  
  ####  ==========================
  ####  Rarefy location data  ####
  #  KEEP locations with 4-hr fix schedule with 2, 6, 10, 14, 18, 22 hr fixes
  #  DISCARD: 
  #    -Extra locations from summer elk collars & collars in mortality-mode   
  #     where collars increase fix schedule to 2-hr or 20-min. interval, respectfully.
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
        Longitude = as.numeric(Longitude)
      ) %>%
      #  Thin data to retain only the 1st location on the hour
      #  Important when collar goes into mortality mode but animal is still alive- 
      #  Fix rate increases but flooring process puts all those times on the hour
      group_by(ID) %>%
      distinct(Floordt, .keep_all = TRUE) %>%   # .keep_all = TRUE saves all columns
      ungroup()
    return(skinny)
  }
  
  md_skinny <- thin_locs(md_clean)
  elk_skinny <- thin_locs(elk_clean)
  wtd_skinny <- thin_locs(wtd_clean)

  #  Gut check- did I drop too many locations?
  nrow(md_clean) - nrow(md_skinny)
  nrow(elk_clean) - nrow(elk_skinny)
  nrow(wtd_clean) - nrow(wtd_skinny)

  
  # #  Is the hour filtering really working? (this ignores really weird times)
  # # wtd_wrong <- with(wtd_clean, wtd_clean[hour(Floordt) == 0 | hour(Floordt) == 4 | 
  # #                            hour(Floordt) == 8 | hour(Floordt) == 12 | 
  # #                            hour(Floordt) == 16 | hour(Floordt) == 20,])
  # not_sched <- with(wtd_clean, wtd_clean[hour(Floordt) != 2 & hour(Floordt) != 6 &
  #                                          hour(Floordt) != 10 & hour(Floordt) != 14 &
  #                                          hour(Floordt) != 18 & hour(Floordt) != 22,])
  # #  Is it all collars or just some?
  # length(unique(not_sched$ID)); length(unique(wtd_clean$ID))  # apparently just some
  # #  Is it just collars that failed or went into mort-mode?
  # drop_locs <- as.data.frame(droplevels(unique(not_sched$ID))) %>%
  #   mutate(Wrong_schedule = "Wrong Sched")  # collars on wrong schedule for any reason
  # colnames(drop_locs) <- "ID"
  # fail_deer <- drop_na(wtd_info, EndCause)  # collars that failed or animal died
  # fail_deer <- droplevels(unique(fail_deer$IndividualIdentifier))
  # fail_deer <- as.data.frame(fail_deer) %>%
  #   mutate(Fail_Dead = "Fail/Dead")
  # colnames(fail_deer) <- c("ID", "Fail/Dead")
  # wtf <- full_join(drop_locs, fail_deer, by = "ID") # not all wrong schedules are from collars that died/failed
  
  
  #  Save locations thinned to correct fix schedule
  write.csv(md_skinny, paste0('md_skinny ', Sys.Date(), '.csv'))
  write.csv(elk_skinny, paste0('elk_skinny ', Sys.Date(), '.csv'))
  write.csv(wtd_skinny, paste0('wtd_skinny ', Sys.Date(), '.csv'))
  
  
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
  md_spdf <- st_as_sf(md_skinny, coords = c("Longitude", "Latitude"), crs = wgs84)
  elk_spdf <- st_as_sf(elk_skinny, coords = c("Longitude", "Latitude"), crs = wgs84)
  wtd_spdf <- st_as_sf(wtd_skinny, coords = c("Longitude", "Latitude"), crs = wgs84)

  #  Plot all locations (takes forever!)
  # ggplot() +
  #   geom_sf(data = md_spdf, aes(colour = ID))
  # ggplot() +
  #   geom_sf(data = elk_spdf, aes(colour = ID))
  # ggplot()+
  #   geom_sf(data = wtd_spdf, aes(colour = ID))

  #  Plot an individual collar
  ind_animal <- group_split(elk_spdf, elk_spdf$ID)
  #  Plot one animal
  ggplot() +
    geom_sf(data = NE_SA, fill = NA) +
    geom_sf(data = ind_animal[[1]], aes(color = ID))


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

  #  Feed elk & wtd through function to map individual telemetry data in NE
  #  Plot an example map
  elk_NE_maps <- plot_telem_NE(elk_spdf)
  plot(elk_NE_maps[[1]])

  wtd_NE_maps <- plot_telem_NE(wtd_spdf)
  plot(wtd_NE_maps[[1]])

  #  Feed mule deer through function to map individual telemetry data in OK
  #  Plot an example map
  md_OK_maps <- plot_telem_OK(md_spdf)
  plot(md_OK_maps[[1]])


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
        geom_sf(data = ind_animal[[i]], aes(color = Floordt)) +
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
  elk_maps <- plot_telem(elk_spdf)
  wtd_maps <- plot_telem(wtd_spdf)
  md_maps <- plot_telem(md_spdf)

  print(elk_maps[[1]])

  # # Save individual plots in a single pdf for each species
  # # With NE or OK study area boundary for context
  # pdf("./Outputs/elk_NE_maps2.pdf")
  # for (i in 1:length(unique(elk_NE_maps))) {
  #   print(elk_NE_maps[[i]])
  # }
  # dev.off()
  # pdf("./Outputs/wtd_NE_maps2.pdf")
  # for (i in 1:length(unique(wtd_NE_maps))) {
  #   print(wtd_NE_maps[[i]])
  # }
  # dev.off()
  # pdf("./Outputs/md_OK_maps2.pdf")
  # for (i in 1:length(unique(md_OK_maps))) {
  #   print(md_OK_maps[[i]])
  # }
  # dev.off()
  # 
  # #  Without study area boundary for context
  # pdf("./Outputs/elk_maps2.pdf")
  # for (i in 1:length(unique(elk_maps))) {
  #   print(elk_maps[[i]])
  # }
  # dev.off()
  # pdf("./Outputs/wtd_maps.pdf")
  # for (i in 1:length(unique(wtd_maps))) {
  #   print(wtd_maps[[i]])
  # }
  # dev.off()
  # pdf("./Outputs/md_maps.pdf")
  # for (i in 1:length(unique(md_maps))) {
  #   print(md_maps[[i]])
  # }
  # dev.off()
  
  #  NOTES
  #  Many mule deer collars still include at least one start point in Winthrop (hotel?)
  #  Handful of wtd, 1 elk, and a couple of mule deer have very few points- probably want to exclude
  
  
  #  Double check mule deer collars for animals that were relocated during collaring
  #  Do these collars need to be truncated differently than others?
  #  Double check individual animals whose collars have a lot of missing data
  
  
  # #  Probably not needed but just in case...
  # #  Take a closer look at animals with odd locations
  # #  Questionable points or walk-abouts
  # elk_review <- filter(elk_spdf, ID == "4830ELK20" | elk_spdf$ID =="3974ELK18" |
  #                      ID == "3712ELK18" | ID == "3686EA17" | ID == "3676EA17")
  # review <- group_split(elk_review, elk_review$ID)
  # 
  # #  Loop through each funky elk and create a plot of its locations
  # names <- c()
  # for(i in 1:length(unique(review))) {
  #   names <- c(names, unique(as.character(review[[i]]$ID)))
  #   plot <- ggplot() +
  #     #geom_sf(data = NE_SA, fill = NA) +
  #     geom_sf(data = review[[i]], aes(color = Floordt)) + #aes(color = ID)
  #     labs(title = paste(names[i], "Locations", sep = " "), x = "Longitude",
  #          y = "Latitude") +
  #     theme(legend.title = element_blank()) 
  #   plot_list[[i]] <- plot
  # }
  # print(plot_list[[2]])
  # 
  # #  Save individual plots in a single pdf
  # pdf("./Outputs/odd_collars_review.pdf")
  # for (i in 1:length(unique(review))) {
  #   print(plot_list[[i]])
  # }
  # dev.off()
  
  
  ####  =============================================
  ####  Summary Stats on Fix Success & Accuracy  ####
  
  #  Function to calculate summary stats on telemetry data
  fix_stats <- function(master, nofix, clean) {
    #  Total number of possible locations (included missed fixes)
    total <- nrow(master)
    #  Number of missed locations
    missed <- nrow(nofix)
    #  Percent of missed fixes out of total possible locations
    perc_missed <- nrow(nofix)/nrow(master)

    #  Total number of successful fixes
    locations <- nrow(clean)
    #  Number of low accuracy fixes (DOP > 10)
    lowacc <- nrow(filter(master, VEC_DOP > 10))
    #  Percent of low accuracy fixes out of all successful locations
    perc_lowacc <- lowacc/nrow(clean)

    #  Combine all summary stats
    stats <- c(total, missed, perc_missed, locations, lowacc, perc_lowacc)

    return(stats)

  }

  #  Run species location data through summary stats function
  md_stats <- fix_stats(md_master, md_nofix, md_clean)
  elk_stats <- fix_stats(elk_master, elk_nofix, elk_clean)
  wtd_stats <- fix_stats(wtd_master, wtd_nofix, wtd_clean)

  telem_stats <- as.data.frame(rbind(md_stats, elk_stats, wtd_stats))
  colnames(telem_stats) <- c("Total_Locs", "Missed_Locs", "Perc_Missed", "Successful_Locs", "LowAccuracy_Locs", "Perc_LowAccuracy")
  print(telem_stats)

  #  Collar-specific issues
  #  Function to calculate missing data and low accuracy fixes per animal
  lame_locs <- function(master) {
    NoFix <- master %>%
      group_by(ID) %>%
      count(VEC_FixType == "No fix") %>%
      ungroup() %>%
      #  Spread data for easier summarizing below (wide format)
      spread(ID, n)
    colnames(NoFix)[1] <- "Problem"

    InaccFix <- master %>%
      group_by(ID) %>%
      na.omit() %>%
      count(VEC_DOP > 10) %>%
      ungroup() %>%
      #  Spread data for easier summarizing below (wide format)
      spread(ID, n)
    colnames(InaccFix)[1] <- "Problem"

    Collar_probs <- rbind(NoFix, InaccFix)
    Problem_Type <- c("Missed Fix", "Missed Fix", "Low Accuracy", "Low Accuracy")
    Collar_probs <- cbind(Problem_Type, Collar_probs)
    #  Gather wide data back into long-format w/ NAs where no missing data occurred
    Collar_problems <- gather(Collar_probs, "Animal_ID", "n", 3:ncol(Collar_probs)) %>%
      #  Reorder so easier to read
      transmute(
        Animal_ID = Animal_ID,
        Nmbr_Locations = n,
        Problem_Type = Problem_Type,
        Problem = Problem
      )

    return(Collar_problems)

  }

  #  Run species-specific data through function and combine
  md_probs <- lame_locs(md_master) %>%
    mutate(Species = "Mule Deer")
  elk_probs <- lame_locs(elk_master) %>%
    mutate(Species = "Elk")
  wtd_probs <- lame_locs(wtd_master) %>%
    mutate(Species = "White-tailed Deer")
  Problem_locs <- rbind(md_probs, elk_probs, wtd_probs) %>%
    relocate(Species)  # default puts this column first

  #  Save summary stats
  write.csv(telem_stats, paste0('Telemetry_stats ', Sys.Date(), '.csv'))
  write.csv(Problem_locs, paste0('Collar_stats ', Sys.Date(), '.csv'))
  
  
  #  Fin
  #  Next step is Collar_Truncating&Filtering.R script

  
  
  # ##  =====================================================  
  # ####  Broken down species by species  ####
  # 
  # #  This version is a little easier for trouble-shooting when the two functions
  # #  aren't working thanks to weird random collar-telemetry mismatches, etc.
  # 
  # 
  # #  Choose an end date if there is no mortality- today
  # lastend <- ymd(Sys.Date())
  # 
  # 
  # ####  Mule deer info ####
  # #  Combine capture and mortality data by unique animal ID
  # md_info <- full_join(md_cap, md_mort, by = "IndividualIdentifier") %>%
  #   right_join(md_capmort, by = c("IndividualIdentifier")) %>%
  #   transmute(
  #     IndividualIdentifier = as.factor(IndividualIdentifier),
  #     IndividualSpecies = IndividualSpecies.x,
  #     IndividualSex = IndividualSex.x,
  #     IndividualID = IndividualID,
  #     CaptureID =  CaptureID,
  #     LifeStage = LifeStage.x,
  #     CaptureDate = mdy(CaptureDate.x),
  #     GPSCollarSerialNumber = GPSCollarSerialNumber,
  #     MortalityID = MortalityID,
  #     MortalityLifeStage = LifeStage.y,
  #     EndDate = mdy(EndDate),
  #     EndCause = EndCause,
  #     MortalityType = MortalityType.x,
  #     MortalitySubType = MortalitySubType.x,
  #     LastTransmission = mdy(LastTransmission),
  #     Notes = Notes
  #   )
  # #  Make sure you got them all
  # length(unique(md_cap$IndividualIdentifier))
  # length(unique(md_info$IndividualIdentifier))
  # 
  # #  Assign mortality date as end date if animal died
  # #md_info$enddate <- mdy(md_info$EstimatedMortalityDate)
  # #  Fill in other end dates with chosen date if animal did not die
  # #md_info$enddate[is.na(md_info$enddate)] <- lastend
  # md_info$EndDate[is.na(md_info$EndDate)] <- lastend
  # 
  # #  Check to make sure each unique individual has location data
  # dat <- as.data.frame(md_info$IndividualIdentifier)
  # dat <- cbind(dat, md_info$GPSCollarSerialNumber)
  # view(dat)  # look for NAs and remove those individuals
  # 
  # #  Remove individuals that don't have a corresponding GPSCollarSerialNumber
  # #  Mule deer: 3969MD17 , 3958MD17
  # which(is.na(md_info$GPSCollarSerialNumber))
  # length(unique(which(is.na(md_info$GPSCollarSerialNumber))))  # number of deer w/o collars
  # md_info <- droplevels(md_info[!is.na(md_info$GPSCollarSerialNumber),])
  # 
  # #  Remove individuals that died on capture date
  # deadcap <- as.character(md_info$IndividualIdentifier[which(md_info$CaptureDate == md_info$EstimatedMortalityDate)])
  # print(deadcap)
  # if(is_empty(deadcap) != TRUE) {
  #   md_info <- md_info[md_info$IndividualIdentifier != deadcap,]
  # } else {
  #   md_info <- md_info
  # }
  # #  Remove this one individual that for some reason ruins the function... bad coding practice
  # #  89MD18 and R89ND18 have the same collar ID b/c 89MD18 died 4 days after capture
  # md_info <- md_info[md_info$IndividualIdentifier != "89MD18",]
  # 
  # #  Identify any GPS collars on captured deer that never generated telemetry data
  # notel <- md_info$GPSCollarSerialNumber[!(md_info$GPSCollarSerialNumber %in% md_tel$CollarID)]
  # print(notel)
  # if(is_empty(notel) != TRUE) {
  #   md_info <- droplevels(md_info[md_info$GPSCollarSerialNumber != notel,])
  # } else {
  #   md_info <- md_info
  # }
  # 
  # 
  # ####  Elk info  ####
  # # Combine capture and mortality data by unique animal ID
  # elk_info <- full_join(elk_cap, elk_mort, by = "IndividualIdentifier") %>%
  #   right_join(elk_capmort, by = c("IndividualIdentifier")) %>%
  #   transmute(
  #     IndividualIdentifier = as.factor(IndividualIdentifier),
  #     IndividualSpecies = IndividualSpecies.x,
  #     IndividualSex = IndividualSex.x,
  #     IndividualID = IndividualID,
  #     CaptureID =  CaptureID,
  #     LifeStage = LifeStage.x,
  #     CaptureDate = mdy(CaptureDate.x),
  #     GPSCollarSerialNumber = GPSCollarSerialNumber,
  #     MortalityID = MortalityID,
  #     MortalityLifeStage = LifeStage.y,
  #     EndDate = mdy(EndDate),
  #     EndCause = EndCause,
  #     MortalityType = MortalityType.x,
  #     MortalitySubType = MortalitySubType.x,
  #     LastTransmission = mdy(LastTransmission),
  #     Notes = Notes
  #   )
  # #  Make sure you got them all
  # length(unique(elk_cap$IndividualIdentifier))
  # length(unique(elk_info$IndividualIdentifier))
  # head(elk_info)
  # 
  # #  Assign mortality date as end date if animal died
  # #elk_info$enddate <- mdy(elk_info$EstimatedMortalityDate)
  # #  Fill in other end dates with chosen date if animal did not die
  # lastend <- ymd(Sys.Date())
  # #elk_info$enddate[is.na(elk_info$enddate)] <- lastend
  # elk_info$EndDate[is.na(elk_info$EndDate)] <- lastend
  # 
  # #  Check to make sure each unique individual has location data
  # dat <- as.data.frame(elk_info$IndividualIdentifier)
  # dat <- cbind(dat, elk_info$GPSCollarSerialNumber)
  # view(dat)
  # 
  # #  Remove individuals that don't have a corresponding GPSCollarSerialNumber
  # #  Elk: 4836ELK20
  # which(is.na(elk_info$GPSCollarSerialNumber))
  # length(unique(which(is.na(elk_info$GPSCollarSerialNumber))))
  # elk_info <- droplevels(elk_info[!is.na(elk_info$GPSCollarSerialNumber),])
  # 
  # #  Remove individuals that died on capture date
  # deadcap <- as.character(elk_info$IndividualIdentifier[which(elk_info$CaptureDate == elk_info$EstimatedMortalityDate)])
  # print(deadcap)
  # if(is_empty(deadcap) != TRUE) {
  #   elk_info <- elk_info[elk_info$IndividualIdentifier != deadcap,]
  # } else {
  #   elk_info <- elk_info
  # }
  # 
  # #  Identify any GPS collars on captured deer that never generated telemetry data
  # notel <- elk_info$GPSCollarSerialNumber[!(elk_info$GPSCollarSerialNumber %in% elk_tel$CollarID)]
  # print(notel)
  # if(is_empty(notel) != TRUE) {
  #   elk_info <- droplevels(elk_info[elk_info$GPSCollarSerialNumber != notel,])
  # } else {
  #   elk_info <- elk_info
  # }
  # 
  # 
  # ####  White-tailed deer info  ####
  # #  Combine capture and mortality data by unique animal ID
  # wtd_info <- full_join(wtd_cap, wtd_mort, by = "IndividualIdentifier") %>%
  #   right_join(wtd_capmort, by = c("IndividualIdentifier")) %>%
  #   transmute(
  #     IndividualIdentifier = as.factor(IndividualIdentifier),
  #     IndividualSpecies = IndividualSpecies.x,
  #     IndividualSex = IndividualSex.x,
  #     IndividualID = IndividualID,
  #     CaptureID =  CaptureID,
  #     LifeStage = LifeStage.x,
  #     CaptureDate = mdy(CaptureDate.x),
  #     GPSCollarSerialNumber = GPSCollarSerialNumber,
  #     MortalityID = MortalityID,
  #     MortalityLifeStage = LifeStage.y,
  #     EndDate = mdy(EndDate),
  #     EndCause = EndCause,
  #     MortalityType = MortalityType.x,
  #     MortalitySubType = MortalitySubType.x,
  #     LastTransmission = mdy(LastTransmission),
  #     Notes = Notes
  #   )
  # #  Make sure you got them all
  # length(unique(wtd_cap$IndividualIdentifier))
  # length(unique(wtd_info$IndividualIdentifier))
  # 
  # #  Assign mortality date as end date if animal died
  # #wtd_info$enddate <- mdy(wtd_info$EstimatedMortalityDate)
  # #  Fill in other end dates with chosen date if animal did not die
  # lastend <- ymd(Sys.Date())
  # #wtd_info$enddate[is.na(wtd_info$enddate)] <- lastend
  # wtd_info$EndDate[is.na(wtd_info$EndDate)] <- lastend
  # 
  # #  Check to make sure each unique individual has location data
  # dat <- as.data.frame(wtd_info$IndividualIdentifier)
  # dat <- cbind(dat, wtd_info$GPSCollarSerialNumber)
  # view(dat)
  # 
  # #  Remove individuals that don't have a corresponding GPSCollarSerialNumber
  # which(is.na(wtd_info$GPSCollarSerialNumber))
  # length(unique(which(is.na(wtd_info$GPSCollarSerialNumber))))  # number of deer w/o collars
  # wtd_info <- droplevels(wtd_info[!is.na(wtd_info$GPSCollarSerialNumber),])
  # 
  # #  Remove individuals that died on capture date
  # #  019WTD17, 70WTD18, etc.
  # deadcap <- as.character(wtd_info$IndividualIdentifier[which(wtd_info$CaptureDate == wtd_info$EndDate)])
  # if(is_empty(deadcap) != TRUE) {
  #   wtd_info <- wtd_info[!(wtd_info$IndividualIdentifier %in% deadcap),]
  #   #wtd_info <- wtd_info[wtd_info$IndividualIdentifier != deadcap,]
  # } else {
  #   wtd_info <- wtd_info
  # }
  # 
  # #  Remove individuals with collars that are not in the telemetry data
  # #  Identify mismatches between GPS collars in wtd_info vs wtd_tel
  # notel <- wtd_info$GPSCollarSerialNumber[!(wtd_info$GPSCollarSerialNumber %in% wtd_tel$CollarID)]
  # print(notel)
  # #  Collars 24833 (90WTD19), 24867 (85WTD19)
  # wtd_info <- droplevels(wtd_info[wtd_info$GPSCollarSerialNumber != notel,])
  # 
  # 
  # #### Merge info & telem for 1 individual  ####
  # # If the IDtelem function is too much, try it for a single individual!
  # 
  # #  Create empty dataframe to fill iteritively
  # MDclean <- data.frame()
  # #  How many individuals are looped over?
  # nrow(elk_info) #md_info
  # #  Loop over every unique individual animal and...
  # for(i in 1:nrow(elk_info)){ #md_info      # DON'T FOR LOOP IT IF ONLY TESTING 1 INDIVIDUAL
  #   #  Take the individual animal ID
  #   mdID <- droplevels(elk_info$IndividualIdentifier[1]) #md_info
  #   #  Take the animal's GPS collar serial number
  #   mdSN <- elk_info$GPSCollarSerialNumber[1] #md_info
  #   #  Buffer capture date to remove locations affected by capture event
  #   #  Suggested to only use data from 2 weeks after the capture data (some papers suggest 1 month)
  #   start <- elk_info$CaptureID[1] #+ 14
  #   #  Exclude locations 1 day before estimated mortality date
  #   end <- elk_info$EndDate[i] #- 1
  # 
  #   #  Subset telemetry data to the specific individual
  #   md <- subset(elk_nofix, CollarID == mdSN) #md_tel
  #   #  Add a new column to the telemetry data with the animal's individual ID
  #   md$ID <- mdID
  #   #  Truncate telemetry data by new start and end dates for that individual
  #   #  Important for collars that are redeployed- ensures locations generated
  #   #  by that specific animal are included, even if collar generates more locations
  #   #  on another animal
  #   mdlive <- subset(md, Finaldt >= start & Finaldt <= end)
  # 
  #   #  Append each unique animal's locations to a clean dataframe
  #   MDclean <- rbind(MDclean, mdlive)
  # }
  # 
  # length(unique(wtd_cap$IndividualIdentifier))
  # length(unique(wtd_info$IndividualIdentifier))
  # length(unique(MDclean$ID))
  # #  FYI 89MD18 & 24MD18 are dropped in when locations are truncated b/c
  # #  animals died within 2 weeks of capture
  # 
  # #  Organize by individual ID and chronological order of locations
  # MDclean <- MDclean %>%
  #   arrange(ID, daytime) #%>%
  #   #  Filter out aberrant locations
  #   # filter(flgLocation != 1) %>%
  #   # filter(flgDate != 1) %>%
  # # flgLocation == 1 indicate inaccurate fixes
  # # flgDate == 1 indicate dates in the future (only issue for Telonics collars)
  # # flgJurisdiction == 1 indicates collar outside WA State jurisdiction (e.g., Tribal land, Canada)
  # # flgActive == 0 indicates locations where the animal that generated those locations is no longer alive (do NOT filter out 0's here)
  # 