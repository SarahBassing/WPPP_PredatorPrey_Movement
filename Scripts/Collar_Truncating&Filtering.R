  ##  Final collar cleaning steps
  ##  Washington Predator-Prey Project
  ##  Nov. 16, 2020
  ##  Sarah Bassing
  ##  =========================================================
  ##  Script takes cleaned master GPS satellite data and takes final  steps to 
  ##  truncate and filter telemetry data for analyses specific to my project. 
  ##  Carnivore data generally excludes large dispersal events that take them
  ##  well outside the extent of the WPPP study areas.
  ##     1. Truncating 3-wks after animal was captured to ensure any movements 
  ##        affected by the capture are excluded from analyses.
  ##     2. Thinning data to only include locations on the WPPP chosen 4-hr fix
  ##        schedule. Location times should be: 2:00, 6:00, 10:00, 14:00, 18:00, 
  ##        & 22:00 for all individuals.
  ##     3. Remove any individuals that have very few locations or are missing
  ##        a lot of data due to collar malfunctions or previous filtering.
  ##     4. Remove all locations associated with mule deer migration times. 
  ##        Migration dates identified by T.Ganz.
  ##  This should produce the final data set to be used with HMMs & RSFs.
  ##
  ##  Cleaned data used below were generated in the Collar_DataCleaning.R and
  ##  Collar_DataCleaning_Mesopredators.R scripts. Cougar & wolf data were cleaned
  ##  by L.Satterfield but were reviewed with Collar_DataCleaning_CougarWolf.R
  ##  =========================================================
  
  #  Clean work space and load libraries
  rm(list = ls())

  library(lubridate)
  library(tidyverse)
  
  #  Turn off scientific notation
  options(scipen = 999) 
  #  Set digits to 15 to ensure GPS coordinates aren't truncated
  options(digits = 15) 
  
  #  Read in data
  #  Make sure specific columns are formatted correctly for data manipulation
  md_info <- read.csv("./Data/md_info 2021-11-09.csv") %>% #2020-11-17
    mutate(
      IndividualIdentifier = as.factor(as.character(IndividualIdentifier)),
      CaptureDate = as_date(CaptureDate)
      ) %>%
    dplyr::select(-X)
  elk_info <- read.csv("./Data/elk_info 2021-11-09.csv") %>%
    mutate(
      IndividualIdentifier = as.factor(as.character(IndividualIdentifier)),
      CaptureDate = as_date(CaptureDate)
    ) %>%
    dplyr::select(-X)
  wtd_info <- read.csv("./Data/wtd_info 2021-11-09.csv") %>%
    mutate(
      IndividualIdentifier = as.factor(as.character(IndividualIdentifier)),
      CaptureDate = as_date(CaptureDate)
    ) %>%
    dplyr::select(-X)
  #  Created based on data provided by L.Satterfield
  cougwolf_info <- read.csv("./Data/cougwolf_info.csv") %>%
    mutate(
      IndividualIdentifier = as.factor(as.character(Animal_ID)),
      GPSCollarSerialNumber = Collar,
      CaptureDate = mdy(AnimalStart, tz = "America/Los_Angeles"),
      CaptureDate = as_date(CaptureDate),
      EndDate = mdy(DataEnd, tz = "America/Los_Angeles"),
      EndDate = as_date(EndDate)
    ) %>%
    dplyr::select(-c(Animal_ID, Collar, AnimalStart, DataEnd)) %>%
    relocate(c(IndividualIdentifier, GPSCollarSerialNumber, CaptureDate, EndDate), .before = StudyArea)
  coug_info <- droplevels(filter(cougwolf_info, Species == "Cougar")) %>%
    #  Remove cougar with no telemetry data
    filter(IndividualIdentifier != "MVC203M")
  wolf_info <- droplevels(filter(cougwolf_info, Species == "Wolf")) %>%
    #  Add "W" to start of AnimalID to match IDs in location data
    mutate(
      IndividualIdentifier = as.factor(as.character(paste0("W", IndividualIdentifier)))
    ) %>%
    #  Remove wolf with no telemetry data
    # filter(IndividualIdentifier != "W110M") %>%
    #'  Remove wolves on very different fix schedule (6hr & 8hr fix schedules)
    #'  L.Satterfield removed these from telemetry data
    filter(IndividualIdentifier != "W51F") %>%
    filter(IndividualIdentifier != "W68M") %>%
    filter(IndividualIdentifier != "W81M") %>%
    filter(IndividualIdentifier != "W85M") %>%
    filter(IndividualIdentifier != "W98F") #%>%
    # filter(GPSCollarSerialNumber != "12438")
  #  Created based on data provided by B.Windell
  meso_info <- read.csv("./Data/meso_info 2021-11-11.csv") %>%
    mutate(
      IndividualIdentifier = as.factor(as.character(IndividualIdentifier)),
      CaptureDate = as_date(CaptureDate),
      EndDate = as_date(EndDate)
    ) %>%
    dplyr::select(-X)
  
  #  Note: WDFW WebApp allowed timezone to shift to PDT so must adjust to only PST
  md_skinny <- read.csv("md_skinny 2021-11-09.csv") %>%
    mutate(StudyArea = "OK",
           daytime = mdy_hms(ObservationDateTimePST, tz = "America/Los_Angeles"),
           UTCdt = with_tz(daytime, "UTC"),
           Finaldt = with_tz(UTCdt, tzone = "Etc/GMT+8"),
           Floordt = floor_date(Finaldt, unit = "hour")) %>%
    dplyr::select("OBJECTID", "ID", "CollarID", "Species", "Sex", "Latitude", "Longitude", 
                  "ObservationDateTimePST", "StudyArea", "daytime", "UTCdt", "Finaldt", "Floordt")
  #  Note: WDFW WebApp allowed timezone to shift to PDT so must adjust to only PST
  elk_skinny <- read.csv("elk_skinny 2021-11-09.csv") %>%
    mutate(StudyArea = "NE",
           daytime = mdy_hms(ObservationDateTimePST, tz = "America/Los_Angeles"),
           UTCdt = with_tz(daytime, "UTC"),
           Finaldt = with_tz(UTCdt, tzone = "Etc/GMT+8"),
           Floordt = floor_date(Finaldt, unit = "hour")) %>%
    dplyr::select("OBJECTID", "ID", "CollarID", "Species", "Sex", "Latitude", "Longitude", 
                  "ObservationDateTimePST", "StudyArea", "daytime", "UTCdt", "Finaldt", "Floordt")
  #  Note: WDFW WebApp allowed timezone to shift to PDT so must adjust to only PST
  wtd_skinny <- read.csv("wtd_skinny 2021-11-09.csv") %>%
    mutate(StudyArea = "NE",
           daytime = mdy_hms(ObservationDateTimePST, tz = "America/Los_Angeles"),
           UTCdt = with_tz(daytime, "UTC"),
           Finaldt = with_tz(UTCdt, tzone = "Etc/GMT+8"),
           Floordt = floor_date(Finaldt, unit = "hour")) %>%
    dplyr::select("OBJECTID", "ID", "CollarID", "Species", "Sex", "Latitude", "Longitude", 
                  "ObservationDateTimePST", "StudyArea", "daytime", "UTCdt", "Finaldt", "Floordt")
  #  Note: data in UTC timezone to begin with
  #  Obvious dispersal events are excluded but large seasonal movements within home range are included
  meso_skinny <- read.csv("meso_skinny 2021-11-12.csv") %>%  
    mutate(
      StudyArea = ifelse(grepl("NE", ID), "NE", "OK"),
      Species = ifelse(Species == "BOB", "Bobcat", Species),
      Species = ifelse(Species == "COY", "Coyote", Species),
      daytime = as.POSIXct(daytime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
      UTCdt = with_tz(daytime, "UTC"),
      Finaldt = with_tz(UTCdt, tzone = "Etc/GMT+8"),
      Floordt = floor_date(Finaldt, unit = "hour")) %>%
    dplyr::select("ID", "CollarID", "Species", "Sex", "StudyArea", "Latitude", "Longitude", 
                  "Acquisition.Time.UTC", "StudyArea", "daytime", "UTCdt", "Finaldt", "Floordt")
  #  Note: L. Satterfield thinned, floored, tz adjusted & AnimalID most of data,
  #  I cleaned Nov. 2020 - on data... LMT_DateTime messed up do to raw data
  #  being formatted differently btwn Satterfield & my data sets- DON'T USE that column!
  #  Large dispersal events are already excluded
  coug_skinny <- read.csv("coug_clean 2021-12-07.csv") %>% #./Data/Cougar_Vectronic_ATS_Spring2021_4hrs.csv
    mutate(Floordt = as.POSIXct(Floordt, format = "%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+8")) %>%
    # mutate(ID = as.factor(as.character(ID)),
    #        CollarID = Collar,
    #        Sex = str_sub(ID, -1),
    #        Latitude = Lat,
    #        Longitude = Long,
    #        StudyArea = ifelse(grepl("NE", ID), "NE", "OK"),
    #        daytime = as.POSIXct(LMT_DateTime, format = "%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+8"),
    #        Finaldt = daytime,
    #        Floordt = daytime) %>%
    dplyr::select("No", "ID", "CollarID", "Sex", "Latitude", "Longitude", "LMT_DateTime", 
                "StudyArea", "daytime", "Finaldt", "Floordt")
  #  Note: L. Satterfield thinned, floored, tz adjusted & AnimalID most of data,
  #  I cleaned Nov. 2020 - on data... LMT_DateTime messed up do to raw data
  #  being formatted differently btwn Satterfield & my data sets- DON'T USE that column!
  #  Large dispersal events are already excluded
  wolf_skinny <- read.csv("wolf_clean 2021-12-07.csv") %>% #./Data/Wolf_Vectronic_Spring2021_4hrs.csv
    mutate(Floordt = as.POSIXct(Floordt, format = "%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+8")) %>%
    # mutate(ID = as.factor(as.character(ID)),
    #        CollarID = Collar,
    #        Sex = str_sub(ID, -1),
    #        Latitude = Lat,
    #        Longitude = Long,
    #        StudyArea = ifelse(grepl("W61M", ID), "OK", "NE"),  
    #        StudyArea = ifelse(grepl("W71F", ID), "OK", StudyArea),
    #        StudyArea = ifelse(grepl("W88M", ID), "OK", StudyArea),
    #        StudyArea = ifelse(grepl("W93M", ID), "OK", StudyArea),
    #        StudyArea = ifelse(grepl("W94M", ID), "OK", StudyArea),
    #        daytime = as.POSIXct(LMT_DateTime, format = "%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+8"),
    #        Finaldt = daytime,
    #        Floordt = daytime) %>%
    dplyr::select("No", "ID", "CollarID", "Sex", "Latitude", "Longitude", "LMT_DateTime", 
            "StudyArea", "daytime", "Finaldt", "Floordt")
  

  #  Save cleaned data
  clean_data <- list(md_skinny, elk_skinny, wtd_skinny, coug_skinny, wolf_skinny, meso_skinny)
  save(clean_data, file = "./Data/Collar_AllSpecies_AllLocations_Clean.RData")

  
  #  Function to truncate telemetry data by excluding first 3 weeks after capture
  #  And removes extra locations that arise when collars go into mortality mode
  Trunk_telem <- function(info, telem) {
    
    #  Truncate telemetry data
    #  Create empty data frame to fill iteratively
    trunk <- data.frame()
    #  Loop over every unique individual animal and...
    for(i in 1:nrow(info)){
      #  Take the individual animal ID
      ID <- droplevels(info$IndividualIdentifier[i])
      #  Take the animal's GPS collar serial number 
      SN <- info$GPSCollarSerialNumber[i]
      #  Buffer capture date to remove locations affected by capture event
      #  Suggested to only use data starting 3 weeks after the capture data
      start <- info$CaptureDate[i] + 20
      #  Identify last day of locations per collar
      end <- info$EndDate[i]
      #  Subset telemetry data to the specific individual
      collar <- subset(telem, CollarID == SN)
      #  Add a new column to the telemetry data with the animal's individual ID
      collar$ID <- ID
      #  truncate telemetry data by new start and end dates for that individual
      collartrunk <- subset(collar, Finaldt >= start & Finaldt <= end) 
      #  Append each unique animal's locations to a clean dataframe
      trunk <- rbind(trunk, collartrunk)
    }
    #  Format data in chronological order by individual & make lat/long numeric
    thin_trunk <- as.data.frame(trunk) %>%
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
    
    return(thin_trunk)
  }
  
  #  Run species-specific ID and telemetry data through the function
  #  Don't forget that all meso data were collected on the "wrong" fix schedule
  md_trunk <- Trunk_telem(md_info, md_skinny)
  elk_trunk <- Trunk_telem(elk_info, elk_skinny)
  wtd_trunk <- Trunk_telem(wtd_info, wtd_skinny) 
  coug_trunk <- Trunk_telem(coug_info, coug_skinny)
  wolf_trunk <- Trunk_telem(wolf_info, wolf_skinny)
  meso_trunk <- Trunk_telem(meso_info, meso_skinny)
  
  # #  Function not working?
  # #  Probably have a mismatched number of unique IDs in telem vs info data
  # t <- as.data.frame(unique(wolf_skinny$ID)) %>%
  #   mutate(d = "telemetry data")
  # colnames(t) <- c("UniqueID", "data source")
  # i <- as.data.frame(unique(wolf_info$IndividualIdentifier)) %>%
  #   mutate(d = "capture data")
  # colnames(i) <- c("UniqueID", "data source")
  # diff <- full_join(t, i, by = "UniqueID")
  
  #  Save truncated (but not thinned or filtered) data
  trunk_data <- list(md_trunk, elk_trunk, wtd_trunk, coug_trunk, wolf_trunk, meso_trunk)
  save(trunk_data, file = "./Data/Collar_AllSpecies_AllLocations_Truncated.RData")
  load("./Data/Collar_AllSpecies_AllLocations_Truncated.RData")
  
  
  #  Thin data to correct fix schedule only
  Thin_telem <- function(trunk) {
    #  Thin truncated data to only include locations on correct fix schedule
    thin_trunk <- with(trunk, trunk[hour(Floordt) == 2 | hour(Floordt) == 6 |
                                      hour(Floordt) == 10 | hour(Floordt) == 14 |
                                      hour(Floordt) == 18 | hour(Floordt) == 22,])
    return(thin_trunk)
  }
  
  #  Run ungulate & large predator data through function to thin locations to "correct" schedule
  #  All meso data collected on "wrong" fix schedule so this drops ALL meso data
  md_thin <- Thin_telem(md_trunk)
  elk_thin <- Thin_telem(elk_trunk)
  wtd_thin <- Thin_telem(wtd_trunk)
  coug_thin <- Thin_telem(coug_trunk)
  wolf_thin <- Thin_telem(wolf_trunk)


  ####  Filter data to desired date ranges  #### 
  #  ----------------------------------------
  #  SUMMER (July 1 - Sept 30): reflects when ungulate offspring are dependent 
  #  but mobile and before most fall migrations begin for mule deer (although
  #  it's kind of odd timing w/ respect to wolf rendezvous site use).
  #  Winter (Dec 1 - Feb 29): reflects when winter conditions are harshest and 
  #  mule deer are on their winter range.
  Seasonal_telem <- function(telem) {
    #  Summer 2018: 07/01/2018 - 09/30/2018 
    telem_summer18 <- telem %>%
      filter(Floordt > "2018-07-01 00:00:00") %>%
      filter(Floordt < "2018-10-01 00:00:00") %>%
      mutate(
        Season = "Summer18",
        Year = "Year1",
        FullID = paste0(ID, "_", Year)
      )
    #  Summer 2019: 07/01/2019 - 09/30/2019 
    telem_summer19 <- telem %>%
      filter(Floordt > "2019-07-01 00:00:00") %>%
      filter(Floordt < "2019-10-01 00:00:00") %>%
      mutate(
        Season = "Summer19",
        Year = "Year2",
        FullID = paste0(ID, "_", Year)
      )
    #  Summer 2020: 07/01/2020 - 09/30/2020
    telem_summer20 <- telem %>%
      filter(Floordt > "2020-07-01 00:00:00") %>%
      filter(Floordt < "2020-10-01 00:00:00") %>%
      mutate(
        Season = "Summer20",
        Year = "Year3",
        FullID = paste0(ID, "_", Year)
      )
    #  Winter 2018-2019: 12/1/2018 - 02/29/2019
    telem_winter1819 <- telem %>%
      filter(Floordt > "2018-12-01 00:00:00") %>%
      filter(Floordt < "2019-03-01 00:00:00") %>%
      mutate(
        Season = "Winter1819",
        Year = "Year1",
        FullID = paste0(ID, "_", Year)
      )
    #  Winter 2019-2020: 12/1/2019 - 02/28/2020 (leap year)
    telem_winter1920 <- telem %>%
      filter(Floordt > "2019-12-01 00:00:00") %>%
      filter(Floordt < "2020-03-01 00:00:00")  %>%
      mutate(
        Season = "Winter1920",
        Year = "Year2",
        FullID = paste0(ID, "_", Year)
      )
    #  Winter 2020-2021: 12/1/2020 - 02/29/2021
    telem_winter2021 <- telem %>%
      filter(Floordt > "2020-12-01 00:00:00") %>%
      filter(Floordt < "2021-03-01 00:00:00")  %>%
      mutate(
        Season = "Winter2021",
        Year = "Year3",
        FullID = paste0(ID, "_", Year)
      )
    #  Combine into single file
    telem_smwtr <- rbind(telem_summer18, telem_winter1819, telem_summer19, telem_winter1920, telem_summer20, telem_winter2021) 
    return(telem_smwtr)
  }

  #  Run species-specific telemetry data through function to filter by data range
  md_season <- Seasonal_telem(md_trunk)
  elk_season <- Seasonal_telem(elk_trunk)
  wtd_season <- Seasonal_telem(wtd_trunk)
  coug_season <- Seasonal_telem(coug_trunk)
  wolf_season <- Seasonal_telem(wolf_trunk)
  meso_season <- Seasonal_telem(meso_trunk)
  
  #  Same thing but with fully thinned ungulate data
  md_season2 <- Seasonal_telem(md_thin)
  elk_season2 <- Seasonal_telem(elk_thin)
  wtd_season2 <- Seasonal_telem(wtd_thin)
  coug_season2 <- Seasonal_telem(coug_thin)
  wolf_season2 <- Seasonal_telem(wolf_thin)

  #  Identify which collars are on incorrect fix schedule
  fix_schedule <- function(trunk) {
    schedj <- trunk %>%
      mutate(hour = as.integer(strftime(Floordt, format = "%H", tz="Etc/GMT+8")),
             schedj = ifelse(hour == 2 | hour == 6 | hour == 10 | hour == 14 | hour == 18 | hour == 22, "schedule_1", "schedule_2"),
             time_gap = as.numeric(difftime(Floordt, lag(Floordt), tz = "Etc/GMT+8", units = "hours")),
             leap_schedj = lead(schedj, n = 2),
             drop_loc = ifelse(schedj == "schedule_2" & time_gap <= 2 & schedj == leap_schedj, "drop", "keep")) %>%
      filter(drop_loc == "keep") %>%
      mutate(new_gap = as.numeric(difftime(Floordt, lag(Floordt), tz = "Etc/GMT+8", units = "hours")),
             new_gap = ifelse(new_gap < 0, "NA", new_gap),
             next_schedj = lead(schedj, n = 1),
             drop_loc1 = ifelse(new_gap == 2 & schedj != next_schedj, "drop", "keep")) %>%
      filter(drop_loc1 == "keep") %>%
      mutate(last_schedj = lead(schedj, n = 1),
             new_track = ifelse(schedj == last_schedj, 0, 1)) %>%
      dplyr::select(-c("schedj", "time_gap", "leap_schedj", "drop_loc", "new_gap", "next_schedj", "drop_loc1", "last_schedj")) #, "new_track"
    return(schedj)
  }
  md_schedj <- fix_schedule(md_season)
  elk_schedj <- fix_schedule(elk_season)
  wtd_schedj <- fix_schedule(wtd_season)
  coug_schedj <- fix_schedule(coug_season)
  wolf_schedj <- fix_schedule(wolf_season)
  meso_schedj <- fix_schedule(meso_season)
  
  # 3709ELK18
  # 32WTD18
  # 3961MD17
  
  
  ####  Fix Success & Gappy Data  ####
  #  ------------------------------
  #  Function to calculate number of locations per individual and season
  #  3 mo study period w/ 6 fixes/day = approx. 540 locations, if no missed fixes
  sum_locs <- function(telem) {
    
    Nmb_locs <- telem %>%
      group_by(ID, Season) %>%
      summarise(count = n())
    
    return(Nmb_locs)
  }
  
  #  Run species-specific data through function based on preferred fix schedule
  md_counts <- sum_locs(md_season) 
  elk_counts <- sum_locs(elk_season) 
  wtd_counts <- sum_locs(wtd_season) 
  coug_counts <- sum_locs(coug_season) 
  wolf_counts <- sum_locs(wolf_season) 
  meso_counts <- sum_locs(meso_season)
  #  How much data am I losing if I stick to only 1 fix schedule for ungulates/large pred?
  md_counts2 <- sum_locs(md_season2)
  elk_counts2 <- sum_locs(elk_season2)
  wtd_counts2 <- sum_locs(wtd_season2)
  coug_counts2 <- sum_locs(coug_season2)
  wolf_counts2 <- sum_locs(wolf_season2)
  #  How much data am I losing if stick to 4 hour fix rate but ignore schedule?
  md_counts3 <- sum_locs(md_schedj)
  elk_counts3 <- sum_locs(elk_schedj)
  wtd_counts3 <- sum_locs(wtd_schedj)
  coug_counts3 <- sum_locs(coug_schedj)
  wolf_counts3 <- sum_locs(wolf_schedj)

  
  #  Visually inspect counts
  #  Looking for time periods with many missing fixes (< 400/season) and collars 
  #  that were on wrong schedule, switched schedules, or took extra locations
  md_cnt <- full_join(md_counts, md_counts3, by = c("ID", "Season"))
  elk_cnt <- full_join(elk_counts, elk_counts3, by = c("ID", "Season"))
  wtd_cnt <- full_join(wtd_counts, wtd_counts3, by = c("ID", "Season"))
  #  Rename predator info for consistency
  coug_cnt <- coug_counts
  wolf_cnt <- wolf_counts
  meso_cnt <- meso_counts

  #  Pull out collars with lots of missing fixes (< 400/season)
  md_400 <- md_cnt[md_cnt$count.x < 401 & md_cnt$count.y < 401 | is.na(md_cnt$count.x) & md_cnt$count.y < 401 | md_cnt$count.x < 401 & is.na(md_cnt$count.y),]
  elk_400 <- elk_cnt[elk_cnt$count.x < 401 & elk_cnt$count.y < 401 | is.na(elk_cnt$count.x) & elk_cnt$count.y < 401 | elk_cnt$count.x < 401 & is.na(elk_cnt$count.y),]
  wtd_400 <- wtd_cnt[wtd_cnt$count.x < 401 & wtd_cnt$count.y < 401 | is.na(wtd_cnt$count.x) & wtd_cnt$count.y < 401 | wtd_cnt$count.x < 401 & is.na(wtd_cnt$count.y),]
  coug_400 <- coug_cnt[coug_cnt$count < 401,]
  wolf_400 <- wolf_cnt[wolf_cnt$count < 401,]
  meso_400 <- meso_cnt[meso_cnt$count < 401,]
  
  #  Are these missing fixes randomly distributed across the season or are there 
  #  large blocks of missing time (e.g., at beginning or end of time period?)
  md_locs <- inner_join(md_schedj, md_400, by = c("ID", "Season"))
  elk_locs <- inner_join(elk_schedj, elk_400, by = c("ID", "Season"))
  wtd_locs <- inner_join(wtd_schedj, wtd_400, by = c("ID", "Season"))
  coug_locs <- inner_join(coug_schedj, coug_400, by = c("ID", "Season"))
  wolf_locs <- inner_join(wolf_schedj, wolf_400, by = c("ID", "Season"))
  meso_locs <- inner_join(meso_schedj, meso_400, by = c("ID", "Season"))
  
  #  Calculate the number of hours between subsequent locations
  #  How often are there big gaps and when are those gaps occurring?
  diftime <- function(locs) {
    time_gap <- locs %>%
      group_by(ID, Season) %>%
      arrange(Floordt) %>%
      mutate(
        time_gap = difftime(Floordt, lag(Floordt), tz = "Etc/GMT+8", units = "hours"),
        time_gap = as.numeric(time_gap),
        gap_4hr = ifelse(time_gap <= 4, 0, time_gap),
        gap_8hr = ifelse(time_gap <= 8, 0, time_gap),
        gap_24hr = ifelse(time_gap <= 24, 0, time_gap),
        date_range = difftime(max(Floordt), min(Floordt), tz = "Etc/GMT+8", units = "days"),
        date_range = round(date_range, digits = 0)
      )
    #  Plot gaps in telemetry locations that are greater than 24hrs in length
    hist(time_gap$gap_24hr, breaks = 50, main = "Histogram of time gaps", xlab = "Hours between locations")
    hist(time_gap$gap_4hr, breaks = 1, main = "Histogram of 2 hour fixs", xlab = "Hours between locations")
    return(time_gap)
  }
  
  #  Run data through function & visually inspect potentially problematic collars
  md_gaps <- as.data.frame(diftime(md_locs))
  #  MULE DEER: 3959MD17, 90MD18, 66MD18, 44MD18, 42MD18, 29MD18, 23MD18, 190MD18, 20MD18 ended early
  #  251MD20, 196MD20, & 200MD20 missing 1 or 2 large chunks fo data in Winter1920
  #  2020 mule deer collars deployed in early Jan so short 1.5 months of Winter1920 data
  elk_gaps <- as.data.frame(diftime(elk_locs))
  #  ELK: 3692EA17, 3705EA17 & 3726ELK18 failed or animal died part way through season
  #  2020 elk collars deployed in early Jan so short 1.5 months of Winter1920 data
  #  3696EA17 & 3697EA17 are missing large chunks of data sporadically thru all seasons
  wtd_gaps <- as.data.frame(diftime(wtd_locs))
  #  WHITE-TAILED DEER: many collars have delayed start in winter season
  #  Some with large sporadic gaps missing
  #  Many collars don't last a full month during a season
  coug_gaps <- as.data.frame(diftime(coug_locs))
  #  COUGAR: Some delayed starts in winter, some sporadic big gaps
  wolf_gaps <- as.data.frame(diftime(wolf_locs))
  #  WOLF: Some delayed starts in winter and summer
  #  Mostly some HUG gaps reduce number of locations to very few sequential points
  #  Not sure how much I'm going to get out of the wolf data... :(
  meso_gaps <- as.data.frame(diftime(meso_locs))
  #  MESO: many bobcat collars have delayed start in winter, coyotes in summer
  #  Some with large gaps missing, many with frequent small gaps


  #  Exclude collars with too few locations in a given season
  #  Too few locations tend to arise from either many HUGE gaps in a full season  
  #  or collars starting/stopping sometime within a season. 
  #  Arbitrary cut-offs to exclude collars with very little data:
  #  Collars must collect 90+ locations/season or be deployed for 20+ days/season
  #  Chose 20 days b/c wanted to include collars that were deployed late in season
  #  or died early but still provided consecutive data for almost a full month
  #  Chose 90 locations b/c that would average to 1 location per day over 90 day season
  #  1. Identify which collars/season need to go
  tf_md <- md_gaps[md_gaps$count.y < 90 | md_gaps$date_range < 20,]
  tf_elk <- elk_gaps[elk_gaps$count.y < 90 | elk_gaps$date_range < 20,]
  tf_wtd <- wtd_gaps[wtd_gaps$count.y < 90 | wtd_gaps$date_range < 20,] # weird NAs introduced, not sure why, but don't seem to be a problem below
  tf_coug <- coug_gaps[coug_gaps$count < 90 | coug_gaps$date_range < 20,]
  tf_wolf <- wolf_gaps[wolf_gaps$count < 90 | wolf_gaps$date_range < 20,]
  tf_meso <- meso_gaps[meso_gaps$count < 90 | meso_gaps$date_range < 20,]
  #  2. Exclude these locations from telemetry data so data are good-to-go for analyses
  md_gtg <- anti_join(md_schedj, tf_md, by = c("ID", "Season"))
  elk_gtg <- anti_join(elk_schedj, tf_elk, by = c("ID", "Season"))
  wtd_gtg <- anti_join(wtd_schedj, tf_wtd, by = c("ID", "Season"))
  coug_gtg <- anti_join(coug_schedj, tf_coug, by = c("ID", "Season"))
  wolf_gtg <- anti_join(wolf_schedj, tf_wolf, by = c("ID", "Season"))
  meso_gtg <- anti_join(meso_schedj, tf_meso, by = c("ID", "Season"))
  #  Don't forget to separate out bobcats from coyotes for actual analyses
  coy_gtg <- filter(meso_gtg, Species == "Coyote")
  bob_gtg <- filter(meso_gtg, Species == "Bobcat")
  
  #  Species_gtg are final data sets for HMM analyses and 3rd-order RSFs
  
  #  Save RData for easy transfer to other computers
  # save.image(paste0("./Data/Collar_Truncating&Filtering_", Sys.Date(), ".RData"))
  save.image(paste0("./Data/Collar_Truncating&Filtering_noDispersal_CleanedFixSchedule_", Sys.Date(), ".RData"))
  load("./Data/Collar_Truncating&Filtering_noDispersal_CleanedFixSchedule_2022-03-08.RData") #2022-02-18 included 2 hr fix schedules (OK for RSFs but bad for HMM)
  
  
  ####  Remove locations associated with mule deer migration tracks  ####
  #  -----------------------------------------------------------------
  #  Migration times for Mule Deer
  #  Use for 3rd-order RSF analyses & HMMs since migration movements are very 
  #  different than within-HR movements that are the focus of this analysis
  #  Migration dates identified by T.Ganz
  md_migtimes <- read.csv("./Data/migtime_trim_062121.csv") %>%
    transmute(
      AnimalID = newUid,
      Year = nsdYear,
      id_yr = id_yr,
      Season = Season,
      Start_mig = lubridate::as_date(start_mig, format = "%m/%d/%Y"),
      End_mig = lubridate::as_date(end_mig, format = "%m/%d/%Y"),
    )
  
  #  Function to identify locations generated during migration 
  #  Note: only identifies locations from dates that overlap summer & winter
  #  seasons of interest. Migrations in Oct - Nov & March - June are ignored & 
  #  automatically excluded via filtering above.
  md_gtg <- mutate(md_gtg, FloorDate = lubridate::as_date(Floordt))  #tz = "Etc/GMT+8"
  #  Create empty data frame to fill with thinned location data
  migtimes <- data.frame()
  #  Loop over every unique individual animal and...
  for(i in 1:nrow(md_migtimes)){
    #  Take the individual animal ID
    AnimalID <- md_migtimes$AnimalID[i]
    #  Identify the start date to exclude
    start <- md_migtimes$Start_mig[i]
    #  Identify last date to exclude
    end <- md_migtimes$End_mig[i]
    #  Subset based on migration data
    migrant <- subset(md_gtg, ID == AnimalID)
    #  Exclude telemetry data outside start & end dates for that individual
    migrations <- filter(migrant, FloorDate >= start & FloorDate <= end) # why is it NOT including the actual start and end dates?!
    #  Append each animal's migration locations to a clean data frame
    migtimes <- rbind(migtimes, migrations)
  }
    
  #  Remove migration times from larger location data frame
  md_gtg_nomig <- anti_join(md_gtg, migtimes) %>%
    #  Format data in chronological order by individual & make lat/long numeric
    arrange(ID, Floordt) %>%
    #  Make sure lat/long are in a numeric format
    mutate(
      Latitude = as.numeric(Latitude),
      Longitude = as.numeric(Longitude)
    )
  
  
  #  Drop random problematic locations from WTD data (for 3rd order RSF only!)
  #  Not sure why these locations are a problem but they prevent kernelUD() from 
  #  estimating a 95% KDE for these individuals. 
  #  ----------------------------------------------
  #  3144WTD19 has one random location that throw off kernelUD function
  #  This location is well away from an otherwise pretty tight home range
  wtd_gtg <- filter(wtd_gtg, wtd_gtg$ID != "3144WTD19" | wtd_gtg$Latitude != "48.04948403")
  #  4806WTD20 has some apparent migrations that throw off kernelUD function
  wtd_gtg <- filter(wtd_gtg, wtd_gtg$ID != "4806WTD20" | wtd_gtg$Longitude != "-117.99087116") 
  wtd_gtg <- filter(wtd_gtg, wtd_gtg$ID != "4806WTD20" | wtd_gtg$Longitude != "-118.03944301")
  wtd_gtg <- filter(wtd_gtg, wtd_gtg$ID != "4806WTD20" | wtd_gtg$Longitude != "-118.03989151")
  #  4828WTD20 has some apparent migrations that throw off kernelUD function
  wtd_gtg <- filter(wtd_gtg, wtd_gtg$ID != "4828WTD20" | wtd_gtg$Latitude != "48.09227101")
  wtd_gtg <- filter(wtd_gtg, wtd_gtg$ID != "4828WTD20" | wtd_gtg$Latitude != "48.02359196")
  wtd_gtg <- filter(wtd_gtg, wtd_gtg$ID != "4828WTD20" | wtd_gtg$Latitude != "48.02237794")
  wtd_gtg <- filter(wtd_gtg, wtd_gtg$ID != "4828WTD20" | wtd_gtg$Latitude != "48.02231165")
  
    
  #  Save RData for easy transfer to other computers
  # save.image(paste0("./Data/Collar_Truncating&Filtering_", Sys.Date(), ".RData"))
  # save.image(paste0("./Data/Collar_Truncating&Filtering_noDispersal_", Sys.Date(), ".RData"))
  save.image(paste0("./Data/Collar_Truncating&Filtering_noDispMig_CleanedFixSchedule_", Sys.Date(), ".RData"))
  
  #  Next step is Collar_Movement_DataPrep.R script

  
