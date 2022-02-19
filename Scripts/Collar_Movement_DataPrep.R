  #'  ============================================
  #'  Movement Data Prep (cam vs collar analysis)
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  April 2021
  #'  ============================================
  #'  Script to format cleaned GPS location data for subsequent HMMM analyses 
  #'  for deer, elk, cougars, wolves, coyotes, and bobcats for summer 2018 
  #'  (7/1/18 - 9/29/18) and winter 2018-2019 (12/1/18 - 3/1/19), respectively. 
  #'  Data were collected & generously provided by WPPP collaborators including
  #'  T.Ganz, T.Roussin, L.Satterfield, B.Windell, and others. Code adapted from
  #'  momentuHMM GitHub, J.Merkel Movement Workshop, L.Satterfield, & R.Emmet.
  #'  Time periods and covariates to match up with single-season occupancy models.
  #'  
  #'  Telemetry data initially cleaned with Collar_DataCleaning.R scripts, 
  #'  Collar_Truncating&Filtering.R script, & by T.Ganz, L.Satterfield & B.Windell.
  #'  Covariate data based on remotely sensed data available from various sources
  #'  (noted in occupancy model scripts in CamTraps_and_Collars repository).
  #'  ============================================
  
  #'  Clear memory
  # rm(list=ls())
  
  #'  Load libraries
  library(momentuHMM)
  library(rgdal)
  library(tidyverse)
  
  #'  Source cleaned telemetry data
  # load("./Data/Collar_Truncating&Filtering_noDispMig_2021-11-16.RData") # includes some deer data with low fix rate, esp. the white-tail data
  # load("./Data/Collar_Truncating&Filtering_noDispMig_2021-12-02.RData") # accidently excludes some coyote and bobcat collars
  load("./Data/Collar_Truncating&Filtering_noDispMig_2022-02-18.RData")
  
  #' I chose to use relocation data that excludes obvious dispersal events that
  #' take carnivores away from extent of study areas and relocation data during
  #' mule deer migrations because these long-distance movements are very different
  #' than the typical within-home range movements that I'm interested in studying.
  #' Dispersals and migrations tend to be very directed and quick (citations?) 
  #' compared to general encamped and exploratory behaviors in an animal's HR.

  
  ####  Data preparation  ####
  #'  Select relevant columns
  #'  Keeping version of datetime that have been floored to beginning of hour
  rawELK <- elk_gtg %>%
    dplyr::select(ID, FullID, Sex, Season, StudyArea, Longitude, Latitude, Floordt) %>%
    arrange(ID, Floordt)
  colnames(rawELK) <- c("ID", "FullID", "Sex", "Season", "StudyArea", "Long", "Lat", "time")
  #'  Only keep first track to practice with
  # rawELK1 <- subset(rawELK, ID == unique(ID)[2])
  rawMD <- md_gtg_nomig %>% #md_gtg  # depends on whether migration movements are excluded
    dplyr::select(ID, FullID, Sex, Season, StudyArea, Longitude, Latitude, Floordt) %>%
    arrange(ID, Floordt)
  colnames(rawMD) <- c("ID", "FullID", "Sex", "Season", "StudyArea", "Long", "Lat", "time")
  rawWTD <- wtd_gtg %>%
    dplyr::select(ID, FullID, Sex, Season, StudyArea, Longitude, Latitude, Floordt)%>%
    arrange(ID, Floordt)
  colnames(rawWTD) <- c("ID", "FullID", "Sex", "Season", "StudyArea", "Long", "Lat", "time")
  rawCOUG <- coug_gtg %>%
    dplyr::select(ID, FullID, Sex, Season, StudyArea, Longitude, Latitude, Floordt)%>%
    arrange(ID, Floordt)
  colnames(rawCOUG) <- c("ID", "FullID", "Sex", "Season", "StudyArea", "Long", "Lat", "time")
  rawWOLF <- wolf_gtg %>%
    dplyr::select(ID, FullID, Sex, Season, StudyArea, Longitude, Latitude, Floordt)%>%
    arrange(ID, Floordt)
  colnames(rawWOLF) <- c("ID", "FullID", "Sex", "Season", "StudyArea", "Long", "Lat", "time")
  rawBOB <- bob_gtg %>%
    dplyr::select(ID, FullID, Sex, Season, StudyArea, Longitude, Latitude, Floordt)%>%
    arrange(ID, Floordt)
  colnames(rawBOB) <- c("ID", "FullID", "Sex", "Season", "StudyArea", "Long", "Lat", "time")
  rawCOY <- coy_gtg %>%
    dplyr::select(ID, FullID, Sex, Season, StudyArea, Longitude, Latitude, Floordt)%>%
    arrange(ID, Floordt)
  colnames(rawCOY) <- c("ID", "FullID", "Sex", "Season", "StudyArea", "Long", "Lat", "time")

  #'  Function to covert times to POSIX & make locations spatial for each species
  prep_raw <- function(raw, plotit = TRUE) {
    raw$time <- as.POSIXct(raw$time, format = "%Y-%m-%d %H:%M:%S", tz = "America/Los_Angeles")

    #'  Make locations spatial and project to UTM coordinates with study area projection
    llcoord <- SpatialPoints(raw[,6:7], proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
    utmcoord <- spTransform(llcoord, CRS("+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs "))
    #'  Add UTM locations to data frame
    raw$x <- attr(utmcoord, "coords")[,1]
    raw$y <- attr(utmcoord, "coords")[,2]

    return(raw)
  }

  #'  Run data from each species through data prep function
  rawMD <- as.data.frame(prep_raw(rawMD))
  rawELK <- as.data.frame(prep_raw(rawELK))
  rawWTD <- as.data.frame(prep_raw(rawWTD))
  rawCOUG <- as.data.frame(prep_raw(rawCOUG))
  rawWOLF <- as.data.frame(prep_raw(rawWOLF))
  rawBOB <- as.data.frame(prep_raw(rawBOB))
  rawCOY <- as.data.frame(prep_raw(rawCOY))

  #'  Quick peak at example data for each species
  plot_collar <- function(raw) {
    #'  Pull out data from first collar only
    first <- raw %>%
      slice_head(n = 1) %>%
      dplyr::select(ID)
    animal1 <- raw[raw$ID == first$ID,]
    #'  Plot the first individual in data set
    ggplot(animal1, aes(x = Long, y = Lat)) +
      geom_point() +
      geom_path() +
      facet_wrap(FullID~Season)
  }
  #'  Plot seasonal locations for one individual
  plot_collar(rawMD)
  plot_collar(rawELK)
  plot_collar(rawWTD)
  plot_collar(rawCOUG)
  plot_collar(rawWOLF)
  plot_collar(rawBOB)
  plot_collar(rawCOY)

  #'  Take a closer look at those seasonal locations, esp. the mule deer
  #'  Are they migrating halfway through a season? Any dispersal?
  #'  Function to plot locations from individual animals
  plot_telem <- function(spdf){
    #'  Split out spatial points df by individual animal ID
    ind_animal <- group_split(spdf, spdf$ID)
    #'  Place holder for unique animal ID
    names <- c()
    #'  Empty list to hold individual maps
    plot_list <- list()
    #'  Loop through all animals one at a time to create maps of their locations
    for(i in 1:length(unique(ind_animal))) {
      names <- c(names, unique(as.character(ind_animal[[i]]$ID)))
      plot <- ggplot(ind_animal[[i]], aes(x = Long, y = Lat, color = time)) +
        geom_point() +
        geom_path() +
        #'  Nifty side-by-side plots
        facet_wrap(FullID~Season)
      plot_list[[i]] <- plot
    }
    return(plot_list)
  }
  #'  Run data from each species through
  #'  Looking for obvious movement from summer to winter range
  #'  within a single season's worth of locations, suggesting migration or dispersal
  rawMD_maps <- plot_telem(rawMD)
  rawELK_maps <- plot_telem(rawELK)
  rawWTD_maps <- plot_telem(rawWTD)
  rawCOUG_maps <- plot_telem(rawCOUG)
  rawWOLF_maps <- plot_telem(rawWOLF)
  rawBOB_maps <- plot_telem(rawBOB)
  rawCOY_maps <- plot_telem(rawCOY)

  #' #'  Save plots as PDF to review and look for evidence of migration or dispersal
  #' pdf("./Outputs/GPSlocs_byseason_maps.pdf")
  #' for (i in 1:length(unique(rawMD_maps))) {
  #'   print(rawMD_maps[[i]])
  #' }
  #' for (i in 1:length(unique(rawELK_maps))) {
  #'   print(rawELK_maps[[i]])
  #' }
  #' for (i in 1:length(unique(rawWTD_maps))) {
  #'   print(rawWTD_maps[[i]])
  #' }
  #' for (i in 1:length(unique(rawCOUG_maps))) {
  #'   print(rawCOUG_maps[[i]])
  #' }
  #' for (i in 1:length(unique(rawWOLF_maps))) {
  #'   print(rawWOLF_maps[[i]])
  #' }
  #' for (i in 1:length(unique(rawBOB_maps))) {
  #'   print(rawBOB_maps[[i]])
  #' }
  #' for (i in 1:length(unique(rawCOY_maps))) {
  #'   print(rawCOY_maps[[i]])
  #' }
  #' dev.off()


  #'  Identify bursts of sequential locations & where there are prolonged gaps
  #'  Source creat.burst.R function to identify bursts
  #'  Requires location data be ordered by ID and time (do this above)
  #'  Script written by J.Merkle & provided at Movement Workshop
  source("./Scripts/creat.burst.R")

  #'  Function to run locations for each species through the creat.burst function
  bursts <- function(rawloc) {
    #'  Tmax = 28.25 hours (in seconds) so that locations are still grouped in a
    #'  single burst if there's a gap of 24hr or less in the data (up to 6 sequential
    #'  fixes missed with 4 hours on each end) but a new burst if gap is >24 hrs.
    rawloc$burst <- creat.burst(data = rawloc, id = TRUE, id_name = "ID", date_name = "time", Tmax = 87300)
    #'  Exclude super short bursts (where burst length is 3 or less locations)
    #'  because need at least 3 points to get a turning angle
    loc_burst <- rawloc[rawloc$burst %in% names(table(rawloc$burst))[table(rawloc$burst) >=3],]
    length(unique(loc_burst$burst))
    head(loc_burst)
    #'  Add burst value ids to Unique ID column to create unique IDs for each track
    loc_burst$UniqueID <- with(loc_burst, paste0(ID, "_", burst))
    loc_burst <- loc_burst %>%
      #'  Change ID from unique animal ID to a combo animal ID & track number identifier
      #'  This is important for crawlWrap function to ensure it does not interpolate
      #'  large chunks of missing data
      mutate(
        AnimalID = ID,
        ID = UniqueID
      ) %>%
      dplyr::select(-UniqueID)
    return(loc_burst)
  }
  #'  Run each species through the function that identifies bursts in the data
  MD_track <- bursts(rawMD)
  ELK_track <- bursts(rawELK)
  WTD_track <- bursts(rawWTD)
  COUG_track_OK <- bursts(rawCOUG[rawCOUG$StudyArea == "OK",])
  COUG_track_NE <- bursts(rawCOUG[rawCOUG$StudyArea == "NE",])
  WOLF_track_OK <- bursts(rawWOLF[rawWOLF$StudyArea == "OK",])
  WOLF_track_NE <- bursts(rawWOLF[rawWOLF$StudyArea == "NE",])
  BOB_track_OK <- bursts(rawBOB[rawBOB$StudyArea == "OK",])
  BOB_track_NE <- bursts(rawBOB[rawBOB$StudyArea == "NE",])
  COY_track_OK <- bursts(rawCOY[rawCOY$StudyArea == "OK",])
  COY_track_NE <- bursts(rawCOY[rawCOY$StudyArea == "NE",])
  
  #'  Save track data sets
  save(MD_track, file = "./Outputs/Telemetry_tracks/MD_track.RData")
  save(ELK_track, file = "./Outputs/Telemetry_tracks/ELK_track.RData")
  save(WTD_track, file = "./Outputs/Telemetry_tracks/WTD_track.RData")
  save(COUG_track_OK, file = "./Outputs/Telemetry_tracks/COUG_track_OK.RData")
  save(COUG_track_NE, file = "./Outputs/Telemetry_tracks/COUG_track_NE.RData")
  save(WOLF_track_OK, file = "./Outputs/Telemetry_tracks/WOLF_track_OK.RData")
  save(WOLF_track_NE, file = "./Outputs/Telemetry_tracks/WOLF_track_NE.RData")
  save(BOB_track_OK, file = "./Outputs/Telemetry_tracks/BOB_track_OK.RData")
  save(BOB_track_NE, file = "./Outputs/Telemetry_tracks/BOB_track_NE.RData")
  save(COY_track_OK, file = "./Outputs/Telemetry_tracks/COY_track_OK.RData")
  save(COY_track_NE, file = "./Outputs/Telemetry_tracks/COY_track_NE.RData")
  
  #'  Create season-specific data sets for crawlWrap function
  MD_smr_track <- MD_track[MD_track$Season == "Summer18" | MD_track$Season == "Summer19", ]
  MD_wtr_track <- MD_track[MD_track$Season == "Winter1819" | MD_track$Season == "Winter1920", ]
  ELK_smr_track <- ELK_track[ELK_track$Season == "Summer18" | ELK_track$Season == "Summer19", ]
  ELK_wtr_track <- ELK_track[ELK_track$Season == "Winter1819" | ELK_track$Season == "Winter1920", ]
  WTD_smr_track <- WTD_track[WTD_track$Season == "Summer18" | WTD_track$Season == "Summer19", ]
  WTD_wtr_track <- WTD_track[WTD_track$Season == "Winter1819" | WTD_track$Season == "Winter1920", ]
  COUG_smr_track_OK <- COUG_track_OK[COUG_track_OK$Season == "Summer18" | COUG_track_OK$Season == "Summer19" | COUG_track_OK$Season == "Summer20", ]
  COUG_smr_track_NE <- COUG_track_NE[COUG_track_NE$Season == "Summer18" | COUG_track_NE$Season == "Summer19" | COUG_track_NE$Season == "Summer20", ]
  COUG_wtr_track_OK <- COUG_track_OK[COUG_track_OK$Season == "Winter1819" | COUG_track_OK$Season == "Winter1920" | COUG_track_OK$Season == "Winter2021", ]
  COUG_wtr_track_NE <- COUG_track_NE[COUG_track_NE$Season == "Winter1819" | COUG_track_NE$Season == "Winter1920" | COUG_track_NE$Season == "Winter2021", ]
  WOLF_smr_track_OK <- WOLF_track_OK[WOLF_track_OK$Season == "Summer18" | WOLF_track_OK$Season == "Summer19" | WOLF_track_OK$Season == "Summer20", ]
  WOLF_smr_track_NE <- WOLF_track_NE[WOLF_track_NE$Season == "Summer18" | WOLF_track_NE$Season == "Summer19" | WOLF_track_NE$Season == "Summer20", ]
  WOLF_wtr_track_OK <- WOLF_track_OK[WOLF_track_OK$Season == "Winter1819" | WOLF_track_OK$Season == "Winter1920" | WOLF_track_OK$Season == "Winter2021", ]
  WOLF_wtr_track_NE <- WOLF_track_NE[WOLF_track_NE$Season == "Winter1819" | WOLF_track_NE$Season == "Winter1920" | WOLF_track_NE$Season == "Winter2021", ]
  BOB_smr_track_OK <- BOB_track_OK[BOB_track_OK$Season == "Summer18" | BOB_track_OK$Season == "Summer19"| BOB_track_OK$Season == "Summer20", ]
  BOB_smr_track_NE <- BOB_track_NE[BOB_track_NE$Season == "Summer18" | BOB_track_NE$Season == "Summer19"| BOB_track_NE$Season == "Summer20", ]
  BOB_wtr_track_OK <- BOB_track_OK[BOB_track_OK$Season == "Winter1819" | BOB_track_OK$Season == "Winter1920" | BOB_track_OK$Season == "Winter2021", ]
  BOB_wtr_track_NE <- BOB_track_NE[BOB_track_NE$Season == "Winter1819" | BOB_track_NE$Season == "Winter1920" | BOB_track_NE$Season == "Winter2021", ]
  COY_smr_track_OK <- COY_track_OK[COY_track_OK$Season == "Summer18" | COY_track_OK$Season == "Summer19" | COY_track_OK$Season == "Summer20", ]
  COY_smr_track_NE <- COY_track_NE[COY_track_NE$Season == "Summer18" | COY_track_NE$Season == "Summer19" | COY_track_NE$Season == "Summer20", ]
  COY_wtr_track_OK <- COY_track_OK[COY_track_OK$Season == "Winter1819" | COY_track_OK$Season == "Winter1920" | COY_track_OK$Season == "Winter2021", ]
  COY_wtr_track_NE <- COY_track_NE[COY_track_NE$Season == "Winter1819" | COY_track_NE$Season == "Winter1920" | COY_track_NE$Season == "Winter2021", ]
  
  #'  Save seasonal tracks
  save(MD_smr_track, file = "./Outputs/Telemetry_tracks/MD_smr_track.RData")
  save(MD_wtr_track, file = "./Outputs/Telemetry_tracks/MD_wtr_track.RData")
  save(ELK_smr_track, file = "./Outputs/Telemetry_tracks/ELK_smr_track.RData")
  save(ELK_wtr_track, file = "./Outputs/Telemetry_tracks/ELK_wtr_track.RData")
  save(WTD_smr_track, file = "./Outputs/Telemetry_tracks/WTD_smr_track.RData")
  save(WTD_wtr_track, file = "./Outputs/Telemetry_tracks/WTD_wtr_track.RData")
  save(COUG_smr_track_OK, file = "./Outputs/Telemetry_tracks/COUG_smr_track_OK.RData")
  save(COUG_smr_track_NE, file = "./Outputs/Telemetry_tracks/COUG_smr_track_NE.RData")
  save(COUG_wtr_track_OK, file = "./Outputs/Telemetry_tracks/COUG_wtr_track_OK.RData")
  save(COUG_wtr_track_NE, file = "./Outputs/Telemetry_tracks/COUG_wtr_track_NE.RData")
  save(WOLF_smr_track_OK, file = "./Outputs/Telemetry_tracks/WOLF_smr_track_OK.RData")
  save(WOLF_smr_track_NE, file = "./Outputs/Telemetry_tracks/WOLF_smr_track_NE.RData")
  save(WOLF_wtr_track_OK, file = "./Outputs/Telemetry_tracks/WOLF_wtr_track_OK.RData")
  save(WOLF_wtr_track_NE, file = "./Outputs/Telemetry_tracks/WOLF_wtr_track_NE.RData")
  save(BOB_smr_track_OK, file = "./Outputs/Telemetry_tracks/BOB_smr_track_OK.RData")
  save(BOB_smr_track_NE, file = "./Outputs/Telemetry_tracks/BOB_smr_track_NE.RData")
  save(BOB_wtr_track_OK, file = "./Outputs/Telemetry_tracks/BOB_wtr_track_OK.RData")
  save(BOB_wtr_track_NE, file = "./Outputs/Telemetry_tracks/BOB_wtr_track_NE.RData")
  save(COY_smr_track_OK, file = "./Outputs/Telemetry_tracks/COY_smr_track_OK.RData")
  save(COY_smr_track_NE, file = "./Outputs/Telemetry_tracks/COY_smr_track_NE.RData")
  save(COY_wtr_track_OK, file = "./Outputs/Telemetry_tracks/COY_wtr_track_OK.RData")
  save(COY_wtr_track_NE, file = "./Outputs/Telemetry_tracks/COY_wtr_track_NE.RData")

  
  #'  Save a giant list of all seasonal tracks
  #'  NOTE the ordering of this list for the predators--- study area then season
  spp_all_tracks <- list(MD_smr_track, MD_wtr_track, ELK_smr_track, ELK_wtr_track,
                     WTD_smr_track, WTD_wtr_track, COUG_smr_track_OK, COUG_wtr_track_OK,
                     COUG_smr_track_NE, COUG_wtr_track_NE, WOLF_smr_track_OK, WOLF_wtr_track_OK, 
                     WOLF_smr_track_NE, WOLF_wtr_track_NE, BOB_smr_track_OK, BOB_wtr_track_OK,
                     BOB_smr_track_NE, BOB_wtr_track_NE, COY_smr_track_OK, COY_wtr_track_OK,
                     COY_smr_track_NE, COY_wtr_track_NE)
  #save(spp_all_tracks, file = "./Outputs/Telemetry_tracks/spp_all_tracks.RData")
  save(spp_all_tracks, file = "./Outputs/Telemetry_tracks/spp_all_tracks_noDis_noMig_SAspecific_updated021822.RData")
  
  
  #'  Load tracks
  load("./Outputs/Telemetry_tracks/spp_all_tracks_noDis_noMig_SAspecific_updated021822.RData")
  
  
  
  #'  Function to interpolate missing fixes based on regular time intervals (4hr)
  #'  Interpolating up to 6 missed fixes (24 hr gaps) within each track/animal
  #'  Use crawlWrap function from momentuHMM where:
  #'   -theta are starting values; crawlWrap defaults to 0 if none are provided
  #'   -fixPar contain all parameter values to be held fixed, if not specified 
  #'    then none are fixed
  #'   -estimated parameters: sigma & beta intercepts... what are these?!
  crwWrp <- function(track) {
    #'  Vector numbering unique tracks
    tracks <- 1:length(unique(track$ID))
    #'  Create list of crawlWrap arguments for each track
    #'  Remember: ID is a unique identifier for animal ID & track number
    theta <- fixPar <- list()
    for(i in unique(track$ID)) {
      theta[[i]] <- c(0, 0)
      fixPar[[i]] <- c(NA, NA)
    }
    #'  Interpolate missing locations within each track
    crwOut <- crawlWrap(obsData = track[which(track$ID %in% unique(track$ID)[tracks]),], 
                        theta = theta, fixPar = fixPar, attempts = 100, 
                        Time.name = "time", timeStep = "4 hours", coord = c("x", "y"))
    
    # return(crwOUT) # LOOK INTO RUNNING THIS IN PARALLEL
  }
  #'  Interpolate missing locations for each species
  #'  KEPP TRACK OF LIST ORDER HERE!!! Predators are ordered by study area then season
  crwOut_MD_smr <- crwWrp(spp_all_tracks[[1]]) #MD_smr_track
  crwOut_MD_wtr <- crwWrp(spp_all_tracks[[2]]) #MD_wtr_track
  crwOut_ELK_smr <- crwWrp(spp_all_tracks[[3]]) #ELK_smr_track
  crwOut_ELK_wtr <- crwWrp(spp_all_tracks[[4]]) #ELK_wtr_track
  crwOut_WTD_smr <- crwWrp(spp_all_tracks[[5]]) #WTD_smr_track
  crwOut_WTD_wtr <- crwWrp(spp_all_tracks[[6]]) #WTD_wtr_track
  crwOut_COUG_smr_OK <- crwWrp(spp_all_tracks[[7]]) #COUG_smr_track_OK
  crwOut_COUG_wtr_OK <- crwWrp(spp_all_tracks[[8]]) # COUG_wtr_track_OK
  crwOut_COUG_smr_NE <- crwWrp(spp_all_tracks[[9]]) #COUG_smr_track_NE
  crwOut_COUG_wtr_NE <- crwWrp(spp_all_tracks[[10]]) # COUG_wtr_track_NE
  crwOut_WOLF_smr_OK <- crwWrp(spp_all_tracks[[11]]) #WOLF_smr_track_OK
  crwOut_WOLF_wtr_OK <- crwWrp(spp_all_tracks[[12]]) #WOLF_wtr_track_OK
  crwOut_WOLF_smr_NE <- crwWrp(spp_all_tracks[[13]]) #WOLF_smr_track_NE
  crwOut_WOLF_wtr_NE <- crwWrp(spp_all_tracks[[14]]) #WOLF_wtr_track_NE
  crwOut_BOB_smr_OK <- crwWrp(spp_all_tracks[[15]]) #BOB_smr_track_OK
  crwOut_BOB_wtr_OK <- crwWrp(spp_all_tracks[[16]]) #BOB_wtr_track_OK
  crwOut_BOB_smr_NE <- crwWrp(spp_all_tracks[[17]]) #BOB_smr_track_NE
  crwOut_BOB_wtr_NE <- crwWrp(spp_all_tracks[[18]]) #BOB_wtr_track_NE
  crwOut_COY_smr_OK <- crwWrp(spp_all_tracks[[19]]) #COY_smr_track_OK
  crwOut_COY_wtr_OK <- crwWrp(spp_all_tracks[[20]]) #COY_wtr_track_OK
  crwOut_COY_smr_NE <- crwWrp(spp_all_tracks[[21]]) #COY_smr_track_NE
  crwOut_COY_wtr_NE <- crwWrp(spp_all_tracks[[22]]) #COY_wtr_track_NE
  
  #'  View interpolated data and new data (step length and turning angle)
  md_move_smr <- crwOut_MD_smr[[2]]
  md_move_wtr <- crwOut_MD_wtr[[2]]
  elk_move_smr <- crwOut_ELK_smr[[2]]
  elk_move_wtr <- crwOut_ELK_wtr[[2]]
  wtd_move_smr <- crwOut_WTD_smr[[2]]
  wtd_move_wtr <- crwOut_WTD_wtr[[2]]
  coug_move_smr_OK <- crwOut_COUG_smr_OK[[2]]
  coug_move_wtr_OK <- crwOut_COUG_wtr_OK[[2]]
  coug_move_smr_NE <- crwOut_COUG_smr_NE[[2]]
  coug_move_wtr_NE <- crwOut_COUG_wtr_NE[[2]]
  wolf_move_smr_OK <- crwOut_WOLF_smr_OK[[2]]
  wolf_move_wtr_OK <- crwOut_WOLF_wtr_OK[[2]]
  wolf_move_smr_NE <- crwOut_WOLF_smr_NE[[2]]
  wolf_move_wtr_NE <- crwOut_WOLF_wtr_NE[[2]]
  bob_move_smr_OK <- crwOut_BOB_smr_OK[[2]]
  bob_move_wtr_OK <- crwOut_BOB_wtr_OK[[2]]
  bob_move_smr_NE <- crwOut_BOB_smr_NE[[2]]
  bob_move_wtr_NE <- crwOut_BOB_wtr_NE[[2]]
  coy_move_smr_OK <- crwOut_COY_smr_OK[[2]]
  coy_move_wtr_OK <- crwOut_COY_wtr_OK[[2]]
  coy_move_smr_NE <- crwOut_COY_smr_NE[[2]]
  coy_move_wtr_NE <- crwOut_COY_wtr_NE[[2]]
  
  #'  Save individual crwOut datasets
  save(crwOut_MD_smr, file = "./Outputs/Telemetry_crwOut/crwOut_MD_smr.RData")
  save(crwOut_MD_wtr, file = "./Outputs/Telemetry_crwOut/crwOut_MD_wtr.RData")
  save(crwOut_ELK_smr, file = "./Outputs/Telemetry_crwOut/crwOut_ELK_smr.RData")
  save(crwOut_ELK_wtr, file = "./Outputs/Telemetry_crwOut/crwOut_ELK_wtr.RData")
  save(crwOut_WTD_smr, file = "./Outputs/Telemetry_crwOut/crwOut_WTD_smr.RData")
  save(crwOut_WTD_wtr, file = "./Outputs/Telemetry_crwOut/crwOut_WTD_wtr.RData")
  save(crwOut_COUG_smr_OK, file = "./Outputs/Telemetry_crwOut/crwOut_COUG_smr_OK.RData")
  save(crwOut_COUG_wtr_OK, file = "./Outputs/Telemetry_crwOut/crwOut_COUG_wtr_OK.RData")
  save(crwOut_COUG_smr_NE, file = "./Outputs/Telemetry_crwOut/crwOut_COUG_smr_NE.RData")
  save(crwOut_COUG_wtr_NE, file = "./Outputs/Telemetry_crwOut/crwOut_COUG_wtr_NE.RData")
  save(crwOut_WOLF_smr_OK, file = "./Outputs/Telemetry_crwOut/crwOut_WOLF_smr_OK.RData")
  save(crwOut_WOLF_wtr_OK, file = "./Outputs/Telemetry_crwOut/crwOut_WOLF_wtr_OK.RData")
  save(crwOut_WOLF_smr_NE, file = "./Outputs/Telemetry_crwOut/crwOut_WOLF_smr_NE.RData")
  save(crwOut_WOLF_wtr_NE, file = "./Outputs/Telemetry_crwOut/crwOut_WOLF_wtr_NE.RData")
  save(crwOut_BOB_smr_OK, file = "./Outputs/Telemetry_crwOut/crwOut_BOB_smr_OK.RData")
  save(crwOut_BOB_wtr_OK, file = "./Outputs/Telemetry_crwOut/crwOut_BOB_wtr_OK.RData")
  save(crwOut_BOB_smr_NE, file = "./Outputs/Telemetry_crwOut/crwOut_BOB_smr_NE.RData")
  save(crwOut_BOB_wtr_NE, file = "./Outputs/Telemetry_crwOut/crwOut_BOB_wtr_NE.RData")
  save(crwOut_COY_smr_OK, file = "./Outputs/Telemetry_crwOut/crwOut_COY_smr_OK.RData")
  save(crwOut_COY_wtr_OK, file = "./Outputs/Telemetry_crwOut/crwOut_COY_wtr_OK.RData")
  save(crwOut_COY_smr_NE, file = "./Outputs/Telemetry_crwOut/crwOut_COY_smr_NE.RData")
  save(crwOut_COY_wtr_NE, file = "./Outputs/Telemetry_crwOut/crwOut_COY_wtr_NE.RData")
  
  #'  List all crwOut datasets together
  crwOut_ALL <- list(crwOut_MD_smr, crwOut_MD_wtr, crwOut_ELK_smr, crwOut_ELK_wtr,
                     crwOut_WTD_smr, crwOut_WTD_wtr, crwOut_COUG_smr_OK, crwOut_COUG_wtr_OK,
                     crwOut_COUG_smr_NE, crwOut_COUG_wtr_NE, crwOut_WOLF_smr_OK, crwOut_WOLF_wtr_OK,
                     crwOut_WOLF_smr_NE, crwOut_WOLF_wtr_NE, crwOut_BOB_smr_OK, crwOut_BOB_wtr_OK,
                     crwOut_BOB_smr_NE, crwOut_BOB_wtr_NE, crwOut_COY_smr_OK, crwOut_COY_wtr_OK,
                     crwOut_COY_smr_NE, crwOut_COY_wtr_NE)
  
  #'  Save list of crwOut datasets
  save(crwOut_ALL, file  = paste0("./Outputs/Telemetry_crwOut/crwOut_ALL_", Sys.Date(), ".RData"))
  
  #'  Save workspace so I never need to rerun crawlWrap function again
  save.image(paste0("./Data/Collar_Movement_DataPrep_", Sys.Date(), ".RData"))
  
  
  ####  Run one species through without function  ####
  
  #' # rawELK1$time <- as.POSIXct(rawELK1$time, tz = "America/Los_Angeles")
  #' rawELK$time <- as.POSIXct(rawELK$time, tz = "America/Los_Angeles")
  #' 
  #' #'  Make locations spatial and project to UTM coordinates with study area projection
  #' # llcoord <- SpatialPoints(rawELK1[,5:6], proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  #' llcoord <- SpatialPoints(rawELK[,5:6], proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  #' utmcoord <- spTransform(llcoord, CRS("+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs "))
  #' #'  Add UTM locations to data frame
  #' # rawELK1$x <- attr(utmcoord, "coords")[,1]
  #' # rawELK1$y <- attr(utmcoord, "coords")[,2]
  #' rawELK$x <- attr(utmcoord, "coords")[,1]
  #' rawELK$y <- attr(utmcoord, "coords")[,2]
  #' 
  #' #'  Quick peak
  #' ggplot(rawELK[rawELK$ID == "1397ELK18",], aes(x = Long, y = Lat)) + 
  #'   geom_point() + 
  #'   geom_path() + 
  #'   facet_wrap(FullID~Season)
  #' 
  #' 
  #' #'  Tmax = 28.25 hours (in seconds) so that locations are still grouped in a 
  #' #'  single burst if there's a gap of 24hr or less in the data (up to 6 sequential   
  #' #'  fixes missed with 4 hours on each end) but a new burst if gap is >24 hrs.
  #' # rawELK1$burst <- creat.burst(data = rawELK1, id = TRUE, id_name = "ID", date_name = "time", Tmax = 87300)
  #' rawELK$burst <- creat.burst(data = rawELK, id = TRUE, id_name = "ID", date_name = "time", Tmax = 87300)   
  #' head(rawELK)
  #' tail(rawELK)
  #' # rawELK <- rawELK1
  #' #'  Summarize length the bursts
  #' res <- rawELK %>% group_by(burst) %>% summarise(Freq = n())
  #' #'  Look at bursts, ordered from shortest to longest
  #' res[order(res$Freq),]
  #' #'  Frequency of counts (5 1-point bursts, 4 2-point bursts, 0 3-point bursts...)
  #' #'  Keep in mind there are multiple bursts per individual animal, some shorter than others
  #' table(res$Freq)
  #' 
  #' #'  Exclude super short bursts (where burst length is 3 or less locations) 
  #' #'  because need at least 3 points to get a turning angle 
  #' ELK_burst <- rawELK[rawELK$burst %in% names(table(rawELK$burst))[table(rawELK$burst) >=3],]
  #' length(unique(ELK_burst$burst)) 
  #' head(ELK_burst)
  #' #'  Add burst value ids to Unique ID column to create unique IDs for each track
  #' ELK_burst$UniqueID <- with(ELK_burst, paste0(ID, "_", burst))
  #' ELK_burst <- ELK_burst %>%
  #'   mutate(
  #'     AnimalID = ID,
  #'     ID = UniqueID
  #'   ) %>%
  #'   dplyr::select(-UniqueID)
  #' #'  Vector for length of unique tracks
  #' tracks <- 1:length(unique(ELK_burst$ID))
  #' 
  #' 
  #' #'  Fit crawl model to interpolate missing fixes
  #' #'  theta are starting values, crawlWrap defaults to 0 if none are provided
  #' #'  fixPar contain all parameter values to be held fixed, if not specified then none are fixed
  #' #'  estimated parameters: sigma & beta intercepts... what are these?!
  #' #'  Vector for length of unique tracks
  #' tracks <- 1:length(unique(ELK_burst$ID))
  #' #'  List of parameters for each individual animal
  #' theta <- fixPar <- list()
  #' for(i in unique(ELK_burst$ID)) {
  #'   theta[[i]] <- c(0, 0)
  #'   fixPar[[i]] <- c(NA, NA)
  #' }
  #' crwOut_ELK <- crawlWrap(ELK_burst, theta = theta, fixPar = fixPar, attempts = 100, Time.name = "time", timeStep = "4 hours", coord = c("x", "y"))
  #' crwOut_ELK_burst <- crawlWrap(obsData = ELK_burst[which(ELK_burst$ID %in% unique(ELK_burst$ID)[tracks]),], theta = theta, fixPar = fixPar, attempts = 100, Time.name = "time", timeStep = "4 hours", coord = c("x", "y"))
  #' 
  #' tst <- crwOut_ELK[[2]]
  #' tstb <- crwOut_ELK_burst[[2]]
  #' 
  #' # crwOut_ELK1 <- crawlWrap(obsData = rawELK1, Time.name = "time", timeStep = "4 hours", coord = c("x", "y"), theta = c(0, 0), fixPar = c(NA, NA), retryFits = 0)
  #' # crwOut_ELK <- crawlWrap(obsData = ELK_burst[which(ELK_burst$ID %in% unique(ELK_burst$ID)[tracks]),], Time.name = "time", timeStep = "4 hours", coord = c("x", "y"), theta = c(0, 0), fixPar = c(NA, NA), retryFits = 100)
  #' # crwOut_ELK <- crawlWrap(obsData = ELK_burst, Time.name = "time", timeStep = "4 hours", coord = c("x", "y"), theta = c(0, 0), fixPar = c(NA, NA), initial.state=list(a=c(0,0),P = diag(c(5000 ^ 2, 10 * 3600 ^ 2))), retryFits = 100)
  #' 
  #' head(crwOut_ELK1[[1]])
  #' cW_ELK1 <- crwOut_ELK1[[2]]
   
  
  #' tracks <- 1:length(unique(WOLF_track$ID))
  #' #'  Create list of crawlWrap arguments for each track
  #' #'  Remember: ID is a unique identifier for animal ID & track number
  #' theta <- fixPar <- list()
  #' for(i in unique(WOLF_track$ID)) {
  #'   theta[[i]] <- c(0, 0)
  #'   fixPar[[i]] <- c(NA, NA)
  #' }
  #' #'  Interpolate missing locations within each track
  #' crwOut <- crawlWrap(obsData = WOLF_track[which(WOLF_track$ID %in% unique(WOLF_track$ID)[tracks]),], 
  #'                     theta = theta, fixPar = fixPar, attempts = 100, 
  #'                     Time.name = "time", timeStep = "4 hours", coord = c("x", "y"))
  #' tst <- crwOut[[2]]  
  
  