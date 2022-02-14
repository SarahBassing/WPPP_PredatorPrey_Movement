  #'  ============================================
  #'  Extract Covariates from Google Earth Engine
  #'  Sarah Bassing
  #'  WPPP
  #'  December 2021
  #'  ============================================
  #'  Script uses the rgee package to call the Google Earth Engine API to extract
  #'  spatially and temporally mapped covariate data for telemtry relocations.
  #'  This requires having an account with Google Earth Engine and Anaconda installed.
  #'  
  #'  Follow these links for instructions to install rgee package (it's complicated).
  #'  https://github.com/ricds/DL_RS_GEE/blob/main/rgee_install_packages.R
  #'  https://www.youtube.com/watch?v=1-k6wNL2hlo
  #'  
  #'  Follow this link as an example of how to extract data using rgee package
  #'  https://smithsonian.github.io/SpatiotemporalMatchingOfAnimalPositionsWithRemotelySensedDataUsingGoogleEarthEngineAndR/
  #'  
  #'  Using movement data prepared through crwOut function from momentuHMM package
  #'  so it includes telemetry relocations and interpolated relocations to be
  #'  used in HMM movement analyses (Collar_Movement_DataPrep.R)
  #'  ============================================
  
  #'  START A NEW SESSION IN R
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load packages and initialize Earth Engine
  library(rgee)
  library(sf)
  library(tidyr)
  library(dplyr)
  rgee::ee_Initialize(drive = T)  # ignore the Welcome to the Earth Engine client blurb

  #'  Load crwOut processed movement data
  load("./Outputs/Telemetry_crwOut/crwOut_ALL_2022-02-03.RData") #2021-12-08

  
  #'  Prepare telemetry data for ee
  data_prep <- function(crwOut_data) {
    #'  Retain dataframe with original & interpolated locations
    full_crwOut <- crwOut_data[[2]]
    #'  Create new column where date/time are formatted correctly
    full_crwOut$Date <- as.POSIXct(full_crwOut$time, format = "%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+8")
    #'  Reformat as a string
    full_crwOut$Date <- as.factor(full_crwOut$Date)
    #'  Reformat so that javascript can read the date (adds a "T" btwn date and time)
    full_crwOut$Date <- sub(" ", "T", full_crwOut$Date)
    #'  Add a unique ID for each observation so data can be joined correctly later
    full_crwOut$ID <- seq(1:nrow(full_crwOut))
    #'  Make spatial using sf (note: use coordinates labeled mu.x & mu.y so that
    #'  the animal relocations & interpolated locations are used)
    crwOut_sf <- st_as_sf(full_crwOut, coords = c('mu.x','mu.y'), crs = "+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs ") 
    #'  Transform projection to WGS84 for Google Earth Engine
    datasf <- st_transform(crwOut_sf, crs = 4326) 
    #'  Retain specific columns (things get wonky if there are additional columns)
    datasf <- dplyr::select(datasf, c("AnimalID", "Season", "Date", "ID", "geometry"))
    
    return(datasf)
  }
  #'  Format all location data for each species and season 
  data_ee_list <- lapply(crwOut_ALL, data_prep)
  
  #'  Check it out
  head(data_ee_list[[5]])
  str(data_ee_list[[5]])

  
  ####  Match & Extract Google Earth Engine image data to telemetry data  ####
  #'  Original functions written by Crego et al. (2021) examples
  #'  These original functions are now embedded in a monster function that requires
  #'  the data set, EE image collection, temporal window, image band, spatial & 
  #'  temporal resolutions of the EE image to be defined.
  match_ee_data <- function(datasf, imagecoll, tempwin, band, sp.res, tmp.res){  
    
    #'  Function to add property with time in milliseconds
    add_date <- function(feature) {
      date <- ee$Date(ee$String(feature$get("Date")))$millis()
      feature$set(list(date_millis = date))
    }
    
    #'  Join EE image & telemetry locations based on a maxDifference Filter within 
    #'  a specified temporal window. Set temporal window in days for filter. This 
    #'  will depend on the remote sensing data used.
    tempwin <- tempwin
    
    #'  Create the filter
    maxDiffFilter <- ee$Filter$maxDifference(
      difference = tempwin*24*60*60*1000, # days * hr * min * sec * milliseconds
      leftField = "date_millis", # Timestamp of the telemetry data
      rightField ="system:time_start" # EE image date
    )
    
    #'  Define how the EE image and telemetry data are joined using the saveBest 
    #'  function, which finds the EE image that best matches the filter (i.e., 
    #'  the EE image closest in time to the particular GPS fix location). 
    #'  This is baller.
    saveBestJoin <- ee$Join$saveBest(
      matchKey = "bestImage",
      measureKey = "timeDiff"
    )
    
    #'  Function to add property with raster pixel value from the matched EE image
    add_value <- function(feature){
      #'  Get the EE image selected by the join
      img1 <- ee$Image(feature$get("bestImage"))$select(band)
      #'  Extract geometry from the feature
      point <- feature$geometry()
      #'  Get pixel value for each point at the desired spatial resolution 
      #'  (adjust this using the scale argument- consider the resolution of the
      #'  EE image!). Also adjust the temporal resolution (tileScale) base on
      #'  temporal res of the EE image????
      pixel_value <- img1$sample(region = point, scale = sp.res, tileScale = tmp.res, dropNulls = F)
      #'  Return the data containing pixel value and image date
      feature$setMulti(list(PixelVal = pixel_value$first()$get(band), DateTimeImage = img1$get('system:index')))
    }
    
    #'  Function to remove image property from features
    removeProperty <- function(feature) {
      #'  Get the properties of the data
      properties = feature$propertyNames()
      #'  Select all items except images
      selectProperties = properties$filter(ee$Filter$neq("item", "bestImage"))
      #'  Return selected features
      feature$select(selectProperties)
    }

    #'  Define EE image collection you want to extract from
    #'  Try to un-hardcode this eventually
    imagecoll <- imagecoll
    
    #'  Name the data band to use (based on EE image)
    #'  Try to un-hardcode this eventually
    band <- band 

    #'  Chunk location data into groups to loop through when extracting pixel values
    #'  This is necessary if using the getInfo argument when converting location
    #'  data to sf object (in ee_as_sf). This is not necessary if using the drive
    #'  or gsc arguments which export data through Google Drive/Cloud. See section
    #'  0.5 in Crego et al. (2021) example for further details
    #'  
    #'  This is for up to 1 million points. To increase the max number of points, 
    #'  increase the value for max repetitions. To change the number of points to 
    #'  run per time, change the value in the argument each (up to 5000).
    datasf$uniq <- rep(1:1000, each = 1000)[1:nrow(datasf)] 
    
    #'  Track amount of time it takes to extract data
    start_time <- Sys.time()
    
    #'  Create empty dataframe to hold extracted pixel values
    dataoutput <- data.frame()
    
    #'  Call functions defined above to extract pixel values from EE images and 
    #'  loop through location data in chunks
    for(x in unique(datasf$uniq)){
      #'  Chunk location data
      data1 <- datasf %>% filter(uniq == x)
      #'  Send sf to GEE
      data <- sf_as_ee(data1)
      #'  Transform day into milliseconds
      data <- data$map(add_date)
      #'  Apply the join
      Data_match <- saveBestJoin$apply(data, imagecoll, maxDiffFilter)
      #'  Add pixel value to the data
      DataFinal <- Data_match$map(add_value)
      #'  Remove image property from the data
      DataFinal <- DataFinal$map(removeProperty)
      #'  Move GEE object into R
      temp <- ee_as_sf(DataFinal, via = 'getInfo')
      #'  Append
      dataoutput <- rbind(dataoutput, temp)
    }
    end_time <- Sys.time()
    #'  How long did that take?
    print(end_time - start_time)
    
    #'  Rename column with pixel values based on the defined band name
    colnames(dataoutput)[colnames(dataoutput) == "PixelVal"] <- band
    
    
    #'  View & return output
    print(dataoutput)
    return(dataoutput)
  }
  #'  Define time window of interest for extracting data
  start <- "2018-06-30"
  end <- "2021-03-01" 
  #'  EE image collections (https://developers.google.com/earth-engine/datasets)
  imageNDVI <- ee$ImageCollection('MODIS/006/MOD13Q1')$filterDate(start,end) # MODIS Terra NDVI/EVI 16day, 250m resolution
  # imageSNOWCOVER <- ee$ImageCollection("MODIS/006/MOD10A1")$filterDate(start,end) # MODIS Terra Snow Cover Daily Global (Normalized Difference Snow Index (NDSI)) daily, 500m resolution, values represent % coverage per pixel?
  # imageSNOWDEPTH <- ee$ImageCollection("NASA/GLDAS/V021/NOAH/G025/T3H")$filterDate(start,end) # NASA Global Land Data Assimilation System (GLDAS) daily, 27830m (or 0.25 degree) resolution, depth measured in meters
  # imageNDSI <- ee$ImageCollection("MODIS/MOD09GA_006_NDSI")$filterDate(start,end) # MODIS Terra Daily Normalized Difference Snow Index (NDSI) daily, 463.313 m resolution, missing data from Jan. 2019
  # imageSWE <- ee$ImageCollection("NASA/ORNL/DAYMET_V4")$filterDate(start,end) # Daymet V4: Daily Surface Weather and Climatological Summaries, Snow-Water Equivalent, daily, 1000m resolution, data up to Dec. 31, 2020
  
  #' #'  Run reformated animal location data through this monster function to
  #' #'  match & extract EE images 
  #' ee_NDVI <- lapply(data_ee_list, match_ee_data, imagecoll = imageNDVI, tempwin = 16, band = "NDVI", sp.res = 250, tmp.res = 16)
  #' ee_SNOWCOVER <- lapply(data_ee_list, match_ee_data, imagecoll = imageSNOWCOVER, tempwin = 1, band = "NDSI_Snow_Cover", sp.res = 500, tmp.res = 1)
  #' ee_SNOWDEPTH <- lapply(data_ee_list, match_ee_data, imagecoll = imageSNOWDEPTH, tempwin = 1, band = "SnowDepth_inst", sp.res = 27830, tmp.res = 1)
  #' ee_NDSI <- lapply(data_ee_list, match_ee_data, imagecoll = imageNDSI, tempwin = 1, band = "NDSI", sp.res = 463.313, tmp.res = 1)
  #' ee_SWE <- lapply(data_ee_list, match_ee_data, imagecoll = imageSWE, tempwin = 1, band = "swe", sp.res = 1000, tmp.res = 1)
  
  #'  Run re-formated summer location data through this monster function to
  #'  match & extract EE image data
  ee_NDVI_md_smr <- match_ee_data(data_ee_list[[1]], imagecoll = imageNDVI, tempwin = 16, band = "NDVI", sp.res = 250, tmp.res = 16)
  ee_NDVI_elk_smr <- match_ee_data(data_ee_list[[3]], imagecoll = imageNDVI, tempwin = 16, band = "NDVI", sp.res = 250, tmp.res = 16)
  ee_NDVI_wtd_smr <- match_ee_data(data_ee_list[[5]], imagecoll = imageNDVI, tempwin = 16, band = "NDVI", sp.res = 250, tmp.res = 16)
  ee_NDVI_coug_smr_OK <- match_ee_data(data_ee_list[[7]], imagecoll = imageNDVI, tempwin = 16, band = "NDVI", sp.res = 250, tmp.res = 16)
  ee_NDVI_coug_smr_NE <- match_ee_data(data_ee_list[[9]], imagecoll = imageNDVI, tempwin = 16, band = "NDVI", sp.res = 250, tmp.res = 16)
  ee_NDVI_wolf_smr_OK <- match_ee_data(data_ee_list[[11]], imagecoll = imageNDVI, tempwin = 16, band = "NDVI", sp.res = 250, tmp.res = 16)
  ee_NDVI_wolf_smr_NE <- match_ee_data(data_ee_list[[13]], imagecoll = imageNDVI, tempwin = 16, band = "NDVI", sp.res = 250, tmp.res = 16)
  ee_NDVI_bob_smr_OK <- match_ee_data(data_ee_list[[15]], imagecoll = imageNDVI, tempwin = 16, band = "NDVI", sp.res = 250, tmp.res = 16)
  ee_NDVI_bob_smr_NE <- match_ee_data(data_ee_list[[17]], imagecoll = imageNDVI, tempwin = 16, band = "NDVI", sp.res = 250, tmp.res = 16)
  ee_NDVI_coy_smr_OK <- match_ee_data(data_ee_list[[19]], imagecoll = imageNDVI, tempwin = 16, band = "NDVI", sp.res = 250, tmp.res = 16)
  ee_NDVI_coy_smr_NE <- match_ee_data(data_ee_list[[21]], imagecoll = imageNDVI, tempwin = 16, band = "NDVI", sp.res = 250, tmp.res = 16)

  #'  List together for safe keeping
  ee_smr_NDVI_list <- list(ee_NDVI_md_smr, ee_NDVI_elk_smr, ee_NDVI_wtd_smr, 
                           ee_NDVI_coug_smr_OK, ee_NDVI_coug_smr_NE, ee_NDVI_wolf_smr_OK, 
                           ee_NDVI_wolf_smr_NE, ee_NDVI_bob_smr_OK, ee_NDVI_bob_smr_NE,
                           ee_NDVI_coy_smr_OK, ee_NDVI_coy_smr_NE)
  
  save(ee_smr_NDVI_list, file = paste0("./Outputs/Telemetry_covs/ee_smr_NDVI_list_", Sys.Date(), ".RData"))
  
  
  ####  Re-scale NDVI values  ####
  #'  ONLY RUN IF EXTRACTING NDVI DATA!
  #'  MODIS data are scaled by a factor of 0.0001 for ease of extraction (see 
  #'  pg. 9 of MODIS User's Guide) but NDVI should range -1 to 1. All extracted 
  #'  NDVI values need to be rescaled by 0.0001 to be in the correct range.
  ndvi_rescale <- function(ee_data, plotit = F) {
    #'  Create new column with re-scaled NDVI values
    dataoutput <- mutate(ee_data, NDVI_scale = NDVI*0.0001,
                         #'  For some reason NAs aren't recognized as NAs here so need to fix that
                         AnimalID = ifelse(AnimalID == "NA", "jerk", AnimalID),
                         AnimalID = ifelse(AnimalID == "jerk", NA, AnimalID),
                         Season = ifelse(Season == "NA", "jerk", Season),
                         Season = ifelse(Season == "jerk", NA, Season)) %>%
                  #'  Fill in missing values with the AnimalID/Season above it
                  fill(AnimalID, .direction = "down") %>%
                  fill(Season, .direction = "down")
    #'  Visualize
    plot(dataoutput$NDVI_scale)
    hist(dataoutput$NDVI_scale)
    
    return(dataoutput)
  }
  ee_NDVI <- lapply(ee_smr_NDVI_list, ndvi_rescale, T)
  
  
  ####  Maximum NDVI  ####
  #'  ======================================================================
  #'  Function to extract MAXIMUM annual NDVI, representing potential forage  
  #'  quality available post-growing season at winter locations.
  find_maxNDVI <- function(crwOut_data, eeimage, start.date, end.date, band.name, band, ee.scale) {
    
    #'  Define EE ImageCollection
    eendvi <- ee$ImageCollection(eeimage) %>%  
      #'  Define date range of interest
      ee$ImageCollection$filterDate(start.date, end.date) %>% 
      ee$ImageCollection$map(
        function(x) {
          date <- ee$Date(x$get("system:time_start"))$format('YYYY_MM_dd')
          name <- ee$String$cat(band.name, date)  
          x$select(band)$rename(name)   
        }
      )
    
    #'  Make telemetry locations spatial using sf (note: use coordinates labeled 
    #'  mu.x & mu.y so that interpolated locations are included)
    crwOut_sf <- st_as_sf(crwOut_data, coords = c('mu.x','mu.y'),  
                          crs = "+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs ") %>%
      #'  Fill in missing values with the AnimalID/Season/StudyArea above it
      fill(AnimalID, .direction = "down") %>%
      fill(Season, .direction = "down") %>%
      fill(StudyArea, .direction = "down")
    #'  Transform projection to WGS84 for Google Earth Engine
    datasf <- st_transform(crwOut_sf, crs = 4326)
    #'  Make a unique ID for every observation 
    datasf$uid <- as.integer(1:nrow(datasf))
    #'  Reformat date/time to a character string
    datasf <- mutate(datasf, times = as.character(time))
    #'  Create a unique ID for each observation based on AnimalID & location time
    #'  Necessary for ee_extract! Use "NAME" for column header
    datasf$NAME <- paste0(datasf$AnimalID, "_", datasf$times)
    # datasf$NAME <- paste0(datasf$AnimalID, "_", datasf$uid)

    #'  Track amount of time it takes to extract data
    start_time <- Sys.time()
    
    EE_allndvi <- ee_extract(
      x = eendvi,
      y = datasf["NAME"],
      scale = ee.scale,
      fun = ee$Reducer$max(),
      via = "drive",
      lazy = TRUE,
      sf = TRUE
    )
    EE_allndvi <- EE_allndvi %>% ee_utils_future_value()
    
    end_time <- Sys.time()
    #'  How long did that take?
    print(end_time - start_time)
    #'  Quick peak
    print(EE_allndvi[1:6,1:4])
    
    #'  Find maximum NDVI values during annual growing season (April - Sept)
    ee_allNDVI <- as.data.frame(EE_allndvi) %>% 
      relocate(geometry, .after = NAME) 
    #'  apply max function row-wise (1 = across each row), excluding 1st 2 columns of df
    ee_allNDVI$maxNDVI <- apply(ee_allNDVI[-c(1:2)], 1, max, na.rm = T)
    ee_allNDVI <- ee_allNDVI %>%
      mutate(
        maxNDVI_scale = maxNDVI*0.0001,
        ID = seq(1:nrow(ee_allNDVI))
      ) %>%
      dplyr::select(c(ID, NAME, maxNDVI_scale, geometry))
    
    #'  Pull out location dates (including interpolated ones)
    # times <- as.Date(crwOut_sf$time, "%Y-%m-%d", tz = "Etc/GMT+8")
    # date.only <- as.data.frame(times)
    TIMES <- as.POSIXct(crwOut_sf$time, format = "%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+8")
    TIMES.only <- as.data.frame(TIMES)
    
    #'  Save only relevant location info
    crwOut_sf$ID2 <- as.integer(1:nrow(crwOut_sf))
    crwOut_sf <- cbind(crwOut_sf, TIMES.only) #date.only
    crwOut_df <- crwOut_sf %>%
      mutate(mu.x = sf::st_coordinates(.)[,1],
             mu.y = sf::st_coordinates(.)[,2]) %>%
      dplyr::select(c(AnimalID, Season, StudyArea, TIMES, mu.x, mu.y, ID2, geometry)) #times
    colnames(crwOut_df) <- c("AnimalID", "Season", "StudyArea", "times", "x", "y", "ID", "geometry") #"Dates"
    
    #'  Join animal location data with maxNDVI values
    maxNDVI_df <- as.data.frame(crwOut_df) %>%
      full_join(ee_allNDVI, by = "ID") %>%
      dplyr::select(-geometry.y)
    colnames(maxNDVI_df) <- c("AnimalID", "Season", "StudyArea", "times", "x", "y", "ID", "geometry", "NAME", "NDVImax") #"Dates"
    
    print(maxNDVI_df[1:6,])
    
    return(maxNDVI_df)
    
  }

  #'  Organize input data a bit more before running through max NDVI function
  #'  Function to drop crwFits lists for each species from crwOut_ALL list of lists 
  drop_list <- function(full.list) {
    short.list <- as.data.frame(full.list[-1]) %>%
      #'  Rename all columns by removing first 11 characters 
      #'  (annoyingly happens when covert each list to a data frame)
      rename_with(~ gsub("^...........", "", .x))
    return(short.list)
  }
  short.list <- lapply(crwOut_ALL, drop_list)
  
  #'  Re-ogranize short list by year, not species
  #'  Note: Summer lists should be all 0; Winter lists that are all 0 mean no 
  #'  collars deployed during that time period (e.g. bobcat wtr 1819, wolf wtr 1920)
  filter_wtr1819 <- function(mylist){
    df <- as.data.frame(mylist) %>%
      filter(time > "2018-11-30 00:00:00" & time < "2019-03-01 00:00:00") %>%
      mutate(Season = "Winter1819") # Pay attention to this!!!!
    return(df)
  }
  wtr1819 <- lapply(short.list, filter_wtr1819)
  filter_wtr1920 <- function(mylist){
    df <- as.data.frame(mylist) %>%
      filter(time > "2019-11-30 00:00:00" & time < "2020-03-01 00:00:00") %>%
      mutate(Season = "Winter1920") # Pay attention to this!!!! 
    return(df)
  }
  wtr1920 <- lapply(short.list, filter_wtr1920)
  filter_wtr2021 <- function(mylist){
    df <- as.data.frame(mylist) %>%
      filter(time > "2020-11-30 00:00:00") %>%
      mutate(Season = "Winter2021") # Pay attention to this!!!!
    return(df)
  }
  wtr2021 <- lapply(short.list, filter_wtr2021)
  
  #'  List data for each species per season (winter observations only)
  wtr_list1819 <- list(wtr1819[[2]], wtr1819[[4]], wtr1819[[6]], wtr1819[[8]], wtr1819[[10]], wtr1819[[12]], wtr1819[[14]], wtr1819[[16]], wtr1819[[18]], wtr1819[[20]], wtr1819[[22]])
  wtr_list1920 <- list(wtr1920[[2]], wtr1920[[4]], wtr1920[[6]], wtr1920[[8]], wtr1920[[10]], wtr1920[[12]], wtr1920[[14]], wtr1920[[16]], wtr1920[[18]], wtr1920[[20]], wtr1920[[22]])
  wtr_list2021 <- list(wtr2021[[2]], wtr2021[[4]], wtr2021[[6]], wtr2021[[8]], wtr2021[[10]], wtr2021[[12]], wtr2021[[14]], wtr2021[[16]], wtr2021[[18]], wtr2021[[20]], wtr2021[[22]])
  #'  New list order: mule deer [1], elk [2], white-tailed deer [3], cougar OK [4], cougar NE [5], wolf OK [6], wolf NE [7], bobcat OK [8], bobcat NE [9], coyote OK [10], coyote NE [11] 
  
  #'  Define start and end dates for EE imageColleciton
  #'  Starting date: beginning of growing season of first year of study (April)
  #'  Ending date: end of growing season of last year of study (Sept)
  start18 <- "2018-04-01"
  end18 <- "2018-10-01"
  start19 <- "2019-04-01"
  end19 <- "2019-10-01"
  start20 <- "2020-04-01"
  end20 <- "2020-10-01"
  
  #'  Define EE imageCollection parameters of interest
  eeimage <- "MODIS/006/MOD13Q1"
  band.name <- "NDVI_"
  band <- "NDVI"
  ee.scale <- 250
  
  #'  Find maximum NDVI for each species in each winter season
  #'  These are spit up by individual winter and species to reduce the amount of 
  #'  EE images that need to be extract from for each species - takes forever 
  #'  otherwise and hits data limits!
  #'  
  #'  List order: mule deer [1], elk [2], white-tailed deer [3], cougar OK [4], 
  #'  cougar NE [5], wolf OK [6], wolf NE [7], bobcat OK [8], bobcat NE [9], 
  #'  coyote OK [10], coyote NE [11]
  #'              
  #'  Split md_wtr1819 data even more because it's apparently too large?!?!
  nrow(wtr_list1819[[1]])
  md_wtr1819a <- wtr_list1819[[1]][1:20000,]
  md_wtr1819b <- wtr_list1819[[1]][20001:45872,]
  md_wtr1819a_NDVImax <- find_maxNDVI(md_wtr1819a, eeimage = eeimage, start.date = start18, end.date = end18, band.name = band.name, band = band, ee.scale = ee.scale)
  md_wtr1819b_NDVImax <- find_maxNDVI(md_wtr1819b, eeimage = eeimage, start.date = start18, end.date = end18, band.name = band.name, band = band, ee.scale = ee.scale)
  md_wtr1819_NDVImax <- rbind(md_wtr1819a_NDVImax, md_wtr1819b_NDVImax)
  md_wtr1920_NDVImax <- find_maxNDVI(wtr_list1920[[1]], eeimage = eeimage, start.date = start19, end.date = end19, band.name = band.name, band = band, ee.scale = ee.scale)
  md_wtr2021_NDVImax <- find_maxNDVI(wtr_list2021[[1]], eeimage = eeimage, start.date = start20, end.date = end20, band.name = band.name, band = band, ee.scale = ee.scale)
  
  elk_wtr1819_NDVImax <- find_maxNDVI(wtr_list1819[[2]], eeimage = eeimage, start.date = start18, end.date = end18, band.name = band.name, band = band, ee.scale = ee.scale)
  elk_wtr1920_NDVImax <- find_maxNDVI(wtr_list1920[[2]], eeimage = eeimage, start.date = start19, end.date = end19, band.name = band.name, band = band, ee.scale = ee.scale)
  elk_wtr2021_NDVImax <- find_maxNDVI(wtr_list2021[[2]], eeimage = eeimage, start.date = start20, end.date = end20, band.name = band.name, band = band, ee.scale = ee.scale)
  
  wtd_wtr1819_NDVImax <- find_maxNDVI(wtr_list1819[[3]], eeimage = eeimage, start.date = start18, end.date = end18, band.name = band.name, band = band, ee.scale = ee.scale)
  wtd_wtr1920_NDVImax <- find_maxNDVI(wtr_list1920[[3]], eeimage = eeimage, start.date = start19, end.date = end19, band.name = band.name, band = band, ee.scale = ee.scale)
  wtd_wtr2021_NDVImax <- find_maxNDVI(wtr_list2021[[3]], eeimage = eeimage, start.date = start20, end.date = end20, band.name = band.name, band = band, ee.scale = ee.scale)
  
  coug_wtr1819_NDVImax_OK <- find_maxNDVI(wtr_list1819[[4]], eeimage = eeimage, start.date = start18, end.date = end18, band.name = band.name, band = band, ee.scale = ee.scale)
  coug_wtr1920_NDVImax_OK <- find_maxNDVI(wtr_list1920[[4]], eeimage = eeimage, start.date = start19, end.date = end19, band.name = band.name, band = band, ee.scale = ee.scale)
  coug_wtr2021_NDVImax_OK <- find_maxNDVI(wtr_list2021[[4]], eeimage = eeimage, start.date = start20, end.date = end20, band.name = band.name, band = band, ee.scale = ee.scale)
  coug_wtr1819_NDVImax_NE <- find_maxNDVI(wtr_list1819[[5]], eeimage = eeimage, start.date = start18, end.date = end18, band.name = band.name, band = band, ee.scale = ee.scale)
  coug_wtr1920_NDVImax_NE <- find_maxNDVI(wtr_list1920[[5]], eeimage = eeimage, start.date = start19, end.date = end19, band.name = band.name, band = band, ee.scale = ee.scale)
  coug_wtr2021_NDVImax_NE <- find_maxNDVI(wtr_list2021[[5]], eeimage = eeimage, start.date = start20, end.date = end20, band.name = band.name, band = band, ee.scale = ee.scale)
  
  wolf_wtr1819_NDVImax_OK <- find_maxNDVI(wtr_list1819[[6]], eeimage = eeimage, start.date = start18, end.date = end18, band.name = band.name, band = band, ee.scale = ee.scale)
  wolf_wtr1920_NDVImax_OK <- find_maxNDVI(wtr_list1920[[6]], eeimage = eeimage, start.date = start19, end.date = end19, band.name = band.name, band = band, ee.scale = ee.scale)
  wolf_wtr2021_NDVImax_OK <- find_maxNDVI(wtr_list2021[[6]], eeimage = eeimage, start.date = start20, end.date = end20, band.name = band.name, band = band, ee.scale = ee.scale)
  wolf_wtr1819_NDVImax_NE <- find_maxNDVI(wtr_list1819[[7]], eeimage = eeimage, start.date = start18, end.date = end18, band.name = band.name, band = band, ee.scale = ee.scale)
  wolf_wtr1920_NDVImax_NE <- find_maxNDVI(wtr_list1920[[7]], eeimage = eeimage, start.date = start19, end.date = end19, band.name = band.name, band = band, ee.scale = ee.scale)
  wolf_wtr2021_NDVImax_NE <- find_maxNDVI(wtr_list2021[[7]], eeimage = eeimage, start.date = start20, end.date = end20, band.name = band.name, band = band, ee.scale = ee.scale)
  
  bob_wtr1819_NDVImax_OK <- find_maxNDVI(wtr_list1819[[8]], eeimage = eeimage, start.date = start18, end.date = end18, band.name = band.name, band = band, ee.scale = ee.scale)
  bob_wtr1920_NDVImax_OK <- find_maxNDVI(wtr_list1920[[8]], eeimage = eeimage, start.date = start19, end.date = end19, band.name = band.name, band = band, ee.scale = ee.scale)
  bob_wtr2021_NDVImax_OK <- find_maxNDVI(wtr_list2021[[8]], eeimage = eeimage, start.date = start20, end.date = end20, band.name = band.name, band = band, ee.scale = ee.scale)
  bob_wtr1819_NDVImax_NE <- find_maxNDVI(wtr_list1819[[9]], eeimage = eeimage, start.date = start18, end.date = end18, band.name = band.name, band = band, ee.scale = ee.scale)
  bob_wtr1920_NDVImax_NE <- find_maxNDVI(wtr_list1920[[9]], eeimage = eeimage, start.date = start19, end.date = end19, band.name = band.name, band = band, ee.scale = ee.scale)
  bob_wtr2021_NDVImax_NE <- find_maxNDVI(wtr_list2021[[9]], eeimage = eeimage, start.date = start20, end.date = end20, band.name = band.name, band = band, ee.scale = ee.scale)
  
  coy_wtr1819_NDVImax_OK <- find_maxNDVI(wtr_list1819[[10]], eeimage = eeimage, start.date = start18, end.date = end18, band.name = band.name, band = band, ee.scale = ee.scale)
  coy_wtr1920_NDVImax_OK <- find_maxNDVI(wtr_list1920[[10]], eeimage = eeimage, start.date = start19, end.date = end19, band.name = band.name, band = band, ee.scale = ee.scale)
  coy_wtr2021_NDVImax_OK <- find_maxNDVI(wtr_list2021[[10]], eeimage = eeimage, start.date = start20, end.date = end20, band.name = band.name, band = band, ee.scale = ee.scale)
  coy_wtr1819_NDVImax_NE <- find_maxNDVI(wtr_list1819[[11]], eeimage = eeimage, start.date = start18, end.date = end18, band.name = band.name, band = band, ee.scale = ee.scale)
  coy_wtr1920_NDVImax_NE <- find_maxNDVI(wtr_list1920[[11]], eeimage = eeimage, start.date = start19, end.date = end19, band.name = band.name, band = band, ee.scale = ee.scale)
  coy_wtr2021_NDVImax_NE <- find_maxNDVI(wtr_list2021[[11]], eeimage = eeimage, start.date = start20, end.date = end20, band.name = band.name, band = band, ee.scale = ee.scale)
  
  #'  Bind annual max NDVI data together for each species
  md_NDVImax <- rbind(md_wtr1819_NDVImax, md_wtr1920_NDVImax, md_wtr2021_NDVImax)
  elk_NDVImax <- rbind(elk_wtr1819_NDVImax, elk_wtr1920_NDVImax, elk_wtr2021_NDVImax)
  wtd_NDVImax <- rbind(wtd_wtr1819_NDVImax, wtd_wtr1920_NDVImax, wtd_wtr2021_NDVImax)
  coug_NDVImax_OK <- rbind(coug_wtr1819_NDVImax_OK, coug_wtr1920_NDVImax_OK, coug_wtr2021_NDVImax_OK)
  coug_NDVImax_NE <- rbind(coug_wtr1819_NDVImax_NE, coug_wtr1920_NDVImax_NE, coug_wtr2021_NDVImax_NE)
  wolf_NDVImax_OK <- rbind(wolf_wtr1819_NDVImax_OK, wolf_wtr1920_NDVImax_OK, wolf_wtr2021_NDVImax_OK)
  wolf_NDVImax_NE <- rbind(wolf_wtr1819_NDVImax_NE, wolf_wtr2021_NDVImax_NE) # pay attention to what happens here with no wolf_wtr1920_NDVImax_NE data
  bob_NDVImax_OK <- rbind(bob_wtr1819_NDVImax_OK, bob_wtr1920_NDVImax_OK, bob_wtr2021_NDVImax_OK)
  bob_NDVImax_NE <- rbind(bob_wtr1920_NDVImax_NE, bob_wtr2021_NDVImax_NE) # pay attention to what happens here with no bob_wtr1819_NDVImax_NE data
  coy_NDVImax_OK <- rbind(coy_wtr1819_NDVImax_OK, coy_wtr1920_NDVImax_OK, coy_wtr2021_NDVImax_OK)
  coy_NDVImax_NE <- rbind(coy_wtr1819_NDVImax_NE, coy_wtr1920_NDVImax_NE, coy_wtr2021_NDVImax_NE)
  
  ee_NDVImax_list <- list(md_NDVImax, elk_NDVImax, wtd_NDVImax, coug_NDVImax_OK, 
                          coug_NDVImax_NE, wolf_NDVImax_OK, wolf_NDVImax_NE, 
                          bob_NDVImax_OK, bob_NDVImax_NE, coy_NDVImax_OK, coy_NDVImax_NE)
  
  save(md_NDVImax, file = paste0("./Outputs/Telemetry_covs/md_NDVImax_", Sys.Date(), ".RData"))
  save(elk_NDVImax, file = paste0("./Outputs/Telemetry_covs/elk_NDVImax_", Sys.Date(), ".RData"))
  save(wtd_NDVImax, file = paste0("./Outputs/Telemetry_covs/wtd_NDVImax_", Sys.Date(), ".RData"))
  save(coug_NDVImax_OK, file = paste0("./Outputs/Telemetry_covs/coug_NDVImax_OK_", Sys.Date(), ".RData"))
  save(coug_NDVImax_NE, file = paste0("./Outputs/Telemetry_covs/coug_NDVImax_NE_", Sys.Date(), ".RData"))
  save(wolf_NDVImax_OK, file = paste0("./Outputs/Telemetry_covs/wolf_NDVImax_OK_", Sys.Date(), ".RData"))
  save(wolf_NDVImax_NE, file = paste0("./Outputs/Telemetry_covs/wolf_NDVImax_NE_", Sys.Date(), ".RData"))
  save(bob_NDVImax_OK, file = paste0("./Outputs/Telemetry_covs/bob_NDVImax_OK_", Sys.Date(), ".RData"))
  save(bob_NDVImax_NE, file = paste0("./Outputs/Telemetry_covs/bob_NDVImax_NE_", Sys.Date(), ".RData"))
  save(coy_NDVImax_OK, file = paste0("./Outputs/Telemetry_covs/coy_NDVImax_OK_", Sys.Date(), ".RData"))
  save(coy_NDVImax_NE, file = paste0("./Outputs/Telemetry_covs/coy_NDVImax_NE_", Sys.Date(), ".RData"))
  
  
  save(ee_NDVImax_list, file = paste0("./Outputs/Telemetry_covs/ee_NDVImax_list_", Sys.Date(), ".RData"))
  
  
  ####  Join Datasets  ####
  # load("G:/My Drive/1_Repositories/WPPP_PredatorPrey_Movement/Outputs/Telemetry_covs/ee_smr_NDVI_list_2022-02-11.RData")
  load("G:/My Drive/1_Repositories/WPPP_PredatorPrey_Movement/Outputs/Telemetry_covs/ee_NDVImax_list_2022-02-14.RData")
  
  #'  Join spatially & temporally matched NDVI data to each summer location
  join_NDVI <- function(crwOut_data, ndvi) { 
    full_crwOut <- crwOut_data[[2]]
    #'  Make sure each observation has the unique Animal ID
    crwOut <- as.data.frame(full_crwOut) %>%
      separate(ID, sep = "_", into = "AnimalID",  extra = "drop") 
    #'  Create a new dataframe (drops the crwPredict structure)
    ee_covs <- as.data.frame(cbind(crwOut$AnimalID, crwOut$Season, crwOut$StudyArea, crwOut$mu.x, crwOut$mu.y))
    #'  Add date/times as a character 
    ee_covs$times <- as.character(crwOut$time)
    colnames(ee_covs) <- c("AnimalID", "Season", "StudyArea", "mu.x", "mu.y", "times")
    #'  Reformat NAs so they are recognized as NAs
    ee_covs <- mutate(ee_covs,
        AnimalID = ifelse(AnimalID == "NA", "jerk", AnimalID),
        AnimalID = ifelse(AnimalID == "jerk", NA, AnimalID),
        Season = ifelse(Season == "NA", "jerk", Season),
        Season = ifelse(Season == "jerk", NA, Season),
        StudyArea = ifelse(StudyArea == "NA", "jerk", StudyArea),
        StudyArea = ifelse(StudyArea == "jerk", NA, StudyArea)) %>%
      #'  Fill in missing values with the AnimalID/Season/StudyArea above it
      fill(AnimalID, .direction = "down") %>%
      fill(Season, .direction = "down") %>%
      fill(StudyArea, .direction = "down")
    #'  Reformat times in NDVI data so they can be matched to original location data
    ndvi$times <- as.character(ndvi$Date)
    #'  Join crwOut location data with NDVI data
    ee_covs <- ee_covs %>%
      full_join(ndvi, by = c("AnimalID", "times")) %>%
      dplyr::select(c(ID, AnimalID, Season.x, StudyArea, mu.x, mu.y, times, NDVI_scale))
    colnames(ee_covs) <- c("ID", "AnimalID", "Season", "StudyArea", "mu.x", "mu.y", "times", "NDVI")
    
    return(ee_covs)
  }
  #'  Join summer observations with spatially & temporally matched NDVI values
  #'  Pay attention to which list is being used for each dataset!
  md_eeNDVI_smr <- join_NDVI(crwOut_ALL[[1]], ee_NDVI[[1]]) 
  elk_eeNDVI_smr <- join_NDVI(crwOut_ALL[[3]], ee_NDVI[[2]]) 
  wtd_eeNDVI_smr <- join_NDVI(crwOut_ALL[[5]], ee_NDVI[[3]]) 
  coug_eeNDVI_smr_OK <- join_NDVI(crwOut_ALL[[7]], ee_NDVI[[4]]) 
  coug_eeNDVI_smr_NE <- join_NDVI(crwOut_ALL[[9]], ee_NDVI[[5]])
  wolf_eeNDVI_smr_OK <- join_NDVI(crwOut_ALL[[11]], ee_NDVI[[6]]) 
  wolf_eeNDVI_smr_NE <- join_NDVI(crwOut_ALL[[13]], ee_NDVI[[7]]) 
  bob_eeNDVI_smr_OK <- join_NDVI(crwOut_ALL[[15]], ee_NDVI[[8]]) 
  bob_eeNDVI_smr_NE <- join_NDVI(crwOut_ALL[[17]], ee_NDVI[[9]])
  coy_eeNDVI_smr_OK <- join_NDVI(crwOut_ALL[[19]], ee_NDVI[[10]])
  coy_eeNDVI_smr_NE <- join_NDVI(crwOut_ALL[[21]], ee_NDVI[[11]])
  
  #'  List all data together
  ee_NDVIsmr_list <- list(md_eeNDVI_smr, elk_eeNDVI_smr, wtd_eeNDVI_smr, 
                          coug_eeNDVI_smr_OK, coug_eeNDVI_smr_NE, wolf_eeNDVI_smr_OK, 
                          wolf_eeNDVI_smr_NE, bob_eeNDVI_smr_OK, bob_eeNDVI_smr_NE, 
                          coy_eeNDVI_smr_OK, coy_eeNDVI_smr_NE)
 
  #'  Join spatially matched MAXIMUM NDVI from previous growing season with winter locations
  join_NDVImax <- function(crwOut_data, ndvimax) { 
    full_crwOut <- crwOut_data[[2]]
    #'  Make sure each observation has the unique Animal ID
    crwOut <- as.data.frame(full_crwOut) %>%
      separate(ID, sep = "_", into = "AnimalID",  extra = "drop") 
    #'  Create a new dataframe (drops the crwPredict structure)
    ee_covs <- as.data.frame(cbind(crwOut$AnimalID, crwOut$Season, crwOut$StudyArea, crwOut$mu.x, crwOut$mu.y))
    #'  Add date/times as a character 
    ee_covs$times <- as.character(crwOut$time)
    colnames(ee_covs) <- c("AnimalID", "Season", "StudyArea", "mu.x", "mu.y", "times") 
    #'  Reformat NAs so they are recognized as NAs
    ee_covs <- mutate(ee_covs,
                      Season = ifelse(Season == "NA", "jerk", Season),
                      Season = ifelse(Season == "jerk", NA, Season),
                      StudyArea = ifelse(StudyArea == "NA", "jerk", StudyArea),
                      StudyArea = ifelse(StudyArea == "jerk", NA, StudyArea)) %>%
      #'  Fill in missing values with the Season and StudyArea above it
      fill(Season, .direction = "down") %>%
      fill(StudyArea, .direction = "down")
    #'  Reformat times in NDVI data so they can be matched to original location data
    ndvimax$times <- as.character(ndvimax$times)
    #'  Join crwOut location data with max NDVI data
    ee_covs <- ee_covs %>%
      #'  That ID is necessary since Dates lack the time and aren't unique
      full_join(ndvimax, by = c("AnimalID", "times")) %>% 
      dplyr::select(c(ID, AnimalID, Season.x, StudyArea.x, mu.x, mu.y, times, NAME, NDVImax))  
    colnames(ee_covs) <- c("ID", "AnimalID", "Season", "StudyArea", "mu.x", "mu.y", "times", "NAME", "maxNDVI") #"Dates"
    
    return(ee_covs)
  }
  #'  Join winter observations with spatially matched max NDVI values from previous growing season
  #'  Pay attention to which list is being used for each dataset!
  md_eeNDVImax_wtr <- join_NDVImax(crwOut_ALL[[2]], ee_NDVImax_list[[1]]) 
  elk_eeNDVImax_wtr <- join_NDVImax(crwOut_ALL[[4]], ee_NDVImax_list[[2]]) 
  wtd_eeNDVImax_wtr <- join_NDVImax(crwOut_ALL[[6]], ee_NDVImax_list[[3]]) 
  coug_eeNDVImax_wtr_OK <- join_NDVImax(crwOut_ALL[[8]], ee_NDVImax_list[[4]]) 
  coug_eeNDVImax_wtr_NE <- join_NDVImax(crwOut_ALL[[10]], ee_NDVImax_list[[5]]) 
  wolf_eeNDVImax_wtr_OK <- join_NDVImax(crwOut_ALL[[12]], ee_NDVImax_list[[6]]) 
  wolf_eeNDVImax_wtr_NE <- join_NDVImax(crwOut_ALL[[14]], ee_NDVImax_list[[7]]) 
  bob_eeNDVImax_wtr_OK <- join_NDVImax(crwOut_ALL[[16]], ee_NDVImax_list[[8]]) 
  bob_eeNDVImax_wtr_NE <- join_NDVImax(crwOut_ALL[[18]], ee_NDVImax_list[[9]]) 
  coy_eeNDVImax_wtr_OK <- join_NDVImax(crwOut_ALL[[20]], ee_NDVImax_list[[10]]) 
  coy_eeNDVImax_wtr_NE <- join_NDVImax(crwOut_ALL[[22]], ee_NDVImax_list[[11]]) 
  
  #'  List all data together
  ee_NDVImax_list <- list(md_eeNDVImax_wtr, elk_eeNDVImax_wtr, wtd_eeNDVImax_wtr, 
                       coug_eeNDVImax_wtr_OK, coug_eeNDVImax_wtr_NE, wolf_eeNDVImax_wtr_OK, 
                       wolf_eeNDVImax_wtr_NE, bob_eeNDVImax_wtr_OK, bob_eeNDVImax_wtr_NE, 
                       coy_eeNDVImax_wtr_OK, coy_eeNDVImax_wtr_NE)
  
  #'  Save extracted EE data
  save(ee_NDVIsmr_list, file = paste0("./Outputs/Telemetry_covs/ee_NDVIsmr_list_", Sys.Date(), ".RData"))
  save(ee_NDVImax_list, file = paste0("./Outputs/Telemetry_covs/ee_NDVImax_list_", Sys.Date(), ".RData"))

  
  #'  Next: add these to other data-extraction script and combine for HMM analyses
  
  
  #'  Citation for MODIS/Terra Vegetation Indices 16-Day L3 Global 250 m SIN Grid:
  #'  See LP DAAC citation policies for citing NASA products (https://lpdaac.usgs.gov/data/data-citation-and-policies/)
  #'  DOI: 10.5067/MODIS/MOD13Q1.006
  #'  Citation for Snow Cover Data
  #'  Hall, D. K., V. V. Salomonson, and G. A. Riggs. 2016. MODIS/Terra Snow Cover Daily L3 Global 500m Grid. Version 6. Boulder, Colorado USA: NASA National Snow and Ice Data Center Distributed Active Archive Center.
  #'  Citation for Snow Depth Data
  #'  Rodell, M., P.R. Houser, U. Jambor, J. Gottschalck, K. Mitchell, C.-J. Meng, K. Arsenault, B. Cosgrove, J. Radakovich, M. Bosilovich, J.K. Entin, J.P. Walker, D. Lohmann, and D. Toll, The Global Land Data Assimilation System, Bull. Amer. Meteor. Soc., 85(3), 381-394, 2004.
  #'  Daymet Snow-Water Equivalent
  #'  Thornton, M.M., R. Shrestha, Y. Wei, P.E. Thornton, S. Kao, and B.E. Wilson. {2021}. Daymet: Daily Surface Weather Data on a 1-km Grid for North America, Version 4. ORNL DAAC, Oak Ridge, Tennessee, USA
    
    
  



