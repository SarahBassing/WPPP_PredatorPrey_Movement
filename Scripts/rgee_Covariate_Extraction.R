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
  load("./Outputs/Telemetry_crwOut/crwOut_ALL_2021-12-08.RData")

  
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
  ee_NDVI_coug_smr <- match_ee_data(data_ee_list[[7]], imagecoll = imageNDVI, tempwin = 16, band = "NDVI", sp.res = 250, tmp.res = 16)
  ee_NDVI_wolf_smr <- match_ee_data(data_ee_list[[9]], imagecoll = imageNDVI, tempwin = 16, band = "NDVI", sp.res = 250, tmp.res = 16)
  ee_NDVI_bob_smr <- match_ee_data(data_ee_list[[11]], imagecoll = imageNDVI, tempwin = 16, band = "NDVI", sp.res = 250, tmp.res = 16)
  ee_NDVI_coy_smr <- match_ee_data(data_ee_list[[13]], imagecoll = imageNDVI, tempwin = 16, band = "NDVI", sp.res = 250, tmp.res = 16)
  
  #'  List together for safe keeping
  ee_smr_NDVI_list <- list(ee_NDVI_md_smr, ee_NDVI_elk_smr, ee_NDVI_wtd_smr, 
                           ee_NDVI_coug_smr, ee_NDVI_wolf_smr, ee_NDVI_bob_smr, 
                           ee_NDVI_coy_smr)
  
  save(ee_smr_NDVI_list, file = paste0("./Outputs/Telemetry_covs/ee_smr_NDVI_list_", Sys.Date(), ".RData"))
  
  
  ####  Re-scale NDVI values  ####
  #'  ONLY RUN IF EXTRACTING NDVI DATA!
  #'  MODIS data are scaled by a factor of 0.0001 for ease of extraction (see 
  #'  pg. 9 of MODIS User's Guide) but NDVI should range -1 to 1. All extracted 
  #'  NDVI values need to be rescaled by 0.0001 to be in the correct range.
  ndvi_rescale <- function(ee_data, plotit = F) {
    #'  Create new column with re-scaled NDVI values
    dataoutput <- mutate(ee_data, NDVI_scale = NDVI*0.0001)
    #'  Visualize
    plot(dataoutput$NDVI_scale)
    hist(dataoutput$NDVI_scale)
    
    return(dataoutput)
  }
  ee_NDVI <- lapply(ee_smr_NDVI_list, ndvi_rescale, T)
  
  #' ####  Clean Snow Cover values  ####
  #' #'  ONLY RUN IF EXTRACTING SNOW COVER DATA!
  #' #'  Any snow cover values >100 indicate a problem in the original ee data so 
  #' #'  need to change these values to NA
  #' snow_cover_na <- function(ee_data) {
  #'   dataoutput <- mutate(ee_data, Snow_Cover = ifelse(NDSI_Snow_Cover >100, NA, NDSI_Snow_Cover))
  #'   return(dataoutput)
  #' }
  #' ee_SNOWCOVER <- lapply(ee_SNOWCOVER, snow_cover_na)
  
  
  
  ####  Maximum NDVI  ####
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
    # full_crwOut <- crwOut_data
    crwOut_sf <- st_as_sf(crwOut_data, coords = c('mu.x','mu.y'),  
                          crs = "+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs ")
    #'  Transform projection to WGS84 for Google Earth Engine
    datasf <- st_transform(crwOut_sf, crs = 4326)
    #'  Make a unique ID for every observation (necessary for ee_extract!)
    datasf$uid <- as.integer(1:nrow(datasf))
    datasf$NAME <- paste0(datasf$AnimalID, "_", datasf$uid)

    #' #'  Chunk location data into groups to loop through with getInfo
    #' datasf$uniq <- rep(1:1000, each = 1000)[1:nrow(datasf)]

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
    ee_allNDVI$maxNDVI <- apply(ee_allNDVI[-c(1:2)], 1, max)
    ee_allNDVI <- ee_allNDVI %>%
      mutate(
        maxNDVI_scale = maxNDVI*0.0001,
        ID = seq(1:nrow(ee_allNDVI))
      ) %>%
      dplyr::select(c(ID, NAME, maxNDVI_scale, geometry))
    
    #'  Pull out location dates (including interpolated ones)
    times <- as.Date(crwOut_sf$time, "%Y-%m-%d", tz = "Etc/GMT+8")
    date.only <- as.data.frame(times)
    
    #'  Save only relevant location info
    crwOut_sf$ID2 <- as.integer(1:nrow(crwOut_sf))
    crwOut_sf <- cbind(crwOut_sf, date.only)
    crwOut_df <- dplyr::select(crwOut_sf, c(AnimalID, Season, StudyArea, times, x, y, ID2, geometry))
    colnames(crwOut_df) <- c("AnimalID", "Season", "StudyArea", "Dates", "x", "y", "ID", "geometry")
    
    #'  Join animal location data with maxNDVI values
    maxNDVI_df <- as.data.frame(crwOut_df) %>%
      full_join(ee_allNDVI, by = "ID") %>%
      dplyr::select(-geometry.y)
    colnames(maxNDVI_df) <- c("AnimalID", "Season", "StudyArea", "Dates", "x", "y", "ID", "geometry", "NAME", "NDVImax")
    
    print(maxNDVI_df[1:6,])
    
    return(maxNDVI_df)
    
  }
  
  
  #' #'  Create empty dataframe to hold extracted pixel values
  #' EE_allNDVI <- data.frame()
  #' 
  #' #'  Loop through location data in chunks to extract pixel values from EE images
  #' for(x in unique(datasf$uniq)) {
  #'   #'  Chunk location data
  #'   data1 <- datasf %>% filter(uniq == x)
  #'   #'  Use getInfo method to extract values from GEE
  #'   dataoutput <- ee_extract(
  #'     x = eendvi,
  #'     y = data1["geometry"],
  #'     scale = ee.scale,
  #'     fun = ee$Reducer$max(),
  #'     via = "getInfo",
  #'     sf = FALSE
  #'   )
  #'   #'  Append each loop
  #'   EE_allNDVI <- rbind(EE_allNDVI, dataoutput)
  #' }
  #' #'  Find maximum NDVI values during annual growing season (April - Sept)
  #' #'  2018 growing season
  #' ee_allNDVI18 <- as.data.frame(EE_allNDVI) %>%
  #'   #'  Rename all columns by removing first 12 characters for easier reading
  #'   rename_with(~ gsub("^............", "", .x)) %>%
  #'   #'  Save only the NDVI values from the growing season (April - Sept)
  #'   dplyr::select(contains(c("2018_04", "2018_05", "2018_06", "2018_07", "2018_08", "2018_09")))
  #' #'  Find max value row-wise and re-scale
  #' ee_allNDVI18$maxNDVI18 <- apply(ee_allNDVI18, 1, max)
  #' ee_allNDVI18 <- ee_allNDVI18 %>%
  #'   mutate(
  #'     maxNDVI18_scale = maxNDVI18*0.0001,
  #'     ID = seq(1:nrow(ee_allNDVI18))
  #'   )
  #' 
  #' #'  2019 growing season
  #' ee_allNDVI19 <- as.data.frame(EE_allNDVI) %>%
  #'   #'  Rename all columns by removing first 12 characters
  #'   rename_with(~ gsub("^............", "", .x)) %>%
  #'   #'  Save only the NDVI values from the growing season (April - Sept)
  #'   dplyr::select(contains(c("2019_04", "2019_05", "2019_06", "2019_07", "2019_09", "2019_09"))) %>%
  #'   mutate(ID = seq(1:nrow(EE_allNDVI)))
  #' #'  Find max value row-wise and re-scale
  #' ee_allNDVI19$maxNDVI19 <- apply(ee_allNDVI19, 1, max)
  #' ee_allNDVI19 <- ee_allNDVI19 %>%
  #'   mutate(
  #'     maxNDVI19_scale = maxNDVI19*0.0001,
  #'     ID = seq(1:nrow(ee_allNDVI19))
  #'   )
  #' 
  #' #'  2020 growing season
  #' ee_allNDVI20 <- as.data.frame(EE_allNDVI) %>%
  #'   #'  Rename all columns by removing first 12 characters
  #'   rename_with(~ gsub("^............", "", .x)) %>%
  #'   #'  Save only the NDVI values from the growing season (April - Sept)
  #'   dplyr::select(contains(c("2020_04", "2020_05", "2020_06", "2020_07", "2020_08", "2020_09"))) %>%
  #'   mutate(ID = seq(1:nrow(EE_allNDVI)))
  #' #'  Find max value row-wise and re-scale
  #' ee_allNDVI20$maxNDVI20 <- apply(ee_allNDVI20, 1, max)
  #' ee_allNDVI20 <- ee_allNDVI20 %>%
  #'   mutate(
  #'     maxNDVI20_scale = maxNDVI20*0.0001,
  #'     ID = seq(1:nrow(ee_allNDVI20))
  #'   )
  #' 
  #' #'  Join max NDVI values across years
  #' maxNDVI <- ee_allNDVI18 %>%
  #'   full_join(ee_allNDVI19, by = "ID") %>%
  #'   full_join(ee_allNDVI20, by = "ID") %>%
  #'   dplyr::select(c(ID, maxNDVI18_scale, maxNDVI19_scale, maxNDVI20_scale))
  #' colnames(maxNDVI) <- c("ID", "max_NDVI18", "max_NDVI19", "max_NDVI20")
  #' 
  #' #'  Pull out location dates (including interpolated ones)
  #' times <- as.Date(crwOut_sf$time, "%Y-%m-%d", tz = "Etc/GMT+8")
  #' date.only <- as.data.frame(times)
  #' 
  #' #'  Save only relevant location info
  #' crwOut_sf$ID2 <- as.integer(1:nrow(crwOut_sf))
  #' crwOut_sf <- cbind(crwOut_sf, date.only)
  #' crwOut_df <- dplyr::select(crwOut_sf, c(AnimalID, Season, StudyArea, times, x, y, ID2, geometry))
  #' #crwOut_df <- as.data.frame(cbind(crwOut_sf$ID, crwOut_sf$Season, crwOut_sf$StudyArea, crwOut_sf$mu.x, crwOut_sf$mu.y, date.only, crwOut_sf$ID2))
  #' colnames(crwOut_df) <- c("AnimalID", "Season", "StudyArea", "Dates", "mu.x", "mu.y", "ID", "geometry") #"mu.x", "mu.y",
  #' 
  #' #'  Join animal location data with maxNDVI values
  #' maxNDVI_df <- as.data.frame(crwOut_df) %>%
  #'   full_join(maxNDVI, by = "ID") %>%
  #'   #'  Retain only the previous growing season's max NDVI value for each location
  #'   #'  Make all other NDVI values "NA"
  #'   mutate(
  #'     newNDVI18 = as.numeric(ifelse(Dates > "2019-03-01", NA, max_NDVI18)),
  #'     newNDVI19 = as.numeric(ifelse(Dates < "2019-12-01" | Dates > "2020-03-01", NA, max_NDVI19)),
  #'     newNDVI20 = as.numeric(ifelse(Dates < "2020-12-01", NA, max_NDVI20))
  #'   ) %>%
  #'   #'  Create a single column of max NDVI values for each location
  #'   rowwise() %>%
  #'   mutate(maxNDVI = max(c(newNDVI18, newNDVI19, newNDVI20), na.rm = TRUE)) %>%
  #'   #'  Drop extra NDVI columns from df
  #'   dplyr::select(-c(max_NDVI18, max_NDVI19, max_NDVI20, newNDVI18, newNDVI19, newNDVI20))
  #' 
  #' return(maxNDVI_df)
  
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
  filter_wtr1819 <- function(mylist){
    df <- as.data.frame(mylist) %>%
      filter(time < "2019-03-01 00:00:00") %>%
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
  wtr_list1819 <- list(wtr1819[[2]], wtr1819[[4]], wtr1819[[6]], wtr1819[[8]], wtr1819[[10]], wtr1819[[12]], wtr1819[[14]])
  wtr_list1920 <- list(wtr1920[[2]], wtr1920[[4]], wtr1920[[6]], wtr1920[[8]], wtr1920[[10]], wtr1920[[12]], wtr1920[[14]])
  wtr_list2021 <- list(wtr2021[[2]], wtr2021[[4]], wtr2021[[6]], wtr2021[[8]], wtr2021[[10]], wtr2021[[12]], wtr2021[[14]])
  #'  New list order: mule deer [1], elk [2], white-tailed deer [3], cougar [4], wolf [5], bobcat [6], coyote [7]
  
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
  #'  List order: mule deer [1], elk [2], white-tailed deer [3], cougar [4], 
  #'              wolf [5], bobcat [6], coyote [7]
  #'              
  #'  Split md_wtr1819 data even more because it's apparently too large?!?!
  nrow(wtr_list1819[[1]])
  md_wtr1819a <- wtr_list1819[[1]][1:20000,]
  md_wtr1819b <- wtr_list1819[[1]][20001:39327,]
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
  
  coug_wtr1819_NDVImax <- find_maxNDVI(wtr_list1819[[4]], eeimage = eeimage, start.date = start18, end.date = end18, band.name = band.name, band = band, ee.scale = ee.scale)
  coug_wtr1920_NDVImax <- find_maxNDVI(wtr_list1920[[4]], eeimage = eeimage, start.date = start19, end.date = end19, band.name = band.name, band = band, ee.scale = ee.scale)
  coug_wtr2021_NDVImax <- find_maxNDVI(wtr_list2021[[4]], eeimage = eeimage, start.date = start20, end.date = end20, band.name = band.name, band = band, ee.scale = ee.scale)
  
  wolf_wtr1819_NDVImax <- find_maxNDVI(wtr_list1819[[5]], eeimage = eeimage, start.date = start18, end.date = end18, band.name = band.name, band = band, ee.scale = ee.scale)
  wolf_wtr1920_NDVImax <- find_maxNDVI(wtr_list1920[[5]], eeimage = eeimage, start.date = start19, end.date = end19, band.name = band.name, band = band, ee.scale = ee.scale)
  wolf_wtr2021_NDVImax <- find_maxNDVI(wtr_list2021[[5]], eeimage = eeimage, start.date = start20, end.date = end20, band.name = band.name, band = band, ee.scale = ee.scale)
  
  bob_wtr1819_NDVImax <- find_maxNDVI(wtr_list1819[[6]], eeimage = eeimage, start.date = start18, end.date = end18, band.name = band.name, band = band, ee.scale = ee.scale)
  bob_wtr1920_NDVImax <- find_maxNDVI(wtr_list1920[[6]], eeimage = eeimage, start.date = start19, end.date = end19, band.name = band.name, band = band, ee.scale = ee.scale)
  bob_wtr2021_NDVImax <- find_maxNDVI(wtr_list2021[[6]], eeimage = eeimage, start.date = start20, end.date = end20, band.name = band.name, band = band, ee.scale = ee.scale)
  
  coy_wtr1819_NDVImax <- find_maxNDVI(wtr_list1819[[7]], eeimage = eeimage, start.date = start18, end.date = end18, band.name = band.name, band = band, ee.scale = ee.scale)
  coy_wtr1920_NDVImax <- find_maxNDVI(wtr_list1920[[7]], eeimage = eeimage, start.date = start19, end.date = end19, band.name = band.name, band = band, ee.scale = ee.scale)
  coy_wtr2021_NDVImax <- find_maxNDVI(wtr_list2021[[7]], eeimage = eeimage, start.date = start20, end.date = end20, band.name = band.name, band = band, ee.scale = ee.scale)
  
  # md_NDVImax_list <- list(md_wtr1819_NDVImax, md_wtr1920_NDVImax, md_wtr2021_NDVImax)
  # elk_NDVImax_list <- list(elk_wtr1819_NDVImax, elk_wtr1920_NDVImax, elk_wtr2021_NDVImax)
  # wtd_NDVImax_list <- list(wtd_wtr1819_NDVImax, wtd_wtr1920_NDVImax, wtd_wtr2021_NDVImax)
  # coug_NDVImax_list <- list(coug_wtr1819_NDVImax, coug_wtr1920_NDVImax, coug_wtr2021_NDVImax)
  # wolf_NDVImax_list <- list(wolf_wtr1819_NDVImax, wolf_wtr1920_NDVImax, wolf_wtr2021_NDVImax)
  # bob_NDVImax_list <- list(bob_wtr1819_NDVImax, bob_wtr1920_NDVImax, bob_wtr2021_NDVImax)
  # coy_NDVImax_list <- list(coy_wtr1819_NDVImax, coy_wtr1920_NDVImax, coy_wtr2021_NDVImax)
  
  ee_NDVImax_list <- list(md_NDVImax_list, elk_NDVImax_list, wtd_NDVImax_list, coug_NDVImax_list, wolf_NDVImax_list, bob_NDVImax_list, coy_NDVImax_list)
  
  md_NDVImax <- rbind(md_wtr1819_NDVImax, md_wtr1920_NDVImax, md_wtr2021_NDVImax)
  elk_NDVImax <- rbind(elk_wtr1819_NDVImax, elk_wtr1920_NDVImax, elk_wtr2021_NDVImax)
  wtd_NDVImax <- rbind(wtd_wtr1819_NDVImax, wtd_wtr1920_NDVImax, wtd_wtr2021_NDVImax)
  coug_NDVImax <- rbind(coug_wtr1819_NDVImax, coug_wtr1920_NDVImax, coug_wtr2021_NDVImax)
  wolf_NDVImax <- rbind(wolf_wtr1819_NDVImax, wolf_wtr1920_NDVImax, wolf_wtr2021_NDVImax)
  bob_NDVImax <- rbind(bob_wtr1819_NDVImax, bob_wtr1920_NDVImax, bob_wtr2021_NDVImax)
  coy_NDVImax <- rbind(coy_wtr1819_NDVImax, coy_wtr1920_NDVImax, coy_wtr2021_NDVImax)
  
  ee_NDVImax_list <- list(md_NDVImax, elk_NDVImax, wtd_NDVImax, coug_NDVImax, wolf_NDVImax, bob_NDVImax, coy_NDVImax)
  
  # ee_NDVImax_wtr1819 <- lapply(wtr_list1819, find_maxNDVI, eeimage = eeimage, start.date = "2018-04-01", end.date = "2018-10-01", band.name = band.name, band = band, ee.scale = ee.scale)
  # ee_NDVImax_wtr1920 <- lapply(wtr_list1920, find_maxNDVI, eeimage = eeimage, start.date = "2019-04-01", end.date = "2019-10-01", band.name = band.name, band = band, ee.scale = ee.scale)
  # ee_NDVImax_wtr2021 <- lapply(wtr_list2021, find_maxNDVI, eeimage = eeimage, start.date = "2020-04-01", end.date = "2020-10-01", band.name = band.name, band = band, ee.scale = ee.scale)
  # 
  # ee_NDVImax_list <- list(ee_NDVImax_wtr1819, ee_NDVImax_wtr1920, ee_NDVImax_wtr2021)
  
  save(ee_NDVImax_list, file = paste0("./Outputs/Telemetry_covs/ee_NDVImax_list_", Sys.Date(), ".RData"))
  
  
  #'  find max NDVI for each species and season in the short list
  # ee_NDVImax <- lapply(short.list, find_maxNDVI, eeimage = eeimage, 
  #                      start.date = start.date, end.date = end.date, 
  #                      band.name = band.name, band = band, ee.scale = ee.scale)
  # ee_NDVImax_md_wtr <- find_maxNDVI(short.list[[2]], eeimage = eeimage, start.date = start.date, end.date = end.date, band.name = band.name, band = band, ee.scale = ee.scale)
  # ee_NDVImax_elk_wtr <- find_maxNDVI(short.list[[4]], eeimage = eeimage, start.date = start.date, end.date = end.date, band.name = band.name, band = band, ee.scale = ee.scale)
  # ee_NDVImax_wtd_wtr <- find_maxNDVI(short.list[[6]], eeimage = eeimage, start.date = start.date, end.date = end.date, band.name = band.name, band = band, ee.scale = ee.scale)
  # ee_NDVImax_coug_wtr <- find_maxNDVI(short.list[[8]], eeimage = eeimage, start.date = start.date, end.date = end.date, band.name = band.name, band = band, ee.scale = ee.scale)
  # ee_NDVImax_wolf_wtr <- find_maxNDVI(short.list[[10]], eeimage = eeimage, start.date = start.date, end.date = end.date, band.name = band.name, band = band, ee.scale = ee.scale)
  # ee_NDVImax_bob_wtr <- find_maxNDVI(short.list[[12]], eeimage = eeimage, start.date = start.date, end.date = end.date, band.name = band.name, band = band, ee.scale = ee.scale)
  # ee_NDVImax_coy_wtr <- find_maxNDVI(short.list[[14]], eeimage = eeimage, start.date = start.date, end.date = end.date, band.name = band.name, band = band, ee.scale = ee.scale)
  # 
  # ee_NDVImax_list <- list(ee_NDVImax_md_wtr, ee_NDVImax_elk_wtr, ee_NDVImax_wtd_wtr,
  #                         ee_NDVImax_coug_wtr, ee_NDVImax_wolf_wtr, ee_NDVImax_bob_wtr,
  #                         ee_NDVImax_coy_wtr)
  
  # load("G:/My Drive/1_Repositories/WPPP_PredatorPrey_Movement/Outputs/Telemetry_covs/ee_covs_list_2021-12-28.RData")

  ####  Join datasets  ####
  join_data <- function(crwOut_data, ndvi, maxndvi, snow_cov, snow_dep) { 
    #'  Make sure each observation has the unique Animal ID
    full_crwOut <- crwOut_data[[2]]
    crwOut <- full_crwOut %>%
      separate(ID, sep = "_", into = "AnimalID",  extra = "drop")
    ee_covs <- as.data.frame(cbind(crwOut$AnimalID, crwOut$Season, crwOut$StudyArea, crwOut$mu.x, crwOut$mu.y))
    #'  Add a unique ID for each observation for easier joining
    # ee_covs$ID <- as.integer(1:nrow(ee_covs))
    ee_covs <- ee_covs %>%
      full_join(ndvi, by = c("AnimalID", "Season")) %>%
        dplyr::select(-c(DateTimeImage, date_millis, uniq, geometry, NDVI)) %>%
      full_join(maxndvi, by = c("ID", "AnimalID", "Season", "StudyArea", "mu.x", "mu.y")) %>%
        dplyr::select(-c(AnimalID2, Dates)) %>%
      # full_join(snow_cov, by = "ID") %>%
      #   dplyr::select(-c(DateTimeImage, date_millis, uniq, geometry, NDSI_Snow_Cover)) %>%
      # full_join(snow_dep, by = "ID") %>%
      #   dplyr::select(-c(Date.y, Date, DateTimeImage, date_millis, uniq, geometry)) %>%
     relocate(ID, .after = (SnowDepth_inst))
    colnames(ee_covs) <- c("AnimalID", "Season", "StudyArea", "mu.x", "mu.y", "Date", "NDVI", "maxNDVI", "ID") #, "Snow.Cover", "Snow.Depth",

    return(ee_covs)
  }
  #'  Combine EE values with species- and season-specific data
  md_ee_smr <- join_data(crwOut_ALL[[1]], ee_NDVI[[1]]) #, ee_NDVImax[[1]], ee_SNOWCOVER[[1]], ee_SNOWDEPTH[[1]])
  md_ee_wtr <- join_data(crwOut_ALL[[2]], ee_NDVImax[[1]]) #, ee_NDVImax[[2]], ee_SNOWCOVER[[2]], ee_SNOWDEPTH[[2]])
  elk_ee_smr <- join_data(crwOut_ALL[[3]], ee_NDVI[[3]]) #, ee_NDVImax[[3]], ee_SNOWCOVER[[3]], ee_SNOWDEPTH[[3]])
  elk_ee_wtr <- join_data(crwOut_ALL[[4]], ee_NDVImax[[2]]) #, ee_NDVImax[[4]], ee_SNOWCOVER[[4]], ee_SNOWDEPTH[[4]])
  wtd_ee_smr <- join_data(crwOut_ALL[[5]], ee_NDVI[[5]]) #, ee_NDVImax[[5]], ee_SNOWCOVER[[5]], ee_SNOWDEPTH[[5]])
  wtd_ee_wtr <- join_data(crwOut_ALL[[6]], ee_NDVImax[[3]]) #, ee_NDVImax[[6]], ee_SNOWCOVER[[6]], ee_SNOWDEPTH[[6]])
  coug_ee_smr <- join_data(crwOut_ALL[[7]], ee_NDVI[[7]]) #, ee_NDVImax[[7]], ee_SNOWCOVER[[7]], ee_SNOWDEPTH[[7]])
  coug_ee_wtr <- join_data(crwOut_ALL[[8]], ee_NDVImax[[4]]) #, ee_NDVImax[[8]], ee_SNOWCOVER[[8]], ee_SNOWDEPTH[[8]])
  wolf_ee_smr <- join_data(crwOut_ALL[[9]], ee_NDVI[[8]]) #, ee_NDVImax[[9]], ee_SNOWCOVER[[9]], ee_SNOWDEPTH[[9]])
  wolf_ee_wtr <- join_data(crwOut_ALL[[10]], ee_NDVImax[[5]]) #, ee_NDVImax[[10]], ee_SNOWCOVER[[10]], ee_SNOWDEPTH[[10]])
  bob_ee_smr <- join_data(crwOut_ALL[[11]], ee_NDVI[[11]]) #, ee_NDVImax[[11]], ee_SNOWCOVER[[11]], ee_SNOWDEPTH[[11]])
  bob_ee_wtr <- join_data(crwOut_ALL[[12]], ee_NDVImax[[6]]) #, ee_NDVImax[[12]], ee_SNOWCOVER[[12]], ee_SNOWDEPTH[[12]])
  coy_ee_smr <- join_data(crwOut_ALL[[13]], ee_NDVI[[13]]) #, ee_NDVImax[[13]], ee_SNOWCOVER[[13]], ee_SNOWDEPTH[[13]])
  coy_ee_wtr <- join_data(crwOut_ALL[[14]], ee_NDVImax[[7]]) #, ee_NDVImax[[14]], ee_SNOWCOVER[[14]], ee_SNOWDEPTH[[14]])

  #'  List all data together
  ee_covs_list <- list(md_ee_smr, md_ee_wtr, elk_ee_smr, elk_ee_wtr, wtd_ee_smr, 
                       wtd_ee_wtr, coug_ee_smr, coug_ee_wtr, wolf_ee_smr, 
                       wolf_ee_wtr, bob_ee_smr, bob_ee_wtr, coy_ee_smr, coy_ee_wtr)
  
  #'  Save extracted EE data
  save(ee_covs_list, file = paste0("./Outputs/Telemetry_covs/ee_covs_list_", Sys.Date(), ".RData"))

  
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
    
    
  



