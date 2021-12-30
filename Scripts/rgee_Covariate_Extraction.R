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
    datasf <- dplyr::select(datasf, c("Date", "ID", "geometry"))
    
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
  imageSNOWCOVER <- ee$ImageCollection("MODIS/006/MOD10A1")$filterDate(start,end) # MODIS Terra Snow Cover Daily Global (Normalized Difference Snow Index (NDSI)) daily, 500m resolution, values represent % coverage per pixel?
  imageSNOWDEPTH <- ee$ImageCollection("NASA/GLDAS/V021/NOAH/G025/T3H")$filterDate(start,end) # NASA Global Land Data Assimilation System (GLDAS) daily, 27830m (or 0.25 degree) resolution, depth measured in meters
  imageNDSI <- ee$ImageCollection("MODIS/MOD09GA_006_NDSI")$filterDate(start,end) # MODIS Terra Daily Normalized Difference Snow Index (NDSI) daily, 463.313 m resolution, missing data from Jan. 2019
  imageSWE <- ee$ImageCollection("NASA/ORNL/DAYMET_V4")$filterDate(start,end) # Daymet V4: Daily Surface Weather and Climatological Summaries, Snow-Water Equivalent, daily, 1000m resolution, data up to Dec. 31, 2020
  
  
  
  #'  Run reformated animal location data through this monster function to
  #'  match & extract EE images 
  # tst <- list(data_ee_list[[13]], data_ee_list[[14]])
  ee_NDVI <- lapply(data_ee_list, match_ee_data, imagecoll = imageNDVI, tempwin = 16, band = "NDVI", sp.res = 250, tmp.res = 16)
  ee_SNOWCOVER <- lapply(data_ee_list, match_ee_data, imagecoll = imageSNOWCOVER, tempwin = 1, band = "NDSI_Snow_Cover", sp.res = 500, tmp.res = 1)
  ee_SNOWDEPTH <- lapply(data_ee_list, match_ee_data, imagecoll = imageSNOWDEPTH, tempwin = 1, band = "SnowDepth_inst", sp.res = 27830, tmp.res = 1)
  ee_NDSI <- lapply(data_ee_list, match_ee_data, imagecoll = imageNDSI, tempwin = 1, band = "NDSI", sp.res = 463.313, tmp.res = 1)
  ee_SWE <- lapply(data_ee_list, match_ee_data, imagecoll = imageSWE, tempwin = 1, band = "swe", sp.res = 1000, tmp.res = 1)
  
  
  
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
  ee_NDVI <- lapply(ee_NDVI, ndvi_rescale, T)
  
  ####  Clean Snow Cover values  ####
  #'  ONLY RUN IF EXTRACTING SNOW COVER DATA!
  #'  Any snow cover values >100 indicate a problem in the original ee data so 
  #'  need to change these values to NA
  snow_cover_na <- function(ee_data) {
    dataoutput <- mutate(ee_data, Snow_Cover = ifelse(NDSI_Snow_Cover >100, NA, NDSI_Snow_Cover))
    return(dataoutput)
  }
  ee_SNOWCOVER <- lapply(ee_SNOWCOVER, snow_cover_na)
  
  
  
  ####  Maximum NDVI  ####
  #'  Function to extract MAXIMUM annual NDVI, representing potential forage  
  #'  quality available post-growing season at winter locations.
  find_maxNDVI <- function(crwOut_data, eeimage, start.date, end.date, band.name, band, ee.scale) {
    
    #'  Define EE ImageCollection
    eendvi <- ee$ImageCollection(eeimage) %>% # "MODIS/006/MOD13Q1" 
      #'  Define date range of interest
      ee$ImageCollection$filterDate(start.date, end.date) %>% # "2018-04-01", "2020-10-01" 
      ee$ImageCollection$map(
        function(x) {
          date <- ee$Date(x$get("system:time_start"))$format('YYYY_MM_dd')
          name <- ee$String$cat(band.name, date) # "NDVI_" 
          x$select(band)$rename(name) #"NDVI"  
        }
      )
    
    #'  Make telemetry locations spatial using sf (note: use coordinates labeled 
    #'  mu.x & mu.y so that interpolated locations are included)
    # full_crwOut <- crwOut_data
    crwOut_sf <- st_as_sf(crwOut_data, coords = c('mu.x','mu.y'),  
                          crs = "+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs ")
    #'  Transform projection to WGS84 for Google Earth Engine
    datasf <- st_transform(crwOut_sf, crs = 4326)
    
    #'  Chunk location data into groups to loop through 
    datasf$uniq <- rep(1:1000, each = 1000)[1:nrow(datasf)] 
    
    #'  Track amount of time it takes to extract data
    start_time <- Sys.time()
    
    #'  Create empty dataframe to hold extracted pixel values
    EE_allNDVI <- data.frame()
    
    #'  Loop through location data in chunks to extract pixel values from EE images 
    for(i in unique(datasf$uniq)) {
      #'  Use getInfo method to extract values from GEE
      dataoutput <- ee_extract(
        x = eendvi,
        y = datasf[i]["geometry"],
        scale = ee.scale,
        fun = ee$Reducer$max(),
        via = "getInfo",
        sf = FALSE
      )
      #'  Append each loop
      EE_allNDVI <- rbind(EE_allNDVI, dataoutput) 
    }
    end_time <- Sys.time()
    #'  How long did that take?
    print(end_time - start_time)
    #'  Quick peak
    print(EE_allNDVI[1:6,1:4])
    
    #'  Find maximum NDVI values during annual growing season (April - Sept)
    #'  2018 growing season
    ee_allNDVI18 <- as.data.frame(EE_allNDVI) %>%
      #'  Rename all columns by removing first 12 characters for easier reading
      rename_with(~ gsub("^............", "", .x)) %>%
      #'  Save only the NDVI values from the growing season (April - Sept)
      dplyr::select(contains(c("2018_04", "2018_05", "2018_06", "2018_07", "2018_08", "2018_09")))
    #'  Find max value row-wise and re-scale
    ee_allNDVI18$maxNDVI18 <- apply(ee_allNDVI18, 1, max)
    ee_allNDVI18 <- ee_allNDVI18 %>%
      mutate(
        maxNDVI18_scale = maxNDVI18*0.0001,
        ID = seq(1:nrow(ee_allNDVI18))
      )
    
    #'  2019 growing season
    ee_allNDVI19 <- as.data.frame(EE_allNDVI) %>%
      #'  Rename all columns by removing first 12 characters
      rename_with(~ gsub("^............", "", .x)) %>%
      #'  Save only the NDVI values from the growing season (April - Sept)
      dplyr::select(contains(c("2019_04", "2019_05", "2019_06", "2019_07", "2019_09", "2019_09"))) %>%
      mutate(ID = seq(1:nrow(EE_allNDVI)))
    #'  Find max value row-wise and re-scale
    ee_allNDVI19$maxNDVI19 <- apply(ee_allNDVI19, 1, max)
    ee_allNDVI19 <- ee_allNDVI19 %>%
      mutate(
        maxNDVI19_scale = maxNDVI19*0.0001,
        ID = seq(1:nrow(ee_allNDVI19))
      )
    
    #'  2020 growing season
    ee_allNDVI20 <- as.data.frame(EE_allNDVI) %>%
      #'  Rename all columns by removing first 12 characters
      rename_with(~ gsub("^............", "", .x)) %>%
      #'  Save only the NDVI values from the growing season (April - Sept)
      dplyr::select(contains(c("2020_04", "2020_05", "2020_06", "2020_07", "2020_08", "2020_09"))) %>%
      mutate(ID = seq(1:nrow(EE_allNDVI)))
    #'  Find max value row-wise and re-scale
    ee_allNDVI20$maxNDVI20 <- apply(ee_allNDVI20, 1, max)
    ee_allNDVI20 <- ee_allNDVI20 %>%
      mutate(
        maxNDVI20_scale = maxNDVI20*0.0001,
        ID = seq(1:nrow(ee_allNDVI20))
      )
    
    #'  Join max NDVI values across years
    maxNDVI <- ee_allNDVI18 %>%
      full_join(ee_allNDVI19, by = "ID") %>%
      full_join(ee_allNDVI20, by = "ID") %>%
      dplyr::select(c(ID, maxNDVI18_scale, maxNDVI19_scale, maxNDVI20_scale))
    colnames(maxNDVI) <- c("ID", "max_NDVI18", "max_NDVI19", "max_NDVI20")
    
    #'  Pull out location dates (including interpolated ones)
    times <- as.Date(full_crwOut$time, "%Y-%m-%d", tz = "Etc/GMT+8")
    date.only <- as.data.frame(times)
    
    #'  Save only relevant location info
    full_crwOut$ID2 <- as.integer(1:nrow(full_crwOut))
    crwOut_df <- as.data.frame(cbind(full_crwOut$ID, full_crwOut$Season, full_crwOut$StudyArea, full_crwOut$mu.x, full_crwOut$mu.y, date.only, full_crwOut$ID2))
    colnames(crwOut_df) <- c("AnimalID2", "Season", "StudyArea", "mu.x", "mu.y", "Dates", "ID")
    
    #'  Join animal location data with maxNDVI values
    maxNDVI_df <- as.data.frame(crwOut_df) %>%
      full_join(maxNDVI, by = "ID") %>%
      #'  Retain only the previous growing season's max NDVI value for each location
      #'  Make all other NDVI values "NA"
      mutate(
        newNDVI18 = as.numeric(ifelse(Dates > "2019-03-01", NA, max_NDVI18)),
        newNDVI19 = as.numeric(ifelse(Dates < "2019-12-01" | Dates > "2020-03-01", NA, max_NDVI19)),
        newNDVI20 = as.numeric(ifelse(Dates < "2020-12-01", NA, max_NDVI20))
      ) %>%
      #'  Create a single column of max NDVI values for each location
      rowwise() %>%
      mutate(maxNDVI = max(c(newNDVI18, newNDVI19, newNDVI20), na.rm = TRUE)) %>%
      #'  Drop extra NDVI columns from df
      dplyr::select(-c(max_NDVI18, max_NDVI19, max_NDVI20, newNDVI18, newNDVI19, newNDVI20))
    
    return(maxNDVI_df)
  }
  #'  Define parameters of interest
  #'  Starting date is beginning of growing season of first year of study (April 2018)
  #'  Ending date is end of growing season of last year of study (Sept 2020) because
  #'  only interested in NDVI values from growing seasons, not winter 2021
  start.date <- "2018-04-01"
  end.date <- "2020-10-01"
  #'  EE imageCollection and relevant details
  eeimage <- "MODIS/006/MOD13Q1"
  band.name <- "NDVI_"
  band <- "NDVI"
  ee.scale <- 250
  
  #'  Function to drop crwFits lists for each species from crwOut_ALL list of lists 
  drop_list <- function(full.list) {
    short.list <- as.data.frame(full.list[-1]) %>%
      #'  Rename all columns by removing first 11 characters 
      #'  (annoyingly happens when covert each list to a data frame)
      rename_with(~ gsub("^...........", "", .x))
    return(short.list)
  }
  short.list <- lapply(crwOut_ALL, drop_list)
  
  #'  find max NDVI for each species and season in the short list
  ee_NDVImax <- lapply(short.list, find_maxNDVI, eeimage = eeimage, 
                       start.date = start.date, end.date = end.date, 
                       band.name = band.name, band = band, ee.scale = ee.scale)

  # load("G:/My Drive/1_Repositories/WPPP_PredatorPrey_Movement/Outputs/Telemetry_covs/ee_covs_list_2021-12-28.RData")

  ####  Join datasets  ####
  join_data <- function(crwOut_data, ndvi, maxndvi, snow_cov, snow_dep) { 
    #'  Make sure each observation has the unique Animal ID
    full_crwOut <- crwOut_data[[2]]
    crwOut <- full_crwOut %>%
      separate(ID, sep = "_", into = "AnimalID",  extra = "drop")
    ee_covs <- as.data.frame(cbind(crwOut$AnimalID, crwOut$Season, crwOut$StudyArea, crwOut$mu.x, crwOut$mu.y))
    #'  Add a unique ID for each observation for easier joining
    ee_covs$ID <- as.integer(1:nrow(ee_covs))
    ee_covs <- ee_covs %>%
      full_join(ndvi, by = "ID") %>%
        dplyr::select(-c(DateTimeImage, date_millis, uniq, geometry, NDVI)) %>%
      full_join(maxndvi, by = c("ID", "Season", "StudyArea", "mu.x", "mu.y")) %>%
        dplyr::select(-c(AnimalID2, Dates)) %>%
      full_join(snow_cov, by = "ID") %>%
        dplyr::select(-c(DateTimeImage, date_millis, uniq, geometry, NDSI_Snow_Cover)) %>%
      full_join(snow_dep, by = "ID") %>%
        dplyr::select(-c(Date.y, Date, DateTimeImage, date_millis, uniq, geometry)) %>%
     relocate(ID, .after = (SnowDepth_inst))
    colnames(ee_covs) <- c("AnimalID", "Season", "StudyArea", "mu.x", "mu.y", "Date", "NDVI", "maxNDVI", "Snow.Cover", "Snow.Depth", "ID")

    return(ee_covs)
  }
  #'  Combine EE values with species- and season-specific data
  md_ee_smr <- join_data(crwOut_ALL[[1]], ee_NDVI[[1]], ee_NDVImax[[1]], ee_SNOWCOVER[[1]], ee_SNOWDEPTH[[1]])
  md_ee_wtr <- join_data(crwOut_ALL[[2]], ee_NDVI[[2]], ee_NDVImax[[2]], ee_SNOWCOVER[[2]], ee_SNOWDEPTH[[2]])
  elk_ee_smr <- join_data(crwOut_ALL[[3]], ee_NDVI[[3]], ee_NDVImax[[3]], ee_SNOWCOVER[[3]], ee_SNOWDEPTH[[3]])
  elk_ee_wtr <- join_data(crwOut_ALL[[4]], ee_NDVI[[4]], ee_NDVImax[[4]], ee_SNOWCOVER[[4]], ee_SNOWDEPTH[[4]])
  wtd_ee_smr <- join_data(crwOut_ALL[[5]], ee_NDVI[[5]], ee_NDVImax[[5]], ee_SNOWCOVER[[5]], ee_SNOWDEPTH[[5]])
  wtd_ee_wtr <- join_data(crwOut_ALL[[6]], ee_NDVI[[6]], ee_NDVImax[[6]], ee_SNOWCOVER[[6]], ee_SNOWDEPTH[[6]])
  coug_ee_smr <- join_data(crwOut_ALL[[7]], ee_NDVI[[7]], ee_NDVImax[[7]], ee_SNOWCOVER[[7]], ee_SNOWDEPTH[[7]])
  coug_ee_wtr <- join_data(crwOut_ALL[[8]], ee_NDVI[[8]], ee_NDVImax[[8]], ee_SNOWCOVER[[8]], ee_SNOWDEPTH[[8]])
  wolf_ee_smr <- join_data(crwOut_ALL[[9]], ee_NDVI[[8]], ee_NDVImax[[9]], ee_SNOWCOVER[[9]], ee_SNOWDEPTH[[9]])
  wolf_ee_wtr <- join_data(crwOut_ALL[[10]], ee_NDVI[[10]], ee_NDVImax[[10]], ee_SNOWCOVER[[10]], ee_SNOWDEPTH[[10]])
  bob_ee_smr <- join_data(crwOut_ALL[[11]], ee_NDVI[[11]], ee_NDVImax[[11]], ee_SNOWCOVER[[11]], ee_SNOWDEPTH[[11]])
  bob_ee_wtr <- join_data(crwOut_ALL[[12]], ee_NDVI[[12]], ee_NDVImax[[12]], ee_SNOWCOVER[[12]], ee_SNOWDEPTH[[12]])
  coy_ee_smr <- join_data(crwOut_ALL[[13]], ee_NDVI[[13]], ee_NDVImax[[13]], ee_SNOWCOVER[[13]], ee_SNOWDEPTH[[13]])
  coy_ee_wtr <- join_data(crwOut_ALL[[14]], ee_NDVI[[14]], ee_NDVImax[[14]], ee_SNOWCOVER[[14]], ee_SNOWDEPTH[[14]])

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
    
    
  



