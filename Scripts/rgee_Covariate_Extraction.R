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
  #'  These original functions are now imbedded in a monster function
  match_ee_data <- function(datasf, imagecoll, tempwin, band, sp.res, tmp.res){ #tempwin, imagecoll, band, sp.res, tmp.res, start.date, end.date, 
    
    #'  Function to add property with time in milliseconds
    add_date <- function(feature) {
      date <- ee$Date(ee$String(feature$get("Date")))$millis()
      feature$set(list(date_millis = date))
    }
    
    #'  Join EE image & telemetry locations based on a maxDifference Filter within 
    #'  a specified temporal window. Set temporal window in days for filter. This 
    #'  will depend on the remote sensing data used.
    tempwin <- tempwin#16 # eventually un-hardcode this
    
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
      pixel_value <- img1$sample(region = point, scale = sp.res, tileScale = tmp.res, dropNulls = F) #scale = 250 #tileScale = 16 #scale = sp.res, tileScale = tmp.res,
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
    
    #' #'  Define start and end date range based on timing of relocation data
    #' start <- "2018-06-30" #start.date#
    #' end <- "2020-11-01" #end.date#
    
    #'  Define EE image collection you want to extract from
    #'  Try to un-hardcode this eventually
    imagecoll <- imagecoll#ee$ImageCollection('MODIS/006/MOD13Q1')$filterDate(start,end) # MODIS NDVI/EVI imagecoll
    
    #'  Name the data band to use (based on EE image)
    #'  Try to un-hardcode this eventually
    band <- band #"NDVI" 
    # band <- "NDSI" # https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD10A1#bands
    
    #'  Chunk location data into groups to loop through when extracting pixel values
    #'  This is necessary if using the getInfo argument in when converting location
    #'  data to sf object (in ee_as_sf). This is not necessary if using the drive
    #'  or gsc arguments which export data through Google Drive/Cloud. See section
    #'  0.5 in Crego et al. (2021) example for further details
    #'  
    #'  This is for up to 1 million points. To increase the max number of points, 
    #'  increase the value for max repetitions. To change the number of points to 
    #'  run per time, change the value in the argument each (up to 5000).
    datasf$uniq <- rep(1:1000, each=1000)[1:nrow(datasf)] 
    
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
    names(dataoutput)[4] <- band
    
    #'  View & return output
    print(dataoutput)
    return(dataoutput)
  }
  #'  Define time window of interest for extracting data
  start <- "2018-06-30"
  end <- "2020-11-01" 
  #'  Define EE image collections
  imageNDVI <- ee$ImageCollection('MODIS/006/MOD13Q1')$filterDate(start,end) # MODIS Terra NDVI/EVI 16day, 250m resolution
  imageSNOW <- ee.ImageCollection("MODIS/006/MOD10A1")$filterDate(start,end) # MODIS Terra Snow Cover Daily Global daily, 500m resolution
  
  
  #'  Run reformated animal location data through this monster function to
  #'  match & extract EE images 
  tst <- list(data_ee_list[[13]], data_ee_list[[14]])  #, data_ee_list[[14]]
  ee_NDVI <- lapply(tst, match_ee_data, imagecoll = imageNDVI, tempwin = 16, band = "NDVI", sp.res = 250, tmp.res = 16)  
  ee_SNOW <- lapply(tst, match_ee_data, imagecoll = imageSNOW, tempwin = 1, band = "SnowCover", sp.res = 500, tmp.res = 1) 
  
  
  ####  Re-scale NDVI values  ####
  #'  ONLY RUN IF WORKING WITH NDVI DATA!
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
  
  ####  Join datasets  ####
  join_data <- function(crwOut_data, ndvi) { #, ndvi, snow
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
    relocate(ID, .after = (NDVI_scale))
    colnames(ee_covs) <- c("AnimalID", "Season", "StudyArea", "mu.x", "mu.y", "Date", "NDVI", "ID")
    
    return(ee_covs)
  }
  coy_ee <- join_data(crwOut_ALL[[13]], ee_NDVI[[13]])

  
  
  #'  Citation for Snow Cover Data
  #'  Hall, D. K., V. V. Salomonson, and G. A. Riggs. 2016. MODIS/Terra Snow Cover Daily L3 Global 500m Grid. Version 6. Boulder, Colorado USA: NASA National Snow and Ice Data Center Distributed Active Archive Center.
  
  
    
    
  



