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
  match_ee_data <- function(datasf){ #tempwin, imagecoll, band, sp.res, tmp.res, start.date, end.date, 
    
    #'  Function to add property with time in milliseconds
    add_date <- function(feature) {
      date <- ee$Date(ee$String(feature$get("Date")))$millis()
      feature$set(list(date_millis=date))
    }
    
    #'  Join EE image & telemetry locations based on a maxDifference Filter within 
    #'  a specified temporal window. Set temporal window in days for filter. This 
    #'  will depend on the remote sensing data used.
    tempwin <- 16#tempwin#16 # eventually un-hardcode this
    
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
      matchKey ="bestImage",
      measureKey ="timeDiff"
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
      pixel_value <- img1$sample(region = point, scale = 250, tileScale = 16, dropNulls = F) #scale = 250 #tileScale = 16 #scale = sp.res, tileScale = tmp.res,
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
    
    #'  Define start and end date range based on timing of relocation data
    start <- "2018-06-30" #start.date#
    end <- "2020-11-01" #end.date#
    
    #'  Define EE image collection you want to extract from
    #'  Try to un-hardcode this eventually
    imagecoll <- ee$ImageCollection('MODIS/006/MOD13Q1')$filterDate(start,end) # MODIS NDVI/EVI imagecoll
    # imagecoll <- ee.ImageCollection("MODIS/006/MOD10A1")$filterDate(start,end) # MODIS Terra Snow Cover Daily Global 500m
    
    #'  Name the data band to use (based on EE image)
    #'  Try to un-hardcode this eventually
    band <- "NDVI"#band #"NDVI" 
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
    dataoutput
    
    #'  NDVI data are extracted on one scale but need to be re-scaled to a more
    #'  meaningful range of values
    dataoutput <- mutate(dataoutput, NDVI_scale = NDVI*0.0001)
    #  NOTE: Scale Factor is 0.0001 for the NDVI values (pg. 9 of MODIS User's Guide)
    #  so all NDVI values need to be rescaled by 0.0001 to be in the right range (-1 to 1)
    
    return(dataoutput)
  }
  #'  Feed reformated animal location data through this monster function to
  #'  match & extract EE images with spatiotemporally varying location data.
  tst <- list(data_ee_list[[13]], data_ee_list[[14]])
  imagecoll <- ee$ImageCollection('MODIS/006/MOD13Q1')$filterDate(start,end) # MODIS NDVI/EVI
  tst2 <- lapply(tst, match_ee_data) #tempwin = 16, imagecoll = imagecoll, band = "NDVI", start.date = "2018-06-30", end.date = "2020-11-01", sp.res = 250, tmp.res = 16, 
  
  
  
  
  
  
  
  
  
  
    
    
  



