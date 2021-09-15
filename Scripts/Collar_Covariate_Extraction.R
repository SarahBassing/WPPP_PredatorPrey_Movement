  #'  ============================================
  #'  Covariate extraction for telemetry locations
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  April 2021
  #'  ============================================
  #'  Script to extract covariate data at each GPS collar location and interpolated
  #'  locations based on the crawlWrap function. Covariates gathered from 
  #'  multiple sources and described in Covariate_Extract.R script from
  #'  WPPP_CameraTrapping repository. Relevant covariates include:
  #'    -Elevation (30m res)
  #'    -Slope (30m res)
  #'    -Human Modified Landscape (1km res)
  #'    -Percent Mixed Forest (within 250m of point)
  #'    -Percent Xeric Grass (within 250m of point)
  #'    -Percent Xeric Shrub (within 205m of point)
  #'    -Distance to nearest road
  #'    -Season
  #'    -Study Area
  #'  ============================================
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(sf)
  library(stars)
  library(rgeos)
  library(raster)
  library(parallel)
  library(doParallel)
  library(future.apply)
  library(tidyverse)
  
  #'  Read in spatial data
  wppp_bound <- st_read("./Shapefiles/WPPP_CovariateBoundary", layer = "WPPP_CovariateBoundary")
  #'  Terrain rasters
  dem <- raster("./Shapefiles/WA DEM rasters/WPPP_DEM_30m.tif")
  Slope <- raster("./Shapefiles/WA DEM rasters/WPPP_slope_aspect.tif", band = 1)
  #'  Cascadia Biodiveristy Watch rasters & shapefiles
  landcov18 <- raster("./Shapefiles/Cascadia_layers/landcover_2018.tif")
  landcov19 <- raster("./Shapefiles/Cascadia_layers/landcover_2019.tif")
  interp_landcov18 <- raster("./Shapefiles/Cascadia_layers/interpolated_landcover_2018.tif")
  interp_landcov19 <- raster("./Shapefiles/Cascadia_layers/interpolated_landcover_2019.tif")
  roads <- st_read("./Shapefiles/Cascadia_layers/roadsForTaylor", layer = "roadsForTaylor")
  #'  Human density and human modified landscapes
  HM <- raster("./Shapefiles/Additional_WPPP_Layers/WPPP_gHM.tif")
  
  #'  Identify projections & resolutions of relevant features
  sa_proj <- projection("+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs ")
  wgs84 <- projection("+proj=longlat +datum=WGS84 +no_defs")
  projection(wppp_bound)
  projection(dem)
  projection(Slope)
  projection(landcov18)
  projection(interp_landcov18)
  projection(HM)
  projection(roads)
  
  res(dem)
  res(landcov18)
  res(interp_landcov18)
  res(HM)
  
  #'  Reproject road shapefile to match animal location projection
  road_reproj <- st_transform(roads, crs = st_crs(sa_proj))
  projection(road_reproj)

  #'  Load animal location data for each species
  load("./Outputs/Telemetry_crwOut/crwOut_ALL_2021-05-03.RData")
  # load("./Outputs/Telemetry_crwOut/crwOut_MD_smr.RData")
  # load("./Outputs/Telemetry_crwOut/crwOut_MD_wtr.RData")
  # load("./Outputs/Telemetry_crwOut/crwOut_ELK_smr.RData")
  # load("./Outputs/Telemetry_crwOut/crwOut_ELK_wtr.RData")
  # load("./Outputs/Telemetry_crwOut/crwOut_WTD_smr.RData")
  # load("./Outputs/Telemetry_crwOut/crwOut_WTD_wtr.RData")
  # load("./Outputs/Telemetry_crwOut/crwOut_COUG_smr.RData")
  # load("./Outputs/Telemetry_crwOut/crwOut_COUG_wtr.RData")
  # load("./Outputs/Telemetry_crwOut/crwOut_WOLF_smr.RData")
  # load("./Outputs/Telemetry_crwOut/crwOut_WOLF_wtr.RData")
  # load("./Outputs/Telemetry_crwOut/crwOut_BOB_smr.RData")
  # load("./Outputs/Telemetry_crwOut/crwOut_BOB_wtr.RData")
  # load("./Outputs/Telemetry_crwOut/crwOut_COY_smr.RData")
  # load("./Outputs/Telemetry_crwOut/crwOut_COY_wtr.RData")
  
  
  #'  Function to make crwOut data a spatial sf object
  spatial_locs <- function(locs) {
    move <- locs[[2]]
    sf_locs <- st_as_sf(move, coords = c("mu.x", "mu.y"), crs = sa_proj)
    return(sf_locs)
  }
  sf_locs <- lapply(crwOut_ALL, spatial_locs)
  
  #'  Extract a few individuals to practice with
  # coy_a <- sf_locs[[1]][sf_locs[[1]]$ID == "NECOY4M_2",]
  # coy_b <- sf_locs[[2]][sf_locs[[2]]$ID == "MVCOY63M_37",]
  # sf_locs <- list(coy_a, coy_b)

  
  ####  COVARIATE EXTRACTION & CALCULATIONS  ####
  #'  Takes forever but running in parallel helps 
  #'  51.49 hrs on lab computer with 30 cores!
  
  #'  Monitor time
  start.time <- Sys.time()

  #'  Setup script to run in parallel
  #'  Extract covariates for each species at once
  #'  Identify how many cores I want to use
  detectCores(logical = FALSE)
  cl <- parallel::makeCluster(4)  # change to 14 when working correctly on lab computer
  #'  Run in parallel on local computer with specified number of cores
  plan(cluster, workers = cl)

  #'  Master function to extract and manipulate covaraite data for each species
  cov_extract <- function(locs) {
    
    #'  1. Extract data from unprojected rasters
    #'  ----------------------------------------
    #'  Reproject location data to match rasters (WGS84)
    reproj_locs <- st_transform(locs, crs = st_crs(wgs84))
    #'  Extract covariates for each location
    elevation <- raster::extract(dem, reproj_locs, df = TRUE)
    slope <- raster::extract(Slope, reproj_locs, df = TRUE)
    modified <- raster::extract(HM, reproj_locs, df = TRUE)
    #'  Merge into a single data frame of covariates
    join_covs <- full_join(elevation, slope, by = "ID") %>%
      full_join(modified, by = "ID") %>%
      transmute(
        obs = ID,
        Elev = WPPP_DEM_30m,
        Slope = round(WPPP_slope_aspect, digits = 2),
        HumanMod = WPPP_gHM
      )
    #'  Pull out unique animal/time information
    animal <- as.data.frame(reproj_locs) %>%
      dplyr::select(c(ID, time))
    #'  Merge animal/time information with covariates
    covs <- as.data.frame(cbind(animal, join_covs))
    
    
    #'  2. Extract data from roads shapefile & calculate distance to nearest road
    #'     for each location
    #'  ------------------------------------------------------------------------
    dist2road <- sapply(1:nrow(locs), function(x) min(st_distance(road_reproj, locs[x, ])))
    dist2road <- as.data.frame(dist2road)
    dist2road$ID <- locs$ID
    dist2road$time <- locs$time
    dist2road$obs <- c(1:nrow(locs))
    #'  Append to covariate data frame
    covs <- covs %>%
      full_join(dist2road, by = c("ID", "time"))
    
    
    #'  3. Extract landcover data from within 250m of each location & calculate
    #'     percent habitat type within that buffer
    #'  ------------------------------------------------------------------------
    #'  Extract landcover value from each pixel within 250m radius of locations
    #'  using interpolated landcover rasters derived from Cascadia landcover
    pixvals18 <- raster::extract(interp_landcov18, locs, factors = TRUE, buffer = 250, df = TRUE)
    pixvals_df18 <- as.data.frame(pixvals18)
    pixvals19 <- raster::extract(interp_landcov19, locs, factors = TRUE, buffer = 250, df = TRUE)
    pixvals_df19 <- as.data.frame(pixvals19)
    #'  Merge together and rename variables
    landcov <- cbind(pixvals_df18, pixvals_df19$interpolated_landcover_2019) 
    colnames(landcov) <- c("obs", "landcover_2018", "landcover_2019")
    #'  Rename categories so they're easier to work with
    landcover <- landcov %>%
      mutate(
        landcover_2018 = ifelse(landcover_2018 == "101", "Water", landcover_2018),
        landcover_2018 = ifelse(landcover_2018 == "121", "Barren", landcover_2018),
        landcover_2018 = ifelse(landcover_2018 == "201", "EmergentWetland", landcover_2018),
        landcover_2018 = ifelse(landcover_2018 == "202", "WoodyWetland", landcover_2018),
        landcover_2018 = ifelse(landcover_2018 == "211", "MesicGrass", landcover_2018),
        landcover_2018 = ifelse(landcover_2018 == "212", "XericGrass", landcover_2018),
        landcover_2018 = ifelse(landcover_2018 == "221", "MesicShrub", landcover_2018),
        landcover_2018 = ifelse(landcover_2018 == "222", "XericShrub", landcover_2018),
        landcover_2018 = ifelse(landcover_2018 == "230", "Forest", landcover_2018),
        landcover_2018 = ifelse(landcover_2018 == "310", "Agriculture", landcover_2018),
        landcover_2018 = ifelse(landcover_2018 == "331", "Commercial", landcover_2018),
        landcover_2018 = ifelse(landcover_2018 == "332", "Residential", landcover_2018),
        landcover_2019 = ifelse(landcover_2019 == "101", "Water", landcover_2019),
        landcover_2019 = ifelse(landcover_2019 == "121", "Barren", landcover_2019),
        landcover_2019 = ifelse(landcover_2019 == "201", "EmergentWetland", landcover_2019),
        landcover_2019 = ifelse(landcover_2019 == "202", "WoodyWetland", landcover_2019),
        landcover_2019 = ifelse(landcover_2019 == "211", "MesicGrass", landcover_2019),
        landcover_2019 = ifelse(landcover_2019 == "212", "XericGrass", landcover_2019),
        landcover_2019 = ifelse(landcover_2019 == "221", "MesicShrub", landcover_2019),
        landcover_2019 = ifelse(landcover_2019 == "222", "XericShrub", landcover_2019),
        landcover_2019 = ifelse(landcover_2019 == "230", "Forest", landcover_2019),
        landcover_2019 = ifelse(landcover_2019 == "310", "Agriculture", landcover_2019),
        landcover_2019 = ifelse(landcover_2019 == "331", "Commercial", landcover_2019),
        landcover_2019 = ifelse(landcover_2019 == "332", "Residential", landcover_2019),
      )
    #'  Add animal ID info to landcover data for further manipulation
    animal <- animal %>% 
      mutate(obs = 1:nrow(.))
    landcover_250m <- full_join(animal, landcover, by = "obs")
    #'  Count the number of cells in each landcover category per location
    #'  2018 landcover and telemetry locations only
    tbl_landcover18 <- as.data.frame(landcover_250m) %>%
      # select(-geometry) %>%
      group_by(obs) %>%
      count(landcover_2018) %>%
      ungroup() %>%
      pivot_wider(names_from = landcover_2018, values_from = n) %>%
      replace(is.na(.), 0) %>% 
      mutate(
        #'  Count number of pixels within 250m of each location
        sumPixels = rowSums(.[2:ncol(.)]),
        #'  Combine similar habitat types
        Forest =  Forest + WoodyWetland + EmergentWetland,
        MesicGrass = MesicGrass,# + Barren,
        MesicMix = MesicShrub + MesicGrass,
        ForestMix = Forest + MesicMix,
        ForestMix2 = Forest + MesicShrub
      ) %>%
      #'  Calculate percent landcover type within 250m of each camera site
      mutate(
        PercForest = round(Forest/sumPixels, 2),
        PercForestMix = round(ForestMix/sumPixels,2),     # Cannot be used in conjunction with Forest or any Mesic landcover types
        PercForestMix2 = round(ForestMix2/sumPixels, 2),
        PercXericShrub = round(XericShrub/sumPixels, 2),
        PercMesicShrub = round(MesicShrub/sumPixels, 2),  # Cannot be used in conjunction with MesicMix
        PercXericGrass = round(XericGrass/sumPixels, 2),
        PercMesicGrass = round(MesicGrass/sumPixels, 2),  # Cannot be used in conjunction with MesicMix
        PercMesicMix = round(MesicMix/sumPixels, 2)      # Cannot be used in conjunction with other Mesic landcover types
      ) %>%
      #'  Join % habitat data to animal location data
      full_join(animal, by = "obs") %>%
      #'  Drop data for year camera was NOT present
      mutate(
        Year = lubridate::year(time),
        Month = lubridate::month(time),
        Season = ifelse(Year == 2018 & Month < 10, "Summer18", NA),
        Season = ifelse(Year == 2018 & Month > 11, "Winter1819", Season),
        Season = ifelse(Year == 2019 & Month < 4, "Winter1819", Season),
        Season = ifelse(Year == 2019 & Month > 6 & Month < 10, "Summer19", Season),
        Season = ifelse(Year == 2019 & Month > 11, "Winter1920", Season),
        Season = ifelse(Year == 2020 & Month < 4, "Winter1920", Season),
        PercForest = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercForest),
        PercForestMix = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercForestMix),
        PercForestMix2 = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercForestMix2),
        PercXericShrub = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercXericShrub),
        PercMesicShrub = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercMesicShrub),
        PercXericGrass = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercXericGrass),
        PercMesicGrass = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercMesicGrass),
        PercMesicMix = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercMesicMix) 
      ) %>%
      #'  Only retain relevant columns
      dplyr::select(obs, sumPixels, PercForest, PercForestMix, PercForestMix2, PercXericShrub,
                    PercMesicShrub, PercXericGrass, PercMesicGrass, PercMesicMix, ID, time, Season) %>%
      #'  Filter out rows with NAs
      filter(!is.na(PercForest))
    #'  Repeat for 2019 landcover data and telemetry locations
    tbl_landcover19 <- as.data.frame(landcover_250m) %>%
      # select(-geometry) %>%
      group_by(obs) %>%
      count(landcover_2019) %>%
      ungroup() %>%
      pivot_wider(names_from = landcover_2019, values_from = n) %>%
      replace(is.na(.), 0) %>% 
      mutate(
        #'  Count number of pixels within 250m of each location
        sumPixels = rowSums(.[2:ncol(.)]),
        #'  Combine similar habitat types
        Forest =  Forest + WoodyWetland + EmergentWetland,
        MesicGrass = MesicGrass, # + Barren,
        MesicMix = MesicShrub + MesicGrass,
        ForestMix = Forest + MesicMix,
        ForestMix2 = Forest + MesicShrub
      ) %>%
      #'  Calculate percent landcover type within 250m of each camera site
      mutate(
        PercForest = round(Forest/sumPixels, 2),
        PercForestMix = round(ForestMix/sumPixels,2),     # Cannot be used in conjunction with Forest or any Mesic landcover types
        PercForestMix2 = round(ForestMix2/sumPixels, 2),
        PercXericShrub = round(XericShrub/sumPixels, 2),
        PercMesicShrub = round(MesicShrub/sumPixels, 2),  # Cannot be used in conjunction with MesicMix
        PercXericGrass = round(XericGrass/sumPixels, 2),
        PercMesicGrass = round(MesicGrass/sumPixels, 2),  # Cannot be used in conjunction with MesicMix
        PercMesicMix = round(MesicMix/sumPixels, 2)      # Cannot be used in conjunction with other Mesic landcover types
      ) %>%
      #'  Join % habitat data to animal location data
      full_join(animal, by = "obs") %>%
      #'  Drop data for year camera was NOT present
      mutate(
        Year = lubridate::year(time),
        Month = lubridate::month(time),
        Season = ifelse(Year == 2018 & Month < 10, "Summer18", NA),
        Season = ifelse(Year == 2018 & Month > 11, "Winter1819", Season),
        Season = ifelse(Year == 2019 & Month < 4, "Winter1819", Season),
        Season = ifelse(Year == 2019 & Month > 6 & Month < 10, "Summer19", Season),
        Season = ifelse(Year == 2019 & Month > 11, "Winter1920", Season),
        Season = ifelse(Year == 2020 & Month < 4, "Winter1920", Season),
        PercForest = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercForest),
        PercForestMix = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercForestMix),
        PercForestMix2 = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercForestMix2),
        PercXericShrub = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercXericShrub),
        PercMesicShrub = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercMesicShrub),
        PercXericGrass = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercXericGrass),
        PercMesicGrass = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercMesicGrass),
        PercMesicMix = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercMesicMix)
      ) %>%
      #'  Only retain relevant columns
      dplyr::select(obs, sumPixels, PercForest, PercForestMix, PercForestMix2, PercXericShrub,
                    PercMesicShrub, PercXericGrass, PercMesicGrass, PercMesicMix, ID, time, Season) %>%
      #'  Filter out rows with NAs
      filter(!is.na(PercForest))
    #'  Merge annual landcover data together so no duplicates
    percHab <- rbind(tbl_landcover18, tbl_landcover19)
    
    #'  4. Join all covatiates together & clean up for inclusion in HMM
    telem_covs <- covs %>%
      full_join(percHab, by = c("ID", "time")) %>%
      transmute(
        ID = ID,
        time = time,
        Season = Season,
        Year = ifelse(Season == "Summer18" | Season == "Winter1819", "Year1", "Year2"),
        Elev = Elev,
        Slope = Slope,
        HumanMod = HumanMod,
        NearestRd = dist2road, 
        PercForMix = PercForestMix2,
        PercXGrass = PercXericGrass,
        PercXShrub = PercXericShrub,
        obs = obs) %>%
      mutate(
        Area = ifelse(grepl("NE", ID), "NE", "OK"),
        Area = ifelse(grepl("MD", ID), "OK", Area),
        Area = ifelse(grepl("EA", ID), "NE", Area),
        Area = ifelse(grepl("ELK", ID), "NE", Area),
        Area = ifelse(grepl("WTD", ID), "NE", Area))
    
    return(telem_covs)
    
  }
  
  #'  Run list of species location data through function in parallel
  #'  This will take AWHILE even in parallel
  # spp_telem_covs <- lapply(sf_locs, cov_extract) # non-parallel approach
  spp_telem_covs <- future_lapply(sf_locs, cov_extract)
  
  
  #'  End time keeping
  end.time <- Sys.time()
  #'  Stop running in parallel
  parallel::stopCluster(cl)
  #'  How long did this take?
  difftime(end.time, start.time, units = "hours")

  #'  Add study area to wolf data
  #'  No easy way of doing this because ID not associated with WPPP study areas
  #'  Double check lists 9 & 10 are wolf summer & winter data
  spp_telem_covs[[9]]$Area <- "NE"   
  spp_telem_covs[[9]] <- mutate(spp_telem_covs[[9]], 
                                Area = ifelse(grepl("W61M", ID), "OK", Area),  #double check no "W71F" in here
                                Area = ifelse(grepl("W88M", ID), "OK", Area),
                                Area = ifelse(grepl("W93M", ID), "OK", Area),
                                Area = ifelse(grepl("W94M", ID), "OK", Area))
  spp_telem_covs[[10]]$Area <- "NE"   
  spp_telem_covs[[10]] <- mutate(spp_telem_covs[[10]], 
                                Area = ifelse(grepl("W61M", ID), "OK", Area),
                                Area = ifelse(grepl("W88M", ID), "OK", Area),
                                Area = ifelse(grepl("W93M", ID), "OK", Area),
                                Area = ifelse(grepl("W94M", ID), "OK", Area))

  #'  Save and hope you never have to run this again!
  save(spp_telem_covs, file = paste0("./Outputs/Telemetry_covs/spp_telem_covs_", Sys.Date(), ".RData"))


  #' #'  Stand-alone functions to extract each covariate
  #' #'  Distance to nearest road
  #' #'  ------------------------------
  #' #'  Calculate distance from each animal location to nearest road and take the 
  #' #'  minimum value -- default units of measurement are in METERS
  #' #'  Make sure features are in either Equidistant or State-specific projection
  #' #'  to preserve distances between features (DO NOT use WGS84)
  #' #'  Albers Equal Area projection good for measuring area
  #' 
  #' 
  #' #'  Function to iterate over all points per species (this takes FOREVER!)
  #' RoadDist <- function(locs) {
  #'   dist2road <- sapply(1:nrow(locs), function(x) min(st_distance(road_reproj, locs[x, ])))
  #'   dist2road <- as.data.frame(dist2road)
  #'   dist2road$ID <- locs$ID
  #'   dist2road$time <- locs$time
  #'   dist2road$obs <- c(1:nrow(locs))
  #'   return(dist2road)
  #' }
  #' 
  #' #'  Run each species through function to calculate distance to nearest road
  #' #'  THIS WILL TAKE FOREVER!!! 
  #' # locs <- list(md_locs, elk_locs, wtd_locs, coug_locs, wolf_locs, bob_locs, coy_locs)
  #' # spp_nearestRd <- lapply(locs, RoadDist)
  #' 
  #' coy_a <- sf_locs[[1]][sf_locs[[1]]$ID == "NECOY4M_2",]
  #' coy_b <- sf_locs[[2]][sf_locs[[2]]$ID == "MVCOY63M_37",]
  #' coy_tst <- rbind(coy_a, coy_b)
  #' coy_rd_tst <- RoadDist(coy_tst)
  #' 
  #' md_rd_smr <- RoadDist(md_locs_smr)
  #' md_rd_wtr <- RoadDist(md_locs_wtr)
  #' elk_rd_smr <- RoadDist(elk_locs_smr)
  #' elk_rd_wtr <- RoadDist(elk_locs_wtr)
  #' wtd_rd_smr <- RoadDist(wtd_locs_smr)
  #' wtd_rd_wtr <- RoadDist(wtd_locs_wtr)
  #' coug_rd_smr <- RoadDist(coug_locs_smr)
  #' coug_rd_wtr <- RoadDist(coug_locs_wtr)
  #' wolf_rd_smr <- RoadDist(wolf_locs_smr)
  #' wolf_rd_wtr <- RoadDist(wolf_locs_wtr)
  #' bob_rd_smr <- RoadDist(bob_locs_smr)
  #' bob_rd_wtr <- RoadDist(bob_locs_wtr)
  #' coy_rd_smr <- RoadDist(coy_locs_smr)
  #' coy_rd_wtr <- RoadDist(coy_locs_wtr)
  #' 
  #' # save(wolf_rd, file = "./Outputs/Telemetry_covs/wolf_rd.RData")
  #' # save(elk_rd_smr, file = "./Outputs/Telemetry_covs/elk_rd_smr.RData")
  #' # save(elk_rd_wtr, file = "./Outputs/Telemetry_covs/elk_rd_wtr.RData")
  #' 
  #' 
  #' 
  #' #'  Terrain & Human Modified covariates
  #' #'  -----------------------------------
  #' 
  #' #'  Function to extract covariates from raster data for each species
  #' rast_extract <- function(locs) {
  #'   #'  Reproject to match rasters projections (WGS84)
  #'   reproj_locs <- st_transform(locs, crs = st_crs(wgs84))
  #'   #'  Extract covariates that share same projection
  #'   elevation <- raster::extract(dem, reproj_locs, df = TRUE)
  #'   slope <- raster::extract(Slope, reproj_locs, df = TRUE)
  #'   modified <- raster::extract(HM, reproj_locs, df = TRUE)
  #'   join_covs <- full_join(elevation, slope, by = "ID") %>%
  #'     full_join(modified, by = "ID") %>%
  #'     transmute(
  #'       obs = ID,
  #'       Elev = WPPP_DEM_30m,
  #'       Slope = round(WPPP_slope_aspect, digits = 2),
  #'       HumanMod = WPPP_gHM
  #'     )
  #'   #'  Pull out unique animal/time information
  #'   animal <- as.data.frame(reproj_locs) %>%
  #'     select(c(ID, time))
  #'   #'  Merge animal/time information with covariates
  #'   covs <- as.data.frame(cbind(animal, join_covs))
  #'   return(covs)
  #' }
  #' 
  #' #'  Run each species through raster extract function
  #' locs <- list(reproj_md_smr, reproj_elk_smr, reproj_wtd_smr, reproj_coug_smr, 
  #'              reproj_wolf_smr, reproj_bob_smr, reproj_coy_smr, reproj_md_wtr, 
  #'              reproj_elk_wtr, reproj_wtd_wtr, reproj_coug_wtr, reproj_wolf_wtr, 
  #'              reproj_bob_wtr, reproj_coy_wtr)
  #' spp_extract <- lapply(locs, rast_extract)
  #' 
  #' #'  Unlist extracted data into separate data frames for each species
  #' names(spp_extract) <- c("md_extract_smr", "elk_extract_smr", "wtd_extract_smr", 
  #'                         "coug_extract_smr", "wolf_extract_smr", "bob_extract_smr", 
  #'                         "coy_extract_smr", "md_extract_wtr", "elk_extract_wtr", 
  #'                         "wtd_extract_wtr", "coug_extract_wtr", "wolf_extract_wtr", 
  #'                         "bob_extract_wtr", "coy_extract_wtr")
  #' invisible(lapply(names(spp_extract),function(x) assign(x,spp_extract[[x]],.GlobalEnv)))
  #' 
  #' #'  Save
  #' save(spp_extract, file = "./Outputs/Telemetry_covs/spp_extract.RData")  
  #' 
  #' # md_tst <- reproj_md[reproj_md$ID == "3959MD17_1",]
  #' # md_cov_tst <- spp_extract(md_tst)
  #' 
  #' 
  #' 
  #' #'  Percent Land Cover
  #' #'  ------------------------
  #' #'  Function to extract landcover value from each pixel within 250m radius of 
  #' #'  animal locations using interpolated landcover rasters derived from 
  #' #'  Cascadia landcover
  #' landcov250 <- function(locs) {
  #'   pixvals18 <- raster::extract(interp_landcov18, locs, factors = TRUE, buffer = 250, df = TRUE)
  #'   pixvals_df18 <- as.data.frame(pixvals18)
  #'   pixvals19 <- raster::extract(interp_landcov19, locs, factors = TRUE, buffer = 250, df = TRUE)
  #'   pixvals_df19 <- as.data.frame(pixvals19)
  #'   #'  Merge together and rename variables
  #'   landcov <- cbind(pixvals_df18, pixvals_df19$interpolated_landcover_2019) 
  #'   colnames(landcov) <- c("obs", "landcover_2018", "landcover_2019")
  #'   landcover <- landcov %>%
  #'     mutate(
  #'       landcover_2018 = ifelse(landcover_2018 == "101", "Water", landcover_2018),
  #'       landcover_2018 = ifelse(landcover_2018 == "121", "Barren", landcover_2018),
  #'       landcover_2018 = ifelse(landcover_2018 == "201", "EmergentWetland", landcover_2018),
  #'       landcover_2018 = ifelse(landcover_2018 == "202", "WoodyWetland", landcover_2018),
  #'       landcover_2018 = ifelse(landcover_2018 == "211", "MesicGrass", landcover_2018),
  #'       landcover_2018 = ifelse(landcover_2018 == "212", "XericGrass", landcover_2018),
  #'       landcover_2018 = ifelse(landcover_2018 == "221", "MesicShrub", landcover_2018),
  #'       landcover_2018 = ifelse(landcover_2018 == "222", "XericShrub", landcover_2018),
  #'       landcover_2018 = ifelse(landcover_2018 == "230", "Forest", landcover_2018),
  #'       landcover_2018 = ifelse(landcover_2018 == "310", "Agriculture", landcover_2018),
  #'       landcover_2018 = ifelse(landcover_2018 == "331", "Commercial", landcover_2018),
  #'       landcover_2018 = ifelse(landcover_2018 == "332", "Residential", landcover_2018),
  #'       landcover_2019 = ifelse(landcover_2019 == "101", "Water", landcover_2019),
  #'       landcover_2019 = ifelse(landcover_2019 == "121", "Barren", landcover_2019),
  #'       landcover_2019 = ifelse(landcover_2019 == "201", "EmergentWetland", landcover_2019),
  #'       landcover_2019 = ifelse(landcover_2019 == "202", "WoodyWetland", landcover_2019),
  #'       landcover_2019 = ifelse(landcover_2019 == "211", "MesicGrass", landcover_2019),
  #'       landcover_2019 = ifelse(landcover_2019 == "212", "XericGrass", landcover_2019),
  #'       landcover_2019 = ifelse(landcover_2019 == "221", "MesicShrub", landcover_2019),
  #'       landcover_2019 = ifelse(landcover_2019 == "222", "XericShrub", landcover_2019),
  #'       landcover_2019 = ifelse(landcover_2019 == "230", "Forest", landcover_2019),
  #'       landcover_2019 = ifelse(landcover_2019 == "310", "Agriculture", landcover_2019),
  #'       landcover_2019 = ifelse(landcover_2019 == "331", "Commercial", landcover_2019),
  #'       landcover_2019 = ifelse(landcover_2019 == "332", "Residential", landcover_2019),
  #'     )
  #'   animal <- locs %>% 
  #'     select(ID, time) %>%
  #'     mutate(obs = 1:nrow(.))
  #'   landcover_250m <- full_join(animal, landcover, by = "obs")
  #'   
  #'   return(landcover_250m)
  #' }
  #' 
  #' #'  Run each species through function to extract land cover values within
  #' #'  250 m of each location. FYI, this takes a long time...
  #' # wolf_tst <- wolf_locs[wolf_locs$ID == "W61M_1" | wolf_locs$ID == "W48F_12",]
  #' # tst <- landcov250(wolf_tst)
  #' # wolf_landcov <- landcov250(wolf_locs)
  #' # save(wolf_landcov, file = "./Outputs/Telemetry_covs/wolf_landcov.RData")
  #' # md_landcov_tst <- landcov250(md_tst)
  #' # unique(md_landcov_tst$landcover_2018); unique(md_landcov_tst$landcover_2019)
  #' 
  #' md_landcov_smr <- landcov250(md_locs_smr)
  #' md_landcov_wtr <- landcov250(md_locs_wtr)
  #' elk_landcov_smr <- landcov250(elk_locs_smr)
  #' elk_landcov_wtr <- landcov250(elk_locs_wtr)
  #' wtd_landcov_smr <- landcov250(wtd_locs_smr)
  #' wtd_landcov_wtr <- landcov250(wtd_locs_wtr)
  #' coug_landcov_smr <- landcov250(coug_locs_smr)
  #' coug_landcov_wtr <- landcov250(coug_locs_wtr)
  #' wolf_landcov_smr <- landcov250(wolf_locs_smr)
  #' wolf_landcov_wtr <- landcov250(wolf_locs_wtr)
  #' bob_landcov_smr <- landcov250(bob_locs_smr)
  #' bob_landcov_wtr <- landcov250(bob_locs_wtr)
  #' coy_landcov_smr <- landcov250(coy_locs_smr)
  #' coy_landcov_wtr <- landcov250(coy_locs_wtr)
  #' 
  #' 
  #' #'  Function to calculate the percent of each landcover type within 250m of 
  #' #'  each telemetry location
  #' perc_landcov <- function(landcover, locs) {
  #'   #'  Create base dataframe from animal ID and location
  #'   animal <- as.data.frame(locs) %>%
  #'     select(c(ID, time)) %>%
  #'     mutate(
  #'       obs = 1:nrow(.)
  #'     )
  #'   #'  Count the number of cells in each landcover category by CameraLocation
  #'   tbl_landcover18 <- as.data.frame(landcover) %>%
  #'     #select(-geometry) %>%
  #'     group_by(obs) %>%
  #'     count(landcover_2018) %>%
  #'     ungroup() %>%
  #'     pivot_wider(names_from = landcover_2018, values_from = n) %>%
  #'     replace(is.na(.), 0) %>% 
  #'     mutate(
  #'       #'  Count number of pixels within 250m of each location
  #'       sumPixels = rowSums(.[2:ncol(.)]),
  #'       #'  Combine similar habitat types
  #'       Forest =  Forest + WoodyWetland + EmergentWetland,
  #'       MesicGrass = MesicGrass + Barren,
  #'       MesicMix = MesicShrub + MesicGrass,
  #'       ForestMix = Forest + MesicMix,
  #'       ForestMix2 = Forest + MesicShrub
  #'     ) %>%
  #'     #'  Calculate percent landcover type within 250m of each camera site
  #'     mutate(
  #'       PercForest = round(Forest/sumPixels, 2),
  #'       PercForestMix = round(ForestMix/sumPixels,2),     # Cannot be used in conjunction with Forest or any Mesic landcover types
  #'       PercForestMix2 = round(ForestMix2/sumPixels, 2),
  #'       PercXericShrub = round(XericShrub/sumPixels, 2),
  #'       PercMesicShrub = round(MesicShrub/sumPixels, 2),  # Cannot be used in conjunction with MesicMix
  #'       PercXericGrass = round(XericGrass/sumPixels, 2),
  #'       PercMesicGrass = round(MesicGrass/sumPixels, 2),  # Cannot be used in conjunction with MesicMix
  #'       PercMesicMix = round(MesicMix/sumPixels, 2)      # Cannot be used in conjunction with other Mesic landcover types
  #'     ) %>%
  #'     #'  Join % habitat data to animal location data
  #'     full_join(animal, by = "obs") %>%
  #'     #'  Drop data for year camera was NOT present
  #'     mutate(
  #'       Year = lubridate::year(time),
  #'       Month = lubridate::month(time),
  #'       Season = ifelse(Year == 2018 & Month < 10, "Summer18", NA),
  #'       Season = ifelse(Year == 2018 & Month > 11, "Winter1819", Season),
  #'       Season = ifelse(Year == 2019 & Month < 4, "Winter1819", Season),
  #'       Season = ifelse(Year == 2019 & Month > 6 & Month < 10, "Summer19", Season),
  #'       Season = ifelse(Year == 2019 & Month > 11, "Winter1920", Season),
  #'       Season = ifelse(Year == 2020 & Month < 4, "Winter1920", Season),
  #'       PercForest = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercForest),
  #'       PercForestMix = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercForestMix),
  #'       PercForestMix2 = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercForestMix2),
  #'       PercXericShrub = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercXericShrub),
  #'       PercMesicShrub = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercMesicShrub),
  #'       PercXericGrass = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercXericGrass),
  #'       PercMesicGrass = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercMesicGrass),
  #'       PercMesicMix = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercMesicMix) 
  #'     ) %>%
  #'     #'  Only retain relevent columns
  #'     select(obs, sumPixels, PercForest, PercForestMix, PercForestMix2, PercXericShrub,
  #'                PercMesicShrub, PercXericGrass, PercMesicGrass, PercMesicMix, ID, time, Season) %>%
  #'     #'  Filter out rows with NAs
  #'     filter(!is.na(PercForest))
  #' 
  #'   
  #'   tbl_landcover19 <- as.data.frame(landcover) %>%
  #'     select(-geometry) %>%
  #'     group_by(obs) %>%
  #'     count(landcover_2019) %>%
  #'     ungroup() %>%
  #'     pivot_wider(names_from = landcover_2019, values_from = n) %>%
  #'     replace(is.na(.), 0) %>% 
  #'     mutate(
  #'       #'  Count number of pixels within 250m of each location
  #'       sumPixels = rowSums(.[2:ncol(.)]),
  #'       #'  Combine similar habitat types
  #'       Forest =  Forest + WoodyWetland + EmergentWetland,
  #'       MesicGrass = MesicGrass + Barren,
  #'       MesicMix = MesicShrub + MesicGrass,
  #'       ForestMix = Forest + MesicMix,
  #'       ForestMix2 = Forest + MesicShrub
  #'     ) %>%
  #'     #'  Calculate percent landcover type within 250m of each camera site
  #'     mutate(
  #'       PercForest = round(Forest/sumPixels, 2),
  #'       PercForestMix = round(ForestMix/sumPixels,2),     # Cannot be used in conjunction with Forest or any Mesic landcover types
  #'       PercForestMix2 = round(ForestMix2/sumPixels, 2),
  #'       PercXericShrub = round(XericShrub/sumPixels, 2),
  #'       PercMesicShrub = round(MesicShrub/sumPixels, 2),  # Cannot be used in conjunction with MesicMix
  #'       PercXericGrass = round(XericGrass/sumPixels, 2),
  #'       PercMesicGrass = round(MesicGrass/sumPixels, 2),  # Cannot be used in conjunction with MesicMix
  #'       PercMesicMix = round(MesicMix/sumPixels, 2)      # Cannot be used in conjunction with other Mesic landcover types
  #'     ) %>%
  #'     #'  Join % habitat data to animal location data
  #'     full_join(animal, by = "obs") %>%
  #'     #'  Drop data for year camera was NOT present
  #'     mutate(
  #'       Year = lubridate::year(time),
  #'       Month = lubridate::month(time),
  #'       Season = ifelse(Year == 2018 & Month < 10, "Summer18", NA),
  #'       Season = ifelse(Year == 2018 & Month > 11, "Winter1819", Season),
  #'       Season = ifelse(Year == 2019 & Month < 4, "Winter1819", Season),
  #'       Season = ifelse(Year == 2019 & Month > 6 & Month < 10, "Summer19", Season),
  #'       Season = ifelse(Year == 2019 & Month > 11, "Winter1920", Season),
  #'       Season = ifelse(Year == 2020 & Month < 4, "Winter1920", Season),
  #'       PercForest = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercForest),
  #'       PercForestMix = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercForestMix),
  #'       PercForestMix2 = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercForestMix2),
  #'       PercXericShrub = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercXericShrub),
  #'       PercMesicShrub = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercMesicShrub),
  #'       PercXericGrass = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercXericGrass),
  #'       PercMesicGrass = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercMesicGrass),
  #'       PercMesicMix = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercMesicMix)
  #'     ) %>%
  #'     #'  Only retain relevant columns
  #'     select(obs, sumPixels, PercForest, PercForestMix, PercForestMix2, PercXericShrub,
  #'            PercMesicShrub, PercXericGrass, PercMesicGrass, PercMesicMix, ID, time, Season) %>%
  #'     #'  Filter out rows with NAs
  #'     filter(!is.na(PercForest))
  #'   
  #'   #'  Merge annual data together so no duplicates
  #'   tbl_landcover <- rbind(tbl_landcover18, tbl_landcover19)
  #' 
  #'   
  #'   return(tbl_landcover)
  #' }
  #' 
  #' #'  Run extracted landcover pixels within 250m of each point AND location data 
  #' #'  through function to calculate percent habitat type at each location
  #' # wolf_percHab <- perc_landcov(wolf_landcov, wolf_locs)
  #' # save(wolf_percHab, file = "./Outputs/Telemetry_covs/wolf_percHab.RData")
  #' # md_percHab_tst <- perc_landcov(md_landcov_tst, md_tst)
  #' 
  #' md_percHab_smr <- perc_landcov(md_landcov_smr, md_locs_smr)
  #' md_percHab_wtr <- perc_landcov(md_landcov_wtr, md_locs_wtr)
  #' elk_percHab_smr <- perc_landcov(elk_landcov_smr, elk_locs_smr)
  #' elk_percHab_wtr <- perc_landcov(elk_landcov_wtr, elk_locs_wtr)
  #' wtd_percHab_smr <- perc_landcov(wtd_landcov_smr, wtd_locs_smr)
  #' wtd_percHab_wtr <- perc_landcov(wtd_landcov_wtr, wtd_locs_wtr)
  #' coug_percHab_smr <- perc_landcov(coug_landcov_smr, coug_locs_smr)
  #' coug_percHab_wtr <- perc_landcov(coug_landcov_wtr, coug_locs_wtr)
  #' wolf_percHab_smr <- perc_landcov(wolf_landcov_smr, wolf_locs_smr)
  #' wolf_percHab_wtr <- perc_landcov(wolf_landcov_wtr, wolf_locs_wtr)
  #' bob_percHab_smr <- perc_landcov(bob_landcov_smr, bob_locs_smr)
  #' bob_percHab_wtr <- perc_landcov(bob_landcov_wtr, bob_locs_wtr)
  #' coy_percHab_smr <- perc_landcov(coy_landcov_smr, coy_locs_smr)
  #' coy_percHab_wtr <- perc_landcov(coy_landcov_wtr, coy_locs_wtr)
  #' 
  #' 
  #' 
  #' #'  Combine all covariate data together and save for HMM analyses
  #' all_covs <- function(locs, elev_etc, rd, percHab) {
  #'   telem_covs <- locs %>%
  #'     select(c(ID, time)) %>%
  #'     full_join(elev_etc, by = c("ID", "time")) %>%
  #'     full_join(rd, by = c("ID", "time")) %>%
  #'     full_join(percHab, by = c("ID", "time")) %>%
  #'     transmute(
  #'       ID = ID,
  #'       time = time,
  #'       Season = Season,
  #'       Year = ifelse(Season == "Summer18" | Season == "Winter1819", "Year1", "Year2"),
  #'       Elev = Elev,
  #'       Slope = Slope,
  #'       HumanMod = HumanMod,
  #'       NearestRd = dist2road, 
  #'       PercForMix = PercForestMix2,
  #'       PercXGrass = PercXericGrass,
  #'       PercXShrub = PercXericShrub,
  #'       obs = obs) %>%
  #'     mutate(
  #'       Area = ifelse(grepl("NE", ID), "NE", "OK"),
  #'       Area = ifelse(grepl("MD", ID), "OK", Area),
  #'       Area = ifelse(grepl("EA", ID), "NE", Area),
  #'       Area = ifelse(grepl("ELK", ID), "NE", Area),
  #'       Area = ifelse(grepl("WTD", ID), "NE", Area))
  #' }
  #' 
  #' #'  Merge all covariates together
  #' md_telem_covs_smr <- all_covs(md_locs_smr, md_extract_smr, md_rd_smr, md_percHab_smr)
  #' md_telem_covs_wtr <- all_covs(md_locs_wtr, md_extract_wtr, md_rd_wtr, md_percHab_wtr)
  #' elk_telem_covs_smr <- all_covs(elk_locs_smr, elk_extract_smr, elk_rd_smr, elk_percHab_smr)
  #' elk_telem_covs_wtr <- all_covs(elk_locs_wtr, elk_extract_wtr, elk_rd_wtr, elk_percHab_wtr)
  #' wtd_telem_covs_smr <- all_covs(wtd_locs_smr, wtd_extract_smr, wtd_rd_smr, wtd_percHab_smr)
  #' wtd_telem_covs_wtr <- all_covs(wtd_locs_wtr, wtd_extract_wtr, wtd_rd_wtr, wtd_percHab_wtr)
  #' coug_telem_covs_smr <- all_covs(coug_locs_smr, coug_extract_smr, coug_rd_smr, coug_percHab_smr)
  #' coug_telem_covs_wtr <- all_covs(coug_locs_wtr, coug_extract_wtr, coug_rd_wtr, coug_percHab_wtr)
  #' wolf_telem_covs_smr <- all_covs(wolf_locs_smr, wolf_extract_smr, wolf_rd_smr, wolf_percHab_smr)
  #' wolf_telem_covs_wtr <- all_covs(wolf_locs_wtr, wolf_extract_wtr, wolf_rd_wtr, wolf_percHab_wtr)
  #' bob_telem_covs_smr <- all_covs(bob_locs_smr, bob_extract_smr, bob_rd_smr, bob_percHab_smr)
  #' bob_telem_covs_wtr <- all_covs(bob_locs_wtr, bob_extract_wtr, bob_rd_wtr, bob_percHab_wtr)
  #' coy_telem_covs_smr <- all_covs(coy_locs_smr, coy_extract_smr, coy_rd_smr, coy_percHab_smr)
  #' coy_telem_covs_wtr <- all_covs(coy_locs_wtr, coy_extract_wtr, coy_rd_wtr, coy_percHab_wtr)
  #' 
  #' #'  Add study area to wolf data (no easy way of doing this b/c ID not associated with WPPP)
  #' wolf_telem_covs_smr <- wolf_telem_covs_smr %>%
  #'   mutate(
  #'     Area = ifelse(grepl("W61M", ID), "OK", "NE"),  #double check no "W71F" in here
  #'     Area = ifelse(grepl("W88M", ID), "OK", Area),  
  #'     Area = ifelse(grepl("W93M", ID), "OK", Area),
  #'     Area = ifelse(grepl("W94M", ID), "OK", Area),
  #'   )
  #' wolf_telem_covs_wtr <- wolf_telem_covs_wtr %>%
  #'   mutate(
  #'     Area = ifelse(grepl("W61M", ID), "OK", "NE"),
  #'     Area = ifelse(grepl("W88M", ID), "OK", Area),
  #'     Area = ifelse(grepl("W93M", ID), "OK", Area),
  #'     Area = ifelse(grepl("W94M", ID), "OK", Area),
  #'   )
  #' 
  #' #'  Save final covariates
  #' all_telem_covs_list <- list(coy_telem_covs_smr, coy_telem_covs_wtr)
  #' save(all_telem_covs_list, file = "./Outputs/Telemetry_covs/all_telem_covs_list.RData")
  #' save(coy_telem_covs_smr, file = "./Outputs/Telemetry_covs/coy_telem_covs_smr.RData")
  #' save(coy_telem_covs_wtr, file = "./Outputs/Telemetry_covs/coy_telem_covs_wtr.RData")


  
  
  
  
  