  #'  =============================================================
  #'  Extract covariate data for projecting RSFs across study areas
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  January 2022
  #'  =============================================================
  #'  Script generated study area-wide grids at varying resolution and extracts
  #'  covariate data for each pixel. Covariate values can be used to project RSF
  #'  results across study area for K-fold cross-validation and to generate maps
  #'  representing probability of use for each species to be used as a covaraite
  #'  for HMM analyses.
  #'  =============================================================
  
  #'  Clear memory
  rm(list=ls())

  #'  Load libraries
  library(adehabitatHR)
  library(sf)
  library(sp)
  library(stars)
  library(rgdal)
  library(rgeos)
  library(raster)
  library(tidyverse)
  library(parallel)
  library(future.apply)
  
  #'  Define desired projections
  sa_proj <- projection("EPSG:2855")  # NAD83(HARN) / Washington North
  wgs84 <- projection("+proj=longlat +datum=WGS84 +no_defs")
  
  #'  Load study area shapefiles
  OK.SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "METHOW_SA") %>%
    st_transform(crs = sa_proj)
  OK.SA <- as(OK.SA, "Spatial")
  NE.SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "NE_SA") %>%
    st_transform(crs = sa_proj)
  NE.SA <- as(NE.SA, "Spatial")
  
  #'  Load telemetry data    
  load("./Outputs/Telemetry_tracks/spp_all_tracks_noDis_noMig.RData")
  
  ####  Create giant MCPs  ####
  #'  Combine all collar data across species for each study area
  md <- rbind(spp_all_tracks[[1]], spp_all_tracks[[2]])
  elk <- rbind(spp_all_tracks[[3]], spp_all_tracks[[4]])
  wtd <- rbind(spp_all_tracks[[5]], spp_all_tracks[[6]])
  coug <- rbind(spp_all_tracks[[7]], spp_all_tracks[[8]])
  wolf <- rbind(spp_all_tracks[[9]], spp_all_tracks[[10]])
  bob <- rbind(spp_all_tracks[[11]], spp_all_tracks[[12]])
  coy <- rbind(spp_all_tracks[[13]], spp_all_tracks[[14]])
  #'  Okanogan
  OK.locs <- rbind(md, coug[coug$StudyArea == "OK",], wolf[wolf$StudyArea == "OK",], bob[bob$StudyArea == "OK",], coy[coy$StudyArea == "OK",])
  #'  Northeast
  NE.locs <- rbind(elk, wtd, coug[coug$StudyArea == "NE",], wolf[wolf$StudyArea == "NE",], bob[bob$StudyArea == "NE",], coy[coy$StudyArea == "NE",])
  
  #'  Make relocations spatial and covert to sp object
  sp_object <- function(locs) {
    sf.ob <- st_as_sf(locs, coords = c("x", "y"), crs = sa_proj)
    sp.ob <- as(sf.ob, "Spatial")
    return(sp.ob)
  }
  NE.locs <- sp_object(NE.locs)
  OK.locs <- sp_object(OK.locs)
  
  #'  Create massive MCP for all collared animals in the study area
  #'  Ignore warnings
  OK.mcp <- mcp(OK.locs, percent = 100)
  NE.mcp <- mcp(NE.locs, percent = 100)
  
  #'  Merge study area and MPCs to extract bbox
  NE.union <- raster::union(NE.SA, NE.mcp); plot(NE.union)
  NE.union.sf <- st_as_sf(NE.union)
  NE.bbox <- bbox(NE.union)
  OK.union <- raster::union(OK.SA, OK.mcp); plot(OK.union)
  OK.union.sf <- st_as_sf(OK.union)
  OK.bbox <- bbox(OK.union)

  
  ####  Create empty rasters  #### 
  #'  WPPP 1km2 grid
  grid_1k <- raster("./Shapefiles/ref_grid_1k.img")
  res(grid_1k); projection(grid_1k)
  bbox_1k <- bbox(grid_1k)
  
  #'  Generate 30m2 grid based on WPPP 1km2 grid
  grid_30m <- raster(extent(bbox_1k), crs = projection(sa_proj), res = 30)
  grid_30m[] <- 1:ncell(grid_30m)
  res(grid_30m); projection(grid_30m)
  
  #'  Crop raster to just the study area unions- will reduce the number of pixels
  #'  that need to extract covariate data for
  NE_1km <- crop(grid_1k, extent(NE.union))
  NE_30m <- crop(grid_30m, extent(NE.union))
  plot(NE_30m); plot(NE.union, add = T)
  
  OK_1km <- crop(grid_1k, extent(OK.union))
  OK_30m <- crop(grid_30m, extent(OK.union))
  plot(OK_30m); plot(OK.union, add = T)

  #'  Convert rasters to pixels (ncessary for extracting from rasters below)
  NE.dots.1km <- as(NE_1km, "SpatialPixelsDataFrame")
  NE.dots.30m <- as(NE_30m, "SpatialPixelsDataFrame")
  OK.dots.1km <- as(OK_1km, "SpatialPixelsDataFrame")
  OK.dots.30m <- as(OK_30m, "SpatialPixelsDataFrame")
  
  #'  Extract coordinates from each pixel (centroid of each cell)
  NE.coord.1km <- coordinates(NE.dots.1km)
  NE.coord.30m <- coordinates(NE.dots.30m)
  OK.coord.1km <- coordinates(OK.dots.1km)
  OK.coord.30m <- coordinates(OK.dots.30m)
  
  #'  Make dots spatial
  NE.centers.1km <- SpatialPoints(NE.dots.1km, proj4string = CRS(sa_proj))
  NE.centers.30m <- SpatialPoints(NE.dots.30m, proj4string = CRS(sa_proj))
  OK.centers.1km <- SpatialPoints(OK.dots.1km, proj4string = CRS(sa_proj))
  OK.centers.30m <- SpatialPoints(OK.dots.30m, proj4string = CRS(sa_proj))

  #'  Reproject to WGS84 (OK.centers take a hot minute to transform)
  NE.centers.1km.wgs84 <- spTransform(NE.centers.1km, crs(wgs84))
  NE.centers.30m.wgs84 <- spTransform(NE.centers.30m, crs(wgs84))
  OK.centers.1km.wgs84 <- spTransform(OK.centers.1km, crs(wgs84))
  OK.centers.30m.wgs84 <- spTransform(OK.centers.30m, crs(wgs84))
  
 
  ####  Extract covariates for each raster pixel  ####
  #'  Read in spatial data (only covs that actually made it into the RSFs)
  #'  Terrain rasters
  dem <- raster("./Shapefiles/WA DEM rasters/WPPP_DEM_30m.tif")
  slope <- raster("./Shapefiles/WA DEM rasters/WPPP_slope_aspect.tif", band = 1)
  #'  Canopy cover change
  canopy18 <- raster("./Shapefiles/Global_Forest_Change/treecov_2018.tif")
  canopy19 <- raster("./Shapefiles/Global_Forest_Change/treecov_2019.tif")
  canopy20 <- raster("./Shapefiles/Global_Forest_Change/treecov_2020.tif")
  #'  Cascadia Biodiveristy Watch rasters & shapefiles
  landcov18 <- raster("./Shapefiles/Cascadia_layers/landcover_2018.tif")
  landcov19 <- raster("./Shapefiles/Cascadia_layers/landcover_2019.tif")
  dist_open18 <- raster("./Shapefiles/Cascadia_layers/Dist2OpenEdge18.tif")
  dist_open19 <- raster("./Shapefiles/Cascadia_layers/Dist2OpenEdge19.tif")
  dist_close18 <- raster("./Shapefiles/Cascadia_layers/Dist2ForestEdge18.tif")
  dist_close19 <- raster("./Shapefiles/Cascadia_layers/Dist2ForestEdge19.tif")
  rdden <- raster("./Shapefiles/Cascadia_layers/roadsForTaylor/RoadDensity_1km.tif")
  #'  Water density
  h2o <- raster("./Shapefiles/WA_DeptEcology_HydroWA/WPPP_dist_to_water_2855.tif")
  #'  Human modified landscapes
  HM <- raster("./Shapefiles/Additional_WPPP_Layers/WPPP_gHM.tif")
  
  #'  Create raster stacks of data with same dimensions and resolution
  terra_stack <- stack(dem, slope)
  dist2edge_stack18 <- stack(dist_close18, dist_open18) # includes open & closed habitats
  dist2edge_stack19 <- stack(dist_close19, dist_open19) # includes open & closed habitats
  
  #'  Check out projections
  projection(dem)
  projection(HM)
  projection(canopy18)
  projection(landcov18)
  projection(dist_open18)
  projection(rdden)
  projection(h2o)
  
  
  #'  COVARIATE EXTRACTION
  #'  =========================
  #'  Monitor time
  start.time <- Sys.time()
  
  #'  Master function to extract and manipulate covaraite data for each raster
  #'  Feed the function raster centroid coordinates in sa_proj & wgs84
  #'  Returns massive dataframes with covaraite values for every pixel
  cov_extract <- function(locs, locs_wgs84) {

        #'  1. Extract data from terrain & anthropogenic rasters
    #'  ----------------------------------------------------
    terrain <- raster::extract(terra_stack, locs_wgs84, df = TRUE) #locs_wgs84
    modified <- raster::extract(HM, locs_wgs84, df = TRUE) #locs_wgs84
    road_den <- raster::extract(rdden, locs, df = TRUE) # note the projection
    water <- raster::extract(h2o, locs, df = TRUE) # note the projection
    #'  Merge into a single data frame of covariates
    join_covs <- terrain %>%
      full_join(modified, by = "ID") %>%
      full_join(road_den, by = "ID") %>%
      full_join(water, by = "ID") %>%
      transmute(
        ID = ID,
        Elev = round(WPPP_DEM_30m, digits = 2),
        Slope = round(WPPP_slope_aspect, digits = 2),
        HumanMod = round(WPPP_gHM, digits = 2),
        RoadDen = round(RoadDensity_1km, digits = 2),
        Dist2Water = round(WPPP_dist_to_water_2855, digits = 2)
      ) %>%
      #'  Need to change NAs to 0 for road density (if NA it means there are no
      #'  roads within that 1km pixel and raster pixel was empty)
      mutate(
        RoadDen = ifelse(is.na(RoadDen), 0, RoadDen),
      )
    
    
    #'  2. Extract annual canopy cover data
    #'  -----------------------------------
    cover18 <- raster::extract(canopy18, locs_wgs84, df = TRUE) %>%
      transmute(
        ID = ID,
        Canopy18 = round(treecov_2018, digits = 2)
      )
    cover19 <- raster::extract(canopy19, locs_wgs84, df = TRUE) %>%
      transmute(
        ID = ID,
        Canopy19 = round(treecov_2019, digits = 2)
      ) 
    cover20 <- raster::extract(canopy20, locs_wgs84, df = TRUE) %>%
      transmute(
        ID = ID,
        Canopy20 = round(treecov_2020, digits = 2)
      ) 
    canopy1820 <- cover18 %>%
      full_join(cover19, by = "ID") %>%
      full_join(cover20, by = "ID")
    
    #'  3. Extract Distance to Edge data, based on Cascadia landcover rasters
    #'  reclassified to Forest (forest + woody wetland) and non-Forest. For each
    #'  pixel then find the minimum distance to open vs closed habitat and retains
    #'  the closest distance to represent distance to edge. Only 2018 & 2019
    #'  landcover rasters available so applying 2019 landcov to 2020 data.
    #'  ------------------------------------------------------------------------
    Dist2Edge18 <- raster::extract(dist2edge_stack18, locs, df = TRUE) %>%
      transmute(
        ID = ID,
        Dist2Forest18 = round(Dist2ForestEdge18, 2),
        Dist2Open18 = round(Dist2OpenEdge18, 2)
      ) %>%
      #'  Combine distance to forest & open edges so one distance to edge value
      #'  for each pixel (essentially drops the NAs and merges into single column)
      rowwise() %>%
      mutate(Dist2Edge18 = min(Dist2Forest18, Dist2Open18, na.rm=TRUE)) 
    Dist2Edge19 <- raster::extract(dist2edge_stack19, locs, df = TRUE) %>%
      transmute(
        ID = ID,
        Dist2Forest19 = round(Dist2ForestEdge19, 2),
        Dist2Open19 = round(Dist2OpenEdge19, 2)
      ) %>%
      rowwise() %>%
      mutate(Dist2Edge19 = min(Dist2Forest19, Dist2Open19, na.rm=TRUE))     
    Dist2Edge <- full_join(Dist2Edge18, Dist2Edge19, by = "ID") %>%
      transmute(
        ID = ID,
        Dist2Edge18 = Dist2Edge18,
        Dist2Edge19 = Dist2Edge19
      )
    
    #'  4. Extract landcover data from Cascadia landcover rasters. Only 2018 & 2019
    #'  landcover rasters available so applying 2019 landcov to 2020 observations.
    #'  ------------------------------------------------------------------------
    landcover18 <- raster::extract(landcov18, locs_wgs84, df = TRUE) %>%
      transmute(
        ID = ID,
        landcov18 = landcover_2018
      ) %>%
      mutate(
        landcover_type18 = ifelse(landcov18 == 101, "Water", landcov18),
        landcover_type18 = ifelse(landcov18 == 111, "Glacier", landcover_type18),
        landcover_type18 = ifelse(landcov18 == 121, "Barren", landcover_type18),
        landcover_type18 = ifelse(landcov18 == 201, "Wetland", landcover_type18),
        landcover_type18 = ifelse(landcov18 == 201, "Wetland", landcover_type18),
        landcover_type18 = ifelse(landcov18 == 202, "Woody Wetland", landcover_type18),
        landcover_type18 = ifelse(landcov18 == 211, "Mesic Grass", landcover_type18),
        landcover_type18 = ifelse(landcov18 == 212, "Xeric Grass", landcover_type18),
        landcover_type18 = ifelse(landcov18 == 221, "Mesic Shrub", landcover_type18),
        landcover_type18 = ifelse(landcov18 == 222, "Xeric Shrub", landcover_type18),
        landcover_type18 = ifelse(landcov18 == 230, "Forest", landcover_type18),
        landcover_type18 = ifelse(landcov18 == 301, "Agriculture", landcover_type18),
        landcover_type18 = ifelse(landcov18 == 310, "Developed", landcover_type18),
        landcover_type18 = ifelse(landcov18 == 331, "Commercial", landcover_type18), 
        landcover_type18 = ifelse(landcov18 == 332, "Developed", landcover_type18)  
      )
    landcover19 <- raster::extract(landcov19, locs_wgs84, df = TRUE) %>%
      transmute(
        ID = ID,
        landcov19 = landcover_2019
      ) %>%
      mutate(
        landcover_type19 = ifelse(landcov19 == 101, "Water", landcov19),
        landcover_type19 = ifelse(landcov19 == 111, "Glacier", landcover_type19),
        landcover_type19 = ifelse(landcov19 == 121, "Barren", landcover_type19),
        landcover_type19 = ifelse(landcov19 == 201, "Wetland", landcover_type19),
        landcover_type19 = ifelse(landcov19 == 201, "Wetland", landcover_type19),
        landcover_type19 = ifelse(landcov19 == 202, "Woody Wetland", landcover_type19),
        landcover_type19 = ifelse(landcov19 == 211, "Mesic Grass", landcover_type19),
        landcover_type19 = ifelse(landcov19 == 212, "Xeric Grass", landcover_type19),
        landcover_type19 = ifelse(landcov19 == 221, "Mesic Shrub", landcover_type19),
        landcover_type19 = ifelse(landcov19 == 222, "Xeric Shrub", landcover_type19),
        landcover_type19 = ifelse(landcov19 == 230, "Forest", landcover_type19),
        landcover_type19 = ifelse(landcov19 == 301, "Agriculture", landcover_type19),
        landcover_type19 = ifelse(landcov19 == 310, "Developed", landcover_type19),
        landcover_type19 = ifelse(landcov19 == 331, "Commercial", landcover_type19), 
        landcover_type19 = ifelse(landcov19 == 332, "Developed", landcover_type19)  
      )
    landcover <- full_join(landcover18, landcover19, by = "ID") %>%
      transmute(
        ID = ID,
        landcover_type18 = landcover_type18,
        landcover_type19 = landcover_type19
      )
      
    #'  5. Join all covariates together & clean up to use for projecting RSFs
    raster_covs <- join_covs %>%
      full_join(canopy1820, by = "ID") %>%
      full_join(Dist2Edge, by = "ID") %>%
      full_join(landcover, by = "ID") %>%
      transmute(
        ID = ID,
        Elev = Elev,
        Slope = Slope,
        RoadDen = RoadDen,
        Dist2Water = Dist2Water,
        HumanMod = HumanMod,
        CanopyCover18 = Canopy18,
        CanopyCover19 = Canopy19,
        CanopyCover20 = Canopy20,
        Dist2Edge18 = Dist2Edge18,
        Dist2Edge19 = Dist2Edge19,
        Landcover_type18 = landcover_type18,
        Landcover_type19 = landcover_type19) 
    
    return(raster_covs)
    
  }
  
  #'  Run list of species used & available location data through function in parallel
  #'  This will take AWHILE even in parallel
  #'  Keep track of which order of selection is being extracted here!!!
  NE.covs.1km <- cov_extract(locs = NE.centers.1km, locs_wgs84 = NE.centers.1km.wgs84)
  NE.covs.30m <- cov_extract(locs = NE.centers.30m, locs_wgs84 = NE.centers.30m.wgs84)
  OK.covs.1km <- cov_extract(locs = OK.centers.1km, locs_wgs84 = OK.centers.1km.wgs84)
  OK.covs.30m <- cov_extract(locs = OK.centers.30m, locs_wgs84 = OK.centers.30m.wgs84)
  

  #'  End time keeping
  end.time <- Sys.time()
  #' #'  Stop running in parallel
  #' parallel::stopCluster(cl)
  #'  How long did this take?
  difftime(end.time, start.time, units = "hours")
  
  
  save(NE.covs.1km, file = paste0("./Outputs/Telemetry_covs/NE_covs_1km_", Sys.Date(), ".RData"))
  save(NE.covs.30m, file = paste0("./Outputs/Telemetry_covs/NE_covs_30m_", Sys.Date(), ".RData"))
  save(OK.covs.1km, file = paste0("./Outputs/Telemetry_covs/OK_covs_1km_", Sys.Date(), ".RData"))
  save(OK.covs.30m, file = paste0("./Outputs/Telemetry_covs/OK_covs_30m_", Sys.Date(), ".RData"))
  
  