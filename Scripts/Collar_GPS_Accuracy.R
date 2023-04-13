  ##  =========================================================
  ##  GPS collar accuracy
  ##  Washington Predator-Prey Project
  ##  April 2023
  ##  Sarah Bassing
  ##  =========================================================
  ##  Calculate proportion of GPS collars that had low accuracy fixes (DOP > 10).
  ##  From the Vectronic website:  "DOP is the abbreviation of “Dilution of Precision” 
  ##  and gives an idea of the location’s precision. High values indicate that the 
  ##  location is likely to be inaccurate while low values indicate a better precision, 
  ##  thus likely a more accurate position. The DOP for stored 2D positions displays 
  ##  HDOP horizontal value with the fix, where 3D locations are stored with PDOP 
  ##  position value in your GPS data file." https://www.vectronic-aerospace.com/faq/
  ##  2D locations are less accurate than 3D locations.
  ##  Plots relating DOP to location error in meters
  ##  https://www.researchgate.net/figure/Relationship-between-location-error-and-PDOP-values-for-all-individual-a-2-D-locations_fig2_229452058
  ##  
  ##  Final data cleaning resulted in no 2D locations and no locations with 
  ##  VEC_Height <0 or >2000 (arbitrary cutoffs). DOP >10 is not good (~20-30 meter 
  ##  accuracy) so used DOP = 8 as cutoff
  ##  =========================================================
  
  
  #'  Libraries
  library(tidyverse)

  #'  Load data
  #'  ALL locations (including low accuracy fixes)
  md_master <- read.csv("./Data/dev_telem_md_11.16.20.csv")
  elk_master <- read.csv("./Data/dev_telem_elk_11.16.20.csv")
  wtd_master <- read.csv("./Data/dev_telem_wtd_11.16.20.csv")
  coug_master <- read.csv("./Data/dev_telem_coug_ACTIVE_COLLARS_through_12.7.2021.csv")
  #'  Cleaned data 
  md_clean <- read.csv("md_skinny 2021-11-09.csv")
  elk_clean <- read.csv("elk_skinny 2021-11-09.csv")
  wtd_clean <- read.csv("wtd_skinny 2021-11-09.csv")
  coug_clean <- read.csv("coug_clean 2021-12-07.csv") %>% 
    mutate(Floordt = as.POSIXct(Floordt, format = "%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+8")) %>%
    dplyr::select("No", "ID", "CollarID", "Sex", "Latitude", "Longitude", "LMT_DateTime", 
                  "StudyArea", "daytime", "Finaldt", "Floordt")

  #  Function to calculate summary stats on telemetry data
  fix_stats <- function(master, clean) {
    #  Total number of successful fixes
    locations <- nrow(clean)
    #  Number of low accuracy fixes (DOP > 8)
    lowacc <- nrow(filter(master, VEC_DOP > 8))
    #  Percent of low accuracy fixes out of all successful locations
    perc_lowacc <- lowacc/nrow(clean)
    
    newclean <- master %>%
      filter(VEC_FixType != "No fix") %>%
      filter(VEC_FixType != "GPS-2D") %>%
      filter(VEC_Height < 2000 & VEC_Height > 0) %>%
      filter(VEC_DOP <= 8)
    
    #'  Calculate location error in meters for each observation
    #'  Using equation published in Lewis et al. 2007, "Effects of habitat on GPS 
    #'  collar performance: Using data screening to reduce location error"
    #'  y = 1.15 + 4.07*x, where x = DOP value for 3D fixes
    y_all <- y_1thru8 <- c()
    for(i in 1:nrow(master)) {
      y_all[i] = 1.15 + (4.07*master$VEC_DOP[i])
    }
    for(i in 1:nrow(newclean)) {
      y_1thru8[i] = 1.15 + (4.07*newclean$VEC_DOP[i])
    }
    
    #  Mean & SE of fix accuracy
    y_all_mean <- mean(y_all, rm.na = T); y_1thru8_mean = mean(y_1thru8, rm.na = T)
    y_all_sd <- sd(y_all); y_1thru8_sd = sd(y_1thru8)
    y_all_se <- y_all_sd/(sqrt(length(y_all))); y_1thru8_se = y_1thru8_sd/(sqrt(length(y_1thru8)))

    loc_error <- c(y_all_mean, y_all_se, y_1thru8_mean, y_1thru8_se)

    #  Combine all summary stats
    stats <- c(locations, lowacc, perc_lowacc, loc_error)
    stats <- round(stats, 4)

    return(stats)

    
  }
  
  #  Run species location data through summary stats function
  md_stats <- fix_stats(md_master, md_clean)
  elk_stats <- fix_stats(elk_master, elk_clean)
  wtd_stats <- fix_stats(wtd_master, wtd_clean)
  coug_stats <- fix_stats(coug_master, coug_clean)
  
  telem_stats <- as.data.frame(rbind(md_stats, elk_stats, wtd_stats, coug_stats))
  colnames(telem_stats) <- c("Successful_Locs", "LowAccuracy_Locs", "Perc_LowAccuracy", 
                             "Mean_Accuracy_ALL_locs", "SE_Accuracy_ALL_locs", 
                             "Mean_Accuracy_DOP8", "SE_Accuracy_DOP8")
  print(telem_stats)
  
  