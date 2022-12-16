
library(mapppdr)
library(sf)
library(foreach)
library(doSNOW)
library(raster)
library(tidyverse)

source("functions.R")

empe_sites <- read.csv("empe_sites_new.csv")

poly_breed <- readRDS("data_poly_500km_empe.rds")
poly_nonbreed <- readRDS("data_poly_850km_empe.rds")


# Transform netcdf to raster ----------------------------------------------

# depth: most variables don't have depth but some include multiple layers of ocean
# depth_count: how many layers of depth to include in the transformation. If the goal is not to 
# integrate across years, this should be a single layer
# Integrate: Integration across multiple years, in which each value in a layer is multiplied with layer length (~ 10 m)
# and all values across layers are then summed.
# Ratio: Mostly for ice concentration which changes between 0 and 100 but transformed to be between 0 and 1.

transform_to_raster <- function(file_name, lon_name, lat_name, var_name, 
                                nyears, depth = F, depth_count, integral = F, ratio = F) {
  
  ncin <- ncdf4::nc_open(file_name)
  
  count2 <- c(-1, 70, -1)
  
  if (depth) {
    
    count1 <- c(-1, 70, depth_count, -1)
    
    z <- ncdf4::ncvar_get(ncin, var_name, count = count1, verbose = FALSE)
    
    if (integral) {
      
      z <- z*10
      z <- apply(z, c(1,2,4), sum)
      
    } 
    
  } else z <- ncdf4::ncvar_get(ncin, var_name, count = count2, verbose = FALSE)
  
  lon <- c(ncdf4::ncvar_get(ncin, lon_name, count = c(-1,70), verbose = FALSE))
  # This conversion of lon ensures MAPPPED sites and rasters overlap
  lon <- sapply(lon, function(x) ifelse(x>180, (x - 360), x))
  
  lat <- c(ncdf4::ncvar_get(ncin, lat_name, count = c(-1, 70), verbose = FALSE))

  
  if (ratio == T) {
    z[which(z<0)] <- NA
    z <- z/100
  } 
  
  grid_coord <- cbind(lon, lat)
  idx <- which(!is.na(grid_coord[,1]))
  grid_coord <- grid_coord[idx,]
  
  r <- list()
  a <- 1
  for (h in 1:(nyears*12)) {
    
    pixels <- SpatialPixelsDataFrame(points = grid_coord,
                                     data = data.frame(z = c(z[,,h])[idx]),
                                     tolerance = 0.0001)
    
    r[[a]] <- raster(pixels[,'z'])
    a <- a + 1
    
  }
  r <- do.call(brick, r)
  
  return(r)
}

r_ice <- transform_to_raster(file_name = "ssmi_monthly_data_gx1v5_1979-2018.nc",
                             lat_name = "Lat",
                             lon_name = "Lon",
                             var_name = "Nasa_Team",
                             nyears = length(1979:2018),
                             ratio = T)

file_name <- "NCAR_forced_NetCDF/g.e21.GOMIPECOIAF_JRA.TL319_g17.CMIP6-omip2.001.cice.h.uatm.030601-036612.nc"
r_uatm <- transform_to_raster(file_name = file_name,
                              lat_name = "TLAT",
                              lon_name = "TLON",
                              var_name = "uatm",
                              nyears = length(1958:2018))

file_name <- "NCAR_forced_NetCDF/g.e21.GOMIPECOIAF_JRA.TL319_g17.CMIP6-omip2.001.cice.h.vatm.030601-036612.nc"
r_vatm <- transform_to_raster(file_name = file_name,
                              lat_name = "TLAT",
                              lon_name = "TLON",
                              var_name = "vatm",
                              nyears = length(1958:2018))

file_name <- "NCAR_forced_NetCDF/g.e21.GOMIPECOIAF_JRA.TL319_g17.CMIP6-omip2.001.pop.h.HMXL.030601-036612.nc"
r_hmxl <- transform_to_raster(file_name = file_name,
                              lat_name = "TLAT",
                              lon_name = "TLONG",
                              var_name = "HMXL",
                              nyears = length(1958:2018))

file_name <- "NCAR_forced_NetCDF/g.e21.GOMIPECOIAF_JRA.TL319_g17.CMIP6-omip2.001.pop.h.photoC_TOT_zint_100m.030601-036612.nc"
r_phtc <- transform_to_raster(file_name = file_name,
                              lat_name = "TLAT",
                              lon_name = "TLONG",
                              var_name = "photoC_TOT_zint_100m",
                              nyears = length(1958:2018))

file_name <- "NCAR_forced_NetCDF/g.e21.GOMIPECOIAF_JRA.TL319_g17.CMIP6-omip2.001.pop.h.TEMP.030601-036612.nc"
r_temp <- transform_to_raster(file_name = file_name,
                              lat_name = "TLAT",
                              lon_name = "TLONG",
                              var_name = "TEMP",
                              nyears = length(1958:2018),
                              depth = T, 
                              integral = F, 
                              depth_count = 1) 

file_name <- "NCAR_forced_NetCDF/g.e21.GOMIPECOIAF_JRA.TL319_g17.CMIP6-omip2.001.pop.h.zooC.030601-036612.nc"
r_zooc <- transform_to_raster(file_name = file_name,
                              lat_name = "TLAT",
                              lon_name = "TLONG",
                              var_name = "zooC",
                              nyears = length(1958:2018),
                              depth = T, 
                              integral = T, 
                              depth_count = 10)

# Grid area raster
# Area raster is used to calculate either total area (or concentration) of a variable within a polygon
# or to take the weighted mean of a variable

# Use aice for getting grid area but all variables have the same grid area
ncin <- ncdf4::nc_open("NCAR_forced_NetCDF/g.e21.GOMIPECOIAF_JRA.TL319_g17.CMIP6-omip2.001.cice.h.aice.030601-036612.nc")
area <- c(ncdf4::ncvar_get(ncin, "tarea", count = c(-1, 70), verbose = FALSE))/1000000

lon <- c(ncdf4::ncvar_get(ncin, "TLON", count = c(-1,70), verbose = FALSE))
# This conversion of lon ensures MAPPPED sites and rasters overlap
lon <- sapply(lon, function(x) ifelse(x>180, (x - 360), x))

lat <- c(ncdf4::ncvar_get(ncin, "TLAT", count = c(-1, 70), verbose = FALSE))

grid_coord <- cbind(lon, lat)
idx <- which(!is.na(grid_coord[,1]))
grid_coord <- grid_coord[idx,]

pixels <- SpatialPixelsDataFrame(points = grid_coord,
                                 data = data.frame(z = area[idx]),
                                 tolerance = 0.0001)
r_area <- raster(pixels[,'z'])

area_sum_breed <- raster::extract(r_area, 
                                  poly_breed, 
                                  weights = F, 
                                  fun = sum, 
                                  na.rm = T,
                                  sp = F)

area_sum_nonbreed <- raster::extract(r_area, 
                                     poly_nonbreed, 
                                     weights = F, 
                                     fun = sum, 
                                     na.rm = T,
                                     sp = F)


# Extract env values at sites ---------------------------------------------

# Extract grid values within a polygon and calculate the area weighted average of a variable across
# months within an ecologically relevant season
# Seasons are:
# nonbreed (Jan to March)
# laying(April and May)
# incubation(June and July)
# rearing (August to December)

extract_env <- function(env, area, poly, area_sum, first_year, last_year, var_name, season_names) {
  
  env_raw <- raster::extract(env*area, 
                             poly, 
                             fun = sum, 
                             na.rm = T,
                             sp = T)
  
  years <- first_year:last_year
  nyears <- length(years)
  nsites <- length(empe_sites$site_id)
  
  for (i in 1:(nyears*12)) {
    env_raw[[2+i]] <- env_raw[[2+i]]/area_sum
  }
  
  years_dat <- rep(rep(years, each = 12), times = nsites)
  months_dat <- rep(rep(1:12, times = nyears), times = nsites)
  
  seasons <- c(rep("nonbreed", 3), rep("laying", 2), rep("incubation", 2), rep("rearing", 5))
  seasons_dat <- rep(rep(seasons, times = nyears), times = nsites)
  
  data_env <- as.data.frame(env_raw) %>%
    pivot_longer(cols = contains("layer."), names_to = "drop") %>%
    dplyr::select(-drop) %>%
    add_column(year = years_dat) %>%
    add_column(month = months_dat) %>%
    add_column(season = seasons_dat)
  
  summarise_season <- function(data, season_name, var_name) {
    z <- filter(data, season == season_name) %>%
      group_by(site_id, new_n, year) %>%
      summarise(n = mean(value))
    
    colnames(z)[4] <- paste(var_name, season_name, sep = "_")
    
    return(z)
  }
  
  if (length(season_names) > 1) {
    
    env_season <- foreach(i = 1:length(season_names)) %do% {
      summarise_season(data = data_env, 
                       season_name = season_names[i], 
                       var_name = var_name)
    }
    
    by_col <- c("site_id", "new_n", "year")
    
    data_env2 <- env_season[[1]]
    for (i in 1:(length(season_names)-1)) {
      data_env2 <- left_join(data_env2, env_season[[i+1]], by = by_col)
    }
    
  } else {
    data_env2 <-  summarise_season(data = data_env, 
                                   season_name = season_names, 
                                   var_name = var_name)
  }
  
  data_env2 <- arrange(data_env2, new_n, year) %>%
    ungroup() %>%
    select(-new_n)
  
  return(data_env2)
  
}

breed <- c("laying", "incubation", "rearing")

ice_breed <- extract_env(env = r_ice,  
                         area = r_area, 
                         poly = poly_breed, 
                         area_sum = area_sum_breed, 
                         first_year = 1979, 
                         last_year = 2018,
                         var_name = "aice",
                         season_names = breed)

ice_nonbreed <- extract_env(env = r_ice,  
                            area = r_area, 
                            poly = poly_nonbreed, 
                            area_sum = area_sum_nonbreed, 
                            first_year = 1979, 
                            last_year = 2018,
                            var_name = "aice",
                            season_names = "nonbreed")

wind_breed <- extract_env(env = sqrt((r_uatm^2) + (r_vatm^2)), 
                          area = r_area, 
                          poly = poly_breed, 
                          area_sum = area_sum_breed, 
                          first_year = 1958, 
                          last_year = 2018, 
                          var_name = "wind",
                          season_names = breed)

wind_nonbreed <- extract_env(env = sqrt((r_uatm^2) + (r_vatm^2)), 
                             area = r_area, 
                             poly = poly_nonbreed, 
                             area_sum = area_sum_nonbreed, 
                             first_year = 1958, 
                             last_year = 2018, 
                             var_name = "wind",
                             season_names = "nonbreed")

hmxl_breed <- extract_env(env = r_hmxl, 
                          area = r_area, 
                          poly = poly_breed, 
                          area_sum = area_sum_breed, 
                          first_year = 1958, 
                          last_year = 2018, 
                          var_name = "hmxl",
                          season_names = breed)

hmxl_nonbreed <- extract_env(env = r_hmxl, 
                             area = r_area, 
                             poly = poly_nonbreed, 
                             area_sum = area_sum_nonbreed, 
                             first_year = 1958, 
                             last_year = 2018, 
                             var_name = "hmxl",
                             season_names = "nonbreed")

phtc_breed <- extract_env(env = r_phtc, 
                          area = r_area, 
                          poly = poly_breed, 
                          area_sum = area_sum_breed, 
                          first_year = 1958, 
                          last_year = 2018, 
                          var_name = "phtc",
                          season_names = breed)

phtc_nonbreed <- extract_env(env = r_phtc, 
                             area = r_area, 
                             poly = poly_nonbreed, 
                             area_sum = area_sum_nonbreed, 
                             first_year = 1958, 
                             last_year = 2018, 
                             var_name = "phtc",
                             season_names = "nonbreed")

zooc_breed <- extract_env(env = r_zooc, 
                          area = r_area, 
                          poly = poly_breed, 
                          area_sum = area_sum_breed, 
                          first_year = 1958, 
                          last_year = 2018, 
                          var_name = "zooc",
                          season_names = breed)

zooc_nonbreed <- extract_env(env = r_zooc, 
                             area = r_area, 
                             poly = poly_nonbreed, 
                             area_sum = area_sum_nonbreed, 
                             first_year = 1958, 
                             last_year = 2018, 
                             var_name = "zooc",
                             season_names = "nonbreed")

temp_breed <- extract_env(env = r_temp, 
                          area = r_area, 
                          poly = poly_breed, 
                          area_sum = area_sum_breed, 
                          first_year = 1958, 
                          last_year = 2018, 
                          var_name = "temp",
                          season_names = breed)

temp_nonbreed <- extract_env(env = r_temp, 
                             area = r_area, 
                             poly = poly_nonbreed, 
                             area_sum = area_sum_nonbreed, 
                             first_year = 1958, 
                             last_year = 2018, 
                             var_name = "temp",
                             season_names = "nonbreed")

# Combine data
by_col <- c("site_id", "year")

env_breed <- left_join(ice_breed, wind_breed, by = by_col) %>%
  left_join(hmxl_breed, by = by_col) %>%
  left_join(zooc_breed, by = by_col) %>%
  left_join(temp_breed, by = by_col) %>%
  left_join(phtc_breed, by = by_col)

env_nonbreed <- left_join(ice_nonbreed, wind_nonbreed, by = by_col) %>%
  left_join(hmxl_nonbreed, by = by_col) %>%
  left_join(zooc_nonbreed, by = by_col) %>%
  left_join(temp_nonbreed, by = by_col) %>%
  left_join(phtc_nonbreed, by = by_col)

data_empe <- left_join(env_breed, env_nonbreed, by = by_col)

saveRDS(data_empe, file = "data_env_empe2.rds")

# Create csv for SIC
data_aice <- select(data_empe, site_id, season, aice_nonbreed, aice_arrival, aice_incubation, aice_rearing)
write.csv(data_empe, "data_aice_empe.csv")

# Create csv for each site
if(!"empe_SIC" %in% list.files()) dir.create("empe_SIC")

for (i in 1:nrow(empe_sites)) {
  
  z <- filter(data_aice, site_id == empe_sites$site_id[i]) %>%
    select(-site_id, -season)
  
  write.csv(z, paste("empe_SIC/SIC", "_col", i, ".csv", sep = ""), row.names = F)
  
}
