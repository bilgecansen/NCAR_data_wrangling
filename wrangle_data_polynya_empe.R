
library(mapppdr)
library(sf)
library(foreach)
library(doSNOW)
library(raster)
library(tidyverse)

source("functions.R")

empe_sites <- read.csv("empe_sites_new.csv")
empe_coord <- rgdal::project(as.matrix(empe_sites[,3:2]), "+init=EPSG:3031")
empe_points <- foreach(i = 1:nrow(empe_coord)) %do% st_point(empe_coord[i,1:2])
empe_sites2 <- st_sf(site_id = empe_sites[,1], new_n = empe_sites[,4], geometry = st_sfc(empe_points))
st_crs(empe_sites2) <- 3031

poly_breed <- readRDS("data_poly_500km_empe.rds")
poly_nonbreed <- readRDS("data_poly_750km_empe.rds")


# Transform netcdf to raster ----------------------------------------------

# Polynya
ncin <- ncdf4::nc_open("SSMI.CDR.85%thresh.polynya2_sh.197901-202012.nc")

z <- ncdf4::ncvar_get(ncin, "polynyas", start =  c(1,7,1), verbose = FALSE)

lon_1d <- ncdf4::ncvar_get(ncin, "tlon1d")
# This conversion of lon ensures MAPPPED sites and rasters overlap
lon_1d <- sapply(lon_1d, function(x) ifelse(x>180, (x - 360), x))

lat_1d <- c(ncdf4::ncvar_get(ncin, "tlat1d", start = 7, verbose = FALSE))

lon_mat <- foreach(i = 1:length(lat_1d), .combine = "cbind") %do% lon_1d
lat_mat <- foreach(i = 1:length(lon_1d), .combine = "rbind") %do% lat_1d

grid_coord <- cbind(c(lon_mat), c(lat_mat))
idx <- which(!is.na(grid_coord[,1]))
grid_coord <- grid_coord[idx,]

r <- list()
a <- 1
for (h in 1:(42*12)) {
  
  pixels <- SpatialPixelsDataFrame(points = grid_coord,
                                   data = data.frame(z = c(z[,,h])[idx]),
                                   tolerance = 0.0001)
  
  r[[a]] <- raster(pixels[,'z'])
  a <- a + 1
  
}

r_pnya <- do.call(brick, r)

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


# Exract polynya area -----------------------------------------------------

# Functions
get_area <- function(env, area, poly, first_year, last_year, var_name, season_names, avg = F) {
  
  env_raw <- raster::extract(env*area, 
                             poly, 
                             fun = sum, 
                             na.rm = T,
                             sp = T)
  
  years <- first_year:last_year
  nyears <- length(years)
  nsites <- length(empe_sites$site_id)
  
  if (avg == T) {
    for (i in 1:(nyears*12)) {
      env_raw[[2+i]] <- env_raw[[2+i]]/area_sum
    }
  }
  
  years_dat <- rep(rep(years, each = 12), times = nsites)
  months_dat <- rep(rep(1:12, times = nyears), times = nsites)
  
  # Seasons are:
  # nonbreed (Jan to March)
  # laying(April and May)
  # incubation(June and July)
  # rearing (August to December)
  
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

# Polynya area
pnya_breed <- get_area(env = r_pnya == 1 | r_pnya == 2,  
                       area = r_area, 
                       poly = poly_breed, 
                       first_year = 1979,
                       last_year = 2020,
                       var_name = "apnya",
                       season_names = breed)

pnya_nonbreed <- get_area(env = r_pnya == 1 | r_pnya == 2,  
                         area = r_area, 
                         poly = poly_nonbreed,
                         first_year = 1979, 
                         last_year = 2020,
                         var_name = "apnya",
                         season_names = "nonbreed")


# Extract food web variables in polynyas ----------------------------------

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

get_sum <- function(env, area, poly, area_sum, first_year, last_year, var_name, season_names) {
  
  env_raw <- raster::extract(env/area, 
                             poly, 
                             fun = sum, 
                             na.rm = T,
                             sp = T)
  
  years <- first_year:last_year
  nyears <- length(years)
  nsites <- length(empe_sites$site_id)
  
  for (i in 1:(nyears*12)) {
    env_raw[[2+i]] <- env_raw[[2+i]]*area_sum
  }
  
  years_dat <- rep(rep(years, each = 12), times = nsites)
  months_dat <- rep(rep(1:12, times = nyears), times = nsites)
  
  # Seasons are:
  # nonbreed (Jan to March)
  # laying(April and May)
  # incubation(June and July)
  # rearing (August to December)
  
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

file_name <- "NCAR_forced_NetCDF/g.e21.GOMIPECOIAF_JRA.TL319_g17.CMIP6-omip2.001.pop.h.photoC_TOT_zint_100m.030601-036612.nc"
r_phtc <- transform_to_raster(file_name = file_name,
                              lat_name = "TLAT",
                              lon_name = "TLONG",
                              var_name = "photoC_TOT_zint_100m",
                              nyears = length(1958:2018))

file_name <- "NCAR_forced_NetCDF/g.e21.GOMIPECOIAF_JRA.TL319_g17.CMIP6-omip2.001.pop.h.zooC.030601-036612.nc"
r_zooc <- transform_to_raster(file_name = file_name,
                              lat_name = "TLAT",
                              lon_name = "TLONG",
                              var_name = "zooC",
                              nyears = length(1958:2018),
                              depth = T, 
                              integral = T, 
                              depth_count = 10)

file_name <- "NCAR_forced_NetCDF/g.e21.GOMIPECOIAF_JRA.TL319_g17.CMIP6-omip2.001.pop.h.HMXL.030601-036612.nc"
r_hmxl <- transform_to_raster(file_name = file_name,
                              lat_name = "TLAT",
                              lon_name = "TLONG",
                              var_name = "HMXL",
                              nyears = length(1958:2018))

# Subset to same years (1979:2018)
r_zooc2 <- subset(r_zooc, 253:732)
r_phtc2 <- subset(r_phtc, 253:732)
r_pnya2 <- subset(r_pnya, 1:480)
r_hmxl2 <- subset(r_hmxl, 253:732)

phtc_breed_pnya <- get_sum(env = r_phtc2*(r_pnya2 == 1 | r_pnya2 == 2), 
                      area = r_area, 
                      poly = poly_breed, 
                      area_sum = area_sum_breed, 
                      first_year = 1979, 
                      last_year = 2018, 
                      var_name = "phtc_sum_pnya",
                      season_names = breed)

phtc_nonbreed_pnya <- get_sum(env = r_phtc2*(r_pnya2 == 1 | r_pnya2 == 2), 
                         area = r_area, 
                         poly = poly_nonbreed, 
                         area_sum = area_sum_nonbreed, 
                         first_year = 1979, 
                         last_year = 2018, 
                         var_name = "phtc_sum_pnya",
                         season_names = "nonbreed")

phtc_breed_poly <- get_sum(env = r_phtc2, 
                           area = r_area, 
                           poly = poly_breed, 
                           area_sum = area_sum_breed, 
                           first_year = 1979, 
                           last_year = 2018, 
                           var_name = "phtc_sum_poly",
                           season_names = breed)

phtc_nonbreed_poly <- get_sum(env = r_phtc2, 
                              area = r_area, 
                              poly = poly_nonbreed, 
                              area_sum = area_sum_nonbreed, 
                              first_year = 1979, 
                              last_year = 2018, 
                              var_name = "phtc_sum_poly",
                              season_names = "nonbreed")

zooc_breed_pnya <- get_sum(env = r_zooc2*(r_pnya2 == 1 | r_pnya2 == 2), 
                      area = r_area, 
                      poly = poly_breed, 
                      area_sum = area_sum_breed, 
                      first_year = 1979, 
                      last_year = 2018, 
                      var_name = "zooc_sum_pnya",
                      season_names = breed)

zooc_nonbreed_pnya <- get_sum(env = r_zooc2*(r_pnya2 == 1 | r_pnya2 == 2), 
                         area = r_area, 
                         poly = poly_nonbreed, 
                         area_sum = area_sum_nonbreed, 
                         first_year = 1979, 
                         last_year = 2018, 
                         var_name = "zooc_sum_pnya",
                         season_names = "nonbreed")

zooc_breed_poly <- get_sum(env = r_zooc2, 
                           area = r_area, 
                           poly = poly_breed, 
                           area_sum = area_sum_breed, 
                           first_year = 1979, 
                           last_year = 2018, 
                           var_name = "zooc_sum_poly",
                           season_names = breed)

zooc_nonbreed_poly <- get_sum(env = r_zooc2, 
                              area = r_area, 
                              poly = poly_nonbreed, 
                              area_sum = area_sum_nonbreed, 
                              first_year = 1979, 
                              last_year = 2018, 
                              var_name = "zooc_sum_poly",
                              season_names = "nonbreed")

hmxl_breed_pnya <- get_area(env = r_hmxl2*(r_pnya2 == 1 | r_pnya2 == 2), 
                            area = r_area, 
                            poly = poly_breed, 
                            first_year = 1979, 
                            last_year = 2018, 
                            var_name = "hmxl_pnya",
                            season_names = breed)

hmxl_nonbreed_pnya <- get_area(env = r_hmxl2*(r_pnya2 == 1 | r_pnya2 == 2), 
                               area = r_area, 
                               poly = poly_nonbreed, 
                               first_year = 1979, 
                               last_year = 2018, 
                               var_name = "hmxl_pnya",
                               season_names = "nonbreed")

data_apnya <- left_join(pnya_nonbreed, pnya_breed, by = c("site_id", "year")) %>%
  left_join(phtc_breed_pnya, by = c("site_id", "year")) %>%
  left_join(phtc_nonbreed_pnya, by = c("site_id", "year")) %>%
  left_join(zooc_breed_pnya, by = c("site_id", "year")) %>%
  left_join(zooc_nonbreed_pnya, by = c("site_id", "year")) %>%
  left_join(phtc_breed_poly, by = c("site_id", "year")) %>%
  left_join(phtc_nonbreed_poly, by = c("site_id", "year")) %>%
  left_join(zooc_breed_poly, by = c("site_id", "year")) %>%
  left_join(zooc_nonbreed_poly, by = c("site_id", "year")) %>%
  left_join(hmxl_breed_pnya, by = c("site_id", "year")) %>%
  left_join(hmxl_nonbreed_pnya, by = c("site_id", "year"))

saveRDS(data_apnya, "data_apnya.rds")
