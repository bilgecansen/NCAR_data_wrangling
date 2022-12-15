
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

extract_env <- function(env, area, poly, area_sum, first_year, last_year, name) {
  
  env_raw <- raster::extract(env*area, 
                             poly, 
                             #weights = T,
                             fun = sum, 
                             na.rm = T,
                             sp = T)
  
  years <- first_year:last_year
  nyears <- length(years)
  nsites <- length(empe_sites$site_id)
  
  for (i in 1:(nyears*12)) {
    env_raw[[2+i]] <- env_raw[[2+i]]/area_sum
  }
  
  seasons <- rep(rep(years, each = 12), times = nsites)
  months <- rep(rep(1:12, times = nyears), times = nsites)
  poly_id <- c("site_id", "new_n")
  
  data_env <- as.data.frame(env_raw) %>%
    dplyr::select(num_range("layer.", 1:(12*nyears)), poly_id) %>%
    pivot_longer(cols = contains("layer."), names_to = "drop") %>%
    dplyr::select(-drop) %>%
    add_column(season = seasons) %>%
    add_column(month = months) %>%
    pivot_wider(names_from = month, 
                values_from = value, 
                names_prefix = paste(name, ".", sep = ""))
  
  return(data_env)
}

ice_breed <- extract_env(env = r_ice,  
                         area = r_area, 
                         poly = poly_breed, 
                         area_sum = area_sum_breed, 
                         first_year = 1979, 
                         last_year = 2018,
                         name = "aice")

ice_nonbreed <- extract_env(env = r_ice,  
                            area = r_area, 
                            poly = poly_nonbreed, 
                            area_sum = area_sum_nonbreed, 
                            first_year = 1979, 
                            last_year = 2018,
                            name = "aice")

wind_breed <- extract_env(env = sqrt((r_uatm^2) + (r_vatm^2)), 
                          area = r_area, 
                          poly = poly_breed, 
                          area_sum = area_sum_breed, 
                          first_year = 1958, 
                          last_year = 2018, 
                          name = "wind")

wind_nonbreed <- extract_env(env = sqrt((r_uatm^2) + (r_vatm^2)), 
                             area = r_area, 
                             poly = poly_nonbreed, 
                             area_sum = area_sum_nonbreed, 
                             first_year = 1958, 
                             last_year = 2018, 
                             name = "wind")

hmxl_breed <- extract_env(env = r_hmxl, 
                          area = r_area, 
                          poly = poly_breed, 
                          area_sum = area_sum_breed, 
                          first_year = 1958, 
                          last_year = 2018, 
                          name = "hmxl")

hmxl_nonbreed <- extract_env(env = r_hmxl, 
                             area = r_area, 
                             poly = poly_nonbreed, 
                             area_sum = area_sum_nonbreed, 
                             first_year = 1958, 
                             last_year = 2018, 
                             name = "hmxl")

phtc_breed <- extract_env(env = r_phtc, 
                          area = r_area, 
                          poly = poly_breed, 
                          area_sum = area_sum_breed, 
                          first_year = 1958, 
                          last_year = 2018, 
                          name = "phtc")

phtc_nonbreed <- extract_env(env = r_phtc, 
                             area = r_area, 
                             poly = poly_nonbreed, 
                             area_sum = area_sum_nonbreed, 
                             first_year = 1958, 
                             last_year = 2018, 
                             name = "phtc")

zooc_breed <- extract_env(env = r_zooc, 
                          area = r_area, 
                          poly = poly_breed, 
                          area_sum = area_sum_breed, 
                          first_year = 1958, 
                          last_year = 2018, 
                          name = "zooc")

zooc_nonbreed <- extract_env(env = r_zooc, 
                             area = r_area, 
                             poly = poly_nonbreed, 
                             area_sum = area_sum_nonbreed, 
                             first_year = 1958, 
                             last_year = 2018, 
                             name = "zooc")

temp_breed <- extract_env(env = r_temp, 
                          area = r_area, 
                          poly = poly_breed, 
                          area_sum = area_sum_breed, 
                          first_year = 1958, 
                          last_year = 2018, 
                          name = "temp")

temp_nonbreed <- extract_env(env = r_temp, 
                             area = r_area, 
                             poly = poly_nonbreed, 
                             area_sum = area_sum_nonbreed, 
                             first_year = 1958, 
                             last_year = 2018, 
                             name = "temp")

# Combine data
by_col = c("site_id", "new_n", "season")

env_breed <- left_join(ice_breed, wind_breed, by = by_col) %>%
  left_join(hmxl_breed) %>%
  left_join(zooc_breed) %>%
  left_join(temp_breed) %>%
  left_join(phtc_breed)

env_nonbreed <- left_join(ice_nonbreed, wind_nonbreed, by = by_col) %>%
  left_join(hmxl_nonbreed) %>%
  left_join(zooc_nonbreed) %>%
  left_join(temp_nonbreed) %>%
  left_join(phtc_nonbreed)


# Summarize data into seasons ---------------------------------------------

summarise_seasons <- function(dat, var_names) {
  
  sites2 <- empe_sites$site_id
  years <- 1979:2018
  
  cl <- makeCluster(6, types = "SOCK")
  registerDoSNOW(cl)
  
  pb <- txtProgressBar(max = length(sites2), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  pack <- c("tidyverse", "foreach")
  
  res <- foreach(i = 1:length(sites2), .packages = pack, .combine = "rbind", .options.snow = opts) %dopar% {
    
    m <-
      foreach(v = 1:length(var_names), .combine = "cbind") %:%
      foreach(t = 1:length(years), .combine = "rbind") %do% {
        lag0 <- filter(dat, site_id == sites2[i], season == years[t]) %>% 
          select(num_range(paste(var_names[v], ".", sep = ""), 1:12)) %>%
          as.matrix()
        
        res <- c(mean(c(lag0[1:3])), mean(lag0[4:5]), mean(lag0[6:7]), mean(lag0[8:12]))
        names(res) <- paste(var_names[v], c("nonbreed", "arrival", "incubation", "rearing"), sep = "_")
        
        return(res)
      }
    
    m2 <- data.frame(site_id = rep(sites2[i], length(years)), season = years, m)
    rownames(m2) <- NULL
    
    return(m2)
    
  }
  
  stopCluster(cl)
  
  return(res)
}

var_names <- c("aice", "wind", "hmxl", "zooc", "temp", "phtc")

env_season_breed <- summarise_seasons(env_breed, var_names = var_names)
env_season_nonbreed <- summarise_seasons(env_nonbreed, var_names = var_names)

# Combine relevant seasons
data_empe <- select(env_season_breed, -contains("_nonbreed")) %>%
  left_join(select(env_season_nonbreed, site_id, season, contains("_nonbreed")), 
            by = c("site_id", "season"))

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
