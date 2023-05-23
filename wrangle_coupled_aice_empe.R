
library(mapppdr)
library(sf)
library(foreach)
library(doSNOW)
library(raster)
library(tidyverse)

empe_sites <- read.csv("empe_sites_new.csv")
poly_nonbreed <- readRDS("data_poly_750km_empe.rds")
poly_breed <- readRDS("data_poly_500km_empe.rds")


# Transform netcdf to raster ----------------------------------------------

transform_to_raster <- function(file_name, time_steps, 
                                netcdf_dir = "NCAR_aice_coupled/NetCDF", 
                                raster_dir = "NCAR_aice_coupled/Raster") {
  
  full_name <- str_split(file_name, ".nc", 2, simplify = T)[1]
  ncin <- ncdf4::nc_open(paste(netcdf_dir, file_name, sep = "/"))
  
  lon <- c(ncdf4::ncvar_get(ncin, "TLON", count = c(-1,70), verbose = FALSE))
  # This conversion of lon ensures MAPPPED sites and rasters overlap
  lon <- sapply(lon, function(x) ifelse(x>180, (x - 360), x))
  
  lat <- c(ncdf4::ncvar_get(ncin, "TLAT", count = c(-1, 70), verbose = FALSE))
  
  z <- ncdf4::ncvar_get(ncin, "aice", count = c(-1, 70, -1, -1), verbose = FALSE)
  z[which(z<0)] <- NA
  
  grid_coord <- cbind(lon, lat)
  idx <- which(!is.na(grid_coord[,1]))
  grid_coord <- grid_coord[idx,]
  
  r <- list()
  a <- 1
  for (h in 1:time_steps) {
    
    pixels <- SpatialPixelsDataFrame(points = grid_coord,
                                     data = data.frame(z = c(z[,,h])[idx]),
                                     tolerance = 0.0001)
    
    r[[a]] <- raster(pixels[,'z'])
    a <- a + 1
    
  }
  r_ice <- do.call(brick, r)
  
  writeRaster(r_ice, filename = paste(raster_dir, paste(full_name, "grd", sep = "."), sep = "/"), overwrite = T)
  
}

forecast_netcdf <- list.files("NCAR_aice_coupled/netcdf")
forecast_netcdf <- forecast_netcdf[str_which(forecast_netcdf, "ssp370")]

historic_netcdf <- list.files("NCAR_aice_coupled/netcdf")
historic_netcdf <- historic_netcdf[str_which(historic_netcdf, "historical")]

for (i in 1:50) transform_to_raster(forecast_netcdf[i], 1032)
for (i in 1:50) transform_to_raster(historic_netcdf[i], 1380)

# Area of the cells
ncin <- ncdf4::nc_open("NCAR_aice_coupled/NetCDF/CESM2LE.ice.historical1900to2014.cice.h.cmip6.1001.r1i1001p1f1.aice.nc")

lon <- c(ncdf4::ncvar_get(ncin, "TLON", count = c(-1,70), verbose = FALSE))
lon <- sapply(lon, function(x) ifelse(x>180, (x - 360), x))
lat <- c(ncdf4::ncvar_get(ncin, "TLAT", count = c(-1, 70), verbose = FALSE))

grid_coord <- cbind(lon, lat)
idx <- which(!is.na(grid_coord[,1]))
grid_coord <- grid_coord[idx,]

area <- c(ncdf4::ncvar_get(ncin, "tarea", count = c(-1, 70), verbose = FALSE))/1000000
pixels <- SpatialPixelsDataFrame(points = grid_coord,
                                 data = data.frame(z = area[idx]),
                                 tolerance = 0.0001)
r_area <- raster(pixels[,'z'])

area_sum <- raster::extract(r_area, 
                            poly_nonbreed, 
                            weights = F, 
                            fun = sum, 
                            na.rm = T,
                            sp = F)

# Extract aice values at colonies -----------------------------------------

extract_aice <- function(raster_dir = "NCAR_aice_coupled/Raster/", raster_names, min_year, max_year, poly) {
    
  r_ice <- brick(paste(raster_dir, raster_names, sep = "/"))
  
  ice_raw <- raster::extract(r_ice*r_area, 
                             poly, 
                             weights = F,
                             fun = sum, 
                             na.rm = T,
                             sp = T)
  
  for (i in 1:(dim(ice_raw)[2]-2)) {
    ice_raw[[2+i]] <- ice_raw[[2+i]]/area_sum
  }
  
  nyears <- max_year - min_year + 1
  years <- min_year:max_year
  nsites <- length(empe_sites$site_id)
  
  seasons <- rep(rep(years, each = 12), times = nsites)
  months <- rep(rep(1:12, times = nyears), times = nsites)
  poly_id <- c("site_id", "new_n")
  
  as.data.frame(ice_raw) %>%
    dplyr::select(num_range("layer.", 1:(12*nyears)), poly_id) %>%
    pivot_longer(cols = contains("layer."), names_to = "drop") %>%
    dplyr::select(-drop) %>%
    add_column(season = seasons) %>%
    add_column(month = months) %>%
    pivot_wider(names_from = month, 
                values_from = value, 
                names_prefix = paste("aice", ".", sep = ""))
}

raster_names_historic_aice <- list.files("NCAR_aice_coupled/Raster/")
raster_names_historic_aice <- raster_names_historic_aice[str_which(raster_names_historic_aice, "historic")]
raster_names_historic_aice <- raster_names_historic_aice[str_which(raster_names_historic_aice, ".grd")]

raster_names_forecast_aice <- list.files("NCAR_aice_coupled/Raster/")
raster_names_forecast_aice <- raster_names_forecast_aice[str_which(raster_names_forecast_aice, "ssp370")]
raster_names_forecast_aice <- raster_names_forecast_aice[str_which(raster_names_forecast_aice, ".grd")]
  
aice_nonbreed_historic <- foreach(i = 1:50) %do% {
  extract_aice(raster_names = raster_names_historic_aice[i],
               min_year = 1900,
               max_year = 2014,
               poly = poly_nonbreed)
} 
  
aice_nonbreed_forecast <- foreach(i = 1:50) %do% { 
  extract_aice(raster_names = raster_names_forecast_aice[i],
               min_year = 2015,
               max_year = 2100,
               poly = poly_nonbreed)
}  

aice_nonbreed <- foreach(i = 1:50) %do% { 
  bind_rows(aice_nonbreed_historic[[i]], aice_nonbreed_forecast[[i]]) %>%
  arrange(new_n, season)
}

aice_breed_historic <- foreach(i = 1:50) %do% {
  extract_aice(raster_names = raster_names_historic_aice[i],
               min_year = 1900,
               max_year = 2014,
               poly = poly_breed)
} 

aice_breed_forecast <- foreach(i = 1:50) %do% { 
  extract_aice(raster_names = raster_names_forecast_aice[i],
               min_year = 2015,
               max_year = 2100,
               poly = poly_breed)
}  

aice_breed <- foreach(i = 1:50) %do% { 
  bind_rows(aice_breed_historic[[i]], aice_breed_forecast[[i]]) %>%
    arrange(new_n, season)
}

summarise_seasons <- function(dat) {
  
  sites2 <- empe_sites$site_id
  years <- 1900:2100
  
  res <- foreach(i = 1:length(sites2), .combine = "rbind") %do% {
    
    m <-
      foreach(t = 1:length(years), .combine = "rbind") %do% {
        lag0 <- filter(dat, site_id == sites2[i], season == years[t]) %>% 
          dplyr::select(num_range(paste("aice", ".", sep = ""), 1:12)) %>%
          as.matrix()
        
        res <- c(mean(c(lag0[1:3])), mean(lag0[4:5]), mean(lag0[6:7]), mean(lag0[8:12]))
        names(res) <- paste("aice", c("nonbreed", "laying", "incubation", "rearing"), sep = "_")
        
        return(res)
      }
    
    m2 <- data.frame(site_id = rep(sites2[i], length(years)), season = years, m)
    rownames(m2) <- NULL
    
    return(m2)
    
  }
  
  return(res)
}

pb <- txtProgressBar(max = 50, style = 3)
data_aice_coupled <- foreach(i = 1:50) %do% {
  z1 <- summarise_seasons(aice_nonbreed[[i]])
  z2 <- summarise_seasons(aice_breed[[i]])

  setTxtProgressBar(pb, i)
      
  left_join(dplyr::select(z1, site_id, season, aice_nonbreed), 
            dplyr::select(z2, -aice_nonbreed), 
            by = c("site_id", "season"))
}  

saveRDS(data_aice_coupled, file = "data_aice_coupled.rds")

# Write csv files
for (i in 1:50) {
  if(!i %in% list.files("empe_SIC/coupled")) dir.create(paste("empe_SIC/coupled", i, sep = "/"))
  
  for (h in 1:66) {
    z <- filter(data_aice_coupled[[i]], site_id == empe_sites$site_id[h]) %>%
      dplyr::select(-site_id)
    write.csv(z, file = paste("empe_SIC/coupled", i, paste("col", h, ".csv", sep = ""), sep = "/"), row.names = F)
  }
}

