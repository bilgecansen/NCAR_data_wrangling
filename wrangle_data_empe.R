
library(mapppdr)
library(sf)
library(foreach)
library(doSNOW)
library(raster)
library(tidyverse)

source("functions.R")

empe_sites <- read.csv("empe_sites_new.csv")


# Create site polygons ----------------------------------------------------

draw_poly <- function(sites, dlon, dlat) {
  
  y1 <- sites$latitude
  y2 <- sites$latitude + dlat
  
  x1 <- sites$longitude - dlon/2
  x2 <- sites$longitude + dlon/2
  
  p1 <- cbind(x1, y1)
  p2 <- cbind(x2, y1)
  p3 <- cbind(x2, y2)
  p4 <- cbind(x1, y2)
  
  site_poly <- foreach(i = 1:nrow(sites)) %do% {
    
    st_polygon(list(rbind(p1[i,], p2[i,], p3[i,], p4[i,], p1[i,])))
    
  }
  
  return(st_sfc(site_poly))
  
}

row.names(empe_sites) <- paste("ID", row.names(empe_sites), sep = "")

poly_breed <- draw_poly(empe_sites, dlat = 20, dlon = 10) %>%
  as_Spatial() %>%
  SpatialPolygonsDataFrame(empe_sites[,-(2:3)])

poly_nonbreed <- draw_poly(empe_sites, dlat = 20, dlon = 20)  %>%
  as_Spatial() %>%
  SpatialPolygonsDataFrame(empe_sites[,-(2:3)])


# Transform netcdf to raster ----------------------------------------------

ncin <- ncdf4::nc_open("ssmi_monthly_data_gx1v5_1979-2018.nc")

lon <- c(ncdf4::ncvar_get(ncin, "Lon", count = c(-1,70), verbose = FALSE))
# This conversion of lon ensures MAPPPED sites and rasters overlap
lon <- sapply(lon, function(x) ifelse(x>180, (x - 360), x))

lat <- c(ncdf4::ncvar_get(ncin, "Lat", count = c(-1, 70), verbose = FALSE))

z <- ncdf4::ncvar_get(ncin, "Nasa_Team", count = c(-1, 70, -1), verbose = FALSE)
z[which(z<0)] <- NA
z <- z/100

grid_coord <- cbind(lon, lat)
idx <- which(!is.na(grid_coord[,1]))
grid_coord <- grid_coord[idx,]

r <- list()
a <- 1
for (h in 1:480) {
  
  pixels <- SpatialPixelsDataFrame(points = grid_coord,
                                   data = data.frame(z = c(z[,,h])[idx]),
                                   tolerance = 0.0001)
  
  r[[a]] <- raster(pixels[,'z'])
  a <- a + 1
  
}
r_ice <- do.call(brick, r)

# Grid area raster
ncin2 <- ncdf4::nc_open("NCAR_forced_NetCDF/g.e21.GOMIPECOIAF_JRA.TL319_g17.CMIP6-omip2.001.cice.h.aice.030601-036612.nc")
area <- c(ncdf4::ncvar_get(ncin2, "tarea", count = c(-1, 70), verbose = FALSE))/1000000
pixels <- SpatialPixelsDataFrame(points = grid_coord,
                                 data = data.frame(z = area[idx]),
                                 tolerance = 0.0001)
r_area <- raster(pixels[,'z'])


# Extract env values at sites ---------------------------------------------

raster_dir <- "NCAR_forced_raster"
raster_names <- list.files(raster_dir)
raster_names <- raster_names[str_which(raster_names, "grd")][c(2,3,7:9)]

cl <- makeCluster(6, types = "SOCK")
registerDoSNOW(cl)

pack <- c("raster", "dplyr", "magrittr", "stringr", "tibble", "tidyr" )

env_breed <- 
  foreach (i=1:length(raster_names), .packages = pack) %dopar%
  extract_env(raster_dir = raster_dir,
              raster_name = raster_names[i],
              data_poly = poly_breed,
              first_year = 1979,
              last_year = 2018,
              tw = NA,
              poly_id = c("site_id", "new_n"))

by_col <- c("site_id", "new_n", "season")
env_breed <- reduce(env_breed, function(x, y) full_join(x, y, by = by_col))

env_nonbreed <- 
  foreach (i=1:length(raster_names), .packages = pack) %dopar%
  extract_env(raster_dir = raster_dir,
              raster_name = raster_names[i],
              data_poly = poly_nonbreed,
              first_year = 1979,
              last_year = 2018,
              tw = NA,
              poly_id = c("site_id", "new_n"))

env_nonbreed <- reduce(env_nonbreed, function(x, y) full_join(x, y, by = by_col))

stopCluster(cl)

# Extract sea ice
ice_raw_breed <- raster::extract(r_ice*r_area, 
                                 poly_breed, 
                                 weights = F,
                                 fun = sum, 
                                 na.rm = T,
                                 sp = T)

ice_raw_nonbreed <- raster::extract(r_ice*r_area, 
                                    poly_nonbreed, 
                                    weights = F, 
                                    fun = sum, 
                                    na.rm = T,
                                    sp = T)

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

for (i in 1:480) {
  ice_raw_breed[[2+i]] <- ice_raw_breed[[2+i]]/area_sum_breed
  ice_raw_nonbreed[[2+i]] <- ice_raw_nonbreed[[2+i]]/area_sum_nonbreed
}

nyears <- 2018 - 1979 + 1
years <- 1979:2018
nsites <- length(empe_sites$site_id)

seasons <- rep(rep(years, each = 12), times = nsites)
months <- rep(rep(1:12, times = nyears), times = nsites)
poly_id <- c("site_id", "new_n")

ice_breed <- as.data.frame(ice_raw_breed) %>%
  dplyr::select(num_range("layer.", 1:(12*nyears)), poly_id) %>%
  pivot_longer(cols = contains("layer."), names_to = "drop") %>%
  dplyr::select(-drop) %>%
  add_column(season = seasons) %>%
  add_column(month = months) %>%
  pivot_wider(names_from = month, 
              values_from = value, 
              names_prefix = paste("aice", ".", sep = ""))

ice_nonbreed <- as.data.frame(ice_raw_nonbreed) %>%
  dplyr::select(num_range("layer.", 1:(12*nyears)), poly_id) %>%
  pivot_longer(cols = contains("layer."), names_to = "drop") %>%
  dplyr::select(-drop) %>%
  add_column(season = seasons) %>%
  add_column(month = months) %>%
  pivot_wider(names_from = month, 
              values_from = value, 
              names_prefix = paste("aice", ".", sep = ""))

env_breed <- left_join(env_breed, ice_breed, by = by_col)
env_nonbreed <- left_join(env_nonbreed, ice_nonbreed, by = by_col)


# Summarize data into seasons ---------------------------------------------

sites2 <- empe_sites$site_id
var_names <- c("aice", "uatm", "vatm", "HMXL")

cl <- makeCluster(6, types = "SOCK")
registerDoSNOW(cl)

pb <- txtProgressBar(max = length(sites), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
pack <- c("tidyverse", "foreach")

env_season_breed <- foreach(i = 1:length(sites2), .packages = pack, .combine = "rbind", .options.snow = opts) %dopar% {
  
  m <-
    foreach(v = 1:length(var_names), .combine = "cbind") %:%
    foreach(t = 1:length(years), .combine = "rbind") %do% {
      lag0 <- filter(env_breed, site_id == sites2[i], season == years[t]) %>% 
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

env_season_nonbreed <- foreach(i = 1:length(sites2), .packages = pack, .combine = "rbind", .options.snow = opts) %dopar% {
  
  m <-
    foreach(v = 1:length(var_names), .combine = "cbind") %:%
    foreach(t = 1:length(years), .combine = "rbind") %do% {
      lag0 <- filter(env_nonbreed, site_id == sites2[i], season == years[t]) %>% 
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

# Combine relevant seasons
data_empe <- select(env_season_breed, -contains("_nonbreed")) %>%
  left_join(select(env_season_nonbreed, site_id, season, contains("_nonbreed")), 
            by = c("site_id", "season"))

saveRDS(data_empe, file = "data_env_empe.rds")

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


# Extract CESM sea ice separately -----------------------------------------

raster_dir <- "NCAR_forced_raster"
raster_names <- list.files(raster_dir)
raster_names <- raster_names[str_which(raster_names, "grd")][1]

cl <- makeCluster(1, types = "SOCK")
registerDoSNOW(cl)

pack <- c("raster", "dplyr", "magrittr", "stringr", "tibble", "tidyr" )

ice_cesm_breed <- 
  foreach (i=1:length(raster_names), .packages = pack) %dopar%
  extract_env(raster_dir = raster_dir,
              raster_name = raster_names[i],
              data_poly = poly_breed,
              first_year = 1958,
              last_year = 2018,
              tw = NA,
              poly_id = c("site_id", "new_n"))

ice_cesm_breed <- ice_cesm_breed[[1]]

ice_cesm_nonbreed <- 
  foreach (i=1:length(raster_names), .packages = pack) %dopar%
  extract_env(raster_dir = raster_dir,
              raster_name = raster_names[i],
              data_poly = poly_nonbreed,
              first_year = 1958,
              last_year = 2018,
              tw = NA,
              poly_id = c("site_id", "new_n"))

ice_cesm_nonbreed <- ice_cesm_nonbreed[[1]]

stopCluster(cl)

var_names <- c("aice")
years2 <- 1958:2018

cl <- makeCluster(6, types = "SOCK")
registerDoSNOW(cl)

pb <- txtProgressBar(max = length(sites), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
pack <- c("tidyverse", "foreach")

ice_cesm_season_breed <- foreach(i = 1:length(sites2), 
                                 .packages = pack, .combine = "rbind", .options.snow = opts) %dopar% {
  
  m <-
    foreach(v = 1:length(var_names), .combine = "cbind") %:%
    foreach(t = 1:length(years2), .combine = "rbind") %do% {
      lag0 <- filter(ice_cesm_breed, site_id == sites2[i], season == years2[t]) %>% 
        select(num_range(paste(var_names[v], ".", sep = ""), 1:12)) %>%
        as.matrix()
      
      res <- c(mean(c(lag0[1:3])), mean(lag0[4:5]), mean(lag0[6:7]), mean(lag0[8:12]))
      names(res) <- paste(var_names[v], c("nonbreed", "arrival", "incubation", "rearing"), sep = "_")
      
      return(res)
    }
  
  m2 <- data.frame(site_id = rep(sites2[i], length(years2)), season = years2, m)
  rownames(m2) <- NULL
  
  return(m2)
  
}

ice_cesm_season_nonbreed <- foreach(i = 1:length(sites2), 
                                    .packages = pack, .combine = "rbind", .options.snow = opts) %dopar% {
  
  m <-
    foreach(v = 1:length(var_names), .combine = "cbind") %:%
    foreach(t = 1:length(years2), .combine = "rbind") %do% {
      lag0 <- filter(ice_cesm_nonbreed, site_id == sites2[i], season == years2[t]) %>% 
        select(num_range(paste(var_names[v], ".", sep = ""), 1:12)) %>%
        as.matrix()
      
      res <- c(mean(c(lag0[1:3])), mean(lag0[4:5]), mean(lag0[6:7]), mean(lag0[8:12]))
      names(res) <- paste(var_names[v], c("nonbreed", "arrival", "incubation", "rearing"), sep = "_")
      
      return(res)
    }
  
  m2 <- data.frame(site_id = rep(sites2[i], length(years2)), season = years2, m)
  rownames(m2) <- NULL
  
  return(m2)
  
}

stopCluster(cl)

# Combine relevant seasons
data_ice_cesm_empe <- select(ice_cesm_season_breed, -contains("_nonbreed")) %>%
  left_join(select(ice_cesm_season_nonbreed, site_id, season, contains("_nonbreed")), 
            by = c("site_id", "season"))

saveRDS(data_ice_cesm_empe, file = "data_aice_cesm_empe.rds")

