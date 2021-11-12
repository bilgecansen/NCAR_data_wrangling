
# Transform netcdf files to raster bricks

library(ncdf4)
library(raster)
library(geodist)
library(sp)
library(magrittr)
library(stringr)
library(purrr)
library(foreach)
library(doSNOW)
library(tidyverse)

netcdf_dir <- "NCAR_forced_NetCDF"
raster_dir <- "NCAR_forced_raster"

source("functions.R")

# Netcdf to raster conversion ---------------------------------------------

# Variables with depth layers to integrate
transform_to_raster(netcdf_dir = netcdf_dir,
                    var_name = "zooC",
                    raster_dir = raster_dir, 
                    long_name = "TLONG", 
                    lat_name = "TLAT", 
                    depth = T, 
                    integral = T, 
                    depth_count = 10, 
                    time_steps = 1:732)

# Variables where only first depth is used
vars_first <- c("TEMP", "aicen")

transform_to_raster(netcdf_dir = netcdf_dir,
                    var_name = "TEMP",
                    raster_dir = raster_dir, 
                    long_name = "TLONG", 
                    lat_name = "TLAT", 
                    depth = T, 
                    integral = F, 
                    depth_count = 1, 
                    time_steps = 1:732)

transform_to_raster(netcdf_dir = netcdf_dir,
                    var_name = "aicen",
                    raster_dir = raster_dir, 
                    long_name = "TLON", 
                    lat_name = "TLAT", 
                    depth = T, 
                    integral = F, 
                    depth_count = 1, 
                    time_steps = 1:732)

# Variables without depth layers
vars_nd1 <- c("HMXL", "photoC_TOT_zint_100m")

for (i in 1:length(vars_nd1)) {
  
  transform_to_raster(netcdf_dir = netcdf_dir,
                      var_name = vars_nd1[i],
                      raster_dir = raster_dir, 
                      long_name = "TLONG", 
                      lat_name = "TLAT", 
                      depth = F, 
                      integral = F, 
                      time_steps = 1:732)
  
}

vars_nd2 <- c("aice", "shear", "divu", "vatm", "uatm", "rain", "ardg")

for (i in 1:length(vars_nd1)) {
  
  transform_to_raster(netcdf_dir = netcdf_dir,
                      var_name = vars_nd2[i],
                      raster_dir = raster_dir, 
                      long_name = "TLON", 
                      lat_name = "TLAT", 
                      depth = F, 
                      integral = F, 
                      time_steps = 1:732)
  
}


# Extract env data from rasters -------------------------------------------

raster_names <- list.files(raster_dir)
raster_names <- raster_names[str_which(raster_names, "grd")]

data_poly_1500 <- readRDS("data_poly_1500km_adpe.rds")
data_poly_500 <- readRDS("data_poly_500km_adpe.rds")
data_poly_100 <- readRDS("data_poly_100km_adpe.rds")

# Calculate the civil twilight for each month
day <- c()
for (i in 1:365) {
  day[i] <- (as.Date(i, origin = "2020-12-31") %>% as.character %>% strsplit("-"))[[1]][3] %>% as.numeric()
}
d <- which(day == 15)
D <- -23.45*cos(((360/365)*(d + 10))*pi/180)
tw <- D-96

# Extract env values
cl <- makeCluster(6, types = "SOCK")
registerDoSNOW(cl)

pack <- c("raster", "dplyr", "magrittr", "stringr", "tibble", "tidyr" )

env_list_1500 <- 
  foreach (i=1:length(raster_names), .packages = pack) %dopar%
    extract_env(raster_dir = raster_dir,
                raster_name = raster_names[i],
                data_poly = data_poly_1500,
                first_year = 1958,
                last_year = 2018,
                tw = tw)

env_list_500 <- 
  foreach (i=1:length(raster_names), .packages = pack) %dopar%
    extract_env(raster_dir = raster_dir,
                raster_name = raster_names[i],
                data_poly = data_poly_500,
                first_year = 1958,
                last_year = 2018,
                tw = tw)

env_list_100 <- 
  foreach (i=1:length(raster_names), .packages = pack) %dopar%
    extract_env(raster_dir = raster_dir,
                raster_name = raster_names[i],
                data_poly = data_poly_100,
                first_year = 1958,
                last_year = 2018,
                tw = tw)

stopCluster(cl)

by_col <- c("site_id", "species_id", "season", "ccamlr_id")
env_dat_1500 <- reduce(env_list_1500, function(x, y) full_join(x, y, by = by_col))
env_dat_500 <- reduce(env_list_500, function(x, y) full_join(x, y, by = by_col))
env_dat_100 <- reduce(env_list_100, function(x, y) full_join(x, y, by = by_col))

saveRDS(env_dat_1500, file = "data_forced_env_1500km_adpe.rds")
saveRDS(env_dat_500, file = "data_forced_env_500km_adpe.rds")
saveRDS(env_dat_100, file = "data_forced_env_100km_adpe.rds")


# Summarize data into seasons ---------------------------------------------

data_rf_1500 <- summarize_env(env_dat_1500, cores = 6)
data_rf_500 <- summarize_env(env_dat_500, cores = 6)
data_rf_100 <- summarize_env(env_dat_100, cores = 6)

data_std_rf_1500 <- group_by(data_rf_1500, site_id) %>%
  filter(season < 2018) %>%
  mutate(across(-season, function(x) (x- mean(x))/sd(x))) %>%
  ungroup()

data_std_rf_500 <- group_by(data_rf_500, site_id) %>%
  filter(season < 2018) %>%
  mutate(across(-season, function(x) (x- mean(x))/sd(x))) %>%
  ungroup()

data_std_rf_100 <- group_by(data_rf_100, site_id) %>%
  filter(season < 2018) %>%
  mutate(across(-season, function(x) (x- mean(x))/sd(x))) %>%
  ungroup()

# Add lags
data_rf_1500_lag <- add_lags(data_rf_1500, cores = 6)
data_rf_500_lag <- add_lags(data_rf_500, cores = 6)
data_rf_100_lag <- add_lags(data_rf_100, cores = 6)

data_std_rf_1500_lag <- add_lags(data_std_rf_1500, cores = 6)
data_std_rf_500_lag <- add_lags(data_std_rf_500, cores = 6)
data_std_rf_100_lag <- add_lags(data_std_rf_100, cores = 6)

saveRDS(data_rf_1500_lag, "data_forced_rf_1500km_adpe.rds")
saveRDS(data_rf_500_lag, "data_forced_rf_500km_adpe.rds")
saveRDS(data_rf_100_lag, "data_forced_rf_100km_adpe.rds")

saveRDS(data_std_rf_1500_lag, "data_forced_std_rf_1500km_adpe.rds")
saveRDS(data_std_rf_500_lag, "data_forced_std_rf_500km_adpe.rds")
saveRDS(data_std_rf_100_lag, "data_forced_std_rf_100km_adpe.rds")

