
library(ncdf4)
library(raster)
library(geodist)
library(sp)
library(magrittr)
library(stringr)
library(purrr)
library(foreach)
library(doSNOW)

source("functions.R")

# Directories setup for seawulf -------------------------------------------

# Main directories where NetCDF files reside
netcdf_fore_dir <- "/gpfs/scratch/bsen/NCAR_coupled_forecast_NetCDF"
netcdf_hist_dir <- "/gpfs/scratch/bsen/NCAR_coupled_historic_NetCDF"

# The main directories where raster files will be written to
raster_fore_dir <- "/gpfs/scratch/bsen/NCAR_coupled_forecast_raster"
raster_hist_dir <- "/gpfs/scratch/bsen/NCAR_coupled_historic_raster"

# Ensemble member directories for NetCDF files, they should be the same for both netcdf directories
member_names <- list.files(netcdf_fore_dir)


# Netcdf to raster conversion ---------------------------------------------

# Variables with depth layers to integrate
for (i in 1:length(member_names)) {
  
  transform_to_raster(netcdf_dir = netcdf_fore_dir,
                      member_name = member_names[i],
                      var_name = "zooC",
                      raster_dir = raster_fore_dir, 
                      long_name = "TLONG", 
                      lat_name = "TLAT", 
                      depth = T, 
                      integral = T, 
                      depth_count = 10, 
                      time_steps = 1:456)
  
  transform_to_raster(netcdf_dir = netcdf_hist_dir,
                      member_name = member_names[i],
                      var_name = "zooC",
                      raster_dir = raster_hist_dir, 
                      long_name = "TLONG", 
                      lat_name = "TLAT", 
                      depth = T, 
                      integral = T, 
                      depth_count = 10, 
                      time_steps = 1:684)
  
}

# Variables where only first depth is used
vars_first <- c("TEMP", "aicen")

for (i in 1:length(member_names)) {

  transform_to_raster(netcdf_dir = netcdf_fore_dir,
                      member_name = member_names[i],
                      var_name = "TEMP",
                      raster_dir = raster_fore_dir, 
                      long_name = "TLONG", 
                      lat_name = "TLAT", 
                      depth = T, 
                      integral = F, 
                      depth_count = 1, 
                      time_steps = 1:456)
  
  transform_to_raster(netcdf_dir = netcdf_hist_dir,
                      member_name = member_names[i],
                      var_name = "TEMP",
                      raster_dir = raster_hist_dir, 
                      long_name = "TLONG", 
                      lat_name = "TLAT", 
                      depth = T, 
                      integral = F, 
                      depth_count = 1, 
                      time_steps = 1:684)
  
}

for (i in 1:length(member_names)) {
  
  transform_to_raster(netcdf_dir = netcdf_fore_dir,
                      member_name = member_names[i],
                      var_name = "aicen",
                      raster_dir = raster_fore_dir, 
                      long_name = "TLON", 
                      lat_name = "TLAT", 
                      depth = T, 
                      integral = F, 
                      depth_count = 1, 
                      time_steps = 1:456)
  
  transform_to_raster(netcdf_dir = netcdf_hist_dir,
                      member_name = member_names[i],
                      var_name = "aicen",
                      raster_dir = raster_hist_dir, 
                      long_name = "TLON", 
                      lat_name = "TLAT", 
                      depth = T, 
                      integral = F, 
                      depth_count = 1, 
                      time_steps = 1:684)
  
}

# Variables without depth layers
vars_nd1 <- c("HMXL", "photoC_TOT_zint_100m")

for (i in 1:length(member_names)) {
  for(h in 1:length(vars_nd1)) {
  
    transform_to_raster(netcdf_dir = netcdf_fore_dir,
                        member_name = member_names[i],
                        var_name = vars_nd1[h],
                        raster_dir = raster_fore_dir, 
                        long_name = "TLONG", 
                        lat_name = "TLAT", 
                        depth = F, 
                        integral = F, 
                        time_steps = 1:456)
    
    transform_to_raster(netcdf_dir = netcdf_hist_dir,
                        member_name = member_names[i],
                        var_name = vars_nd1[h],
                        raster_dir = raster_hist_dir, 
                        long_name = "TLONG", 
                        lat_name = "TLAT", 
                        depth = F, 
                        integral = F,  
                        time_steps = 1:684)
    
  }
}

vars_nd2 <- c("aice", "shear", "divu", "vatm", "uatm", "rain", "ardg")

for (i in 1:length(member_names)) {
  for(h in 1:length(vars_nd2)) {
    
    transform_to_raster(netcdf_dir = netcdf_fore_dir,
                        member_name = member_names[i],
                        var_name = vars_nd2[h],
                        raster_dir = raster_fore_dir, 
                        long_name = "TLON", 
                        lat_name = "TLAT", 
                        depth = F, 
                        integral = F, 
                        time_steps = 1:456)
    
    transform_to_raster(netcdf_dir = netcdf_hist_dir,
                        member_name = member_names[i],
                        var_name = vars_nd2[h],
                        raster_dir = raster_hist_dir, 
                        long_name = "TLON", 
                        lat_name = "TLAT", 
                        depth = F, 
                        integral = F,  
                        time_steps = 1:684)
  }
}


# Extract env data from rasters -------------------------------------------

raster_fore_names <- foreach(i = 1:length(member_names)) %do% {
  z <- list.files(paste(raster_fore_dir, member_names[i], sep = "/"))
  z[which(str_detect(z, ".grd"))]
}

raster_hist_names <- foreach(i = 1:length(member_names)) %do% {
  z <- list.files(paste(raster_hist_dir, member_names[i], sep = "/"))
  z[which(str_detect(z, ".grd"))]
}

data_poly_1500 <- readRDS("data_poly_1500km_adpe.rds")
data_poly_100 <- readRDS("data_poly_100km_adpe.rds")

# Extract env values from rasters -----------------------------------------

# Calculate the civil twilight for each month
day <- c()
for (i in 1:365) {
  day[i] <- (as.Date(i, origin = "2020-12-31") %>% as.character %>% strsplit("-"))[[1]][3] %>% as.numeric()
}
d <- which(day == 15)
D <- -23.45*cos(((360/365)*(d + 10))*pi/180)
tw <- D-96

cl <- makeCluster(6, types = "SOCK")
registerDoSNOW(cl)

pack <- c("raster", "dplyr", "magrittr", "stringr", "tibble", "tidyr" )

nras <- length(raster_fore_names[[1]])

env_list_fore_1500 <- 
  foreach (i = 1:length(member_names), .packages = pack) %:%
    foreach(h = 1:nras) %dopar% {
      extract_env(raster_dir = raster_fore_dir,
                  raster_name = raster_fore_names[[i]][h],
                  member_name = member_names[i],
                  data_poly = data_poly_1500,
                  first_year = 2015,
                  last_year = 2052,
                  tw = tw)
    }

env_list_fore_100 <- 
  foreach (i = 1:length(member_names), .packages = pack) %:%
    foreach(h = 1:nras) %dopar% {
      extract_env(raster_dir = raster_fore_dir,
                  raster_name = raster_fore_names[[i]][h],
                  member_name = member_names[i],
                  data_poly = data_poly_100,
                  first_year = 2015,
                  last_year = 2052,
                  tw = tw)
    }
    
env_list_hist_1500 <- 
  foreach (i = 1:length(member_names), .packages = pack) %:%
    foreach(h = 1:nras) %dopar% {
      extract_env(raster_dir = raster_hist_dir,
                  raster_name = raster_hist_names[[i]][h],
                  member_name = member_names[i],
                  data_poly = data_poly_1500,
                  first_year = 1958,
                  last_year = 2014,
                  tw = tw)
    }

env_list_hist_100 <- 
  foreach (i = 1:length(member_names), .packages = pack) %:%
    foreach(h = 1:nras) %dopar% {
      extract_env(raster_dir = raster_hist_dir,
                  raster_name = raster_hist_names[[i]][h],
                  member_name = member_names[i],
                  data_poly = data_poly_100,
                  first_year = 1958,
                  last_year = 2014,
                  tw = tw)
      }

stopCluster(cl)

by_col <- c("site_id", "species_id", "season", "ccamlr_id")

env_dat_fore_1500 <- map(env_list_fore_1500, function(x) reduce(x, function(x, y) full_join(x, y, by = by_col))) 
env_dat_fore_100 <- map(env_list_fore_100, function(x) reduce(x, function(x, y) full_join(x, y, by = by_col))) 

env_dat_hist_1500 <- map(env_list_hist_1500, function(x) reduce(x, function(x, y) full_join(x, y, by = by_col))) 
env_dat_hist_100 <- map(env_list_hist_100, function(x) reduce(x, function(x, y) full_join(x, y, by = by_col)))

env_dat_1500 <- foreach(i = 1:length(env_dat_fore_1500)) %do% {
  rbind(env_dat_hist_1500[[i]], env_dat_fore_1500[[i]]) %>%
    arrange(site_id, season)
}
  
env_dat_100 <- foreach(i = 1:length(env_dat_fore_100)) %do% {
  rbind(env_dat_hist_100[[i]], env_dat_fore_100[[i]]) %>%
    arrange(site_id, season)
}

saveRDS(env_dat_1500, file = "data_coupled_env_1500km_adpe.rds")
saveRDS(env_dat_100, file = "data_coupled_env_100km_adpe.rds")

# Summarize env data for forecasts ----------------------------------------

data_rf_1500 <- foreach(i = 1:length(env_dat_1500)) %do% summarize_env(env_dat_1500[[i]], cores = 6)
data_rf_100 <- foreach(i = 1:length(env_dat_100)) %do% summarize_env(env_dat_100[[i]], cores = 6)

# Standardize data with the mean and sd of 1958-2017 period
mean_1500 <- foreach(i = 1:length(data_rf_1500)) %do% {
  filter(data_rf_1500[[i]], season <2018) %>%
    group_by(site_id) %>% 
    summarise(across(-season, mean))
}

sd_1500 <- foreach(i = 1:length(data_rf_1500)) %do% {
  filter(data_rf_1500[[i]], season <2018) %>%
    group_by(site_id) %>% 
    summarise(across(-season, sd))
}

mean_100 <- foreach(i = 1:length(data_rf_100)) %do% {
  filter(data_rf_100[[i]], season <2018) %>%
    group_by(site_id) %>% 
    summarise(across(-season, mean))
}

sd_100 <- foreach(i = 1:length(data_rf_100)) %do% {
  filter(data_rf_100[[i]], season <2018) %>%
    group_by(site_id) %>% 
    summarise(across(-season, sd))
}

sites <- unique(data_rf_1500[[1]]$site_id)

data_std_rf_1500 <- 
  foreach(i = length(data_rf_1500)) %:%
    foreach(h = 1:length(sites), .combine = "rbind") %do% {
      z <- filter(data_rf_1500[[i]], site_id == sites[h])
      z1 <- z[,1:2]
      z2 <- as.matrix(z[,3:ncol(z)])
      
      m <- filter(mean_1500[[i]], site_id == sites[h]) %>%
        select(-site_id) %>%
        as.matrix()
      
      sdev <- filter(sd_1500[[i]], site_id == sites[h]) %>%
        select(-site_id) %>%
        as.matrix()
      
      z3 <- t(apply(z2, 1, function(x) (x - m)/sdev))
      colnames(z3) <- colnames(m)
      
      data.frame(z1, z3)
  }

data_std_rf_100 <- 
  foreach(i = length(data_rf_100)) %:%
    foreach(h = 1:length(sites), .combine = "rbind") %do% {
      z <- filter(data_rf_100[[i]], site_id == sites[h])
      z1 <- z[,1:2]
      z2 <- as.matrix(z[,3:ncol(z)])
      
      m <- filter(mean_100[[i]], site_id == sites[h]) %>%
        select(-site_id) %>%
        as.matrix()
      
      sdev <- filter(sd_100[[i]], site_id == sites[h]) %>%
        select(-site_id) %>%
        as.matrix()
      
      z3 <- t(apply(z2, 1, function(x) (x - m)/sdev))
      colnames(z3) <- colnames(m)
      
      data.frame(z1, z3)
  }

# Add lags
data_rf_1500_lag <- foreach(i = 1:length(data_rf_1500)) %do% add_lags(data_rf_1500[[i]], cores = 6)
data_rf_100_lag <- foreach(i = 1:length(data_rf_100)) %do% add_lags(data_rf_100[[i]], cores = 6)

data_std_rf_1500_lag <- foreach(i = 1:length(data_std_rf_1500)) %do% add_lags(data_std_rf_1500[[i]], cores = 6)
data_std_rf_100_lag <- foreach(i = 1:length(data_std_rf_100)) %do% add_lags(data_std_rf_100[[i]], cores = 6)

saveRDS(data_rf_1500_lag, "data_coupled_rf_1500km_adpe.rds")
saveRDS(data_rf_100_lag, "data_coupled_rf_100km_adpe.rds")

saveRDS(data_std_rf_1500_lag, "data_coupled_std_rf_1500km_adpe.rds")
saveRDS(data_std_rf_100_lag, "data_coupled_std_rf_100km_adpe.rds")

