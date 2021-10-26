
# Transform netcdf files to raster bricks

library(ncdf4)
library(raster)
library(geodist)
library(sp)
library(magrittr)
library(stringr)
library(purrr)
library(foreach)
library(tidyverse)
library(doSNOW)

# Directories setup for seawulf -------------------------------------------

netcdf_dir <- "NCAR_historic"
member_names <- list.files(netcdf_dir)

member_dir <- foreach(i = 1:length(member_names), .combine = "c") %do% 
  paste(netcdf_dir, member_names[i], sep = "/")

netcdf_names <- foreach(i = 1:length(member_names)) %do% {
  z <- list.files(member_dir[i])
  z[which(str_detect(z, ".nc"))]
} 


# Netcdf to raster conversion ---------------------------------------------

transform_to_raster <- function(netcdf_names, var_name, member_dir, long_name, lat_name, 
                                depth, integral, count, time_steps) {
  
  member_name <- str_split(member_dir, "/")[[1]][2]
  x <- paste("h.cmip6", member_name, sep = ".")
  
  idx <- str_which(netcdf_names, var_name)
  netcdf <- netcdf_names[idx]
  
  #var_name <- str_split(netcdf_name, x, 2, simplify = T)[2] %>%
  #str_split(., "\\.", 2, simplify = T) %>%
  #.[1]
  
  ncin <- ncdf4::nc_open(paste(member_dir, netcdf, sep = "/"))
  
  lon <- c(ncdf4::ncvar_get(ncin, long_name, count = c(-1,60), verbose = FALSE))
  # This conversion of lon ensures MAPPPED sites and rasters overlap
  lon <- sapply(lon, function(x) ifelse(x>180, (x - 360), x))
  
  lat <- c(ncdf4::ncvar_get(ncin, lat_name, count = c(-1, 60), verbose = FALSE))
  
  if (depth) {
    
    z <- ncdf4::ncvar_get(ncin, var_name, count = c(-1, 60, count, -1, -1), verbose = FALSE)
    
    if (integral) {
      
      z <- z*10
      z <- apply(z, c(1,2,4), sum)
      
    } 
    
  } else z <- ncdf4::ncvar_get(ncin, var_name, count = c(-1, 60, -1, -1), verbose = FALSE)
  
  grid_coord <- cbind(lon, lat)
  idx <- which(!is.na(grid_coord[,1]))
  grid_coord <- grid_coord[idx,]
  
  r <- list()
  a <- 1
  for (h in time_steps) {
    
    pixels <- SpatialPixelsDataFrame(points = grid_coord,
                                     data = data.frame(z = c(z[,,h])[idx]),
                                     tolerance = 0.0001)
    
    r[[a]] <- raster(pixels[,'z'])
    a <- a + 1
    
  }
  r_brick <- do.call(brick, r)
  
  full_name <- str_split(netcdf, ".nc", 2, simplify = T)[1]
  file_name <- paste(member_dir, full_name, sep = "/") %>%
    paste(., "grd", sep = ".")
  
  writeRaster(r_brick, filename = file_name, overwrite = T)
  
}

# Variables with depth layers to integrate
for (i in 1:length(member_dir)) {
  
  transform_to_raster(netcdf_names = netcdf_names[[i]], 
                      member_dir = member_dir[i],
                      var_name = "zooC",
                      long_name = "TLONG", 
                      lat_name = "TLAT", 
                      depth = T, 
                      integral = T, 
                      count = 10, 
                      time_steps = 1297:1980)
  
}

# Variables where only first depth is used
vars_first <- c("TEMP", "aicen")

for (i in 1:length(member_dir)) {
  
  transform_to_raster(netcdf_names = netcdf_names[[i]], 
                      member_dir = member_dir[i],
                      var_name = "TEMP",
                      long_name = "TLONG", 
                      lat_name = "TLAT", 
                      depth = T, 
                      integral = F, 
                      count = 1, 
                      time_steps = 1297:1980)
  
}

for (i in 1:length(member_dir)) {
  
  transform_to_raster(netcdf_names = netcdf_names[[i]], 
                      member_dir = member_dir[i],
                      var_name = "aicen",
                      long_name = "TLON", 
                      lat_name = "TLAT", 
                      depth = T, 
                      integral = F, 
                      count = 1, 
                      time_steps = 1297:1980)
  
}

# Variables without depth layers
vars_nd1 <- c("HMXL", "photoC_TOT_zint_100m")

for (i in 1:length(member_dir)) {
  for(h in 1:length(vars_nd1))
    
    transform_to_raster(netcdf_names = netcdf_names[[i]], 
                        member_dir = member_dir[i],
                        var_name = vars_nd1[h],
                        long_name = "TLONG", 
                        lat_name = "TLAT", 
                        depth = F, 
                        integral = F, 
                        time_steps = 1297:1980)
  
}

vars_nd2 <- c("aice", "shear", "divu", "vatm", "uatm", "rain", "ardg")

for (i in 1:length(member_dir)) {
  for(h in 1:length(vars_nd2))
    
    transform_to_raster(netcdf_names = netcdf_names[[i]], 
                        member_dir = member_dir[i],
                        var_name = vars_nd2[h],
                        long_name = "TLON", 
                        lat_name = "TLAT", 
                        depth = F, 
                        integral = F, 
                        time_steps = 1297:1980)
  
}


# Extract data from rasters -----------------------------------------------

raster_names <- foreach(i = 1:length(member_names)) %do% {
  z <- list.files(member_dir[i])
  z[which(str_detect(z, ".grd"))]
} 

data_poly_1500 <- readRDS("data_poly_1500km_adpe.rds")
data_poly_100 <- readRDS("data_poly_100km_adpe.rds")

# Extract env values from rasters -----------------------------------------

extract_env <- function(member_dir, raster_name, member_name, data_poly, first_year, last_year) {
  
  x <- paste("h.cmip6.", member_name, ".", sep = "")
  
  var_name <- str_split(raster_name, x, 2, simplify = T)[2] %>%
    str_split(., "\\.", 2, simplify = T) %>%
    .[1]
  
  r_brick <- raster::brick(paste(member_dir, raster_name, sep = "/"))
  
  env_raw <- raster::extract(r_brick, 
                             data_poly, 
                             weights = T, 
                             fun = mean, 
                             na.rm = T,
                             sp = T)
  
  nyears <- last_year - first_year + 1
  years <- first_year:last_year
  nsites <- as.data.frame(env_raw)$site_id %>%
    unique() %>%
    length()
  
  # When transformed to data frame env_raw has columns names with z.
  # Each z corresponds to a month in year 
  # For example, z.1 to z.12 is 12 months of the first year
  
  seasons <- rep(rep(years, each = 12), times = nsites)
  months <- rep(rep(1:12, times = nyears), times = nsites)
  
  env_dat <- as.data.frame(env_raw) %>%
    dplyr::select(-site_name, -region) %>%
    pivot_longer(cols = contains("z."), names_to = "drop") %>%
    dplyr::select(-drop) %>%
    add_column(season = seasons) %>%
    add_column(month = months) %>%
    pivot_wider(names_from = month, 
                values_from = value, 
                names_prefix = paste(var_name, ".", sep = ""))
  
  return(env_dat)
  
}

cl <- makeCluster(6, types = "SOCK")
registerDoSNOW(cl)

pack <- c("raster", "dplyr", "magrittr", "stringr", "tibble", "tidyr" )

nras <- length(raster_names[[1]])

env_list_1500 <- 
  foreach (i = 1:length(member_dir), .packages = pack) %:%
  foreach(h = 1:nras) %dopar%
  extract_env(member_dir = member_dir[i],
              raster_name = raster_names[[i]][h],
              member_name = member_names[i],
              data_poly = data_poly_1500,
              first_year = 1958,
              last_year = 2014)

env_list_100 <- 
  foreach (i = 1:length(member_dir), .packages = pack) %:%
  foreach(h = 1:nras) %dopar% {
    extract_env(member_dir = member_dir[i],
                raster_name = raster_names[[i]][h],
                member_name = member_names[i],
                data_poly = data_poly_100,
                first_year = 1958,
                last_year = 2014)
  }

stopCluster(cl)

by_col <- c("site_id", "species_id", "season", "ccamlr_id")

env_dat_1500 <- map(env_list_1500, function(x) reduce(x, function(x, y) full_join(x, y, by = by_col)))
env_dat_100 <- map(env_list_100, function(x) reduce(x, function(x, y) full_join(x, y, by = by_col)))

# Add data until 2017
env_dat_1500_fore <- readRDS("data_env_s") %>%
  filter(season %in% c(2017, 2016, 2015, 2014, 2013, 2012))

env_dat_100_hist <- readRDS("data_env_50km_adpe.rds") %>%
  filter(season %in% c(2017, 2016, 2015, 2014, 2013, 2012))

env_dat_500 <- map(env_dat_500, function(x) rbind(env_dat_500_hist, x)) %>%
  map(function(x) arrange(x, site_id, season))

env_dat_50 <- map(env_dat_50, function(x) rbind(env_dat_50_hist, x)) %>%
  map(function(x) arrange(x, site_id, season))

# Summarize env data for forecasts ----------------------------------------

sites <- unique(env_dat_500[[1]]$site_id)
seasons_500 <- unique(env_dat_500[[1]]$season)
seasons_50 <- unique(env_dat_50[[1]]$season)

var_names <- str_split(colnames(env_dat_500[[1]])[5:ncol(env_dat_500[[1]])], "[.]") %>% 
  map_chr(function(x) x[1]) %>% 
  unique()

cl <- makeCluster(6, types = "SOCK")
registerDoSNOW(cl)

pb <- txtProgressBar(max = length(sites)*length(member_dir), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
pack <- c("tidyverse", "foreach")

env_fore_500 <- foreach(h = 1:length(env_dat_500)) %:%
  foreach(i = 1:length(sites), .packages = pack, .combine = "rbind", .options.snow = opts) %dopar% {
    
    m <- 
      foreach(v = 1:length(var_names), .combine = "cbind") %:%
      foreach(t = 7:length(seasons_500), .combine = "rbind") %do% {
        
        lag0 <- filter(env_dat_500[[h]], site_id == sites[i], season == seasons_500[t]) %>% 
          select(num_range(paste(var_names[v], ".", sep = ""), 1:12)) %>%
          as.matrix()
        
        # lags in the breeding season
        lag3 <- filter(env_dat_500[[h]], site_id == sites[i], season == seasons_500[t-3]) %>% 
          select(num_range(paste(var_names[v], ".", sep = ""), 1:12)) %>%
          as.matrix()
        
        lag4 <- filter(env_dat_500[[h]], site_id == sites[i], season == seasons_500[t-4]) %>% 
          select(num_range(paste(var_names[v], ".", sep = ""), 1:12)) %>%
          as.matrix()
        
        lag5 <- filter(env_dat_500[[h]], site_id == sites[i], season == seasons_500[t-5]) %>% 
          select(num_range(paste(var_names[v], ".", sep = ""), 1:12)) %>%
          as.matrix()
        
        lag6 <- filter(env_dat_500[[h]], site_id == sites[i], season == seasons_500[t-6]) %>% 
          select(num_range(paste(var_names[v], ".", sep = ""), 1:12)) %>%
          as.matrix()
        
        # Winter is April to September
        res <- c(mean(lag0[4:9]), mean(lag3[4:9]), mean(lag4[4:9]), mean(lag5[4:9]), mean(lag6[4:9]))
        names(res) <- paste(var_names[v], c("winter", "winter_lag3", "winter_lag4", 
                                            "winter_lag5", "winter_lag6"), sep = "_")
        
        return(res)
        
      }
    
    m2 <- data.frame(site_id = rep(sites[i], length(seasons_500)-6), season = seasons_500[-c(1:6)], m)
    rownames(m2) <- NULL
    
    return(m2)
    
  }

env_fore_50 <- foreach(h = 1:length(env_dat_50)) %:%
  foreach(i = 1:length(sites), .packages = pack, .combine = "rbind", .options.snow = opts) %dopar% {
    
    m <- 
      foreach(v = 1:length(var_names), .combine = "cbind") %:%
      foreach(t = 7:length(seasons_50), .combine = "rbind") %do% {
        
        lag0 <- filter(env_dat_50[[h]] , site_id == sites[i], season == seasons_50[t]) %>% 
          select(num_range(paste(var_names[v], ".", sep = ""), 1:12)) %>%
          as.matrix()
        
        # lags in the breeding season
        lag3a <- filter(env_dat_50[[h]] , site_id == sites[i], season == seasons_50[t-2]) %>% 
          select(num_range(paste(var_names[v], ".", sep = ""), 12)) %>%
          as.matrix()
        
        lag3b <- filter(env_dat_50[[h]], site_id == sites[i], season == seasons_50[t-3]) %>% 
          select(num_range(paste(var_names[v], ".", sep = ""), 1:3)) %>%
          as.matrix()
        
        lag4a <- filter(env_dat_50[[h]], site_id == sites[i], season == seasons_50[t-3]) %>% 
          select(num_range(paste(var_names[v], ".", sep = ""), 12)) %>%
          as.matrix()
        
        lag4b <- filter(env_dat_50[[h]], site_id == sites[i], season == seasons_50[t-4]) %>% 
          select(num_range(paste(var_names[v], ".", sep = ""), 1:3)) %>%
          as.matrix()
        
        lag5a <- filter(env_dat_50[[h]], site_id == sites[i], season == seasons_50[t-4]) %>% 
          select(num_range(paste(var_names[v], ".", sep = ""), 12)) %>%
          as.matrix()
        
        lag5b <- filter(env_dat_50[[h]], site_id == sites[i], season == seasons_50[t-5]) %>% 
          select(num_range(paste(var_names[v], ".", sep = ""), 1:3)) %>%
          as.matrix()
        
        lag6a <- filter(env_dat_50[[h]], site_id == sites[i], season == seasons_50[t-5]) %>% 
          select(num_range(paste(var_names[v], ".", sep = ""), 12)) %>%
          as.matrix()
        
        lag6b <- filter(env_dat_50[[h]], site_id == sites[i], season == seasons_50[t-6]) %>% 
          select(num_range(paste(var_names[v], ".", sep = ""), 1:3)) %>%
          as.matrix()
        
        # March is transition, April to September is Winter, 
        # October and November are breeding months 
        res <- c(lag0[10:11], mean(c(lag3a, lag3b[1:2])), 
                 mean(c(lag4a, lag4b[1:2])), mean(c(lag5a, lag5b[1:2])), 
                 mean(c(lag6a, lag6b[1:2])), lag3b[3], lag4b[3], lag5b[3], lag6b[3])
        
        names(res) <- paste(var_names[v], 
                            c("oct", "nov", "summer_lag3", "summer_lag4", 
                              "summer_lag5", "summer_lag6", "march_lag3",  "march_lag4",
                              "march_lag5",  "march_lag6"), sep = "_")
        
        return(res)
        
      }
    
    m2 <- data.frame(site_id = rep(sites[i], length(seasons_50)-6), season = seasons_50[-c(1:6)], m)
    rownames(m2) <- NULL
    
    return(m2)
    
  }

stopCluster(cl)


# Save finalized data -----------------------------------------------------

saveRDS(env_fore_500, "data_fore_500km_adpe.rds")
saveRDS(env_fore_50, "data_fore_50km_adpe.rds")


