
# Functions for data wrangling

# Netcdf to raster transformation -----------------------------------------

# Raster transformation for data with atmospheric forcing
transform_to_raster <- function(netcdf_dir, var_name, member_name = NULL, raster_dir, long_name, lat_name, 
                                depth, integral, depth_count = NA, time_steps, skip = T) {
  
  if (is.null(member_name)) {
    
    netcdf_names <- list.files(netcdf_dir)
    
    idx <- str_which(netcdf_names, var_name)
    netcdf <- netcdf_names[idx]
    
    full_name <- str_split(netcdf, ".nc", 2, simplify = T)[1]
    file_name <- paste(raster_dir, full_name, sep = "/") %>%
      paste(., "grd", sep = ".")
    
    if (skip == T) {
      if (file.exists(file_name)) {
        print(paste(file_name, "already exists", sep = " "))
        return(NULL)
      }
    } 
    
    ncin <- ncdf4::nc_open(paste(netcdf_dir, netcdf, sep = "/"))
    
  } else {
    
    netcdf_names <- list.files(paste(netcdf_dir, member_name, sep = "/"))
    
    idx <- str_which(netcdf_names, var_name)
    netcdf <- netcdf_names[idx]
    
    full_name <- str_split(netcdf, ".nc", 2, simplify = T)[1]
    folder_name <- paste(raster_dir, member_name, sep = "/")
    if (!member_name %in% list.files(raster_dir)) dir.create(folder_name)
    file_name <- paste(folder_name, full_name, sep = "/") %>%
      paste(., "grd", sep = ".")
    
    if (skip == T) {
      if (file.exists(file_name)) {
        print(paste(file_name, "already exists", sep = " "))
        return(NULL)
      }
    }
    
    ncin <- ncdf4::nc_open(paste(netcdf_dir, member_name, netcdf, sep = "/"))
    
  }

  lon <- c(ncdf4::ncvar_get(ncin, long_name, count = c(-1,60), verbose = FALSE))
  # This conversion of lon ensures MAPPPED sites and rasters overlap
  lon <- sapply(lon, function(x) ifelse(x>180, (x - 360), x))
  
  lat <- c(ncdf4::ncvar_get(ncin, lat_name, count = c(-1, 60), verbose = FALSE))
  
  if (is.null(member_name)) {
    count1 <- c(-1, 60, depth_count, -1)
    count2 <- c(-1, 60, -1)
  
  } else {
    count1 <- c(-1, 60, depth_count, -1, -1)
    count2 <- c(-1, 60, -1, -1)
  }
  
  if (depth) {
    
    z <- ncdf4::ncvar_get(ncin, var_name, count = count1, verbose = FALSE)
    
    if (integral) {
      
      z <- z*10
      z <- apply(z, c(1,2,4), sum)
      
    } 
    
  } else z <- ncdf4::ncvar_get(ncin, var_name, count = count2, verbose = FALSE)
  
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
  
  writeRaster(r_brick, filename = file_name, overwrite = T)
  
}


# Extract values from rasters ---------------------------------------------

extract_env <- function(raster_dir, raster_name, member_name = NULL, data_poly, first_year, last_year, tw) {
  
  nyears <- last_year - first_year + 1
  years <- first_year:last_year
  
  if (is.null(member_name)) {
    
    var_name <- str_split(raster_name, "h.", 2, simplify = T)[2] %>%
      str_split(., "\\.", 2, simplify = T) %>%
      .[1]
    
    r_brick <- brick(paste(raster_dir, raster_name, sep = "/"))
  
  } else {
    
    if (str_detect(member_name, "_")) {
      member_name2 <- str_split(member_name2, "_")[[1]][1]
      x <- paste("h.cmip6.", member_name2, ".", sep = "")
    } else {
      x <- paste("h.cmip6.", member_name, ".", sep = "") 
    }
    
    var_name <- str_split(raster_name, x, 2, simplify = T)[2] %>%
      str_split(., "\\.", 2, simplify = T) %>%
      .[1]
    
    r_brick <- raster::brick(paste(raster_dir, member_name, raster_name, sep = "/"))
  
  }
  
  # Assign NAs to lats under civil twilight 
  for (i in 1:12) {
    
    z <- coordinates(r_brick[[i]])
    idx <- which(z[,2]<tw[i])
    
    if (length(idx) == 0) next
    
    r_brick[[i]][idx] <- NA
    
    for (h in 1:(nyears-1)) {
      
      z <- coordinates(r_brick[[i+h*12]])
      idx <- which(z[,2]<tw[i])
      
      r_brick[[i+h*12]][idx] <- NA
      
    }
  }
  
  env_raw <- raster::extract(r_brick, 
                             data_poly, 
                             weights = T, 
                             fun = mean, 
                             na.rm = T,
                             sp = T)
  
  
  nsites <- as.data.frame(env_raw)$site_id %>%
    unique() %>%
    length()
  
  # When transformed to data frame env_raw has columns z.1 to z.732 (or some other value)
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


# Summarize data into seasons ---------------------------------------------

summarize_env <- function(env_dat, cores) {
  
  var_names <- str_split(colnames(env_dat)[5:ncol(env_dat)], "[.]") %>% 
    map_chr(function(x) x[1]) %>% 
    unique()
  
  sites <- unique(env_dat$site_id)
  seasons <- unique(env_dat$season)
  
  cl <- makeCluster(cores, types = "SOCK")
  registerDoSNOW(cl)
  
  pb <- txtProgressBar(max = length(sites), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  pack <- c("tidyverse", "foreach")
  
  dat <- foreach(i = 1:length(sites), .packages = pack, .combine = "rbind", .options.snow = opts) %dopar% {
    
    m <- 
      foreach(v = 1:length(var_names), .combine = "cbind") %:%
        foreach(t = 2:length(seasons), .combine = "rbind") %do% {
          
          lag0a <- filter(env_dat, site_id == sites[i], season == seasons[t-1]) %>% 
            select(num_range(paste(var_names[v], ".", sep = ""), 12)) %>%
            as.matrix()
          
          lag0b <- filter(env_dat, site_id == sites[i], season == seasons[t]) %>% 
            select(num_range(paste(var_names[v], ".", sep = ""), 1:11)) %>%
            as.matrix()
          
          # Winter is April to September, Summer is December to February
          res <- c(mean(c(lag0a, lag0b[1:2])), lag0b[3], mean(lag0b[4:9]), lag0b[10:11])
          names(res) <- paste(var_names[v], c("summer", "winter", "march", "oct", "nov"), sep = "_")
          
          return(res)
        
      }
    
    m2 <- data.frame(site_id = rep(sites[i], length(seasons)-1), season = seasons[-1], m)
    rownames(m2) <- NULL
    
    return(m2)
    
  }
  
  stopCluster(cl)
  
  return(dat)
  
}


# Add lags to summarized data ---------------------------------------------

add_lags <- function(env_dat, cores) {
  
  sites <- unique(env_dat$site_id)
  seasons <- unique(env_dat$season)
  
  cl <- makeCluster(cores, types = "SOCK")
  registerDoSNOW(cl)
  
  pb <- txtProgressBar(max = length(sites), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  pack <- c("tidyverse", "foreach")
  
  dat_lag <- foreach(i = 1:length(sites), .packages = pack, .combine = "rbind", .options.snow = opts) %dopar% {
    
    m <- 
      foreach(t = 7:length(seasons), .combine = "rbind") %do% {
        
        lag3 <- filter(env_dat, site_id == sites[i], season == seasons[t-2]) %>% 
          select(-site_id, -season, -contains("oct"), -contains("nov")) %>%
          as.matrix()
        
        lag4 <- filter(env_dat, site_id == sites[i], season == seasons[t-3]) %>% 
          select(-site_id, -season, -contains("oct"), -contains("nov")) %>%
          as.matrix()

        lag5 <- filter(env_dat, site_id == sites[i], season == seasons[t-4]) %>% 
          select(-site_id, -season, -contains("oct"), -contains("nov")) %>%
          as.matrix()

        lag6 <- filter(env_dat, site_id == sites[i], season == seasons[t-5]) %>% 
          select(-site_id, -season, -contains("oct"), -contains("nov")) %>%
          as.matrix()

        res <- c(lag3, lag4, lag5, lag6)
        
        names(res) <- c(paste(colnames(lag3), "lag3", sep = "_"), 
                        paste(colnames(lag4), "lag4", sep = "_"),
                        paste(colnames(lag5), "lag5", sep = "_"),
                        paste(colnames(lag6), "lag6", sep = "_"))
        
        return(res)
        
      }
    
    m2 <- data.frame(site_id = rep(sites[i], length(seasons)-6), season = seasons[-c(1:6)], m)
    rownames(m2) <- NULL
    
    return(m2)
    
  }
  
  stopCluster(cl)
  
  dat <- left_join(dat_lag, env_dat, by = c("site_id", "season"))
  
  return(dat)
  
}
