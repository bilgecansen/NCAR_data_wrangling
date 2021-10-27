
library(ncdf4)

# Main directories where NetCDF files reside
netcdf_fore_dir <- "/gpfs/scratch/bsen/NCAR_coupled_forecast_NetCDF"
netcdf_hist_dir <- "/gpfs/scratch/bsen/NCAR_coupled_historic_NetCDF"

# Ensemble member directories for NetCDF files, they should be the same for both netcdf directories
member_names <- list.files(netcdf_fore_dir)

m_hist <- matrix(nrow = length(member_names), ncol = 12)
for (i in 1:length(member_names)) {
  netcdf_names <- list.files(paste(netcdf_hist_dir, member_names[i], sep = "/"))
  
  for (h in 1:length(netcdf_names)) {
    ncin <- ncdf4::nc_open(netcdf_names[h])
    m_hist[i,h] <- ncin$var$time_bound$varsize[2] 
  }
}

rownames(m_hist) <- member_names

m_fore <- matrix(nrow = length(member_names), ncol = 12)
for (i in 1:length(member_names)) {
  netcdf_names <- list.files(paste(netcdf_fore_dir, member_names[i], sep = "/"))
  
  for (h in 1:length(netcdf_names)) {
    ncin <- ncdf4::nc_open(netcdf_names[h])
    m_fore[i,h] <- ncin$var$time_bound$varsize[2] 
  }
}

rownames(m_fore) <- member_names

res <- list(m_hist, m_fore)
saveRDS(res, "ts_length.rds")



