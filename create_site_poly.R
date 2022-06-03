
library(rspice)
library(mapppdr)
library(sf)
library(foreach)


# Adelie Penguins ---------------------------------------------------------

adpe_sites <- sites_sf %>%
dplyr::right_join(site_species %>% dplyr::filter(species_id == "ADPE"), by = "site_id")

adpe_poly_50 <- rangeField(sites = adpe_sites, max_range = 50)[[1]] %>%
  st_transform(crs = 4326) %>%
  st_wrap_dateline() %>%
  as_Spatial()

adpe_poly_100 <- rangeField(sites = adpe_sites, max_range = 100)[[1]] %>%
  st_transform(crs = 4326) %>%
  st_wrap_dateline() %>%
  as_Spatial()

adpe_poly_250 <- rangeField(sites = adpe_sites, max_range = 250)[[1]] %>%
  st_transform(crs = 4326) %>%
  st_wrap_dateline() %>%
  as_Spatial()

adpe_poly_500 <- rangeField(sites = adpe_sites, max_range = 500)[[1]] %>%
  st_transform(crs = 4326) %>%
  st_wrap_dateline() %>%
  as_Spatial()
  
adpe_poly_1000 <- rangeField(sites = adpe_sites, max_range = 1000)[[1]] %>%
  st_transform(crs = 4326) %>%
  st_wrap_dateline() %>%
  as_Spatial()

adpe_poly_1500 <- rangeField(sites = adpe_sites, max_range = 1500)[[1]] %>%
  st_transform(crs = 4326) %>%
  st_wrap_dateline() %>%
  as_Spatial()

saveRDS(adpe_poly_50, "data_poly_50km_adpe.rds")
saveRDS(adpe_poly_100, "data_poly_100km_adpe.rds")
saveRDS(adpe_poly_250, "data_poly_250km_adpe.rds")
saveRDS(adpe_poly_500, "data_poly_500km_adpe.rds")
saveRDS(adpe_poly_1000, "data_poly_1000km_adpe.rds")
saveRDS(adpe_poly_1500, "data_poly_1500km_adpe.rds")


# Emperor Penguins --------------------------------------------------------

empe_sites <- read.csv("empe_sites_new.csv")
empe_points <- foreach(i = 1:nrow(empe_sites)) %do% st_point(as.matrix(empe_sites[i,3:2]))
empe_sites2 <- st_sf(site_id = empe_sites[,1], new_n = empe_sites[,4], geometry = st_sfc(empe_points))
st_crs(empe_sites2) <- 4326

empe_poly_100 <- rangeField(sites = empe_sites2, max_range = 100)[[1]] %>%
  st_transform(crs = 4326) %>%
  st_wrap_dateline() %>%
  as_Spatial()

empe_poly_500 <- rangeField(sites = empe_sites2, max_range = 500)[[1]] %>%
  st_transform(crs = 4326) %>%
  st_wrap_dateline() %>%
  as_Spatial()

empe_poly_600 <- rangeField(sites = empe_sites2, max_range = 500)[[1]] %>%
  st_transform(crs = 4326) %>%
  st_wrap_dateline() %>%
  as_Spatial()

empe_poly_750 <- rangeField(sites = empe_sites2, max_range = 750)[[1]] %>%
  st_transform(crs = 4326) %>%
  st_wrap_dateline() %>%
  as_Spatial()

empe_poly_850 <- rangeField(sites = empe_sites2, max_range = 750)[[1]] %>%
  st_transform(crs = 4326) %>%
  st_wrap_dateline() %>%
  as_Spatial()

empe_poly_1000 <- rangeField(sites = empe_sites2, max_range = 1000)[[1]] %>%
  st_transform(crs = 4326) %>%
  st_wrap_dateline() %>%
  as_Spatial()

empe_poly_1500 <- rangeField(sites = empe_sites2, max_range = 1500)[[1]] %>%
  st_transform(crs = 4326) %>%
  st_wrap_dateline() %>%
  as_Spatial()

empe_poly_2000 <- rangeField(sites = empe_sites2, max_range = 2000)[[1]] %>%
  st_transform(crs = 4326) %>%
  st_wrap_dateline() %>%
  as_Spatial()

saveRDS(empe_poly_100, "data_poly_100km_empe.rds")
saveRDS(empe_poly_500, "data_poly_500km_empe.rds")
saveRDS(empe_poly_600, "data_poly_600km_empe.rds")
saveRDS(empe_poly_750, "data_poly_750km_empe.rds")
saveRDS(empe_poly_850, "data_poly_850km_empe.rds")
saveRDS(empe_poly_1000, "data_poly_1000km_empe.rds")
saveRDS(empe_poly_1500, "data_poly_1500km_empe.rds")
saveRDS(empe_poly_2000, "data_poly_2000km_empe.rds")

