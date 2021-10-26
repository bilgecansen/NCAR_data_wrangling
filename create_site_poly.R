
library(rspice)
library(mapppdr)
library(sf)

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
