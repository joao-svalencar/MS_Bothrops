# predictors --------------------------------------------------------------
# raw
folders <- dir(here::here("data","raw","predictors"))

#original
preds.pres <- terra::rast(here::here("data","processed","predictors","current.tif"))

preds.fut11 <- terra::rast(here::here("data","processed","predictors","future_2041-2070_ssp370.tif"))
preds.fut12 <- terra::rast(here::here("data","processed","predictors","future_2041-2070_ssp585.tif"))

preds.fut21 <- terra::rast(here::here("data","processed","predictors","future_2071-2100_ssp370.tif"))
preds.fut22 <- terra::rast(here::here("data","processed","predictors","future_2071-2100_ssp585.tif"))

# processed 10 min
preds.pres_10m <- terra::rast(here::here("data","processed","predictors","current_10m.tif"))

preds.fut11_10m <- terra::rast(here::here("data","processed","predictors","2041-2070_ssp370_10m.tif"))
preds.fut12_10m <- terra::rast(here::here("data","processed","predictors","2041-2070_ssp585_10m.tif"))

preds.fut21_10m <- terra::rast(here::here("data","processed","predictors","2071-2100_ssp370_10m.tif"))
preds.fut22_10m <- terra::rast(here::here("data","processed","predictors","2071-2100_ssp585_10m.tif"))

# study area --------------------------------------------------------------
sa <- terra::vect(here::here("data","area","south_america.shp"))
topo <- terra::rast(here::here("data","area","southamerica_topo.tif"))

# species occurrences -----------------------------------------------------
dfs <- read.csv(here::here("data","raw","occurrences", "group_alternatus.csv"), h=T)

# outputs for processing --------------------------------------------------
files <- list.files(here::here("outputs","rasters","vars_10m","binary"), pattern='tif$', full.names = TRUE)
files_or <- list.files(here::here("outputs","rasters","original vars","binary"), pattern = ".tif", full.names = TRUE)

# analyses ----------------------------------------------------------------
centroids <- read.csv(here::here("outputs","tables","centroids_bin.csv"), h=T)
areas <- read.csv(here::here("outputs","tables","area_scenarios.csv"), h=T)
altitude <- read.csv(here::here("outputs","tables","elevation_scenarios.csv"), h=T)
