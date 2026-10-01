# Loading packages --------------------------------------------------------

library(here)
library(spThin)
library(terra)

# Processing predictor variables ------------------------------------------

# Attention: read object 'folders' in R_data.R line 3
# Attention: read object 'sa' in R_data.R line 16

for(i in seq_along(folders)){
  preds.Chelsa <- list.files(here::here("data","raw","predictors", folders[i]), pattern='tif$', full.names = T)
  preds <- terra::rast(preds.Chelsa)
  preds.sa <- terra::crop(preds, sa)
  preds.m <- terra::mask(preds.sa, sa)
  preds.10m <- terra::aggregate(preds.m,fact=20, fun=median)
  
  #terra::writeRaster(preds.m, filename=here::here("data","processed","predictors", paste0(folders[i], ".tif")), overwrite=TRUE)
  terra::writeRaster(preds.10m, filename=here::here("data","processed","predictors", paste0(folders[i], "_10m.tif")), overwrite=TRUE)
}

# Processing species occurrences ------------------------------------------

sds <- split(dfs, dfs$species)

for(t in seq_along(sds)){
  df <- sds[[t]]
  df1 <- df[,c("species", "longitude", "latitude")]
  
  suppressWarnings(
    thin(loc.data = df1, # = data.frame com todos os pontos de ocorrencia
         lat.col = "latitude", long.col = "longitude", # nomes das colunas de coordenadas
         spec.col = "species", # nome da coluna de taxon
         thin.par = 20, reps = 10, # distancia em km de um ponto a outro; reps = numero de replicas de reamostragem
         locs.thinned.list.return = FALSE, 
         write.files = TRUE, 
         max.files = 1, 
         out.dir = here::here("data","processed","occurrences_thin"), # diretorio de saida
         out.base = df1[1,1], # nome do arquivo de saida seguindo de _thin1
         write.log.file = TRUE,
         log.file = here::here("data","processed","occurrences_thin", paste0(df1[1,1],"_log.txt")),
         verbose=FALSE))
}


# Processing model outputs ------------------------------------------------

# extracting and exporting centroids coordinates from binary rasters
#files <- list.files(here::here("outputs","rasters","original vars","binary"), pattern='tif$', full.names = TRUE) # does not work, too much data
files <- list.files(here::here("outputs","rasters","vars_10m","binary"), pattern='tif$', full.names = TRUE)

all_coords <- list()  # lista geral

for (i in seq_along(names(sds))) {
  spp_files <- grep(names(sds)[i], files, value = TRUE)
  btps_bin <- terra::rast(spp_files)
  
  for (j in 1:terra::nlyr(btps_bin)) {
    lyr <- btps_bin[[j]]
    lyr_name <- terra::sources(lyr)
    
    idx <- which(terra::values(lyr) == 1) # presence only
    
    if (length(idx) > 0) {
      coords <- terra::xyFromCell(lyr, idx) # get coordinates
      
      all_coords[[length(all_coords) + 1]] <- data.frame(
        species = names(sds)[i],
        x = coords[, 1],
        y = coords[, 2],
        layer = lyr_name
      )
    }
  }
}

centroids <- do.call(rbind, all_coords)

centroids$layer <- ifelse(grepl("cur.tif$", centroids$layer), "present", centroids$layer)
centroids$layer <- ifelse(grepl("fut11.tif$", centroids$layer), "future_11", centroids$layer)
centroids$layer <- ifelse(grepl("fut12.tif$", centroids$layer), "future_12", centroids$layer)
centroids$layer <- ifelse(grepl("fut21.tif$", centroids$layer), "future_21", centroids$layer)
centroids$layer <- ifelse(grepl("fut22.tif$", centroids$layer), "future_22", centroids$layer)

names(centroids) <- c("species", "longitude", "latitude", "scenario")
centroids$scenario <- factor(centroids$scenario, levels = c("present", "future_11", "future_12", "future_21", "future_22"))

write.csv(centroids, here::here("outputs","tables","centroids_bin.csv"), row.names = FALSE) # object used for modelling

# processing for modelling ------------------------------------------------
centroids$scenario <- factor(centroids$scenario,
                             levels = c(
                               "present",
                               "future_11",
                               "future_12",
                               "future_21",
                               "future_22")
)

centroids$spGroup <- ifelse(
  grepl("alternatus|itapetiningae", centroids$species),
  "Open Areas & Plateau", "Forest & Mountain")

centroids$spGroup <- factor(centroids$spGroup, levels = c("Forest & Mountain", "Open Areas & Plateau"))
centroids$species <- factor(centroids$species, levels = c("Bothrops fonsecai", "Bothrops cotiara", "Bothrops itapetiningae", "Bothrops alternatus"))

#ready to model (D_analyses line 6)

# processing to create figure latitude ~ scenarios ------------------------
# calculating mean latitude
mean_df <- aggregate(
  latitude ~ species + scenario,
  data = centroids,
  FUN = mean
)
names(mean_df)[names(mean_df) == "latitude"] <- "mean_lat"

# calculating standard deviation
sd_df <- aggregate(
  latitude ~ species + scenario,
  data = centroids,
  FUN = sd
)
names(sd_df)[names(sd_df) == "latitude"] <- "sd_lat"


# calculating sample size
n_df <- aggregate(
  latitude ~ species + scenario,
  data = centroids,
  FUN = length
)
names(n_df)[names(n_df) == "latitude"] <- "n_lat"


lat_summary <- merge(mean_df, sd_df,
                     by = c("species","scenario"))
lat_summary <- merge(lat_summary, n_df,
                     by = c("species","scenario"))

lat_summary$lwr <- lat_summary$mean_lat - lat_summary$sd_lat / sqrt(lat_summary$n_lat)
lat_summary$upr <- lat_summary$mean_lat + lat_summary$sd_lat / sqrt(lat_summary$n_lat)

present_rows <- lat_summary[lat_summary$scenario == "present", ]

present_future1 <- present_rows
present_future1$futureScenario <- "SSP 370"

present_future2 <- present_rows
present_future2$futureScenario <- "SSP 585"

presents <- rbind(
  present_future1,
  present_future1,
  present_future2,
  present_future2
)

presents$year <- rep(c("2050", "2090"), each = 4, times = 2)


lat_summary <- lat_summary[lat_summary$scenario !="present",] # remove os presentes

lat_summary$futureScenario <- ifelse(grepl("1$", lat_summary$scenario), "SSP 370", "SSP 585") #cria futureScenario
lat_summary$year <- ifelse(grepl("_1", lat_summary$scenario), "2050", "2090") #cria year

lat_summary <- rbind(presents, lat_summary)

lat_summary <- lat_summary[order(lat_summary$species),]

lat_summary

lat_summary$spGroup <- ifelse(
  grepl("alt|ita", lat_summary$species),
  "Open Areas & Plateau", "Forest & Mountain")
lat_summary$spGroup <- factor(lat_summary$spGroup, levels = c("Forest & Mountain", "Open Areas & Plateau"))

lat_summary$horizon <- ifelse(lat_summary$scenario == "present", "Present", "Future")
lat_summary$horizon <- factor(lat_summary$horizon, levels = c("Present", "Future"))

lat_summary$shape <- ifelse(
  lat_summary$horizon == "Present",
  "Present",
  as.character(lat_summary$year)
)

lat_summary$futureScenario <- factor(lat_summary$futureScenario, levels = c("SSP 370", "SSP 585"))
lat_summary$year <- factor(lat_summary$year, levels = c("2050", "2090"))
lat_summary$species <- factor(lat_summary$species, levels = c("Bothrops fonsecai", "Bothrops cotiara", "Bothrops itapetiningae", "Bothrops alternatus"))

write.csv(lat_summary, here::here("outputs", "tables", "fig_lat.csv"), row.names = FALSE)


# calculating area of binary rasters --------------------------------------
crs_aea_wgs84 <- "+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

files <- list.files(here::here("outputs","rasters","vars_10m","binary"), pattern = '.tif', full.names = T)

dfs <- read.csv(here::here("data","raw","occurrences","group_alternatus.csv"), h=T)

sds <- split(dfs, dfs$species)

species <- names(sds)

rasters <- lapply(files, terra::rast) # lendo os arquivos como raster

areas <- data.frame(
  species = character(),
  scenario = character(),
  area_km2 = numeric()
) # criando uma lista para armazenar resultados

for (i in 1:length(rasters)) {
  try({
    r <- rasters[[i]]
    ### Raster to Vector -------------------------------
    
    rast <- raster::stack(r)  
    spp <- raster::rasterToPolygons(rast)
    
    ### Separando presenca de ausencia -----------------------
    sp <- sf::st_as_sf(spp)

    colnames(sp)[1] <- "presence"
    sp1 <- sp[which(sp$presence==1),] # separei em um objeto apenas as presencas

    ### Criando vetores sem grid e calculando area ------------------
    
    sp2 <- aggregate(sp1, by=list(sp1$presence), FUN = sum)
    
    sp2 <- sf::st_transform(sp2, crs=crs_aea_wgs84)
    area <- sf::st_area(sp2)
    area_km2 <- units::set_units(area, "km^2")
    
    areas[i,1] <- terra::sources(r)
    areas[i,2] <- files[i]
    areas[i,3] <- area_km2
    
  })
}

areas$scenario <- ifelse(grepl("cur", areas$scenario), "present", 
                         ifelse(grepl("fut11", areas$scenario), "future_11", 
                         ifelse(grepl("fut12", areas$scenario), "future_12",
                         ifelse(grepl("fut21", areas$scenario), "future_21",
                         ifelse(grepl("fut22", areas$scenario), "future_22", NA)))))


areas$species <- ifelse(grepl("alt", areas$species), "Bothrops alternatus",
                        ifelse(grepl("coti", areas$species), "Bothrops cotiara",
                        ifelse(grepl("fons", areas$species), "Bothrops fonsecai",
                        ifelse(grepl("itape", areas$species), "Bothrops itapetiningae", NA))))



# processing for figure 2 -------------------------------------------------

present_rows <- areas[areas$scenario == "present", ]

present_future1 <- present_rows
present_future1$futureScenario <- "SSP 370"

present_future2 <- present_rows
present_future2$futureScenario <- "SSP 585"

presents <- rbind(
  present_future1,
  present_future1,
  present_future2,
  present_future2
)

presents$year <- rep(c("2050", "2090"), each = 4, times = 2)


areas <- areas[areas$scenario !="present",] # remove os presentes

areas$futureScenario <- ifelse(grepl("1$", areas$scenario), "SSP 370", "SSP 585") #cria futureScenario
areas$year <- ifelse(grepl("_1", areas$scenario), "2050", "2090") #cria year

areas <- rbind(presents, areas)

areas <- areas[order(areas$species),]

areas

areas$spGroup <- ifelse(grepl("alt|ita", areas$species), "Open Areas & Plateau", "Forest & Mountain")
areas$horizon <- ifelse(areas$scenario == "present", "Present", "Future")

areas$spGroup <- factor(areas$spGroup, levels = c("Forest & Mountain", "Open Areas & Plateau"))
areas$species <- factor(areas$species, levels = c("Bothrops fonsecai", "Bothrops cotiara", "Bothrops itapetiningae", "Bothrops alternatus"))

areas$horizon <- factor(areas$horizon, levels = c("Present", "Future"))

areas$shape <- ifelse(
  areas$horizon == "Present",
  "Present",
  as.character(areas$year)
)


write.csv(areas, here::here("outputs","tables","area_scenarios.csv"), row.names = FALSE)




# processing for refugia maps ---------------------------------------------
sds <- split(dfs, dfs$species)

for(i in seq_along(names(sds))){
  bin <- terra::rast(files[grepl(names(sds[i]), files)])  
  refugia <- sum(bin, na.rm = TRUE)
  refugia <- terra::ifel(refugia == 0, NA, refugia)
  terra::writeRaster(refugia, here::here("outputs", "rasters", "refugia", paste0(names(sds[i]), "_refugia.tif")), overwrite = TRUE)
}


# extracting elevation data -----------------------------------------------
alt <- list()

#for 10m resolution
for(i in seq_along(files)){
  rast <- terra::rast(files[i])
  vect <- terra::as.polygons(rast, aggregate=FALSE)
  
  if(names(vect)=="layer"){
    ext <- terra::zonal(topo, vect[vect$layer==1,], fun = 'median')  
  }else{
    ext <- terra::zonal(topo, vect[vect$ensemble_weighted==1,], fun = 'median')  
  }
  
  name <- basename(files[i])
  n <- tools::file_path_sans_ext(name)
  
  alt[[length(alt) + 1]] <- data.frame(species = n, elevation = ext)
  
}
altitude <- do.call(rbind, alt)

names(altitude)[2] <- "elevation"

# altitude$scenario <- ifelse(grepl("cur", altitude$scenario), "present", 
#                             ifelse(grepl("fut11", altitude$scenario), "future_11", 
#                                    ifelse(grepl("fut12", altitude$scenario), "future_12",
#                                           ifelse(grepl("fut21", altitude$scenario), "future_21",
#                                                  ifelse(grepl("fut22", altitude$scenario), "future_22", NA)))))

altitude$scenario <- factor(altitude$scenario, levels = c("present", "future_11", "future_12", "future_21", "future_22"))

altitude$species <- ifelse(grepl("alt", altitude$species), "Bothrops alternatus",
                           ifelse(grepl("coti", altitude$species), "Bothrops cotiara",
                                  ifelse(grepl("fons", altitude$species), "Bothrops fonsecai",
                                         ifelse(grepl("itape", altitude$species), "Bothrops itapetiningae", NA))))

altitude$species <- factor(altitude$species, levels = c("Bothrops fonsecai", "Bothrops cotiara", "Bothrops itapetiningae", "Bothrops alternatus"))

write.csv(altitude, here::here("outputs","tables","elevation_scenarios.csv"), row.names = FALSE)


# Elevation summary - FIGURE ----------------------------------------------

mean_df <- aggregate(
  elevation ~ species + scenario,
  data = altitude,
  FUN = mean
)
names(mean_df)[names(mean_df) == "elevation"] <- "mean_elevation"

sd_df <- aggregate(
  elevation ~ species + scenario,
  data = altitude,
  FUN = sd
)
names(sd_df)[names(sd_df) == "elevation"] <- "sd_elevation"

n_df <- aggregate(
  elevation ~ species + scenario,
  data = altitude,
  FUN = length
)
names(n_df)[names(n_df) == "elevation"] <- "n_elevation"

elev_summary <- merge(mean_df, sd_df,
                      by = c("species","scenario"))

elev_summary <- merge(elev_summary, n_df,
                      by = c("species","scenario"))

elev_summary$lwr <- elev_summary$mean_elevation - elev_summary$sd_elevation / sqrt(elev_summary$n_elevation)
elev_summary$upr <- elev_summary$mean_elevation + elev_summary$sd_elevation / sqrt(elev_summary$n_elevation)

# present rows duplication
present_rows <- elev_summary[elev_summary$scenario == "present", ]

present_future1 <- present_rows
present_future1$futureScenario <- "SSP 370"

present_future2 <- present_rows
present_future2$futureScenario <- "SSP 585"

presents <- rbind(
  present_future1,
  present_future1,
  present_future2,
  present_future2
)

presents$year <- rep(c("2050", "2090"), each = 4, times = 2)

elev_summary <- elev_summary[elev_summary$scenario !="present",] # remove os presentes

elev_summary$futureScenario <- ifelse(grepl("1$", elev_summary$scenario), "SSP 370", "SSP 585") #cria futureScenario
elev_summary$year <- ifelse(grepl("_1", elev_summary$scenario), "2050", "2090") #cria year

elev_summary <- rbind(presents, elev_summary)

elev_summary <- elev_summary[order(elev_summary$species),]

elev_summary$spGroup <- ifelse(
  grepl("alt|ita", elev_summary$species),
  "Open Areas & Plateau", "Forest & Mountain")
elev_summary$spGroup <- factor(elev_summary$spGroup, levels = c("Forest & Mountain", "Open Areas & Plateau"))

elev_summary$futureScenario <- factor(elev_summary$futureScenario, levels = c("SSP 370", "SSP 585"))
elev_summary$year <- factor(elev_summary$year, levels = c("2050", "2090"))
elev_summary$species <- factor(elev_summary$species, levels = c("Bothrops fonsecai", "Bothrops cotiara", "Bothrops itapetiningae", "Bothrops alternatus"))

elev_summary$horizon <- ifelse(elev_summary$scenario == "present", "Present", "Future")
elev_summary$horizon <- factor(elev_summary$horizon, levels = c("Present", "Future"))

elev_summary$shape <- ifelse(elev_summary$horizon == "Present", "Present", as.character(elev_summary$year))

#go to figures with elev_summary
# analyses ----------------------------------------------------------------
names(altitude)
str(altitude)

altitude$species <- ifelse(grepl("alt", altitude$species), "Bothrops alternatus", altitude$species)
altitude$species <- ifelse(grepl("coti", altitude$species), "Bothrops cotiara", altitude$species)
altitude$species <- ifelse(grepl("fons", altitude$species), "Bothrops fonsecai", altitude$species)
altitude$species <- ifelse(grepl("itape", altitude$species), "Bothrops itapetiningae", altitude$species)

altitude$species <- factor(altitude$species, levels = c("Bothrops fonsecai", "Bothrops cotiara", "Bothrops itapetiningae", "Bothrops alternatus"))

altitude$spGroup <- ifelse(
  grepl("alt|ita", altitude$species),
  "Open Areas & Plateau", "Forest & Mountain")

altitude$spGroup <- factor(altitude$spGroup, levels = c("Forest & Mountain", "Open Areas & Plateau"))

altitude$species <- factor(altitude$species, levels = c("Bothrops fonsecai", "Bothrops cotiara", "Bothrops itapetiningae", "Bothrops alternatus"))

elev_summary

library(lmerTest)
m <- lmerTest::lmer(elevation ~ scenario * spGroup + (1|species), data = altitude)
anova(m)

m2 <- lmerTest::lmer(elevation ~ scenario + spGroup + (1|species), data = altitude)

anova(m, m2)
summary(m)
