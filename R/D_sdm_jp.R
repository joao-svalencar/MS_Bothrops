# loading packages --------------------------------------------------------

library(dplyr)
library(here)
library(terra)
library(sf)
library(geosphere)
library(flexsdm) # calibration area and colinearity vars
library(enmSdmX) # maybe delete (see line 80)
library(sdm)
library(ecospat)
library(magick)

# Models loop --------------------------------------------------------
sds <- split(dfs, dfs$species) # read dfs in R_data.R line 22

# Attention: read processed predictor variables in R_data.R lines 6 to 16
# Attention: read object 'sa' in R_data.R line 19

for(i in seq_along(sds)){
  
  df <- sds[[i]]
  sp <- df[,c("species", "longitude", "latitude")]
  
# buffer background -------------------------------------------------------

  mdist <- geosphere::distm(sp[, c('longitude','latitude')])
  d <- max(mdist)
  w <- d/2
  
  area <- flexsdm::calib_area(sp,'longitude','latitude', 
                              method = c("bmcp", width = w),
                              crs = terra::crs(sa))

# variables selection and background cropping -----------------------------

  preds.buffer <- terra::crop(preds.pres_10m, area)
  preds.buffer.mask <- terra::mask(preds.buffer, area)

  vif_var <- flexsdm::correct_colinvar(preds.buffer.mask, method = c("vif", th = "10")) #testar outro pacote
  vars <- terra::subset(preds.buffer.mask, vif_var$vif_table$Variables)

# verifying spatial autocorrelation ---------------------------------------
  
  sp_th <- read.csv(here::here("data","processed","occurrences_thin",paste0(sp[1,1],"_thin1.csv")))
  
  extract_pts <- terra::extract(vars, sp_th[,c(2,3)])
  
  matriz.autocorrelacao <- cbind(sp_th, extract_pts)
  matriz.autocorrelacao <- na.omit(matriz.autocorrelacao)
  
  ## testando autocorrelacao para cada grau
  png(file=here::here("outputs","autocorrelation", paste0(sp[1,1],"_10m.png")))
  
  ecospat::ecospat.mantel.correlogram(matriz.autocorrelacao, 
                                      colxy=2:3,n=100, 
                                      colvar=5:dim(matriz.autocorrelacao)[2],
                                      max=8, nclass=80, nperm=100) 
  
  dev.off()
  
# creating pseudo-absences ------------------------------------------------

  sp1 <- matriz.autocorrelacao[,2:3]
  sp1$Occurrence <- 1 
  sp::coordinates(sp1) <- c('longitude','latitude')
  pseudo <- terra::spatSample(preds.buffer.mask[[1]], nrow(sp1)*10, 
                              method="random", 
                              xy=T, cells=T, lonlat=T)
  pseudo$Occurrence <- 0
  pseudo <- na.omit(pseudo)
  sp::coordinates(pseudo) <- c('x','y')
  terra::crs(pseudo) <- terra::crs(sa)
  
  pseudo@data <- pseudo@data[,"Occurrence", drop=F]
  
  buffer_pts <- flexsdm::calib_area(sp1@coords,'longitude','latitude',
                                    method = c("buffer", width = 50000),
                                    crs = terra::crs(sa))
  
  buffer_pts <- enmSdmX::spatVectorToSpatial(buffer_pts)
  
  pts_outside <- pseudo[is.na(sp::over(pseudo, buffer_pts)),]
  
  terra::crs(sp1) <- terra::crs(sa)
  
  pts_outside@data$Occurrence <- 0 
  pts_outside@data <- pts_outside@data[,"Occurrence",drop=F]
  
  speciesN <- rbind(sp1, pts_outside)
  
  ### 9 - ajuste e avaliacao dos modelos 
  ## criando objeto smdata para implementar algoritmos no pacote sdm
  vars1 <- raster::stack(vars) # transforma para RasterStack
  d <- sdm::sdmData(Occurrence~., train=speciesN, predictors=vars1)
 
  ### 10 - implementando o modelo, separando 30% dos pontos para teste
  m <- sdm::sdm(Occurrence~., d, methods = c('glm','rf','bioclim'), test.p=30, replication='boot', n=10)
  var_imp <- getVarImp(m, id=1, wtest = 'test.dep')
  #sdm::gui(m)
  capture.output(var_imp, file = here::here("outputs","var_imp", paste0(sp[1,1],"_10m_vars.txt")))
  capture.output(m, file = here::here("outputs","model_performance", paste0(sp[1,1],"_10m_sdm.txt")))
  
  ### 11 - projetando modelos no espaco e fazendo consenso
  # fazendo um modelo consenso pesado por AUC e TSS:
  
  # PRESENTE
  pres_consenso_10m <- sdm::ensemble(m, vars1, filename = here::here("outputs","rasters","vars_10m","adequability",paste0(sp[1,1],"_10m.tif")), setting = list(method='weighted', stat='TSS', expr='auc > 0.7 & tss > 0.5', opt=1), overwrite = T)
  
  # FUT1.1 - 2041-2070 SSP370
  vars_f11_10m <- terra::subset(preds.fut11_10m,vif_var$vif_table$Variables)
  buf <- terra::crop(vars_f11_10m, area)
  vars2 <- terra::mask(buf, area)
  
  fut1.1_consenso_10m <- sdm::ensemble(m, vars2, filename = here::here("outputs","rasters","vars_10m","adequability", paste0(sp[1,1],"_fut11_10m.tif")), setting = list(method='weighted', stat='TSS', expr='auc > 0.7 & tss > 0.5', opt=1), overwrite = T)
  
  # FUT1.2 - 2041-2070 SSP585
  vars_f12_10m <- terra::subset(preds.fut12_10m,vif_var$vif_table$Variables)
  buf <- terra::crop(vars_f12_10m, area)
  vars3 <- terra::mask(buf, area)
  
  fut1.2_consenso_10m <- sdm::ensemble(m, vars3, filename = here::here("outputs","rasters","vars_10m","adequability", paste0(sp[1,1],"_fut12_10m.tif")), setting = list(method='weighted', stat='TSS', expr='auc > 0.7 & tss > 0.5', opt=1), overwrite = T)
  
  # FUT2.1 - 2071-2100 SSP370
  vars_f21_10m <- terra::subset(preds.fut21_10m,vif_var$vif_table$Variables)
  buf <- terra::crop(vars_f21_10m, area)
  vars4 <- terra::mask(buf, area)
  
  fut2.1_consenso_10m <- sdm::ensemble(m, vars4, filename = here::here("outputs","rasters","vars_10m","adequability", paste0(sp[1,1],"_fut21_10m.tif")), setting = list(method='weighted', stat='TSS', expr='auc > 0.7 & tss > 0.5', opt=1), overwrite = T)
  
  # FUT2.2 - 2071-2100 SSP585
  vars_f22_10m <- terra::subset(preds.fut22_10m,vif_var$vif_table$Variables)
  buf <- terra::crop(vars_f22_10m, area)
  vars5 <- terra::mask(buf, area)
  
  fut2.2_consenso_10m <- sdm::ensemble(m, vars5, filename = here::here("outputs","rasters","vars_10m","adequability", paste0(sp[1,1],"_fut22_10m.tif")), setting = list(method='weighted', stat='TSS', expr='auc > 0.7 & tss > 0.5', opt=1), overwrite = T)
  
  ###  12 - definindo limiar binario (threshold)
  speciesn <- sf::st_as_sf(speciesN)
  area1 <- sf::st_as_sf(area)
  suppressWarnings(
    speciesM <- sf::st_intersection(speciesn,area1)
  )
  
  obs <- speciesM$Occurrence
  pred <- terra::extract(pres_consenso_10m, speciesM) 
  ev <- sdm::evaluates(obs, pred)
  ev1 <- ev@threshold_based$threshold
  ev2 <- ev1[1] #se=sp
  
  binario <- pres_consenso_10m>=ev2
  binario <- terra::rast(binario)
  writeRaster(binario, filename=here::here("outputs","rasters","vars_10m","binary", paste0(sp[1,1],"_10m_bin_cur.tif")), overwrite=TRUE)
  
  binario_fut1.1 <- fut1.1_consenso_10m>=ev2
  writeRaster(binario_fut1.1, filename=here::here("outputs","rasters","vars_10m","binary", paste0(sp[1,1],"_10m_bin_fut11.tif")), overwrite=TRUE)
 
  binario_fut1.2 <- fut1.2_consenso_10m>=ev2
  writeRaster(binario_fut1.2, filename=here::here("outputs","rasters","vars_10m","binary", paste0(sp[1,1],"_10m_bin_fut12.tif")), overwrite=TRUE)
  
  binario_fut2.1 <- fut2.1_consenso_10m>=ev2
  writeRaster(binario_fut2.1, filename=here::here("outputs","rasters","vars_10m","binary", paste0(sp[1,1],"_10m_bin_fut21.tif")), overwrite=TRUE)
  
  binario_fut2.2 <- fut2.2_consenso_10m>=ev2
  writeRaster(binario_fut2.2, filename=here::here("outputs","rasters","vars_10m","binary", paste0(sp[1,1],"_10m_bin_fut22.tif")), overwrite=TRUE)
  
  # gif adequability:
  rasters <- list(pres_consenso_10m,
                  fut1.1_consenso_10m,
                  fut1.2_consenso_10m,
                  fut2.1_consenso_10m,
                  fut2.2_consenso_10m)
  
  titles <- paste(sp[1,1], c("Present", "Futuro 1.1", "Futuro 1.2", "Futuro 2.1", "Futuro 2.2"), sep=" - ")
  
  imgs <- mapply(function(r, title) {
    tmp <- tempfile(fileext = ".png")
    png(tmp, width = 800, height = 600)
    
    # Optional styling
    par(cex.main = 1.5, font.main = 2)
    
    # Each plot gets its own scale automatically
    plot(r, main = title)
    
    dev.off()
    image_read(tmp)
  }, rasters, titles, SIMPLIFY = FALSE) #trocar para rasters_bin para binarios
  
  gif <- image_animate(image_join(imgs), fps = 1)
  image_write(gif, here::here("outputs","gifs", paste0(sp[1,1],"_10m.gif")))
  
  #gifs binarios:
  rasters_bin <- list(binario,
                      binario_fut1.1,
                      binario_fut1.2,
                      binario_fut2.1,
                      binario_fut2.2)
  
  imgs_bin <- mapply(function(r, title) {
    tmp <- tempfile(fileext = ".png")
    png(tmp, width = 800, height = 600)
    
    # Optional styling
    par(cex.main = 1.5, font.main = 2)
    
    # Each plot gets its own scale automatically
    plot(r, main = title)
    
    dev.off()
    image_read(tmp)
  }, rasters_bin, titles, SIMPLIFY = FALSE) #trocar para rasters_bin para binarios
  
  gif_bin <- image_animate(image_join(imgs_bin), fps = 1)
  image_write(gif_bin, here::here("outputs","gifs", paste0(sp[1,1],"_10m_bin.gif")))
  
  
  ### 13 - criando máscara com binarios para recortar variaveis originais
  #binario[binario==0] <- NA
  
  #binario_fut1.1[binario_fut1.1==0] <- NA
  
  #binario_fut1.2[binario_fut1.2==0] <- NA
  
  #binario_fut2.1[binario_fut2.1==0] <- NA
  
  #binario_fut2.2[binario_fut2.2==0] <- NA
  
  ### 14 - projetando nas variáveis originais
  # Present
  vars_pres <- terra::subset(preds.pres,vif_var$vif_table$Variables)
  bin_vec <- terra::as.polygons(binario)
  buf <- terra::crop(vars_pres, bin_vec)
  
  #vars_pres_or <- terra::mask(buf, bin_vec)
  vars_pres_or <- terra::mask(buf, bin_vec[bin_vec$layer==1,])
  
  pres_consenso <- sdm::ensemble(m, vars_pres_or, filename = here::here("outputs","rasters","original vars","adequability",paste0(sp[1,1],"_cur.tif")),setting = list(method='weighted', stat='TSS', expr='auc > 0.7 & tss > 0.5', opt=1), overwrite=T)
  
  # Fut1.1 - 2041-2070 ssp370
  vars_f11 <- terra::subset(preds.fut11,vif_var$vif_table$Variables)
  bin_vec <- terra::as.polygons(binario_fut1.1)
  buf <- terra::crop(vars_f11, bin_vec)
  
  #vars_f11_or <- terra::mask(buf, bin_vec)
  vars_f11_or <- terra::mask(buf, bin_vec[bin_vec$ensemble_weighted==1,])
  
  fut1.1_consenso <- sdm::ensemble(m, vars_f11_or, filename = here::here("outputs","rasters","original vars","adequability",paste0(sp[1,1],"_fut11.tif")), setting = list(method='weighted', stat='TSS', expr='auc > 0.7 & tss > 0.5', opt=1), overwrite = T)
  
  # Fut1.2 - 2041-2070 ssp585
  vars_f12 <- terra::subset(preds.fut12,vif_var$vif_table$Variables)
  bin_vec <- terra::as.polygons(binario_fut1.2)
  buf <- terra::crop(vars_f12, bin_vec)
  
  #vars_f12_or <- terra::mask(buf, bin_vec)
  vars_f12_or <- terra::mask(buf, bin_vec[bin_vec$ensemble_weighted==1,])
  
  fut1.2_consenso <- sdm::ensemble(m, vars_f12_or, filename = here::here("outputs","rasters","original vars","adequability",paste0(sp[1,1],"_fut12.tif")), setting = list(method='weighted', stat='TSS', expr='auc > 0.7 & tss > 0.5', opt=1), overwrite = T)
  
  # Fut2.1 - 2071-2100 ssp370
  vars_f21 <- terra::subset(preds.fut21,vif_var$vif_table$Variables)
  bin_vec <- terra::as.polygons(binario_fut2.1)
  buf <- terra::crop(vars_f21, bin_vec)
  
  #vars_f21_or <- terra::mask(buf, bin_vec)
  vars_f21_or <- terra::mask(buf, bin_vec[bin_vec$ensemble_weighted==1,])
  
  fut2.1_consenso <- sdm::ensemble(m, vars_f21_or, filename = here::here("outputs","rasters","original vars","adequability",paste0(sp[1,1],"_fut21.tif")), setting = list(method='weighted', stat='TSS', expr='auc > 0.7 & tss > 0.5', opt=1), overwrite = T)
  
  # Fut2.2 - 2071-2100 ssp585
  vars_f22 <- terra::subset(preds.fut22,vif_var$vif_table$Variables)
  bin_vec <- terra::as.polygons(binario_fut2.2)
  buf <- terra::crop(vars_f22, bin_vec)
  
  #vars_f22_or <- terra::mask(buf, bin_vec)
  vars_f22_or <- terra::mask(buf, bin_vec[bin_vec$ensemble_weighted==1,])
  
  fut2.2_consenso <- sdm::ensemble(m, vars_f22_or, filename = here::here("outputs","rasters","original vars","adequability", paste0(sp[1,1],"_fut22.tif")), setting = list(method='weighted', stat='TSS', expr='auc > 0.7 & tss > 0.5', opt=1), overwrite = T)
  
  # 15 - criando binários das projeções refinadas
  
  binario_or <- pres_consenso>=ev2
  writeRaster(binario_or, filename=here::here("outputs","rasters","original vars","binary", paste0(sp[1,1],"_bin_cur.tif")), overwrite=TRUE)
  
  binario_fut1.1_or <- fut1.1_consenso>=ev2
  writeRaster(binario_fut1.1_or, filename=here::here("outputs","rasters","original vars","binary", paste0(sp[1,1],"_bin_fut11.tif")), overwrite=TRUE)
  
  binario_fut1.2_or <- fut1.2_consenso>=ev2
  writeRaster(binario_fut1.2_or, filename=here::here("outputs","rasters","original vars","binary", paste0(sp[1,1],"_bin_fut12.tif")), overwrite=TRUE)
  
  binario_fut2.1_or <- fut2.1_consenso>=ev2
  writeRaster(binario_fut2.1_or, filename=here::here("outputs","rasters","original vars","binary", paste0(sp[1,1],"_bin_fut21.tif")), overwrite=TRUE)
  
  binario_fut2.2_or <- fut2.2_consenso>=ev2
  writeRaster(binario_fut2.2_or, filename=here::here("outputs","rasters","original vars","binary", paste0(sp[1,1],"_bin_fut22.tif")), overwrite=TRUE)
 
  print(paste(sp[1,1], "model", "done!", sep=" "))
}

# end of script -----------------------------------------------------------
