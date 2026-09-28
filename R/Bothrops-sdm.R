########################################################
########################################################
##       MODELAGEM DE DISTRIBUICAO DE ESPECIES        ##
########################################################
########################################################

### 1 - Instalando e carregando pacotes ----------
library(dplyr)
library(here)
library(terra)
library(sf)
library(spThin)
library(geosphere)
library(flexsdm)
library(enmSdmX)
library(sdm)
library(ecospat)

# Foi criada uma pasta com o script, a planilha com os pontos e uma pasta com as variaveis
# Em seguida foi criado um projeto no R dentro dessa pasta para que tudo seja feito dentro da pasta, sem necessidade de estabelecer um setwd()

### 2 - Criando DF para cada especie ----------
dfs <- read.csv(here::here("data","raw","group_alternatus.csv"), h=T)

sds <- split(dfs, dfs$species)

for(t in 1:length(sds)){
  df <- sds[[t]]
  df1 <- df[,c("species", "longitude", "latitude")] # garantindo que as colunas de todos os df estejam na mesma ordem
  
  write.table(df1,paste(here::here("data","raw"),"/",df1[1,1],".txt",sep=""),sep="\t", row.names = F)
}

### 3 - Preparando variaveis bioclimaticas na area de plotagem ----------
# Selecionei um subconjunto de 12 variaveis excluindo variais que misturam precipitação e temperatura e as médias. Deixando apenas variaveis extremas.
# Adicionei as variaveis de todos os cenarios nas respectivas pastas dentro de "2_Predictors".

## Manipulando variaveis
eco <- terra::vect(here::here("data","area","ecoregions_2017.gpkg"))

# presente
preds.Chelsa <- list.files(here::here("data","predictors","current"), pattern='tif$', full.names = T)
preds <- terra::rast(preds.Chelsa)
preds.eco <- terra::crop(preds, eco)
preds.m <- terra::mask(preds.eco, eco)
preds.10m <- terra::aggregate(preds.m,fact=20, fun=mean)

# FUT1.1 (2041-2070 ssp370)
preds.Chelsa1 <- list.files(here::here("data","predictors","2041-2070_ssp370"), pattern='tif$', full.names = T)
preds1 <- terra::rast(preds.Chelsa1)
preds.eco1 <- terra::crop(preds1, eco)
preds.m1 <- terra::mask(preds.eco1, eco)
predsf11.10m <- terra::aggregate(preds.m1,fact=20, fun=mean)

# FUT1.2 (2041-2070 ssp585)
preds.Chelsa2 <- list.files(here::here("data","predictors","2041-2070_ssp585"), pattern='tif$', full.names = T)
preds2 <- terra::rast(preds.Chelsa2)
preds.eco2 <- terra::crop(preds2, eco)
preds.m2 <- terra::mask(preds.eco2, eco)
predsf12.10m <- terra::aggregate(preds.m2,fact=20, fun=mean)

# FUT2.1 (2071-2100 ssp370)
preds.Chelsa3 <- list.files(here::here("data","predictors","2071-2100_ssp370"), pattern='tif$', full.names = T)
preds3 <- terra::rast(preds.Chelsa3)
preds.eco3 <- terra::crop(preds3, eco)
preds.m3 <- terra::mask(preds.eco3, eco)
predsf21.10m <- terra::aggregate(preds.m3,fact=20, fun=mean)

# FUT2.2 (2071-2100 ssp585)
preds.Chelsa4 <- list.files(here::here("data","predictors","2071-2100_ssp585"), pattern='tif$', full.names = T)
preds4 <- terra::rast(preds.Chelsa4)
preds.eco4 <- terra::crop(preds4, eco)
preds.m4 <- terra::mask(preds.eco4, eco)
predsf22.10m <- terra::aggregate(preds.m4,fact=20, fun=mean)

# salvar variáveis agregadas e cortadas para SA
#salvas em "data" > "processed" > "predictors"
terra::writeRaster(preds.10m, filename=here::here("data","processed","predictor","current", "bio.tif"))
terra::writeRaster(predsf11.10m, filename=here::here("data","processed","predictor","f11", "bio.tif"))


# Models loop --------------------------------------------------------

spp <- names(sds)

for(i in 1:length(spp)){
  sp <- read.table(paste(here::here("data","raw"),"/",spp[i],".txt",sep=""), header = T)
  species <- sp::SpatialPointsDataFrame(sp[,2:3], sp)
  
  ### 4 - criando buffer para Background de modelos
  mdist <- geosphere::distm(sp[, c('longitude','latitude')]) # matriz de distancia entre pontos de ocorrência
  d <- max(mdist) # maior distância
  w <- d/2 # dividindo por dois para usar como area de calibração dos modelos
  
  area <- flexsdm::calib_area(sp,'longitude','latitude', 
                              method = c("bmcp",width = w),
                              crs = terra::crs(eco))
  
  terra::plot(area)
  terra::plot(eco, add=T)
  terra::plot(species, col="red", pch=16, add=T)
  
  ### 5 - tratando tendencia de coleta
  spThin::thin(loc.data = sp, # = data.frame com todos os pontos de ocorrencia
      lat.col = "latitude", long.col = "longitude", # nomes das colunas de coordenadas
      spec.col = "species", # nome da coluna de taxon
      # escolha uma distancia inicial minima entre um ponto e outro, digamos 50 km
      thin.par = 20, reps = 10, # thin.par = distancia em km de um ponto a outro; reps = numero de replicas de reamostragem
      locs.thinned.list.return = FALSE, 
      write.files = TRUE, 
      max.files = 1, 
      out.dir = here::here("data","processed"), # cada um cria sua pasta de output nesse passo ou usa getwd() para salvar os arquivos gerados no diretorio de trabalho ja definido
      out.base = sp[1,1], # out.base  nome do arquivo de saida seguindo de _thin1
      write.log.file = TRUE,
      log.file = here::here("data","processed","thin",paste0(sp[1,1],".txt")),
      verbose=F)
  
  ### 6 - selecionando e preparando variaveis preditoras recortadas por background
  preds.buffer <- terra::crop(preds.10m, area)
  preds.buffer.mask <- terra::mask(preds.buffer, area)
  #terra::plot(preds.buffer.mask)
  
  preds.buffer <- terra::crop(preds.m, area) #teste
  preds.buffer.mask <- terra::mask(preds.buffer, area) #teste
  
  vif_var <- flexsdm::correct_colinvar(preds.buffer.mask, method = c("vif", th = "10"))
  vif_var
  #vif_var$removed_variables
  vars <- terra::subset(preds.buffer.mask,vif_var$vif_table$Variables)
  #terra::plot(vars)
  
  ### 7 - analisando autocorrelacao espacial das variaveis
  sp_th <- read.csv(here::here("data","processed",paste0(sp[1,1],"_thin1.csv")))
  extract_pts <- terra::extract(vars, sp_th[,2:3])
  
  matriz.autocorrelacao <- cbind(sp_th, extract_pts)
  matriz.autocorrelacao <- na.omit(matriz.autocorrelacao)
  #class(matriz.autocorrelacao) #dataframe
  
  ## testando autocorrelacao para cada grau
  png(file=here::here("outputs","autocorrelation",paste(sp[1,1],".png",sep="")))
  
  ecospat::ecospat.mantel.correlogram(matriz.autocorrelacao, 
                                      colxy=2:3,n=100, 
                                      colvar=5:dim(matriz.autocorrelacao)[2],
                                      max=8, nclass=80, nperm=100) 
  
  dev.off()
  
  ### 8 - gerando pontos de pseudoausencia
  sp1 <- matriz.autocorrelacao[,2:3]
  sp1$Occurrence <- 1 
  sp::coordinates(sp1) <- c('longitude','latitude')
  pseudo <- terra::spatSample(preds.buffer.mask[[1]], nrow(sp1)*4, 
                              method="random", 
                              xy=T, cells=T, lonlat=T)
  pseudo$Occurrence <- 0
  pseudo <- na.omit(pseudo)
  sp::coordinates(pseudo) <- c('x','y')
  terra::crs(pseudo) <- "EPSG:4326"
  terra::plot(pseudo, cex=0.2, pch=16)
  
  pseudo@data <- pseudo@data[,"Occurrence", drop=F]
  buffer_pts <- flexsdm::calib_area(sp1@coords,'longitude','latitude',
                                    method = c("bmcp", 
                                               width = 50000),
                                    crs = terra::crs(eco))
  
  terra::crs(buffer_pts) <- "EPSG:4326"
  
  buffer_pts <- enmSdmX::spatVectorToSpatial(buffer_pts)
  pts_outside <- pseudo[is.na(sp::over(pseudo, buffer_pts)),]
  
  sp::proj4string(sp1) <- sp::proj4string(pseudo)
  sp::proj4string(species) <- sp::proj4string(pseudo)
  sp::proj4string(pts_outside) <- sp::proj4string(pseudo)
  pts_outside@data$Occurrence <- 0 
  pts_outside@data <- pts_outside@data[,"Occurrence",drop=F]
  
  # plotando objetos para verificar se tudo foi criado corretamente
  #terra::plot(eco) # plotando ecorregioes
  #terra::plot(pseudo, pch=16, col="grey", add=T) # plotando pseudoausencias (geradas apos thin)
  #terra::plot(buffer_pts, add=T, col="red") # plotando buffer em torno de cada presenca
  terra::plot(sp1, cex=0.2, pch=16, col="blue", add=T) # plotando pontos tratados (apos thin)
  #terra::plot(area, add=T) # plotando buffer em torno do conjunto de presencas
  terra::plot(pts_outside, cex=0.2, pch=16, col="black", add=T) # plotando pseudo-ausencias fora do buffer das presencas
  
  speciesN <- rbind(sp1, pts_outside)
  
  ### 9 - ajuste e avaliacao dos modelos 
  ## criando objeto smdata para implementar algoritmos no pacote sdm
  vars1 <- raster::stack(vars)
  d <- sdm::sdmData(Occurrence~., train=speciesN, predictors=vars1)
  d
  
  ### 10 - implementando o modelo, separando 30% dos pontos para teste
  m <- sdm::sdm(Occurrence~., d, methods = c('bioclim','glm','rf'), test.p=30, replication='boot', n=10)
  var_imp <- getVarImp(m, id=1, wtest = 'test.dep')
  #sdm::gui(m)
  
  capture.output(var_imp, file = here::here("outputs","var_imp", paste(sp[1,1],"_vars.txt", sep="")))
  capture.output(m, file = here::here("outputs","model_performance", paste(sp[1,1],"_sdm.txt", sep="")))
  
  ### 11 - projetando modelos no espaco e fazendo consenso
  # fazendo um modelo consenso pesado por AUC e TSS:
  # PRESENTE
  pres_consenso <- sdm::ensemble(m, vars1, filename = here::here("outputs","rasters","adequability",paste(sp[1,1],".tif",sep="")), 
                                 setting = list(method='weighted', stat='TSS', expr='auc > 0.7 & tss > 0.5', opt=1), overwrite = T)
  terra::plot(pres_consenso)
  
  p <- predict(m,vars1,filename=here::here("outputs","rasters","adequability",paste(sp[1,1],".tif",sep="")),mean=T,overwrite=T)
  terra::plot(p)
  
  # FUT1.1 - 2041-2070 SSP370
  # preparando variaveis
  vars_f11 <- terra::subset(predsf11.10m,vif_var$vif_table$Variables)
  buf <- terra::crop(vars_f11, area)
  buffer <- terra::mask(buf, area)
  vars2 <- stack(buffer)
  # projetando modelo - consenso
  fut1.1_consenso <- sdm::ensemble(m, vars2, filename = here::here("outputs","rasters","adequability", paste(sp[1,1],"_fut1-1.tif", sep="")), 
                                   setting = list(method='weighted', stat='TSS', expr='auc > 0.7 & tss > 0.5', opt=1), overwrite = T)
  terra::plot(fut1.1_consenso)
  
  p <- predict(m,vars2,filename=here::here("outputs","rasters","adequability",paste(sp[1,1],"_fut11.tif",sep="")),mean=T,overwrite=T)
  terra::plot(p)
  
  # FUT1.2 - 2041-2070 SSP585
  vars_f12 <- terra::subset(predsf12.10m,vif_var$vif_table$Variables)
  buf <- terra::crop(vars_f12, area)
  buffer <- terra::mask(buf, area)
  vars3 <- stack(buffer)
  # projetando modelo - consenso
  fut1.2_consenso <- sdm::ensemble(m, vars3, filename = here::here("outputs","rasters","adequability", paste(sp[1,1],"_fut1-2.tif", sep="")), 
                                   setting = list(method='weighted', stat='TSS', expr='auc > 0.7 & tss > 0.5', opt=1), overwrite = T)
  terra::plot(fut1.2_consenso)
  
  # FUT2.1 - 2071-2100 SSP370
  vars_f21 <- terra::subset(predsf21.10m,vif_var$vif_table$Variables)
  buf <- terra::crop(vars_f21, area)
  buffer <- terra::mask(buf, area)
  vars4 <- stack(buffer)
  # projetando modelo - consenso
  fut2.1_consenso <- sdm::ensemble(m, vars4, filename = here::here("outputs","rasters","adequability", paste(sp[1,1],"_fut2-1.tif", sep="")), 
                                   setting = list(method='weighted', stat='TSS', expr='auc > 0.7 & tss > 0.5', opt=1), overwrite = T)
  terra::plot(fut2.1_consenso)
  
  # FUT2.2 - 2071-2100 SSP585
  vars_f22 <- terra::subset(predsf22.10m,vif_var$vif_table$Variables)
  buf <- terra::crop(vars_f22, area)
  buffer <- terra::mask(buf, area)
  vars5 <- stack(buffer)
  # projetando modelo - consenso
  fut2.2_consenso <- sdm::ensemble(m, vars5, filename = here::here("outputs","rasters","adequability", paste(sp[1,1],"_fut2-2.tif", sep="")), 
                                   setting = list(method='weighted', stat='TSS', expr='auc > 0.7 & tss > 0.5', opt=1), overwrite = T)
  terra::plot(fut2.2_consenso)
  
  ###  12 - definindo limiar binario (threshold)
  speciesn <- sf::st_as_sf(speciesN)
  area1 <- sf::st_as_sf(area)
  speciesM <- sf::st_intersection(speciesn,area1)
  terra::plot(speciesM)
  obs <- speciesM$Occurrence
  pred <- raster::extract(pres_consenso, speciesM) 
  ev <- sdm::evaluates(obs, pred)
  ev1 <- ev@threshold_based$threshold
  ev2 <- ev1[1]
  
  binario <- pres_consenso>=ev2
  terra::plot(binario)
  writeRaster(binario, filename=here::here("outputs","rasters","binary", paste(sp[1,1],"_bin_cur.tif", sep="")), overwrite=TRUE)
  
  binario_fut1.1 <- fut1.1_consenso>=ev2
  terra::plot(binario_fut1.1)
  writeRaster(binario_fut1.1, filename=here::here("outputs","rasters","binary", paste(sp[1,1],"_bin_fut11.tif", sep="")), overwrite=TRUE)
  
  binario_fut1.2 <- fut1.2_consenso>=ev2
  terra::plot(binario_fut1.2)
  writeRaster(binario_fut1.2, filename=here::here("outputs","rasters","binary", paste(sp[1,1],"_bin_fut12.tif", sep="")), overwrite=TRUE)
  
  binario_fut2.1 <- fut2.1_consenso>=ev2
  terra::plot(binario_fut2.1)
  writeRaster(binario_fut2.1, filename=here::here("outputs","rasters","binary", paste(sp[1,1],"_bin_fut21.tif", sep="")), overwrite=TRUE)
  
  binario_fut2.2 <- fut2.2_consenso>=ev2
  terra::plot(binario_fut2.2)
  writeRaster(binario_fut2.2, filename=here::here("outputs","rasters","binary", paste(sp[1,1],"_bin_fut22.tif", sep="")), overwrite=TRUE)
  
  print(paste(sp[1,1], "model", "done!", sep=" "))
}
