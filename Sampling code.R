library(rgdal)
library(raster)
library(rgeos)
library(spsann)
library(sp)
library(ICSNP)
library(velox)
library(spcosa)
library(spatstat)
library(dplyr)
library(tibble)
library(tidyr)



area_studio <- readOGR(dsn = 'c:/Users/MARCO/Dropbox/Dottorato/Cartografia/Layers', layer = 'Pinetamacchia')
ndvi_clip <- raster("C:/Users/MARCO/Dropbox/Dottorato/Cartografia/Raster remote sensing/NDVImap_SITEB.tif")
sr <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 " 
ndvi_clip <- projectRaster(ndvi_clip, crs = sr)



##### FUNZIONE

sampleboost <- function(x, ignorance, boundary, nplot, perm, quant, approach)
  {
  ndvi.vx <-velox(x)
  igno.vx <- velox(ignorance)
  result<-list()
  distanze<-matrix(ncol=1, nrow = perm)
  pb <- txtProgressBar(min = 0, max = perm, style = 3)
  for (i in 1:perm){
    punti_random <- spsample(boundary, n=nplot, type='random')
    sampling_points <- as(punti_random, "data.frame")
    xy <- sampling_points[,c(1,2)]
    
    spdf <- SpatialPointsDataFrame(coords = xy, data = sampling_points,
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    
    
    estratti <- data.frame(coordinates(spdf),ndvi.vx$extract_points(sp = spdf))
    names(estratti) <- c("x", "y", "ndvi")
    estratti$ignorance <- igno.vx$extract_points(sp = spdf)
    result[[i]]<-data.frame(estratti)
    dataset_points <- cbind(xy, ID = 1:NROW(xy))
    
    pairwise_distances <- distm(dataset_points[,1:2])
    distanze[[i]] <- total.dist(pairwise_distances, method = "euclidean", square = FALSE, p = 0)
    setTxtProgressBar(pb, i)
  }}


