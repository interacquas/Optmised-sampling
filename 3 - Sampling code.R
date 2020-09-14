if(!require(rgdal)){install.packages("rgdal"); library(rgdal)} 
if(!require(raster)){install.packages("raster"); library(raster)} 
if(!require(rgeos)){install.packages("rgeos"); library(rgeos)} 
if(!require(spsann)){install.packages("spsann"); library(spsann)} 
if(!require(sp)){install.packages("sp"); library(sp)} 
if(!require(ICSNP)){install.packages("ICSNP"); library(ICSNP)} 
if(!require(spatstat)){install.packages("spatstat"); library(spatstat)} 
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)} 
if(!require(tibble)){install.packages("tibble"); library(tibble)} 
if(!require(tidyr)){install.packages("tidyr"); library(tidyr)} 
if(!require(geosphere)){install.packages("geosphere"); library(geosphere)} 
if(!require(Rfast)){install.packages("Rfast"); library(Rfast)} 
library(ggplot2)

### CARICO I LAYER

# area studio
site <- readOGR(dsn = 'gis zone campionamento', layer = 'campionamento meno esclusione')


# ndvi
igno_map <- raster("MATERIALE PER LA FUNZIONE/Mappa Ignoranza 5 Km.tif")

# ignoranza

ndvi_map <- raster("MATERIALE PER LA FUNZIONE//MAPPA NDVI area studio_28m.tif")

newproj <- '+init=EPSG:3035'
ndvi_map <- projectRaster(ndvi_map, crs=newproj) # riproietto al sistema 3035
plot(ndvi_map)


#####################################
### DEFINISCO LA FUNZIONE ###########

sampleboost <- function(ndvi, ignorance, boundary, samp_strategy, nplot, areaplot, perm, ndvi.weight, igno.weight, dist.weight){
  normalize <- function(x) {
    return ((x - min(x)) / (max(x) - min(x)))
  } # funzione per normalizzare
  
  result<-list()
  distanze<-matrix(ncol=1, nrow = perm)
  check <- c()
  boundary <- spTransform(boundary, crs(ignorance))
  
  pb <- txtProgressBar(min = 0, max = perm, style = 3)
  for (i in 1:perm){
    punti_random <- spsample(boundary, n=nplot, type= samp_strategy, iter = 5)
    sampling_points <- as(punti_random, "data.frame")
    xy <- sampling_points[,c(1,2)]
    
    spdf <- SpatialPointsDataFrame(coords = xy, data = sampling_points,
                                   proj4string = crs(ignorance))
    
    
    spdf_buffer <- gBuffer(spdf, width=sqrt(areaplot/pi), byid=TRUE )

    # Test self intersection
    combos <- combn(nrow(spdf_buffer@data),2)
    int <- c()
    
    for(k in 1:ncol(combos)){
      ii <- combos[1,k]
      j <- combos[2,k]
      
      int[k] <- gIntersects(spdf_buffer[ii, ], spdf_buffer[j,]) # questa riga salta quando i poligoni non si intersecano
      
      }
      
      if (any(int) == TRUE) {check[[i]] <- TRUE} else {check[[i]] <- FALSE}
      
    spectral_values <- raster::extract(ndvi, spdf) # campiono i valori del raster di NDVI
    igno_values <- raster::extract(ignorance, spdf)  # campiono i valori del rater di ignoranza
    
    ## Calcolare distanze con CRS metrico
    # 1. obtain a ppp object from imported data
    m <- ppp(xy$x, xy$y, range(xy$x), range(xy$y))
    # 2. calculate Euclidean distance matrix
    pairwise_distances <- pairdist.ppp(m)
    distanze[[i]] <- sum(pairwise_distances)
    distance_values <- rep(distanze[[i]], nplot)
    
    
    estratti <- data.frame(coordinates(spdf),spectral_values, igno_values,  distance_values)
    names(estratti) <- c("x", "y", "ndvi", "ignorance", "distances")
    
    estratti$INTERSECTION <- check[[i]]
    
    result[[i]]<-data.frame(estratti)
    
    setTxtProgressBar(pb, i)
    
  }
  
  new_mat<-plyr::ldply(result, data.frame)
  new_mat$try<-as.factor(rep(1:perm, each= nplot))
  
  agg1<-aggregate(new_mat$ndvi,by=list(new_mat$try),FUN=var)
  agg_igno<-aggregate(new_mat$ignorance,by=list(new_mat$try),FUN=mean)
  
  
  agg2<-data.frame(agg1, distanze, agg_igno[[2]], unlist(check))
  colnames(agg2)<-c('Try','Variance','Mean Dist', 'Mean Ignorance', "INTERSECTION")
  agg2 <- na.omit(agg2)
  
  agg2$ndvi_score <- agg2$Variance * ndvi.weight
  agg2$igno_score <- agg2$`Mean Ignorance` * igno.weight
  agg2$spatial_score <- agg2$`Mean Dist` * dist.weight
  
  agg2$ndvi_norm <- normalize(agg2$ndvi_score)
  agg2$igno_norm <- normalize(agg2$igno_score)
  agg2$spatial_norm <- normalize(agg2$spatial_score)

  agg2$FINAL_SCORE <- agg2$ndvi_norm * agg2$igno_norm * agg2$spatial_norm
  
  agg2 <- agg2[agg2$INTERSECTION=="FALSE",] ## elimino le configurazioni dove c'è intersezione
  
  
  ordered_solutions <- agg2[order(agg2[,'FINAL_SCORE'], decreasing = TRUE),]
  Index <- as.numeric(ordered_solutions[1,1])
  sol <- subset(new_mat[new_mat$try %in% Index,])
  sol2 <- subset(agg2[agg2$Try %in% Index,])
  
  ## Plot best solution
  
  xy_out1 <- sol$Best[,c(1,2)]
  out1_points <- SpatialPointsDataFrame(coords = xy_out1, data = sol$Best, proj4string = crs(boundary))
  out1_buffers <- gBuffer(out1_points, width=sqrt(areaplot/pi), byid=TRUE )
  
  site2 <- spTransform(site, newproj)
  
   plot(site2)
   plot(out1_buffers, add=TRUE)
  
   p <- rasterVis::levelplot(ndvi, layers=1, margin = list(FUN = median))+
        latticeExtra::layer(sp.points(out1_points, lwd= 1.5, col='black'))
  
   p1 <- rasterVis::levelplot(ignorance, layers=1, margin = list(FUN = median))+
         latticeExtra::layer(sp.points(out1_points, lwd= 0.8, col='darkgray'))
  
  
   p2 <- ggplot(out1$`Full matrix`, aes(x = ndvi, group = try)) +
         geom_density(colour = "lightgrey")+
         theme(legend.position = "none")+
         geom_density(data = sol$Best, aes(x = ndvi, colour = "red"))
  
  
   
 
  return(list("Full matrix"=new_mat, "Aggregated matrix"=agg2, "Best"= sol, "Variance of sampling points"=sol2[,'Variance'],
              "Mean Ignorance" = sol2[,'igno_score'],
              "Spatial Median of Distance"= sol2[,'Mean Dist'], "Final score"= sol2[,'FINAL_SCORE'], p, p1, p2))
  
  
}

out1 <- sampleboost(ndvi = ndvi_map, ignorance = igno_map, samp_strategy='random', nplot= 10,  areaplot = 10^6, perm = 10, boundary=site,
                    ndvi.weight = 1, igno.weight=1, dist.weight=1)

out1






##### PLOTTO LA SOLUZIONE OUT1, BEST SOLUTION
xy_out1 <- out1$Best[,c(1,2)]

out1_points <- SpatialPointsDataFrame(coords = xy_out1, data = out1$Best,
                                      proj4string = CRS("+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))



plot(site)
plot(out1_points, add=TRUE)

sp::plot(site)
sp::plot(mfi2, add=TRUE)
sp::plot(out1_points, add=TRUE)

plot(mfi2)


plot(mystratification@centroids, add=TRUE, col="red")


#### Ordino per valori decrescenti di varianza con il quantile 0.99 #########
out1 <- tent6


out_filter <- na.omit(out1$`Aggregated matrix`)

out_new <- out_filter[out_filter$Variance > quantile(out_filter$Variance, quantile_threshold),]

ordered_solutions <- out_new[order(out_new[,2], decreasing = TRUE),]
ordered_solutions_2 <- ordered_solutions[order(ordered_solutions[,3], decreasing = TRUE),]

head(ordered_solutions_2, 10)


########################################################################
####### Plotto la soluzione scelta #####################################
prova <- out1$`Full matrix`[is.element(out1$`Full matrix`$try, 3313),]

xy_prova <- prova[,c(1,2)]

prova_points <- SpatialPointsDataFrame(coords = xy_prova, data = prova,
                                       proj4string = CRS("+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

plot(ndvi_clip)
#plotRGB(rgb_crop, r= 1, g= 2, b = 3, stretch = "lin")
plot(area_studio, add=TRUE)
plot(prova_points, add=TRUE, col="black")

##############################


REFERENCE_SYSTEM <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
file_export <- spTransform(prova_points, crs(REFERENCE_SYSTEM))


#### salvo il csv ---- rinominare il file

write.csv(file_export, "Sampling points_a.csv", row.names = TRUE)

saveRDS(out1, file = "Sampling points_a.rds")






a <- df %>% 
  column_to_rownames("ID") %>% #make the ID the rownames. dist will use these> NB will not work on a tibble
  dist() %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ID.x") %>% #capture the row IDs
  gather(key = ID.y, value = dist, -ID.x) %>% 
  filter(ID.x < ID.y) %>% 
  as_tibble()


a <-as.data.frame(a)
distance <- sum(a$dist)
