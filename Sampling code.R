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
library(geosphere)
library(Rfast)


sampleboost <- function(x, ignorance, boundary, nplot, perm, quant, approach){
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
    
  }
  
  new_mat<-plyr::ldply(result, data.frame)
  new_mat$try<-as.factor(rep(1:perm, each= nplot))
  
  agg1<-aggregate(new_mat$ndvi,by=list(new_mat$try),FUN=var)
  agg_igno<-aggregate(new_mat$ignorance,by=list(new_mat$try),FUN=mean)
  agg2<-data.frame(agg1,distanze,agg_igno[[2]])
  colnames(agg2)<-c('Try','Variance','Mean Dist', 'Mean Ignorance')
  agg2 <- na.omit(agg2)
  agg3 <- agg2[agg2$Variance > quantile(agg2$Variance, quant),]
  ordered_solutions <- agg3[order(agg3[,'Mean Dist'], decreasing = TRUE),]
  best <- ordered_solutions[order(ordered_solutions[,3], decreasing = TRUE),]
  Index <- as.numeric(best[1,1])
  sol <- subset(new_mat[new_mat$try %in% Index,])
  sol2 <- subset(agg2[agg2$Try %in% Index,])
  return(list("Full matrix"=new_mat, "Aggregated matrix"=agg2, "Best"= sol, "Variance of sampling points"=sol2[,'Variance'],
              "Spatial Median of Distance"= sol2[,'Mean Dist']))
  
  ## Plot best solution
  
  xy_out1 <- out1$Best[,c(1,2)]
  out1_points <- SpatialPointsDataFrame(coords = xy_out1, data = out1$Best,
                                        proj4string = crs(boundary))
  p <- rasterVis::levelplot(x, layers=1, margin = list(FUN = median))+
    latticeExtra::layer(sp.points(out1_points, lwd= 0.8, col='darkgray'))
  
  p2 <- rasterVis::levelplot(mfi2, layers=1, margin = list(FUN = median))+
    latticeExtra::layer(sp.points(out1_points, lwd= 0.8, col='darkgray'))
  
  
  p3 <- ggplot(out1$`Full matrix`, aes(x = ndvi, group = try)) +
    geom_density(colour = "lightgrey")+
    theme(legend.position = "none")+
    geom_density(data = out1$Best, aes(x = ndvi, colour = "red"))
  
  print(p)
  print(p2)
  print(p3)
  
}

out1 <- sampleboost(x=ndvi_clip, ignorance = ndvi_clip, nplot= 9, quant = 0.99, perm = 1000, boundary=area_studio)

out1


##### PLOTTO LA SOLUZIONE OUT1, BEST SOLUTION
xy_out1 <- out1$Best[,c(1,2)]

out1_points <- SpatialPointsDataFrame(coords = xy_out1, data = out1$Best,
                                      proj4string = CRS("+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))


plot(ndvi_clip)
plot(site, add=TRUE)
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
prova <- out1$`Full matrix`[is.element(out1$`Full matrix`$try,  3313),]

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
