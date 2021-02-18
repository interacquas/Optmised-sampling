library(tidyverse)
library(reshape2)
library(readxl)
library(splitstackshape)
library(vegan)
library(ggplot2)
library(raster)
library(rgdal)
library(sp)
library(maptools)
library(rgeos)

### Dati spettrali

zone <- raster("gis zone campionamento/MAPPA NDVI_28m.tif") #raster ndvi


##### carico il dataframe presenza/assenza specie con nome df

df <- read_excel("specie-per-plot_coord.xlsx", sheet = "plots")
df <- as.data.frame(df)


##### Unità di campionamento #####

##### carico il dataframe coordinate con nome df_coord

df_coord <- read_excel("specie-per-plot_coord.xlsx", 
                       sheet = "coordinates")
df_coord <- as.data.frame(df_coord)
df_coord$X <- as.numeric(df_coord$X)
df_coord$Y <- as.numeric(df_coord$Y)


#1  Con i buffer

xy <- df_coord[,3:2]
df_points <- SpatialPointsDataFrame(coords = xy, data = df_coord, proj4string = CRS("+init=epsg:4326"))
df_buffer <- spTransform(df_points,CRS("+init=epsg:3035")) #df_buffer?


sampling <- gBuffer(df_buffer, byid=TRUE, width=564)

#2 Con i punti 

sampling <- df_points

#3 Con le tracce di campionamento

track1 <- readOGR('lavoro su ndvi/01/buffer01.shp')
track2 <- readOGR('lavoro su ndvi/02/buffer02.shp')
track3 <- readOGR('lavoro su ndvi/03/buffer03.shp')
track4 <- readOGR('lavoro su ndvi/04/buffer04.shp')
track5 <- readOGR('lavoro su ndvi/05/buffer05.shp')
track6 <- readOGR('lavoro su ndvi/06/buffer06.shp')
track7 <- readOGR('lavoro su ndvi/07/buffer07.shp')
track8 <- readOGR('lavoro su ndvi/08/buffer08.shp')
track9 <- readOGR('lavoro su ndvi/09/buffer09.shp')
track10 <- readOGR('lavoro su ndvi/10/buffer10.shp')
track11 <- readOGR('lavoro su ndvi/11/buffer11.shp')
track12 <- readOGR('lavoro su ndvi/12/buffer12.shp')

sampling <- rbind(track1, track2, track3, track4, track5, track6, track7,
                  track8,track9,track10,track11, track12)

###### 
# Estraggo i valori

values <- extract(zone, sampling, na.rm =TRUE, fun = median) # faccio il sampling del raster

df_spectra <- as.data.frame(cbind(df_coord$ID, values))
names(df_spectra) <- c("ID", "ndvi")
df_spectra$ndvi<- as.numeric(df_spectra$ndvi)


#####  creo un dataframe 'presenza-assenza' (cioè 0/1)

df2 <- df %>% gather(Plot, Species) %>%
  filter(!is.na(Species)) %>%
  mutate(value = 1) %>%
  dcast(Plot~Species, value.var = "value", fill = 0)

# create a vector with letters in the desired order
x <- c("plot1", "plot2", "plot3", "plot4", "plot5", "plot6", "plot7", "plot8", "plot9", "plot10", "plot11", "plot12")

df2 <- df2 %>% slice(match(x, Plot))

plot_list <- df2$Plot
df_matrix <- df2[,2:ncol(df2)]

df_matrix[df_matrix >=1] <- 1
df2 <- cbind(plot_list, df_matrix)


### Applico la funzione centroidpattern (non servono)

centroid_pattern <- function (z, longlat) {
  
  centroid <- function (x) {
    return(c(X = mean(x$X), Y =mean(x$Y)))
  }
  
  
  dist_1 <- function (p1, p2) {
    return(sqrt((p2['X'] - p1['X'])^2+ (p2['Y'] - p1['Y'])^2))
  }
  
  require(sp)
  if (longlat== TRUE) {  
    xy <- z[,c(2,3)]
    colnames(xy) <- c("X", "Y")
    
    centroiddata <- SpatialPointsDataFrame(coords = xy, data = z, proj4string = CRS("+init=epsg:4326"))
    centroiddata <- spTransform(centroiddata,CRS("+init=epsg:3035"))
    
    centroiddata <- cbind(as.character(centroiddata$ID), centroiddata@coords)
    centroiddata <- as.data.frame(centroiddata)
    
    centroiddata$X <- as.numeric(as.numeric(as.character(centroiddata$X)))
    centroiddata$Y <- as.numeric(as.numeric(as.character(centroiddata$Y)))
    colnames(centroiddata) <- c("ID", "X", "Y")
    z <- centroiddata
  }
  
  supersequence <- c()
  for(j in 1:nrow(z)) {
    df_filter <- c(j)
    df_nofilter <- (1:nrow(z))[-df_filter]
    df_centre <- centroid(z[df_filter,])
    distance_vector <- c()
    sequence_vector <- as.character(z$ID[df_filter])
    while(length(df_nofilter>0)) {
      for (i in 1:length(df_nofilter)) {
        distance_vector <- c(distance_vector, dist_1(df_centre, z[df_nofilter[i],]))
        distance_vector <- unlist(distance_vector, use.names=FALSE)
        
      }
      print(paste("Minimum distance",min(distance_vector)))
      md <- which(distance_vector==min(distance_vector))
      print(paste("Index of distance_vector:",md))
      print(paste("Corresponding to pointdata Index",df_nofilter[md]))
      print(paste("Nearest point:",z$ID[df_nofilter[md]]))
      sequence_vector <- c(sequence_vector,as.character(z$ID[df_nofilter[md]]))
      distance_vector <- c()
      df_filter <- c(df_filter,df_nofilter[md])
      df_nofilter <- (1:nrow(z))[-df_filter]
      df_centre <- centroid(z[df_filter,])
      
    }
    supersequence <- c(supersequence,sequence_vector)
  }
  seqmatrix <- matrix(supersequence,ncol = nrow(z),byrow = TRUE)
  seqmatrix <- as.data.frame(seqmatrix)
  colnames(seqmatrix) <- c("Plot", 1:(nrow(z)-1))
  return(seqmatrix)
}
mxp_all <- centroid_pattern(df_coord, longlat = TRUE)


##### Applico la funzione Spatially Explicit rarefaction (non serve)
SCR <- function(community, spatial_order) {
  library(vegan)
  f <- nrow(spatial_order)
  n <- ncol(spatial_order)
  change <- names(community)
  change[1] <- "com"
  colnames(community) <- change
  result <- array(dim = c(f , n))
  for(i in 1:n) {
    frame<-data.frame(spatial_order[,i ])
    colnames(frame) <- "ordered"
    agg <- merge(frame, community, by.x="ordered", by.y="com", sort=FALSE)
    c <- specaccum( agg[,2 :ncol(agg)], method="collector")
    result[,i ] <- c$richness
  }
  average <- rowMeans(result)
  IC_plus <- average + (1.96*(sd(t(result))/sqrt(n )))
  IC_neg <- average - (1.96*(sd(t(result))/sqrt(n )))
  SCR <- data.frame(as.matrix(average), IC_neg,IC_plus)
  names(SCR) <- c("SCR ", "95 %IC Negative", "95%IC Positive")
  return(SCR)
}

classic <-specaccum(df2[2:ncol(df2 )], method="exact")
explicit_curve <-SCR(df2, t(mxp_all))


############################################################################################
##### ORDINAMENTO PLOTS PER UN CRITERIO (SPETTRALI O SPAZIALI), metodo mio

library(rdist)

centroid_pattern_max <- function (z, longlat) {
  
  centroid <- function (x) {
    return(c(X = mean(x$X), Y =mean(x$Y)))
  }
  
  
  dist_1 <- function (p1, p2) {
    return(sqrt((p2['X'] - p1['X'])^2+ (p2['Y'] - p1['Y'])^2))
  }
  
  require(sp)
  if (longlat== TRUE) {  
    xy <- z[,c(2,3)]
    colnames(xy) <- c("X", "Y")
    
    centroiddata <- SpatialPointsDataFrame(coords = xy, data = z, proj4string = CRS("+init=epsg:4326"))
    centroiddata <- spTransform(centroiddata,CRS("+init=epsg:3035"))
    
    centroiddata <- cbind(as.character(centroiddata$ID), centroiddata@coords)
    centroiddata <- as.data.frame(centroiddata)
    
    centroiddata$X <- as.numeric(as.numeric(as.character(centroiddata$X)))
    centroiddata$Y <- as.numeric(as.numeric(as.character(centroiddata$Y)))
    colnames(centroiddata) <- c("ID", "X", "Y")
    z <- centroiddata
  }
  
  supersequence <- c()
  for(j in 1:nrow(z)) {
    df_filter <- c(j)
    df_nofilter <- (1:nrow(z))[-df_filter]
    df_centre <- centroid(z[df_filter,])
    distance_vector <- c()
    sequence_vector <- as.character(z$ID[df_filter])
    while(length(df_nofilter>0)) {
      for (i in 1:length(df_nofilter)) {
        distance_vector <- c(distance_vector, dist_1(df_centre, z[df_nofilter[i],]))
        distance_vector <- unlist(distance_vector, use.names=FALSE)
        
      }
      print(paste("Maximum distance",max(distance_vector)))
      md <- which(distance_vector==max(distance_vector))
      print(paste("Index of distance_vector:",md))
      print(paste("Corresponding to pointdata Index",df_nofilter[md]))
      print(paste("Nearest point:",z$ID[df_nofilter[md]]))
      sequence_vector <- c(sequence_vector,as.character(z$ID[df_nofilter[md]]))
      distance_vector <- c()
      df_filter <- c(df_filter,df_nofilter[md])
      df_nofilter <- (1:nrow(z))[-df_filter]
      df_centre <- centroid(z[df_filter,])
      
    }
    supersequence <- c(supersequence,sequence_vector)
  }
  seqmatrix <- matrix(supersequence,ncol = nrow(z),byrow = TRUE)
  seqmatrix <- as.data.frame(seqmatrix)
  colnames(seqmatrix) <- c("Plot", 1:(nrow(z)-1))
  return(seqmatrix)
}

mxp_all_2 <- centroid_pattern_max(df_coord, longlat = TRUE)
explicit_curve_max <- SCR(df2, t(mxp_all_2))


spectral_dist<-vegdist(df_spectra[,2], upper=TRUE, method="euclidean",diag=TRUE)
yy <- as.matrix(spectral_dist)
sampling_order <- data.frame()

for(j in 1:ncol(yy)) {
  
  VEC <- yy[,j]
  
  
for(i in 1:(nrow(yy)-1)) { 
 
  SAMP <- max(VEC, na.rm = TRUE)
  check <- match(SAMP, VEC)
  sampling_order[j, i] <- check
  VEC[sampling_order[j, i]] <- NA # questa riga deve rimanere uguale
 
 
} 
  print(j)
} ### faccio la matrice per la rarefazione spettralmente esplicita!

sampling_order<- cbind(Plot = 1:12, sampling_order)
sampling_order[] <- lapply(sampling_order, function(x) paste("plot", x, sep=""))
names(sampling_order) <- c("Plot", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")

spectral_curve_max <-SCR(df2, t(sampling_order))
explicit_curve_max <-SCR(df2, t(mxp_all_2))


#Faccio il grafico
plot(seq(1:12), spectral_curve_max[,1], xlim= c(1, 12), ylim=c(0,max(explicit_curve[,1 ])+1),
     xlab="Number of plots", ylab="Species Richness", type="p", pch=10, ann=TRUE,col ="white", main="Rarefaction (maximised distances)")
lines(seq(1:12),spectral_curve_max[,1], lwd=3, col="red")
axis(side=1, at=c(1:12))
lines(seq(1:12),explicit_curve_max[,1], lwd =3, col ="black")
#lines(seq(1:12),classic$richness, lwd =3, col ="dark gray")
lines(seq(1:12),classic$richness + classic$sd, lwd =1, lty=2, col ="dark gray")
lines(seq(1:12),classic$richness - classic$sd, lwd =1, lty=2, col ="dark gray")
legend(6, 40, c("Spatial (12 comb.)", "Spectral (12 comb.)", "Random (479.001.600 comb.)"), lty=c(1,1),lwd=c(2,2),col =
         c("black","red", "dark grey"), bty="n", text.col = "black", merge = TRUE, cex=1.3, pt.cex=1)



##########  Beta diversity Baselga

library(betapart)


# get betapart objects

ceram.s.core <- betapart.core(df2[,2:ncol(df2)])

# multiple site measures

ceram.s.multi <- beta.multi(ceram.s.core)


# sampling across equal sites

ceram.s.samp <- beta.sample(ceram.s.core, sites=12, samples=100)

# plotting the distributions of components

dist.s <- ceram.s.samp$sampled.values

plot(density(dist.s$beta.SOR), xlim=c(0,0.8), ylim=c(0, 27), xlab='Beta diversity', main='', lwd=3)

lines(density(dist.s$beta.SNE), lty=1, lwd=2)

lines(density(dist.s$beta.SIM), lty=2, lwd=2)


# pairwise 

pair.s <- beta.pair(df2[,2:ncol(df2)])

# plotting clusters

dist.s <- ceram.s.samp$sampled.values

plot(hclust(pair.s$beta.sim, method='average'), hang=-1, main='', sub='', xlab='')
title(xlab=expression(beta[sim]), line=0.3)

plot(hclust(pair.s$beta.sne, method='average'), hang=-1, main='', sub='', xlab='')
title(xlab=expression(beta[sne]), line=0.3)


shapiro.test(spatialdist)
shapiro.test(pair.s$beta.sim)
shapiro.test(spectral_dist2)

library(ade4)
mt1 <- mantel.randtest(spatialdist,pair.s$beta.sim,nrepet=10^6)
plot(mt1)

mt2 <- mantel.randtest(spectral_dist, pair.s$beta.sim, nrepet=10^6)
plot(mt2)


### Plotto tutte le singole 12 curve di accumulo dello spettrale

spec1 <- specaccum(comm_matrix[gsub("plot", "", sampling_order[1,]),], method="collector")
spec2 <- specaccum(comm_matrix[gsub("plot", "", sampling_order[2,]),], method="collector")
spec3 <- specaccum(comm_matrix[gsub("plot", "", sampling_order[3,]),], method="collector")
spec4 <- specaccum(comm_matrix[gsub("plot", "", sampling_order[4,]),], method="collector")
spec5 <- specaccum(comm_matrix[gsub("plot", "", sampling_order[5,]),], method="collector")
spec6 <- specaccum(comm_matrix[gsub("plot", "", sampling_order[6,]),], method="collector")
spec7 <- specaccum(comm_matrix[gsub("plot", "", sampling_order[7,]),], method="collector")
spec8 <- specaccum(comm_matrix[gsub("plot", "", sampling_order[8,]),], method="collector")
spec9 <- specaccum(comm_matrix[gsub("plot", "", sampling_order[9,]),], method="collector")
spec10 <- specaccum(comm_matrix[gsub("plot", "", sampling_order[10,]),], method="collector")
spec11 <- specaccum(comm_matrix[gsub("plot", "", sampling_order[11,]),], method="collector")
spec12 <- specaccum(comm_matrix[gsub("plot", "", sampling_order[12,]),], method="collector")


SCR_values <- as.data.frame(cbind(spec1$richness, spec2$richness, spec3$richness, spec4$richness, spec5$richness, 
                    spec6$richness, spec7$richness, spec8$richness, spec9$richness,
                    spec10$richness, spec11$richness,spec11$richness ))
SCR_values$mean <- apply(SCR_values, 1, mean)



plot(seq(1:12), spec1$richness, xlim= c(1, 12), ylim=c(0,max(explicit_curve[,1 ])+1),
     xlab="Number of plots", ylab="Species Richness", type="p", pch=10, ann=TRUE,col ="white", main="Spectral Rarefaction")
lines(seq(1:12),spec1$richness, lwd=1, col="blue")
lines(seq(1:12),spec2$richness, lwd=1, col="blue")
lines(seq(1:12),spec3$richness, lwd=1, col="blue")
lines(seq(1:12),spec4$richness, lwd=1, col="blue")
lines(seq(1:12),spec5$richness, lwd=1, col="blue")
lines(seq(1:12),spec6$richness, lwd=1, col="blue")
lines(seq(1:12),spec7$richness, lwd=1, col="blue")
lines(seq(1:12),spec8$richness, lwd=1, col="blue")
lines(seq(1:12),spec9$richness, lwd=1, col="blue")
lines(seq(1:12),spec10$richness, lwd=1, col="blue")
lines(seq(1:12),spec11$richness, lwd=1, col="blue")
lines(seq(1:12),spec12$richness, lwd=1, col="blue")
lines(seq(1:12),spectral_curve_max[,1], lwd=2, lty=2, col="red")
axis(side=1, at=c(1:12))
#legend(8, 40, c("Classic", "SER", "Spectral"), lty=c(1,1),lwd=c(2,2),col =
         #c("black","blue", "red"), bty="n", text.col = "black", merge = TRUE, cex=1.3, pt.cex=1)


### Plotto tutte le singole 12 curve di accumulo dello spaziale

spaz1 <- specaccum(comm_matrix[gsub("plot", "", mxp_all_2[1,]),], method="collector")
spaz2 <- specaccum(comm_matrix[gsub("plot", "", mxp_all_2[2,]),], method="collector")
spaz3 <- specaccum(comm_matrix[gsub("plot", "", mxp_all_2[3,]),], method="collector")
spaz4 <- specaccum(comm_matrix[gsub("plot", "", mxp_all_2[4,]),], method="collector")
spaz5 <- specaccum(comm_matrix[gsub("plot", "", mxp_all_2[5,]),], method="collector")
spaz6 <- specaccum(comm_matrix[gsub("plot", "", mxp_all_2[6,]),], method="collector")
spaz7 <- specaccum(comm_matrix[gsub("plot", "", mxp_all_2[7,]),], method="collector")
spaz8 <- specaccum(comm_matrix[gsub("plot", "", mxp_all_2[8,]),], method="collector")
spaz9 <- specaccum(comm_matrix[gsub("plot", "", mxp_all_2[9,]),], method="collector")
spaz10 <- specaccum(comm_matrix[gsub("plot", "", mxp_all_2[10,]),], method="collector")
spaz11 <- specaccum(comm_matrix[gsub("plot", "", mxp_all_2[11,]),], method="collector")
spaz12 <- specaccum(comm_matrix[gsub("plot", "", mxp_all_2[12,]),], method="collector")


SCR_values_spaz <- as.data.frame(cbind(spaz1$richness, spaz2$richness, spaz3$richness, spaz4$richness, spaz5$richness, 
                                  spaz6$richness, spaz7$richness, spaz8$richness, spaz9$richness,
                                  spaz10$richness, spaz11$richness, spaz12$richness))
SCR_values_spaz$mean <- apply(SCR_values_spaz, 1, mean)



plot(seq(1:12), spec1$richness, xlim= c(1, 12), ylim=c(0,max(explicit_curve[,1 ])+1),
     xlab="Number of plots", ylab="Species Richness", type="p", pch=10, ann=TRUE,col ="white", main="Spatial Rarefaction")
lines(seq(1:12),spaz1$richness, lwd=1, col="blue")
lines(seq(1:12),spaz2$richness, lwd=1, col="blue")
lines(seq(1:12),spaz3$richness, lwd=1, col="blue")
lines(seq(1:12),spaz4$richness, lwd=1, col="blue")
lines(seq(1:12),spaz5$richness, lwd=1, col="blue")
lines(seq(1:12),spaz6$richness, lwd=1, col="blue")
lines(seq(1:12),spaz7$richness, lwd=1, col="blue")
lines(seq(1:12),spaz8$richness, lwd=1, col="blue")
lines(seq(1:12),spaz9$richness, lwd=1, col="blue")
lines(seq(1:12),spaz10$richness, lwd=1, col="blue")
lines(seq(1:12),spaz11$richness, lwd=1, col="blue")
lines(seq(1:12),spaz12$richness, lwd=1, col="blue")
lines(seq(1:12),explicit_curve_max[,1], lwd=2, lty=2, col="red")
axis(side=1, at=c(1:12))
#legend(8, 40, c("Classic", "SER", "Spectral"), lty=c(1,1),lwd=c(2,2),col =
#c("black","blue", "red"), bty="n", text.col = "black", merge = TRUE, cex=1.3, pt.cex=1)

