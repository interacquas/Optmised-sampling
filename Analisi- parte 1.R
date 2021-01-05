library(tidyverse)
library(reshape2)
library(readxl)
library(splitstackshape)
library(vegan)
library(ggplot2)

##### carico il dataframe specie con nome df

df <- read_excel("C:/Users/MARCO/Documents/Tesi laureandi/Giuseppe Antonelli/specie-per-plot_coord.xlsx", sheet = "plots")
df <- as.data.frame(df)

##### carico il dataframe coordinate con nome df_coord

df_coord <- read_excel("C:/Users/MARCO/Documents/Tesi laureandi/Giuseppe Antonelli/specie-per-plot_coord.xlsx", 
                       sheet = "coordinates")
df_coord <- as.data.frame(df_coord)
df_coord$X <- as.numeric(df_coord$X)
df_coord$Y <- as.numeric(df_coord$Y)


df_spectra <- read_excel("C:/Users/MARCO/Documents/Tesi laureandi/Giuseppe Antonelli/specie-per-plot_coord.xlsx", sheet = "spectral")
df_spectra <- as.data.frame(df_spectra) 

#####  creo un dataframe 'presenza-assenza'

df2 <- df %>% gather(Plot, Species) %>%
  filter(!is.na(Species)) %>%
  mutate(value = 1) %>%
  dcast(Plot~Species, value.var = "value", fill = 0)


### Applico la funzione centroidpattern

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


##### Applico la funzione Spatially Explicit rarefaction
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



### PLOT THE CURVES ###

plot(seq(1:12), explicit_curve[,1], xlim= c(1, 12), ylim=c(0,max(explicit_curve[,1 ])+5),
     xlab="Number of plots", ylab="Species Richness", type="p", pch=10, ann=TRUE,col ="white", main="Rarefaction")
lines(seq(1:12),explicit_curve[,1], lwd=3, col="dark gray")
axis(side=1, at=c(1:12))
lines(seq(1:12),classic$richness, lwd =3, col ="black")
legend(8, 40, c("Classic", "SER"), lty=c(1,1),lwd=c(2,2),col =
         c("black","darkgrey"), bty="n", text.col = "black", merge = TRUE, cex=1.3, pt.cex=1)



############## diretional SAC ################

require(Rarefy)
require(vegan)

comm_matrix <- df2[,2:ncol(df2)]
mite.xy <- df_coord[,2:3]

# Spatially-explicit curves can be obtained as follows
spatialdist <- dist(mite.xy) # to calculate the geographic
spectraldist <- dist(df_spectra, method = "euclidean")

# distance between plots, i.e. the Euclidean distance # between the coordinates of the plots)
betas <- directionalSAC(comm_matrix, spatialdist) # to calculate directional
betas_spectral <- directionalSAC(comm_matrix, spectraldist) # to calculate directional spectral distance


# and non directional beta diversity
plot(1:12, betas$N_Exact, xlab="M", ylab="Species richness", ylim=range(c(betas$N_Exact,
                                                                          betas$N_SCR, betas$Alpha, mean(apply(comm_matrix, 1, function(x) length(x[x>0]))))))
points(1:12,rep( mean(apply(comm_matrix, 1, function(x) length(x[x>0]))), 12), pch=2)
points(1:12, betas$N_SCR, pch=3)
points(1:12, betas$Alpha_dir, pch=4)
legend("right", legend=c("Non-directional SAC",
                         "Non-directional alpha diversity", "Directional SAC",
                         "Directional alpha diversity"), pch=1:4)

# click on the figure, on an empty area of the figure, to place the legend.
# M is the number of plots
plot(1:12, betas$Beta_M, xlab="M", ylab="Beta diversity",
     ylim=range(c(betas$Beta_M_dir, betas$Beta_M)))
points(1:12, betas$Beta_M_dir, pch=2)
legend("right", legend=c("Non-directional beta", "Directional beta"), pch=1:2)

# click on the figure, on an empty area of the figure, to place the legend.
plot(2:12, betas$Beta_N[2:12], xlab="M", ylab="Normalized beta diversity",
     ylim=range(c(betas$Beta_N_dir[2:12], betas$Beta_N[2:12])))
points(2:12, betas$Beta_N_dir[2:12], pch=2)
legend("right", legend=c("Non-directional beta", "Directional beta"), pch=1:2)

# click on the figure, on an empty area of the figure, to place the legend.
plot(2:12, betas$Beta_Autocor[2:12], xlab="M",
     ylab="Normalized measure of autocorrelation")




##### Plot spaziale vs. spettrale


plot(seq(1:12), betas[,1], xlim= c(1, 12), ylim=c(0,max(explicit_curve[,1 ])+5),
     xlab="Number of plots", ylab="Species Richness", type="p", pch=10, ann=TRUE,col ="white", main="Rarefaction")
lines(seq(1:12),betas[,1], lwd=3, col="dark gray")
axis(side=1, at=c(1:12))
lines(seq(1:12),betas_spectral[,1], lwd =3, col ="black")
legend(8, 40, c("Spectral", "Spatial"), lty=c(1,1),lwd=c(2,2),col =
         c("black","darkgrey"), bty="n", text.col = "black", merge = TRUE, cex=1.3, pt.cex=1)










##### BASELGA


