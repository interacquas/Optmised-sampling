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
df_buffer <- spTransform(df_points,CRS("+init=epsg:3035"))


sampling <- gBuffer(df_buffer, byid=TRUE, width=564)

#2 Con i punti 

sampling <- df_points

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
df2 <- df2[-c(1)]

ord1 <- specaccum(df2[sample(nrow(df2)),], "collector")
ord2 <- specaccum(df2[sample(nrow(df2)),], "collector")
ord3 <- specaccum(df2[sample(nrow(df2)),], "collector")
ord4 <- specaccum(df2[sample(nrow(df2)),], "collector")
ord5 <- specaccum(df2[sample(nrow(df2)),], "collector")
ord6 <- specaccum(df2[sample(nrow(df2)),], "collector")
ord7 <- specaccum(df2[sample(nrow(df2)),], "collector")
ord8 <- specaccum(df2[sample(nrow(df2)),], "collector")
ord9 <- specaccum(df2[sample(nrow(df2)),], "collector")
ord10 <- specaccum(df2[sample(nrow(df2)),], "collector")
ord11 <- specaccum(df2[sample(nrow(df2)),], "collector")
ord12 <- specaccum(df2[sample(nrow(df2)),], "collector")
ord13 <- specaccum(df2[sample(nrow(df2)),], "collector")
ord14 <- specaccum(df2[sample(nrow(df2)),], "collector")
ord15 <- specaccum(df2[sample(nrow(df2)),], "collector")
ord16 <- specaccum(df2[sample(nrow(df2)),], "collector")
ord17 <- specaccum(df2[sample(nrow(df2)),], "collector")
ord18 <- specaccum(df2[sample(nrow(df2)),], "collector")
ord19 <- specaccum(df2[sample(nrow(df2)),], "collector")
ord20 <- specaccum(df2[sample(nrow(df2)),], "collector")
ord21 <- specaccum(df2[sample(nrow(df2)),], "collector")
rar <- specaccum(df2, "random")

plot(ord1, ylim=c(0,110), col="red", lwd=4, main=("Accumulation"), xlab=("Number of plots"), ylab=("Species Richness"))
axis(side=1, at=c(1:12))

plot(ord1, ylim=c(0,110), col="black", lwd=1, main=("Accumulation"), xlab=("Number of plots"), ylab=("Species Richness"))
lines(ord2, col="red", lwd=4)
axis(side=1, at=c(1:12))

plot(ord1, ylim=c(0,110), col="black", lwd=1, main=("Accumulation"), xlab=("Number of plots"), ylab=("Species Richness"))
lines(ord2, col="black", lwd=1)
lines(ord3, col="red", lwd=4)
axis(side=1, at=c(1:12))

plot(ord1, ylim=c(0,110), col="black", lwd=1, main=("Accumulation"), xlab=("Number of plots"), ylab=("Species Richness"))
lines(ord2, col="black", lwd=1)
lines(ord3, col="black", lwd=1)
lines(ord4, col="red", lwd=4)
axis(side=1, at=c(1:12))

plot(ord1, ylim=c(0,110), col="black", lwd=1, main=("Accumulation"), xlab=("Number of plots"), ylab=("Species Richness"))
lines(ord2, col="black", lwd=1)
lines(ord3, col="black", lwd=1)
lines(ord4, col="black", lwd=1)
lines(ord5, col="red", lwd=4)
axis(side=1, at=c(1:12))

plot(ord1, ylim=c(0,110), col="black", lwd=1, main=("Accumulation"), xlab=("Number of plots"), ylab=("Species Richness"))
lines(ord2, col="black", lwd=1)
lines(ord3, col="black", lwd=1)
lines(ord4, col="black", lwd=1)
lines(ord5, col="black", lwd=1)
lines(ord6, col="black", lwd=1)
lines(ord7, col="black", lwd=1)
lines(ord8, col="black", lwd=1)
lines(ord9, col="black", lwd=1)
lines(ord10, col="black", lwd=1)
lines(ord11, col="black", lwd=1)
axis(side=1, at=c(1:12))

plot(ord1, ylim=c(0,110), col="black", lwd=1, main=("Accumulation"), xlab=("Number of plots"), ylab=("Species Richness"))
lines(ord2, col="black", lwd=1)
lines(ord3, col="black", lwd=1)
lines(ord4, col="black", lwd=1)
lines(ord5, col="black", lwd=1)
lines(ord6, col="black", lwd=1)
lines(ord7, col="black", lwd=1)
lines(ord8, col="black", lwd=1)
lines(ord9, col="black", lwd=1)
lines(ord10, col="black", lwd=1)
lines(ord11, col="black", lwd=1)
lines(ord12, col="black", lwd=1)
lines(ord13, col="black", lwd=1)
lines(ord14, col="black", lwd=1)
lines(ord15, col="black", lwd=1)
lines(ord16, col="black", lwd=1)
axis(side=1, at=c(1:12))

plot(ord1, ylim=c(0,110), col="black", lwd=1, main=("Accumulation"), xlab=("Number of plots"), ylab=("Species Richness"))
lines(ord2, col="black", lwd=1)
lines(ord3, col="black", lwd=1)
lines(ord4, col="black", lwd=1)
lines(ord5, col="black", lwd=1)
lines(ord6, col="black", lwd=1)
lines(ord7, col="black", lwd=1)
lines(ord8, col="black", lwd=1)
lines(ord9, col="black", lwd=1)
lines(ord10, col="black", lwd=1)
lines(ord11, col="black", lwd=1)
lines(ord12, col="black", lwd=1)
lines(ord13, col="black", lwd=1)
lines(ord14, col="black", lwd=1)
lines(ord15, col="black", lwd=1)
lines(ord16, col="black", lwd=1)
lines(ord17, col="black", lwd=1)
lines(ord18, col="black", lwd=1)
lines(ord19, col="black", lwd=1)
lines(ord20, col="black", lwd=1)
lines(ord21, col="black", lwd=1)
axis(side=1, at=c(1:12))

plot(ord1, ylim=c(0,110), col="grey", lwd=1, main=("Accumulation"), xlab=("Number of plots"), ylab=("Species Richness"))
lines(ord2, col="grey", lwd=1)
lines(ord3, col="grey", lwd=1)
lines(ord4, col="grey", lwd=1)
lines(ord5, col="grey", lwd=1)
lines(ord6, col="grey", lwd=1)
lines(ord7, col="grey", lwd=1)
lines(ord8, col="grey", lwd=1)
lines(ord9, col="grey", lwd=1)
lines(ord10, col="grey", lwd=1)
lines(ord11, col="grey", lwd=1)
lines(ord12, col="grey", lwd=1)
lines(ord13, col="grey", lwd=1)
lines(ord14, col="grey", lwd=1)
lines(ord15, col="grey", lwd=1)
lines(ord16, col="grey", lwd=1)
lines(ord17, col="grey", lwd=1)
lines(ord18, col="grey", lwd=1)
lines(ord19, col="grey", lwd=1)
lines(ord20, col="grey", lwd=1)
lines(ord21, col="grey", lwd=1)
axis(side=1, at=c(1:12))

plot(ord1, ylim=c(0,110), col="grey", lwd=1, main=("Rarefaction"), xlab=("Number of plots"), ylab=("Species Richness"))
lines(ord2, col="grey", lwd=1)
lines(ord3, col="grey", lwd=1)
lines(ord4, col="grey", lwd=1)
lines(ord5, col="grey", lwd=1)
lines(ord6, col="grey", lwd=1)
lines(ord7, col="grey", lwd=1)
lines(ord8, col="grey", lwd=1)
lines(ord9, col="grey", lwd=1)
lines(ord10, col="grey", lwd=1)
lines(ord11, col="grey", lwd=1)
lines(ord12, col="grey", lwd=1)
lines(ord13, col="grey", lwd=1)
lines(ord14, col="grey", lwd=1)
lines(ord15, col="grey", lwd=1)
lines(ord16, col="grey", lwd=1)
lines(ord17, col="grey", lwd=1)
lines(ord18, col="grey", lwd=1)
lines(ord19, col="grey", lwd=1)
lines(ord20, col="grey", lwd=1)
lines(ord21, col="grey", lwd=1)
lines(rar, col="red", lwd=4, ci=0)
axis(side=1, at=c(1:12))

