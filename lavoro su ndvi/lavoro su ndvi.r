library(raster)
library(rgdal)
library(sp)
library(maptools)

buffer_plot <- raster("C:/Users/Giuseppe Antonelli/Desktop/tirocinio/Optmised-sampling/lavoro su ndvi/NDVI buffer plot.tif")
zone_buffer <- raster("C:/Users/Giuseppe Antonelli/Desktop/tirocinio/Optmised-sampling/lavoro su ndvi/NDVI zone campionate buffer.tif")
zone <- raster("C:/Users/Giuseppe Antonelli/Desktop/tirocinio/Optmised-sampling/lavoro su ndvi/NDVI zone campionate.tif")

plot(buffer_plot)
