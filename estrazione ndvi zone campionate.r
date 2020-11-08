library(raster)
library(rgdal)
library(sp)
library(maptools)

zone <- readOGR(dsn = 'gis zone campionamento', layer = 'zone campionate 3035')

red <- raster("gis zone campionamento/REMOTE SENSING/LC08_L1TP_192030_20170611_20170627_01_T1_B4.tif")
green <- raster("gis zone campionamento/REMOTE SENSING/LC08_L1TP_192030_20170611_20170627_01_T1_B3.tif")
blue <- raster("gis zone campionamento/REMOTE SENSING/LC08_L1TP_192030_20170611_20170627_01_T1_B2.tif")
nir <- raster("gis zone campionamento/REMOTE SENSING/LC08_L1TP_192030_20170611_20170627_01_T1_B5.tif")

multiband <- stack(red, green, blue) # la banda green e blue in realtà serve solo per fare una mappa tipo ortofoto
plot(multiband)

plotRGB(multiband, r= 1, g= 2, b = 3, stretch = "lin")

zone_repr  <- spTransform(zone, crs(multiband))
rgb_crop <- crop(multiband, zone_repr)

ndvi <- (nir - red)/(nir+red) # per la mappa NDVI servono solo le bande rosso e infrarosso
ndvi_crop <- crop(ndvi, extent(zone_repr))
ndvi_zone <- mask(ndvi_crop, zone_repr)
plot(ndvi_zone)

writeRaster(ndvi_zone, "NDVI zone campionate", format = "GTiff")
