library(raster)
library(rgdal)
library(sp)
library(maptools)

library(ignobioR)
data("unsuitablezone")

site <- readOGR(dsn = 'gis zone campionamento', layer = 'campionamento meno esclusione')



###########################################################################
################### CREO LA MAPPA DI IGNORANZA ############################

wiki_final <- read.csv("Wikiplantbase/wiki_final.csv")
wiki_final <- wiki_new[, 2:ncol(wiki_new)]

igno_map <- ignorance_map(wiki_final, site= site, year_study=2020, excl_areas = unsuitablezone, CRS.new = 3035, tau =20, cellsize= 5000) 


##########################################################################
################## CREO LA MAPPA SPETTRALE ###############################

red <- raster("gis zone campionamento/REMOTE SENSING/LC08_L1TP_192030_20170611_20170627_01_T1_B4.tif")
green <- raster("gis zone campionamento/REMOTE SENSING/LC08_L1TP_192030_20170611_20170627_01_T1_B3.tif")
blue <- raster("gis zone campionamento/REMOTE SENSING/LC08_L1TP_192030_20170611_20170627_01_T1_B2.tif")
nir <- raster("gis zone campionamento/REMOTE SENSING/LC08_L1TP_192030_20170611_20170627_01_T1_B5.tif")

multiband <- stack(red, green, blue) # la banda green e blue in realtà serve solo per fare una mappa tipo ortofoto
plot(multiband)

plotRGB(multiband, r= 1, g= 2, b = 3, stretch = "lin")

site_repr  <- spTransform(site, crs(multiband))
rgb_crop <- crop(multiband, site_repr)



plotRGB(rgb_crop, r= 1, g= 2, b = 3, stretch = "lin")
plot(site_repr, border="red", lines =3, add=TRUE)     

ndvi <- (nir - red)/(nir+red) # per la mappa NDVI servono solo le bande rosso e infrarosso
ndvi_crop <- crop(ndvi, extent(site_repr))
ndvi_crop2 <- mask(ndvi_crop, site_repr)
plot(ndvi_crop2)
plot(site_repr, border="red", lines =3, add=TRUE)


##### ESPORTA FILE RASTER ######

#writeRaster(ndvi_crop2, "MAPPA NDVI_28m", format = "GTiff")










