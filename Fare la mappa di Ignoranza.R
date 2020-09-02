# Creare la mappa di ignoranza delle Toscana


library(raster)
library(sp)
library(ignobioR)

italia <- raster::getData('GADM', country='ITA', level=2) # creo lo shapefile dell'Italia
plot(italia)

toscana <- italia[italia$NAME_1 == "Toscana",] # subsetto la Toscana
plot(toscana)

pisalucca <- italia[italia$NAME_2 == "Pisa" | italia$NAME_2 == "Lucca",] # subsetto la Toscana
plot(pisalucca)

map <- ignorance_map(wiki_final[1:1000,], site=pisalucca, year_study=2020, excl_areas = NULL, CRS.new = 3035, tau =20, cellsize= 5000) 

