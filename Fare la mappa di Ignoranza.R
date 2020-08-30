# Creare la mappa di ignoranza delle Toscana


library(raster)
library(sp)
library(ignobioR)

italia <- raster::getData('GADM', country='ITA', level=1) # creo lo shapefile dell'Italia
plot(italia)

toscana <- toscana[toscana$NAME_1 == "Toscana",] # subsetto la Toscana
plot(toscana)


map <- ignorance_map(wiki_final, site=toscana, year_study=2020, excl_areas = NULL, CRS.new = 3035, tau =20, cellsize= 1000)

