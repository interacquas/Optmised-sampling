# Creare la mappa di ignoranza delle Toscana


library(raster)
library(sp)
library(ignobioR)

italia <- raster::getData('GADM', country='ITA', level=2) # creo lo shapefile dell'Italia
plot(italia)

toscana <- italia[italia$NAME_1 == "Toscana",] # subsetto la Toscana
plot(toscana)

livorno <- italia[italia$NAME_2 == "Livorno",] # subsetto la Toscana
plot(livorno)

wiki_final <- wiki_final(-c(1))
map <- ignorance_map(wiki_final, site=livorno, year_study=2020, excl_areas = NULL, CRS.new = 3035, tau =20, cellsize= 2000) 

