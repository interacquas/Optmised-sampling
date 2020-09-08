# Creare la mappa di ignoranza delle Toscana


library(raster)
library(sp)
library(ignobioR)

italia_reg <- raster::getData('GADM', country='ITA', level=1) # creo lo shapefile dell'Italia
italia_prov <- raster::getData('GADM', country='ITA', level=2)
plot(italia_prov)

toscana <- italia_reg[italia_reg$NAME_1 == "Toscana",] # subsetto la Toscana
plot(toscana)

pisalucca<- italia_prov[italia_prov$NAME_2 == "Pisa" | italia_prov$NAME_2 == "Lucca",] # subsetto la provincia
plot(pisalucca)

wiki_final <- wiki_final[-c(1)] #tolgo la prima colonna a wiki_final
map <- ignorance_map(wiki_final, site=toscana, year_study=2020, excl_areas = NULL, CRS.new = 3035, tau =20, cellsize= 1000) 

