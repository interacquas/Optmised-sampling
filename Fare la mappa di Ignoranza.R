library(raster)
library(sp)

italia <- raster::getData('GADM', country='ITA', level=1) # creo lo shapefile dell'Italia
plot(italia)

toscana <- toscana[toscana$NAME_1 == "Toscana",] # subsetto la Toscana
plot(toscana)


ignorance_map(wiki_final, )