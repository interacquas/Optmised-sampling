library(raster)
library(ignobioR)
library(sp)
library(readr)

data(floratus)
data(unsuitablezone)  #layer mare

wiki_final <- read_csv("wiki_final.csv")
wiki_final <- wiki_final[-c(1)]
head(wiki_final)

area_studio <- readOGR(dsn = 'CELLA TARGET (materiale extra)', layer = 'Cella per campionamento_WGS84')
plot(area_studio)

italia_prov <- raster::getData('GADM', country='ITA', level=2)
plot(italia_prov)

livorno <- italia_prov[italia_prov$NAME_2 == 'Livorno',]
plot(livorno)

vfl <- virtual_list(data_flor=wiki_final, excl_areas = livorno, site=area_studio, tau=20)
