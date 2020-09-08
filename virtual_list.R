library(raster)
library(ignobioR)
library(sp)
library(readr)
library(rgdal)

data(floratus)
data(unsuitablezone)  #layer mare

wiki_final <- read_csv("wiki_final.csv")
wiki_final <- wiki_final[-c(1)]
is.data.frame(wiki_final)
head(wiki_final)

area_studio <- readOGR(dsn = 'CELLA TARGET (materiale extra)', layer = 'Cella per campionamento_WGS84')
plot(area_studio)

italia_prov <- raster::getData('GADM', country='ITA', level=2)
plot(italia_prov)

livorno <- italia_prov[italia_prov$NAME_2 == 'Livorno',]
plot(livorno)

wiki_final <- as.data.frame(wiki_final)
vfl <- virtual_list(data_flor=wiki_final,  site=livorno, excl_areas = unsuitablezone, tau=20)