# Caricare gli oggetti (ma attenzione a cambiare i percorsi dei files)
# su GITHUB i layers li trovi nella cartella 'OBJECTS'


area_studio <- readOGR(dsn = 'c:/Users/MARCO/Dropbox/Dottorato/Cartografia/Layers', layer = 'Pinetamacchia')
ignorance_map <- raster("C:/Users/MARCO/Documents/Ignorance Map.tif")
ndvi_clip <- raster("C:/Users/MARCO/Dropbox/Dottorato/Cartografia/Raster remote sensing/NDVImap_SITEB.tif")

sr <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 " 
ndvi_clip <- projectRaster(ndvi_clip, crs = sr)
ignorance_map <- projectRaster(ignorance_map, crs = sr)
