# Caricare gli oggetti (MARCO)


area_studio <- readOGR(dsn = 'c:/Users/MARCO/Dropbox/Dottorato/Cartografia/Layers', layer = 'Pinetamacchia')
ndvi_clip <- raster("C:/Users/MARCO/Dropbox/Dottorato/Cartografia/Raster remote sensing/NDVImap_SITEB.tif")
sr <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 " 
ndvi_clip <- projectRaster(ndvi_clip, crs = sr)