#'Clip GEDI Level1Bdt data
#'
#'@description Clip GEDI Level1Bdt data within a given bounding coordinates
#'
#'
#'@param level1Bdt level1Bdt; S4 object of class "gedi.level1b.dt"
#'@param extent Extent object of a Spatial object .
#'@return Returns An object of class "gedi.level1b.dt"; subset of GEDI Level1B data
#'@examples
#'
#'#' GEDI level1B file path
#'level1_filepath = system.file("extdata", "lvis_level1_clip.h5", package="rLVIS")
#'
#'#' Reading LVIS level 2 file
#'level1_waveform = readLevel1b(level1Bfilepath)
#'
#'# Rectangle area for cliping
#'xmin = -116.4683
#'xmax = -116.3177
#'ymin = 46.75208
#'ymax = 46.84229
#'
#'# creating an exent object
#'ext<-extent(c(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax))
#'
#'clipped_level1Bdt = clipLevel1(level1Bdt, extent=ext)
#'
#'library(leaflet)
#'leaflet() %>%
#'  addCircleMarkers(clipped_level1Bdt@dt$longitude_bin0,
#'                   clipped_level1Bdt@dt$latitude_bin0,
#'                   radius = 1,
#'                   opacity = 1,
#'                   color = "red")  %>%
#'  addScaleBar(options = list(imperial = FALSE)) %>%
#'  addProviderTiles(providers$Esri.WorldImagery)
#'@export
clipLevel1Bdt = function(level1bdt,extent){
  level1bdt<-level1bdt@dt
  # xleft ybottom xright ytop
  mask =
    level1bdt$longitude_bin0 >= extent[1] &
    level1bdt$longitude_bin0 <= extent[2] &
    level1bdt$latitude_bin0 >= extent[3] &
    level1bdt$latitude_bin0 <= extent[4] &
    level1bdt$longitude_lastbin >= extent[1] &
    level1bdt$longitude_lastbin <= extent[2] &
    level1bdt$latitude_lastbin >= extent[3] &
    level1bdt$latitude_lastbin <= extent[4]

  mask = (1:length(level1bdt$longitude_bin0))[mask]
  newFile<- new("gedi.level1b.dt", dt = level1bdt[mask,])
  return (newFile)
}

#'Clip GEDI Level1Bdt data by geometry
#'
#'@description Clip GEDI Level1Bdt data within a given geometry area
#'
#'
#'@param level1Bdt level1Bdt; S4 object of class "gedi.level1b.dt"
#'@param polygon_spdf SpatialDataFrame. A polygon dataset for clipping the waveform
#'@return Returns An object of class "gedi.level1b.dt"; subset of GEDI Level1B data
#'@examples
#'
#'#' GEDI level1B file path
#'level1Bfilepath = system.file("extdata", "lvis_level1_clip.h5", package="rGEDI")
#'
#'#' Reading GEDI level1B file
#'level1b = readLevel1b(level1Bfilepath)
#'
#'#' Creating GEDI level1Bdt object
#'level1Bdt = level1Bdt(level1b)
#'
#'# Polgons file path
#'polygon_filepath <- system.file("extdata", "clip_polygon.shp", package="rGEDI")
#'
#'# Reading GEDI level2B file
#'polygon_spdf<-raster::shapefile(polygons_filepath)
#'
#'clipped_level1Bdt = clipLevel1Geometry(level1Bdt, polygon_spdf)
#'
#'library(leaflet)
#'leaflet() %>%
#'  addCircleMarkers(clipped_level1Bdt@dt$longitude_bin0,
#'                   clipped_level1Bdt@dt$latitude_bin0,
#'                   radius = 1,
#'                   opacity = 1,
#'                   color = "red")  %>%
#'  addScaleBar(options = list(imperial = FALSE)) %>%
#'  addPolygons(data=polygon_spdf,weight=1,col = 'white',
#'              opacity = 1, fillOpacity = 0) %>%
#'  addProviderTiles(providers$Esri.WorldImagery)
#'@export
clipLevel1BdtGeometry = function(level1bdt, polygon_spdf) {
  exshp<-extent(polygon_spdf)
  level1bdt2<-clipLevel1Bdt(level1bdt, extent=exshp)
  points = sp::SpatialPointsDataFrame(coords=matrix(c(level1bdt2@dt$longitude_bin0, level1bdt2@dt$latitude_bin0), ncol=2),
                                      data=data.frame(id=1:length(level1bdt2@dt$longitude_bin0)), proj4string = polygon_spdf@proj4string)
  pts = raster::intersect(points, polygon_spdf)
  mask = as.integer(pts@data$id)
  newFile<- new("gedi.level1b.dt", dt = level1bdt2@dt[mask,])
  return (newFile)
}



