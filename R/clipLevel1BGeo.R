#'Clip GEDI x data
#'
#'@description Clip GEDI x data within a given bounding coordinates
#'
#'
#'@param x x; S4 object of class "gedi.level1b.dt"
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
#'clipped_x = clipLevel1(x, extent=ext)
#'
#'library(leaflet)
#'leaflet() %>%
#'  addCircleMarkers(clipped_x@dt$longitude_bin0,
#'                   clipped_x@dt$latitude_bin0,
#'                   radius = 1,
#'                   opacity = 1,
#'                   color = "red")  %>%
#'  addScaleBar(options = list(imperial = FALSE)) %>%
#'  addProviderTiles(providers$Esri.WorldImagery)
#'@export
clipx = function(x,xleft, xright, ybottom, ytop){
  # xleft ybottom xright ytop
  mask =
    x$longitude_bin0 >= xleft &
    x$longitude_bin0 <= xright &
    x$latitude_bin0 >= ybottom &
    x$latitude_bin0 <=  ytop &
    x$longitude_lastbin >= xleft &
    x$longitude_lastbin <= xright &
    x$latitude_lastbin >= ybottom &
    x$latitude_lastbin <=  ytop

  mask = (1:length(x$longitude_bin0))[mask]
  newFile<-x[mask,]
  #newFile<- new("gedi.level1b.dt", dt = x[mask,])
  if (nrow(newFile) == 0) {print("The polygon does not overlap the GEDI data")} else {
    return (newFile)
  }

}

#'Clip GEDI x data by geometry
#'
#'@description Clip GEDI x data within a given geometry area
#'
#'
#'@param x x; S4 object of class "gedi.level1b.dt"
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
#'#' Creating GEDI x object
#'x = x(level1b)
#'
#'# Polgons file path
#'polygon_filepath <- system.file("extdata", "clip_polygon.shp", package="rGEDI")
#'
#'# Reading GEDI level2B file
#'polygon_spdf<-raster::shapefile(polygons_filepath)
#'
#'clipped_x = clipLevel1Geometry(x, polygon_spdf)
#'
#'library(leaflet)
#'leaflet() %>%
#'  addCircleMarkers(clipped_x@dt$longitude_bin0,
#'                   clipped_x@dt$latitude_bin0,
#'                   radius = 1,
#'                   opacity = 1,
#'                   color = "red")  %>%
#'  addScaleBar(options = list(imperial = FALSE)) %>%
#'  addPolygons(data=polygon_spdf,weight=1,col = 'white',
#'              opacity = 1, fillOpacity = 0) %>%
#'  addProviderTiles(providers$Esri.WorldImagery)
#'@export
clipxGeometry = function(x, polygon_spdf, split_by="id") {
  exshp<-raster::extent(polygon_spdf)
  x<-clipLevel2BVPM(x, xleft=exshp[1], xright=exshp[2], ybottom=exshp[3], ytop=exshp[4])
  if (nrow(x) == 0) {print("The polygon does not overlap the GEDI data")} else {
    points = sp::SpatialPointsDataFrame(coords=matrix(c(x$lon_lowestmode, x$lat_lowestmode), ncol=2),
                                        data=data.frame(id=1:length(x$lon_lowestmode)), proj4string = polygon_spdf@proj4string)
    points(points, col="red")
    pts = raster::intersect(points, polygon_spdf)
    if (!is.null(split_by)){

      if ( any(names(polygon_spdf)==split_by)){
        mask = as.integer(pts@data$id)
        newFile<-x[mask,]
        newFile$poly_id<-pts@data[,split_by]
      } else {stop(paste("The",split_by,"is not included in the attribute table.
                       Please check the names in the attribute table"))}

    } else {
      mask = as.integer(pts@data$id)
      newFile<-x[mask,]
    }
    #newFile<- new("gedi.level1b.dt", dt = x2@dt[mask,])
    return (newFile)}
}



