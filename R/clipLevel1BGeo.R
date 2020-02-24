#'Clip Level1BGeo data by Coordinates
#'
#'@description This function clips GEDI level1B extracted geolocation (level1BGeo)
#' data within given bounding coordinates
#'
#'@usage clipLevel1BGeo(level1BGeo, xleft, xright, ybottom, ytop)
#'
#'@param level1BGeo A GEDI Level1b object (output of \code{\link[rGEDI:readLevel1B]{readLevel1B}} function). A S4 object of class "gedi.level1b".
#'@param xleft Numeric. West longitude (x) coordinate of bounding rectangle, in decimal degrees.
#'@param xright Numeric. East longitude (x) coordinate of bounding rectangle, in decimal degrees.
#'@param ybottom Numeric. South latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'@param ytop Numeric. North latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'
#'@return An S4 object of class \code{\link[data.table:data.table]{data.table-class}}.
#'
#'@examples
#'# specify the path to GEDI Level 1B data
#'level1bpath <- system.file("extdata", "GEDIexample_level01B.h5", package="rGEDI")
#'
#'# Reading GEDI level1B data
#'level1b <- readLevel1B(level1bpath)
#'
#'# Get GEDI level1B geolocations
#'level1BGeo<-getLevel1BGeo(level1b)
#'
#'# Bounding rectangle coordinates
#'xleft = -116.4683
#'xright = -116.5583
#'ybottom = 46.75208
#'ytop = 46.84229
#'
#'# clip by boundary box coordinates
#'level1bGeo_clip <- clipLevel1BGeo(level1BGeo,xleft, xright, ybottom, ytop)
#'
#'library(leaflet)
#'leaflet() %>%
#'  addCircleMarkers(level1bGeo_clip@dt$longitude_bin0,
#'                   level1bGeo_clip@dt$latitude_bin0,
#'                   radius = 1,
#'                   opacity = 1,
#'                   color = "red")  %>%
#'  addScaleBar(options = list(imperial = FALSE)) %>%
#'  addProviderTiles(providers$Esri.WorldImagery)
#'@export
clipLevel1BGeo = function(level1BGeo,xleft, xright, ybottom, ytop){
  # xleft ybottom xright ytop
  mask =
    level1BGeo$longitude_bin0 >= xleft &
    level1BGeo$longitude_bin0 <= xright &
    level1BGeo$latitude_bin0 >= ybottom &
    level1BGeo$latitude_bin0 <=  ytop &
    level1BGeo$longitude_lastbin >= xleft &
    level1BGeo$longitude_lastbin <= xright &
    level1BGeo$latitude_lastbin >= ybottom &
    level1BGeo$latitude_lastbin <=  ytop

  mask = (1:length(level1BGeo$longitude_bin0))[mask]
  newFile<-level1BGeo[mask,]
  #newFile<- new("gedi.level1b.dt", dt = x[mask,])
  if (nrow(newFile) == 0) {print("The polygon does not overlap the GEDI data")} else {
    return (newFile)
  }

}

#'Clip Level1BGeo data by geometry
#'
#'@description This function clips GEDI level1B extracted geolocation (level1BGeo)
#' data within given geometry
#'
#'@usage clipLevel1BGeo(level1BGeo, polygon_spdf, output)
#'
#'@param level1BGeo A GEDI Level1b object (output of \code{\link[rGEDI:readLevel1B]{readLevel1B}} function). A S4 object of class "gedi.level1b".
#'@param polygon_spdf Polygon. An object of class \code{\link[sp]{SpatialPolygonsDataFrame-class}},
#'which can be loaded as an ESRI shapefile using \code{\link[rgdal:readOGR]{readOGR}} function in the \emph{rgdal} package.
#'@param split_by Polygon id. If defined, GEDI data will be clipped by each polygon using the polygon id from table of attribute defined by the user
#'@return A S4 object of class \code{\link[data.table:data.table]{data.table-class}}.
#'
#'@examples
#'
#'# specify the path to GEDI Level 1B data
#'level1bpath <- system.file("extdata", "GEDIexample_level01B.h5", package="rGEDI")
#'
#'# Reading GEDI level1B data
#'level1b <- readLevel1B(level1bpath)
#'
#'# Get GEDI level1B geolocations
#'level1BGeo<-getLevel1BGeo(level1b)
#'
#'# specify the path to shapefile
#'polygon_filepath <- system.file("extdata", "clip_polygon.shp", package="rGEDI")
#'
#'# Reading shapefile as SpatialPolygonsDataFrame object
#'library(rgdal)
#'polygon_spdf<-readOGR(polygons_filepath)
#'
#'level1bGeo_clip = clipLevel1BGeometry(level1bGeo, polygon_spdf, split_by="id")
#'
#'library(leaflet)
#'leaflet() %>%
#'  addCircleMarkers(level1b_clip@dt$longitude_bin0,
#'                   level1b_clip@dt$latitude_bin0,
#'                   radius = 1,
#'                   opacity = 1,
#'                   color = "red")  %>%
#'  addScaleBar(options = list(imperial = FALSE)) %>%
#'  addPolygons(data=polygon_spdf,weight=1,col = 'white',
#'              opacity = 1, fillOpacity = 0) %>%
#'  addProviderTiles(providers$Esri.WorldImagery)
#'@export
clipxGeometry = function(level1BGeo, polygon_spdf, split_by="id") {
  exshp<-raster::extent(polygon_spdf)
  level1BGeo<-clipLevel2BVPM(level1BGeo, xleft=exshp[1], xright=exshp[2], ybottom=exshp[3], ytop=exshp[4])
  if (nrow(level1BGeo) == 0) {print("The polygon does not overlap the GEDI data")} else {
    points = sp::SpatialPointsDataFrame(coords=matrix(c(level1BGeo$lon_lowestmode, level1BGeo$lat_lowestmode), ncol=2),
                                        data=data.frame(id=1:length(level1BGeo$lon_lowestmode)), proj4string = polygon_spdf@proj4string)
    points(points, col="red")
    pts = raster::intersect(points, polygon_spdf)
    if (!is.null(split_by)){

      if ( any(names(polygon_spdf)==split_by)){
        mask = as.integer(pts@data$id)
        newFile<-level1BGeo[mask,]
        newFile$poly_id<-pts@data[,split_by]
      } else {stop(paste("The",split_by,"is not included in the attribute table.
                       Please check the names in the attribute table"))}

    } else {
      mask = as.integer(pts@data$id)
      newFile<-level1BGeo[mask,]
    }
    #newFile<- new("gedi.level1b.dt", dt = x2@dt[mask,])
    return (newFile)}
}



