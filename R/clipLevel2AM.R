#'Clip Level2AM data by Coordinates
#'
#'@description This function clips GEDI Level2A extracted Elevation and Height Metrics (Level2AM)
#' within given bounding coordinates
#'
#'@usage cliplevel2AM(level2a, xleft, xright, ybottom, ytop, output)
#'
#'@param level2AM A GEDI Level2A object (output of \code{\link[rGEDI:readLevel2A]{readLevel2A}} function). A S4 object of class "gedi.level2a".
#'@param xleft Numeric. West longitude (x) coordinate of bounding rectangle, in decimal degrees.
#'@param xright Numeric. East longitude (x) coordinate of bounding rectangle, in decimal degrees.
#'@param ybottom Numeric. South latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'@param ytop Numeric. North latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'@param output Optional character path where to save the new hdf5file. The default stores a temporary file only.
#'
#'@return An S4 object of class "gedi.level2a".
#'
#'@examples
#'#' GEDI level2A file path
#'level2afilepath = system.file("extdata", "lvis_level1_clip.h5", package="rGEDI")
#'
#'#' Reading GEDI level2A file
#'level2a = readLevel2A(level2afilepath)
#'
#'#' Extracting GEDI Elevation and Height Metrics
#'level2AM = getLevel2AM(level2a)
#'
#'# Bounding rectangle coordinates
#'xleft = -116.4683
#'xright = -116.5583
#'ybottom = 46.75208
#'ytop = 46.84229
#'
#'# clip by extent boundary box
#'level2AM_clip <- clipLevel2AM(level2AM,xleft, xright, ybottom, ytop)
#'
#'@import hdf5r
#'@export
cliplevel2AM = function(level2AM,xleft, xright, ybottom, ytop){
  # xleft ybottom xright ytop
  mask =
    level2AM$lon_lowestmode >= xleft &
    level2AM$lon_lowestmode <= xright &
    level2AM$lat_lowestmode >= ybottom &
    level2AM$lat_lowestmode <=  ytop &
    level2AM$lon_lowestmode >= xleft &
    level2AM$lon_lowestmode <= xright &
    level2AM$lat_lowestmode >= ybottom &
    level2AM$lat_lowestmode <=  ytop

  mask = (1:length(level2AM$longitude_bin0))[mask]
  newFile<-level2AM[mask,]
  #newFile<- new("gedi.level1b.dt", dt = level1bdt[mask,])
  return (newFile)
}

#'Clip Level2AM data by Coordinates
#'
#'@description This function clips GEDI Level2A extracted Elevation and Height Metrics (Level2AM)
#' within given bounding coordinates
#'
#'@usage cliplevel2AM(level2a, polygon_spdf, split_by)
#'
#'@param level2AM A GEDI Level2A object (output of \code{\link[rGEDI:readLevel2A]{readLevel2A}} function). A S4 object of class "gedi.level2a".
#'@param polygon_spdf Polygon. An object of class \code{\link[sp]{SpatialPolygonsDataFrame-class}},
#'which can be loaded as an ESRI shapefile using \code{\link[rgdal:readOGR]{readOGR}} function in the \emph{rgdal} package.
#'@param split_by Polygon id. If defined, GEDI data will be clipped by each polygon using the polygon id from table of attribute defined by the user
#'
#'@return A S4 object of class \code{\link[data.table:data.table]{data.table-class}}.
#'
#'@examples
#'
#'#' GEDI level2A file path
#'level2afilepath = system.file("extdata", "lvis_level1_clip.h5", package="rGEDI")
#'
#'#' Reading GEDI level2A file
#'level2a = readLevel2A(level2afilepath)
#'
#'#' Extracting GEDI Elevation and Height Metrics
#'level2AM = getLevel2AM(level2a)
#'
#'# specify the path to shapefile
#'polygon_filepath <- system.file("extdata", "clip_polygon.shp", package="rGEDI")
#'
#'# Reading shapefile as SpatialPolygonsDataFrame object
#'library(rgdal)
#'polygon_spdf<-readOGR(polygons_filepath)
#'
#'level2AM_clip = clipLevel2AMGeometry(level2AM, polygon_spdf, split_by="id")
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
clipLevel2AMGeometry = function(level2AM, polygon_spdf, split_by="id") {
  exshp<-raster::extent(polygon_spdf)
  level2adt<-clipLevel2AM(level2AM, xleft=exshp[1], xright=exshp[2], ybottom=exshp[3], ytop=exshp[4])
  if (nrow(level2adt) == 0) {print("The polygon does not overlap the GEDI data")} else {
    points = sp::SpatialPointsDataFrame(coords=matrix(c(level2adt$lon_lowestmode, level2adt$lat_lowestmode), ncol=2),
                                        data=data.frame(id=1:length(level2adt$lon_lowestmode)), proj4string = polygon_spdf@proj4string)
    pts = raster::intersect(points, polygon_spdf)
    if (!is.null(split_by)){

      if ( any(names(polygon_spdf)==split_by)){
        mask = as.integer(pts@data$id)
        newFile<-level2adt[mask,]
        newFile$poly_id<-pts@data[,split_by]
      } else {stop(paste("The",split_by,"is not included in the attribute table.
                       Please check the names in the attribute table"))}

    } else {
      mask = as.integer(pts@data$id)
      newFile<-level2adt[mask,]
    }
    #newFile<- new("gedi.level1b.dt", dt = level2adt2@dt[mask,])
    return (newFile)}
}


