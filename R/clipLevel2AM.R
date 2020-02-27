#'Clip Level2AM data by Coordinates
#'
#'@description This function clips GEDI Level2A extracted Elevation and Height Metrics (Level2AM)
#' within given bounding coordinates
#'
#'@usage clipLevel2AM(level2a, xmin, xmax, ymin, ymax, output)
#'
#'@param level2AM A GEDI Level2A object (output of \code{\link[rGEDI:readLevel2A]{readLevel2A}} function). A S4 object of class "gedi.level2a".
#'@param xmin Numeric. West longitude (x) coordinate of bounding rectangle, in decimal degrees.
#'@param xmax Numeric. East longitude (x) coordinate of bounding rectangle, in decimal degrees.
#'@param ymin Numeric. South latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'@param ymax Numeric. North latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'
#'@return An S4 object of class \code{\link[data.table:data.table]{data.table-class}}
#'containing the elevation and relative heights.
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi02_av001/
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
#'xmin = -44.15036
#'xmax = -44.10066
#'ymin = -13.75831
#'ymax = -13.71244
#'
#'# clip by extent boundary box
#'level2AM_clip <- clipLevel2AM(level2AM,xmin, xmax, ymin, ymax)
#'
#'@import hdf5r
#'@export
clipLevel2AM = function(level2AM,xmin, xmax, ymin, ymax){
  # xmin ymin xmax ymax
  mask =
    level2AM$lon_lowestmode >= xmin &
    level2AM$lon_lowestmode <= xmax &
    level2AM$lat_lowestmode >= ymin &
    level2AM$lat_lowestmode <=  ymax &
    level2AM$lon_lowestmode >= xmin &
    level2AM$lon_lowestmode <= xmax &
    level2AM$lat_lowestmode >= ymin &
    level2AM$lat_lowestmode <=  ymax

  mask = (1:length(level2AM$lat_lowestmode))[mask]
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
#'polygon_filepath <- system.file("extdata", "stands_cerrado.shp", package="rGEDI")
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
  level2adt<-clipLevel2AM(level2AM, xmin=exshp[1], xmax=exshp[2], ymin=exshp[3], ymax=exshp[4])
  if (nrow(level2adt) == 0) {print("The polygon does not overlap the GEDI data")} else {
    points = sp::SpatialPointsDataFrame(coords=matrix(c(level2adt$lon_lowestmode, level2adt$lat_lowestmode), ncol=2),
                                        data=data.frame(id=1:length(level2adt$lon_lowestmode)), proj4string = polygon_spdf@proj4string)
    pts = raster::intersect(points, polygon_spdf)
    colnames(pts@data)<-c("rowids",names(polygon_spdf))

    if (!is.null(split_by)){

      if ( any(names(polygon_spdf)==split_by)){

        mask = as.integer(pts@data$rowids)
        newFile<-level2adt[mask,]
        newFile$poly_id<-pts@data[,split_by]
      } else {stop(paste("The",split_by,"is not included in the attribute table.
                       Please check the names in the attribute table"))}

    } else {
      mask = as.integer(pts@data$rowids)
      newFile<-level2adt[mask,]
    }
    #newFile<- new("gedi.level1b.dt", dt = level2adt2@dt[mask,])
    return (newFile)}
}


