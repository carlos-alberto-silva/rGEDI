#'Clip level2BVPM data by Coordinates
#'
#'@description This function clips GEDI level2B-derived
#'Canopy Cover and Vertical Profile metrics within given bounding coordinates
#'
#'
#'@usage clipLevel2BVPM(level2BVPM, xleft, xright, ybottom, ytop, output="")
#'
#'
#'@param level2BVPM A GEDI Level2B object (output of \code{\link[rGEDI:readLevel1B]{readLevel1B}} function). A S4 object of class "gedi.level2b".
#'@param xleft Numeric. West longitude (x) coordinate of bounding rectangle, in decimal degrees.
#'@param xright Numeric. East longitude (x) coordinate of bounding rectangle, in decimal degrees.
#'@param ybottom Numeric. South latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'@param ytopNumeric. North latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'@param output Optional character path where to save the new hdf5file. The default stores a temporary file only.
#'
#'@return An S4 object of class "gedi.level2b".
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi02_bv001/
#'
#'@examples
#'# specify the path and data file and read it
#'level2bpath <- system.file("extdata", "GEDIexample_level02B.h5", package="rGEDI")
#'level2b <- readLevel2B(level2bpath)
#'
#'# Get canopy cover and vertical profile metrics
#'level2BVPM<-getlevel2BVPM(level2b)
#'
#'# Bounding rectangle coordinates
#'xleft = -116.4683
#'xright = -116.5583
#'ybottom = 46.75208
#'ytop = 46.84229
#'
#'# clip level2BVPM by extent boundary box
#'level2b_clip <- level2BVPM(level2BVPM,xleft, xright, ybottom, ytop)
#'
#'library(leaflet)
#'leaflet() %>%
#'  addCircleMarkers(level2b_clip@dt$longitude_bin0,
#'                   level2b_clip@dt$latitude_bin0,
#'                   radius = 1,
#'                   opacity = 1,
#'                   color = "red")  %>%
#'  addScaleBar(options = list(imperial = FALSE)) %>%
#'  addPolygons(data=polygon_spdf,weight=1,col = 'white',
#'              opacity = 1, fillOpacity = 0) %>%
#'  addProviderTiles(providers$Esri.WorldImagery)
#'@export
clipLevel2BVPM = function(level2BVPM,xleft, xright, ybottom, ytop){
  # xleft ybottom xright ytop
  mask =
    level2BVPM$lon_lowestmode >= xleft &
    level2BVPM$lon_lowestmode <= xright &
    level2BVPM$lat_lowestmode >= ybottom &
    level2BVPM$lat_lowestmode <=  ytop &
    level2BVPM$lon_lowestmode >= xleft &
    level2BVPM$lon_lowestmode <= xright &
    level2BVPM$lat_lowestmode >= ybottom &
    level2BVPM$lat_lowestmode <=  ytop

  mask = (1:length(level2BVPM$longitude_bin0))[mask]
  newFile<-level2BVPM[mask,]
  #newFile<- new("gedi.level1b.dt", dt = level1bdt[mask,])
  return (newFile)
}

#'Clip level2BVPM data by geometry
#'
#'@description This function clips GEDI level2B-derived
#'Canopy Cover and Vertical Profile metrics within a given geometry
#'
#'@usage clipLevel2BVPM(level2BVPM, polygon_spdf,split_by)
#'
#'
#'@param level2BVPM A GEDI Level2B object (output of \code{\link[rGEDI:readLevel1B]{readLevel1B}} function). A S4 object of class "gedi.level2b".
#'@param polygon_spdf Polygon. An object of class \code{\link[sp]{SpatialPolygonsDataFrame-class}},
#'which can be loaded as an ESRI shapefile using \code{\link[rgdal:readOGR]{readOGR}} function in the \emph{rgdal} package.
#'@param split_by Polygon id. If defined, GEDI data will be clipped by each polygon using the polygon id from table of attribute defined by the user
#'
#'@return An S4 object of class \code{\link[data.table:data.table]{data.table-class}}.
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi02_bv001/
#'
#'@examples
#'
#'# specify the path and data file and read it
#'level2bpath <- system.file("extdata", "GEDIexample_level02B.h5", package="rGEDI")
#'
#'#'# reading GEDI level2B data
#'level2b <- readLevel2B(level2bpath)
#'
#'# Get canopy cover and vertical profile metrics
#'level2BVPM<-getlevel2BVPM(level2b)
#'
#'# specify the path to shapefile
#'polygon_filepath <- system.file("extdata", "clip_polygon.shp", package="rGEDI")
#'
#'# Reading shapefile as SpatialPolygonsDataFrame object
#'library(rgdal)
#'polygon_spdf<-readOGR(polygons_filepath)
#'
#'# clip level2BVPM by geometry
#'level2b_clip_geometry <- clipLevel2BVPMGeometry(level2BVPM,polygon_spdf,split_by="id")
#'
#'library(leaflet)
#'leaflet() %>%
#'  addCircleMarkers(level2b_clip_geometry@dt$longitude_bin0,
#'                   level2b_clip_geometry@dt$latitude_bin0,
#'                   radius = 1,
#'                   opacity = 1,
#'                   color = "red")  %>%
#'  addScaleBar(options = list(imperial = FALSE)) %>%
#'  addPolygons(data=polygon_spdf,weight=1,col = 'white',
#'              opacity = 1, fillOpacity = 0) %>%
#'  addProviderTiles(providers$Esri.WorldImagery)
#'@export
clipLevel2BVPMGeometry = function(level2BVPM, polygon_spdf, split_by=NULL) {
  exshp<-raster::extent(polygon_spdf)
  level2bdt<-clipLevel2BVPM(level2BVPM, xleft=exshp[1], xright=exshp[2], ybottom=exshp[3], ytop=exshp[4])

  if (nrow(level2bdt) == 0) {print("The polygon does not overlap the GEDI data")} else {
  points = sp::SpatialPointsDataFrame(coords=matrix(c(level2bdt$lon_lowestmode, level2bdt$lat_lowestmode), ncol=2),
                                      data=data.frame(id=1:length(level2bdt$lon_lowestmode)), proj4string = polygon_spdf@proj4string)
  pts = raster::intersect(points, polygon_spdf)

  if (!is.null(split_by)){

    if ( any(names(polygon_spdf)==split_by)){
      mask = as.integer(pts@data$id)
      newFile<-level2bdt[mask,]
      newFile$poly_id<-pts@data[,split_by]
    } else {stop(paste("The",split_by,"is not included in the attribute table.
                       Please check the names in the attribute table"))}

  } else {
  mask = as.integer(pts@data$id)
  newFile<-level2bdt[mask,]
  }
  #newFile<- new("gedi.level1b.dt", dt = level2bdt2@dt[mask,])
  return (newFile)}
}



