#'Clip level2BPAIProfile data by Coordinates
#'
#'@description This function clips GEDI level2B-derived
#'Plant Area Index profile within given bounding coordinates
#'
#'
#'@usage clipLevel2BPAIProfile(level2BPAIProfile, xleft, xright, ybottom, ytop, output="")
#'
#'
#'@param level2BPAIProfile A GEDI Level2B object (output of \code{\link[rGEDI:getLevel2BPAIProfile]{getLevel2BPAIProfile}} function). A S4 object of class "gedi.level2b".
#'@param xleft Numeric. West longitude (x) coordinate of bounding rectangle, in decimal degrees.
#'@param xright Numeric. East longitude (x) coordinate of bounding rectangle, in decimal degrees.
#'@param ybottom Numeric. South latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'@param ytopNumeric. North latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'@param output Optional character path where to save the new hdf5file. The default stores a temporary file only.
#'
#'@return An S4 object of class "gedi.level2b".
#'
#'@examples
#'# specify the path and data file and read it
#'level2bpath <- system.file("extdata", "GEDIexample_level02B.h5", package="rGEDI")
#'
#'# Reading GEDI level2B data
#'level2b <- readLevel2B(level2bpath)
#'
#'# Get Plant Area Index profile
#'level2BPAIProfile<-getLevel2BPAIProfile(level2b)
#'
#'# Bounding rectangle coordinates
#'xleft = -116.4683
#'xright = -116.5583
#'ybottom = 46.75208
#'ytop = 46.84229
#'
#'# clip level2BVPM by extent boundary box
#'level2b_clip <- clipLevel2BPAIProfile(level2BPAIProfile,xleft, xright, ybottom, ytop)
#'
#'@export
clipLevel2BPAIProfile = function(x,xleft, xright, ybottom, ytop){
  # xleft ybottom xright ytop
  mask =
    x$lon_lowestmode >= xleft &
    x$lon_lowestmode <= xright &
    x$lat_lowestmode >= ybottom &
    x$lat_lowestmode <=  ytop

  mask = (1:length(x$lon_lowestmode))[mask]
  newFile<-x[mask,]
  #newFile<- new("gedi.level1b.dt", dt = level2bdt[mask,])
  if (nrow(newFile) == 0) {print("The polygon does not overlap the GEDI data")} else {
    return (newFile)
  }

}

#'Clip level2BPAIProfile data by geometry
#'
#'@description This function clips GEDI level2B-derived
#'Plant Area Index profile within given geometry
#'
#'@usage clipLevel2BPAIProfileGeometry(level2BPAIProfile, polygon_spdf, output)
#'
#'@param level2BPAIProfile A GEDI Level2B object (output of \code{\link[rGEDI:getLevel2BPAIProfile]{getLevel2BPAIProfile}} function). A S4 object of class "gedi.level2b".
#'@param polygon_spdf Polygon. An object of class \code{\link[sp]{SpatialPolygonsDataFrame-class}},
#'which can be loaded as an ESRI shapefile using \code{\link[rgdal:readOGR]{readOGR}} function in the \emph{rgdal} package.
#'@param split_by Polygon id. If defined, GEDI data will be clipped by each polygon using the polygon id from table of attribute defined by the user
#'
#'@return An S4 object of class \code{\link[data.table:data.table]{data.table-class}}.
#'
#'@examples
#'
#'# specify the path and data file and read it
#'level2bpath <- system.file("extdata", "GEDIexample_level02B.h5", package="rGEDI")
#'
#'#'# reading GEDI level2B data
#'level2b <- readLevel2B(level2bpath)
#'
#'# Get Plant Area Index profile
#'level2BPAIProfile<-getLevel2BPAIProfile(level2b)
#'
#'# specify the path to shapefile
#'polygon_filepath <- system.file("extdata", "clip_polygon.shp", package="rGEDI")
#'
#'# Reading shapefile as SpatialPolygonsDataFrame object
#'library(rgdal)
#'polygon_spdf<-readOGR(polygons_filepath)
#'
#'# clip level2BPAIProfile by geometry
#'level2b_clip_geometry <- clipLevel2BPAIGeometry(level2BPAIProfile,polygon_spdf, split_by="id")
#'
#'@export
clipLevel2BPAIProfileGeometry = function(x, polygon_spdf, split_by=NULL) {
  exshp<-raster::extent(polygon_spdf)
  level2bdt<-clipLevel2BPAIProfile(x, xleft=exshp[1], xright=exshp[2], ybottom=exshp[3], ytop=exshp[4])

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



