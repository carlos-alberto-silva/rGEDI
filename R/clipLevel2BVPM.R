#'Clip GEDI Canopy Cover and Vertical Profile Metrics by Coordinates
#'
#'@description This function clips GEDI level2B derived
#'Canopy Cover and Vertical Profile metrics a within given bounding coordinates
#'
#'
#'@usage clipLevel2BVPM(level2BVPM, xmin, xmax, ymin, ymax)
#'
#'
#'@param level2BVPM A GEDI Level2B object (output of \code{\link[rGEDI:readLevel1B]{readLevel1B}} function).
#'An S4 object of class "data.table".
#'@param xmin Numeric. West longitude (x) coordinate of the bounding rectangle, in decimal degrees.
#'@param xmax Numeric. East longitude (x) coordinate of the bounding rectangle, in decimal degrees.
#'@param ymin Numeric. South latitude (y) coordinate of the bounding rectangle, in decimal degrees.
#'@param ymax Numeric. North latitude (y) coordinate of the bounding rectangle, in decimal degrees.
#'
#'@return An S4 object of class \code{\link[data.table:data.table]{data.table-class}}
#'containing the Canopy Cover and Vertical Profile metrics.
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi02_bv001/
#'
#'@examples
#'# Specifying the path to GEDI level2B data (zip file)
#'level2B_fp_zip <- system.file("extdata",
#'                   "GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'                   package="rGEDI")
#'
#'# Unzipping GEDI level2A data
#'level2Bpath <- unzip(level2B_fp_zip,exdir = dirname(level2B_fp_zip))
#'
#'# Reading GEDI level2B data (h5 file)
#'level2b<-readLevel2B(level2Bpath=level2Bpath)
#'
#'# Extracting canopy cover and vertical profile metrics
#'level2BVPM<-getLevel2BVPM(level2b)
#'
#'# Bounding rectangle coordinates
#'xmin = -44.15036
#'xmax = -44.10066
#'ymin = -13.75831
#'ymax = -13.71244
#'
#'# Clipping level2BVPM by extent boundary box
#'level2b_clip <- clipLevel2BVPM(level2BVPM,xmin, xmax, ymin, ymax)
#'
#'\donttest{
#'library(leaflet)
#'leaflet() %>%
#'  addCircleMarkers(level2b_clip$longitude_bin0,
#'                   level2b_clip$latitude_bin0,
#'                   radius = 1,
#'                   opacity = 1,
#'                   color = "red")  %>%
#'  addScaleBar(options = list(imperial = FALSE)) %>%
#'  addProviderTiles(providers$Esri.WorldImagery)
#'}
#'
#'close(level2b)
#'@export
clipLevel2BVPM = function(level2BVPM,xmin, xmax, ymin, ymax){
  # xmin ymin xmax ymax
  mask =
    level2BVPM$longitude_bin0 >= xmin &
    level2BVPM$longitude_bin0 <= xmax &
    level2BVPM$latitude_bin0 >= ymin &
    level2BVPM$latitude_bin0 <=  ymax &
    level2BVPM$longitude_lastbin >= xmin &
    level2BVPM$longitude_lastbin <= xmax &
    level2BVPM$latitude_lastbin >= ymin &
    level2BVPM$latitude_lastbin <=  ymax

  mask = (1:length(level2BVPM$longitude_bin0))[mask]
  newFile<-level2BVPM[mask,]
  return (newFile)
}

#'Clip GEDI Canopy Cover and Vertical Profile Metrics by geometry
#'
#'@description This function clips GEDI level2B derived
#'Canopy Cover and Vertical Profile metrics within a given geometry
#'
#'@usage clipLevel2BVPMGeometry(level2BVPM, polygon_spdf, split_by)
#'
#'
#'@param level2BVPM A GEDI Level2B object (output of \code{\link[rGEDI:readLevel1B]{readLevel1B}} function).
#'An S4 object of class "gedi.level2b".
#'@param polygon_spdf Polygon. An object of class \code{\link[sp]{SpatialPolygonsDataFrame-class}},
#'which can be loaded as an ESRI shapefile using \code{\link[raster:shapefile]{raster::shapefile()}} function in the \emph{raster} package.
#'@param split_by Polygon id. If defined, GEDI data will be clipped by each polygon using the attribute specified by \code{split_by} from the attribute table.
#'
#'@return An S4 object of class \code{\link[data.table:data.table]{data.table-class}}
#'containing the Canopy Cover and Vertical Profile metrics.
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi02_bv001/
#'
#'@examples
#'# Specifying the path to GEDI level2B data (zip file)
#'level2B_fp_zip <- system.file("extdata",
#'                   "GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'                   package="rGEDI")
#'
#'# Unzipping GEDI level2A data
#'level2Bpath <- unzip(level2B_fp_zip,exdir = dirname(level2B_fp_zip))
#'
#'# Reading GEDI level2B data (h5 file)
#'level2b<-readLevel2B(level2Bpath=level2Bpath)
#'
#'# Extracting canopy cover and vertical profile metrics
#'level2BVPM<-getLevel2BVPM(level2b)
#'
#'# Specifying the path to shapefile
#'polygon_filepath <- system.file("extdata", "stands_cerrado.shp", package="rGEDI")
#'
#'# Reading shapefile as SpatialPolygonsDataFrame object
#'library(raster)
#'polygon_spdf<-shapefile(polygon_filepath)
#'
#'# Clipping level2BVPM by geometry
#'level2b_clip_geometry <- clipLevel2BVPMGeometry(level2BVPM,polygon_spdf,split_by="id")
#'
#'\donttest{
#'library(leaflet)
#'leaflet() %>%
#'  addCircleMarkers(level2b_clip_geometry$longitude_bin0,
#'                   level2b_clip_geometry$latitude_bin0,
#'                   radius = 1,
#'                   opacity = 1,
#'                   color = "red")  %>%
#'  addScaleBar(options = list(imperial = FALSE)) %>%
#'  addPolygons(data=polygon_spdf,weight=1,col = 'white',
#'              opacity = 1, fillOpacity = 0) %>%
#'  addProviderTiles(providers$Esri.WorldImagery)
#'}
#'
#'close(level2b)
#'@export
clipLevel2BVPMGeometry = function(level2BVPM, polygon_spdf, split_by=NULL) {
  exshp<-raster::extent(polygon_spdf)
  level2bdt<-clipLevel2BVPM(level2BVPM, xmin=exshp[1], xmax=exshp[2], ymin=exshp[3], ymax=exshp[4])

  if (nrow(level2bdt) == 0) {print("The polygon does not overlap the GEDI data")} else {
  points = sp::SpatialPointsDataFrame(coords=matrix(c(level2bdt$longitude_lastbin, level2bdt$latitude_lastbin), ncol=2),
                                      data=data.frame(id=1:length(level2bdt$longitude_lastbin)), proj4string = polygon_spdf@proj4string)
  pts = raster::intersect(points, polygon_spdf)
  colnames(pts@data)<-c("rowids",names(polygon_spdf))

  if (!is.null(split_by)){

    if ( any(names(polygon_spdf)==split_by)){

      mask = as.integer(pts@data$rowids)
      newFile<-level2bdt[mask,]
      newFile$poly_id<-pts@data[,split_by]
    } else {stop(paste("The",split_by,"is not included in the attribute table.
                       Please check the names in the attribute table"))}

  } else {
  mask = as.integer(pts@data$rowids)
  newFile<-level2bdt[mask,]
  }
  return (newFile)}
}



