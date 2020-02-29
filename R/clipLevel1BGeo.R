#'Clip Level1BGeo data by Coordinates
#'
#'@description This function clips GEDI level1B extracted geolocation (level1BGeo)
#' data within given bounding coordinates
#'
#'@usage clipLevel1BGeo(level1BGeo, xmin, xmax, ymin, ymax)
#'
#'@param level1BGeo A GEDI Level1b object (output of \code{\link[rGEDI:readLevel1B]{readLevel1B}} function). A S4 object of class "gedi.level1b".
#'@param xmin Numeric. West longitude (x) coordinate of bounding rectangle, in decimal degrees.
#'@param xmax Numeric. East longitude (x) coordinate of bounding rectangle, in decimal degrees.
#'@param ymin Numeric. South latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'@param ymax Numeric. North latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'
#'@return An S4 object of class \code{\link[data.table:data.table]{data.table-class}}.
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi01_bv001/
#'
#'@examples
#'\dontrun{
#'# specify the path to download GEDI example dataset
#'outdir<-getwd()
#'
#'# downloading GEDI example dataset (zip file)
#'download.file("https://github.com/carlos-alberto-silva/rGEDI/releases/download/examples/examples.zip",destfile=outdir)
#'
#'# unzip the file
#'unzip(paste0(outdir,"\\examples.zip"))
#'
#'# specify the path to GEDI level1B data
#'level1bpath = paste0(outdir,"\\GEDI01_B_2019108080338_O01964_T05337_02_003_01_sub.h5"))
#'
#'# Reading GEDI level1B file
#'level1b<-readLevel1b(gedilevel1b)
#'
#'# Bounding rectangle coordinates
#'xmin = -44.15036
#'xmax = -44.10066
#'ymin = -13.75831
#'ymax = -13.71244
#'
#'# clip by boundary box coordinates
#'level1bGeo_clip <- clipLevel1BGeo(level1BGeo,xmin, xmax, ymin, ymax)
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
#'}
#'@export
clipLevel1BGeo = function(level1BGeo,xmin, xmax, ymin, ymax){
  # xmin ymin xmax ymax
  mask =
    level1BGeo$longitude_bin0 >= xmin &
    level1BGeo$longitude_bin0 <= xmax &
    level1BGeo$latitude_bin0 >= ymin &
    level1BGeo$latitude_bin0 <=  ymax &
    level1BGeo$longitude_lastbin >= xmin &
    level1BGeo$longitude_lastbin <= xmax &
    level1BGeo$latitude_lastbin >= ymin &
    level1BGeo$latitude_lastbin <=  ymax

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
#'@seealso https://lpdaac.usgs.gov/products/gedi01_bv001/
#'
#'@examples
#'\dontrun{
#'# specify the path to download GEDI example dataset
#'outdir<-getwd()
#'
#'# downloading GEDI example dataset (zip file)
#'download.file("https://github.com/carlos-alberto-silva/rGEDI/releases/download/examples/examples.zip",destfile=outdir)
#'
#'# unzip the file
#'unzip(paste0(outdir,"\\examples.zip"))
#'
#'# specify the path to GEDI level1B data
#'level1bpath = paste0(outdir,"\\GEDI01_B_2019108080338_O01964_T05337_02_003_01_sub.h5"))
#'
#'# Reading GEDI level1B file
#'level1b<-readLevel1b(gedilevel1b)
#'
#'# Get GEDI level1B geolocations
#'level1BGeo<-getLevel1BGeo(level1b)
#'
#'# specify the path to shapefile
#'polygon_filepath <- system.file("extdata", "stands_cerrado.shp", package="rGEDI")
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
#'}
#'@export
clipLevel1BGeoGeometry = function(level1BGeo, polygon_spdf, split_by="id") {
  exshp<-raster::extent(polygon_spdf)
  level1BGeo<-clipLevel1BGeo(level1BGeo,xmin=exshp[1], xmax=exshp[2], ymin=exshp[3], ymax=exshp[4])

  if (nrow(level1BGeo) == 0) {print("The polygon does not overlap the GEDI data")} else {
    points = sp::SpatialPointsDataFrame(coords=matrix(c(level1BGeo$longitude_bin0, level1BGeo$latitude_lastbin), ncol=2),
                                        data=data.frame(id=1:length(level1BGeo$latitude_bin0)), proj4string = polygon_spdf@proj4string)

    pts = raster::intersect(points, polygon_spdf)
    colnames(pts@data)<-c("rowids",names(polygon_spdf))

    if (!is.null(split_by)){

      if ( any(names(polygon_spdf)==split_by)){
        mask = as.integer(pts@data$rowids)
        newFile<-level1BGeo[mask,]
        newFile$poly_id<-pts@data[,split_by]
      } else {stop(paste("The",split_by,"is not included in the attribute table.
                       Please check the names in the attribute table"))}

    } else {
      mask = as.integer(pts@data$rowids)
      newFile<-level1BGeo[mask,]
    }
    #newFile<- new("gedi.level1b.dt", dt = x2@dt[mask,])
    return (newFile)}
}
