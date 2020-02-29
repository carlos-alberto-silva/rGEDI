#'Compute descriptive statistics of GEDI Level2A-derived Metrics
#'
#'@description Computes a Series of Statistics from GEDI-derived Elevation and Height Metrics (Level2A)
#'within defined polygon ids or entire area
#'
#'@usage polyStatsLevel2AM(level2AM, func, id)
#'
#'@param level2AM A GEDI Level2AM object (output of \code{\link[rGEDI:getLevel2AM]{getLevel2AM}} function). A S4 object of class "data.table".
#'@param func the function to be applied for computing the defined statistics
#'@param id a vector contatining the polygon id for each GEDI observation. Defaut is NULL
#'
#'@return A S4 object of class \code{\link[data.table:data.table]{data.table-class}}
#'containting Statistics of GEDI level2A defined metrics
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi02_av001/
#'
#'@examples
#'\dontrun{
#'# specify the path to download GEDI example dataset
#'outdir<-getwd()
#'
#'# downloading GEDI example dataset (zip file)
#'download.file("https://github.com/carlos-alberto-silva/rGEDI/releases/download/examples/examples.zip",destfile=paste0(outdir,"/examples.zip"))
#'
#'# unzip the file
#'unzip(paste0(outdir,"\\examples.zip"))
#'
#'# specify the path to GEDI level2A data
#'level2apath = paste0(outdir,"\\GEDI02_A_2019108080338_O01964_T05337_02_001_01_sub.h5"))
#'
#'# Reading GEDI level2A data
#'level2a<-readLevel2A(level2apath)
#'
#'# specify the path to shapefile
#'polygon_filepath <- system.file("extdata", "clip_polygon.shp", package="rGEDI")
#'
#'# Reading shapefile as SpatialPolygonsDataFrame object
#'library(rgdal)
#'polygon_spdf<-readOGR(polygons_filepath)
#'
#'# Get GEDI Eleveation and Relative Metrics (level2A)
#'level2AM<-getLevel2AM(level2a)
#'head(level2AM)
#'
#'# clipping level2a by geometry
#'level2AM_clip = clipLevel2AMGeometry(level2AM, polygon_spdf, split_by="id")
#'
#'#' Define your own function
#'mySetOfMetrics = function(x)
#'{
#'metrics = list(
#'    min =min(x), # Min of x
#'    max = max(x), # Max of x
#'    mean = mean(x), # Mean of x
#'    sd = sd(x)# Sd of x
#'  )
#'  return(metrics)
#'}
#'
#'# Computing the maximum of RH100
#'RH100max<-polyStatsLevel2AM(level2AM_clip,func=max(RH100), id=NULL)
#'
#'# Computing the maximum of RH100 stratified by polygon
#'RH100max_poly<-polyStatsLevel2AM(level2AM_clip,func=max(RH100), id=NULL)
#'
#'# Computing a serie statistics for GEDI metrics stratified by polygon
#'RH100metrics<-polyStatsLevel2AM(level2AM_clip,func=mySetOfMetrics(RH100),
#'                      id=level2AM_clip@data$id)
#'}
#'@import data.table lazyeval
#'@export
polyStatsLevel2AM = function(level2AM, func=mean(rh100), id = NULL)
{

  # this code has been adapted from the grid_metrics function in lidR package (Roussel et al. 2019)
  # https://github.com/Jean-Romain/lidR/blob/master/R/grid_metrics.r

  requireNamespace("data.table")

  is_formula <- tryCatch(lazyeval::is_formula(func), error = function(e) FALSE)
  if (!is_formula) func <- lazyeval::f_capture(func)

  # Add data.table operator
  `:=` <- data.table::`:=`

  func<- lazyeval::f_interp(func)
  call<- lazyeval::as_call(func)

  if ( is.null(id)) {
    metrics   <- with(level2AM, level2AM[, c(eval(call))])
    metrics<-data.table::data.table(metrics)
    if (ncol(metrics) < 2) {
      colnames(metrics)<-paste0(call)[1]
    }
  } else {
    metrics   <- with(level2AM[, c(eval(call)), by = id])
    if (ncol(metrics) < 3) {
      colnames(metrics)[2]<-paste0(call)[1]
    }
  }

  return(metrics)
}
