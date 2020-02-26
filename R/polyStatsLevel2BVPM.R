#'Compute descriptive statistics of GEDI Level2BVPM-derived Metrics
#'
#'@description Computes a Series of Statistics of GEDI-derived Canopy Cover and Vertical Profile metrics (Level2BVPM)
#'for all obsercation or only those defined within a giving polygon
#'
#'@usage polyStatsLevel2BVPM(level2AM, func, id)
#'
#'@param level2BVPM A GEDI Level2BVPM object (output of \code{\link[rGEDI:getLevel2BVPM]{getLevel2BVPM}} function). A S4 object of class "data.table".
#'@param func the function to be applied for computing the defined statistics
#'@param id a vector contatining the polygon id for each GEDI observation. Defaut is NULL
#'
#'@return A S4 object of class \code{\link[data.table:data.table]{data.table-class}}
#'containting Statistics of GEDI level2BVPM defined metrics
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi02_bv001/
#'
#'@examples
#'\dontrun{
#'# specify the path to GEDI Level 2B data
#'level2bpath <- system.file("extdata", "GEDIexample_level01B.h5", package="rGEDI")
#'
#'# Reading GEDI level2B data
#'level2b <- readLevel2B(level2bpath)
#'
#'# specify the path to shapefile
#'polygon_filepath <- system.file("extdata", "clip_polygon.shp", package="rGEDI")
#'
#'# Reading shapefile as SpatialPolygonsDataFrame object
#'library(rgdal)
#'polygon_spdf<-readOGR(polygons_filepath)
#'
#'# Get GEDI Canopy Cover and Vertical Profile Metrics
#'level2BVPM<-getLevel2BVPM(level2b)
#'head(level2BVPM)
#'
#'# clipping level2BVPM by geometry
#'level2BVPM_clip = clipLevel2BVPMGeometry(level2BVPM, polygon_spdf, split_by="id")
#'
#'# Define your own function
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
#'# Computing the max of the Total Plant Area Index
#'pai_max<-polyStatsLevel2BVPM(level2BVPM_clip,func=max(pai), id=NULL)
#'pai_max
#'
#'# Computing the max of the Total Plant Area Index stratified by polygon
#'pai_max_poly<-polyStatsLevel2BVPM(level2BVPM_clip,func=max(pai), id="id")
#'head(pai_max_poly)
#'
#'# Computing the serie of statistics of Foliage Clumping Index stratified by polygon
#'omega_metrics<-polyStatsLevel2BVPM(level2BVPM_clip,func=mySetOfMetrics(omega),
#'                      id=level2BM_clip@data$id)
#'head(omega_metrics)
#'}
#'@export
polyStatsLevel2BVPM = function(level2BVPM, func=mean(pai), id = NULL)
{
  # this code has been adapted from the grid_metrics function in lidR package (Roussel et al. 2019)
  # https://github.com/Jean-Romain/lidR/blob/master/R/grid_metrics.r

  is_formula <- tryCatch(lazyeval::is_formula(func), error = function(e) FALSE)
  if (!is_formula) func <- lazyeval::f_capture(func)

  # Add data.table operator
  `:=` <- data.table::`:=`

  func<- lazyeval::f_interp(func)
  call<- lazyeval::as_call(func)

  if ( is.null(id)) {
    metrics   <- level2BVPM[, c(eval(call))]
    metrics<-data.table::data.table(metrics)
    if (ncol(metrics) < 2) {
      colnames(metrics)<-paste0(call)[1]
    }
  } else {
    metrics   <- level2BVPM[, c(eval(call)), by = id]
    if (ncol(metrics) < 3) {
      colnames(metrics)[2]<-paste0(call)[1]
    }
  }

  return(metrics)
}