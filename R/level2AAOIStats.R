#'Statistics of LVIS Level2-derived Metrics
#'
#'@description Computes a Series of Statistics from LVIS Level2-derived Metrics
#'
#'@usage L2Stats(level2_spdf, func, id)
#'
#'@param level2_spdf LVIS l2 dataset; object of class \code{SpatialPointsDataFrame}
#'@param func the function to be applied to each cell
#'@param id a vector contatining the id for each LVIS Level2 observation. Defaut is NULL
#'@return Returns a \code{data.table} object containting Statistics of LVIS Level2-derived Metrics
#'@author Carlos Alberto Silva (This function has been adapted from the grid_metrics function in lidR package. All credits to Roussel et al. 2019).
#'@examples
#'\dontrun{
#'#' LVIS level 2 file path
#'level2_filepath_zip <- system.file("extdata", "LVIS_Mondah_level2.zip", package="rLVIS")
#'unzip(level2_filepath_zip, exdir = tempdir())
#'level2_filepath <- file.path(tempdir(), "LVIS_Mondah_level2.txt")
#'
#'#' Polygons file path
#'polygons_filepath <- system.file("extdata", "LVIS_Mondah_polygons.shp", package="rLVIS")
#'
#'#' Reading LVIS level 2 file
#'level2_spdf<-readLevel2(level2path=level2_filepath,spdf=TRUE,glatlon=TRUE)
#'
#'#' Plot LVIS Level2 data
#'plotLevel2(level2_spdf=level2_spdf, color = "RH100",
#'           colorPalette = c("blue","green","yellow","red"))
#'
#'#' Reading Polygons
#'library(rgdal)
#'plots<-readOGR(polygons_filepath)
#'proj4string(plots) <- CRS("+proj=longlat +datum=WGS84")
#'plot(plots, add=TRUE, border="black", lwd=2)
#'
#'#'Clipping LVIS Level2 data
#'level2_spdf_sub<-clipLevel2(level2_spdf=level2_spdf,polygon_spdf=plots)
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
#'#'Computing single LVIS metrics
#'RH100max<-level2AAOIStats(x=level2_spdf,func=~max(RH100), id=NULL)
#'
#'#'Computing LVIS metrics by id
#'RH100metrics<-level2AAOIStats(level2_spdf=level2_spdf_sub,func=~mySetOfMetrics(RH100),
#'                      id=level2_spdf_sub@data$CLIPID)
#'}
#'@export
level2AAOIStats = function(x, func=mean(rh100), id = NULL)
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
    metrics   <- x[, c(eval(call))]
    metrics<-data.table::data.table(metrics)
    if (ncol(metrics) < 2) {
      colnames(metrics)<-paste0(call)[1]
    }
  } else {
    metrics   <- x[, c(eval(call)), by = id]
    if (ncol(metrics) < 3) {
      colnames(metrics)[2]<-paste0(call)[1]
    }
  }

  return(metrics)
}
