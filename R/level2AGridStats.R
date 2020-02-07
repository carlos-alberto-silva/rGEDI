#'Compute a Series of Grid Metrics
#'
#'@description This function computes a series of user-defined descriptive statistics for a LVIS dataset within
#'each grid cell
#'
#'@param level2a.dt LVIS l2 dataset; object of class \code{SpatialPointsDataFrame}
#'@param func the function to be applied to each cell
#'@param res spatial resolution for the output raster layer
#'@return Returns raster layer (s) of selected LVIS metric (s)
#'@author Carlos Alberto Silva (This function has been adapted from the grid_metrics function in lidR package. All credits to Roussel et al. 2019).
#'@examples
#'\dontrun{
#'#' LVIS level 2 file path
#'level2_filepath_zip <- system.file("extdata", "LVIS_Mondah_level2.zip", package="rLVIS")
#'unzip(level2_filepath_zip, exdir = tempdir())
#'level2_filepath <- file.path(tempdir(), "LVIS_Mondah_level2.txt")
#'
#'#' Reading LVIS level 2 file
#'level2a.dt<-readLevel2(level2path=level2_filepath,spdf=TRUE,glatlon=TRUE)
#'
#'#' Define your own function
mySetOfMetrics = function(x)
{
metrics = list(
    min =min(x), # Min of z
    max = max(x), # Max of z
    mean = mean(x), # Mean of z
    sd = sd(x)# Sd of z
  )
  return(metrics)
}
#'
#'#' Computing a serie of LVIS metrics
#'mlvis<-GridMetrics(level2a.dt=level2a.dt,func=mySetOfMetrics(ZT), res=0.0005)
#'
#'#' Computing the maximum of RH100 only
#'maxRH100<-GridMetrics(level2a.dt=level2a.dt,func=max(RH100), res=0.0005)
#'
#'#' Computing the mean of ZG only
#'ZGmean<-GridMetrics(level2a.dt=level2a.dt,func=mean(ZG), res=0.0005)
#'rasterVis::plot3D(ZGmean, col="gray")
#'}
#'@importFrom utils read.table
#'@import data.table
#'@export
#'

level2AGridStats = function(x, func=max(rh100), res = 0.5)
{
  requireNamespace("data.table")
  # this code has been adapted from the grid_metrics function in lidR package (Roussel et al. 2019)
  # https://github.com/Jean-Romain/lidR/blob/master/R/grid_metrics.r

  is_formula <- tryCatch(lazyeval::is_formula(func), error = function(e) FALSE)
  if (!is_formula) func <- lazyeval::f_capture(func)

    # Add data.table operator
  `:=` <- data.table::`:=`
  func<-lazyeval::f_interp(func)
  vars<-all.names(func)[3:length(all.names(func))]
  level2a.dt <- x[,names(x) %in% c("lon_lowestmode","lat_lowestmode",vars), with=FALSE]
  level2a.dt<-setNames(level2a.dt,c("y","x",vars))
  layout    <- raster::raster(raster::extent(level2a.dt), res=res)
  call      <- lazyeval::as_call(func)
  cells     <- raster::cellFromXY(layout, na.omit(level2a.dt[,2:1]))
  metrics   <- level2a.dt[,eval(call), by = cells]
  xy_coords <- raster::xyFromCell(layout, metrics[[1]])
  metrics[, cells := NULL]
  output.dt<-na.omit(cbind(xy_coords,metrics))
  output <- sp::SpatialPixelsDataFrame(output.dt[,1:2], output.dt[,-c(1:2)])#, proj4string = level2a.dt@proj4string)
  if (names(metrics)[1]=="V1") {
    names(output)<-all.names(func)[2]
  } else {names(output) <- names(metrics)}
  if (length(names(metrics)) > 1 ) {output<-raster::brick(output)} else {output<-raster::raster(output)}
  rm(level2a.dt)
  return(output)
}
