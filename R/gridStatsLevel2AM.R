#'Compute Grids with descriptive statistics of GEDI-derived Elevation and Height Metrics (Level2A)
#'
#'@description This function computes a series of user-defined descriptive statistics within
#'each grid cell for GEDI-derived Elevation and Height Metrics (Level2A)
#'
#'@usage gridStatsLevel2AM(level2AM, func, res)
#'
#'@param level2AM A GEDI Level2AM object (output of \code{\link[rGEDI:getLevel2AM]{getLevel2AM}} function). A S4 object of class "data.table".
#'@param func the function(s) to be applied to each cell
#'@param res spatial resolution in decimal degrees for the output raster layer
#'
#'@return Returns raster layer(s) of selected GEDI Elevation and Height Metric(s)
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi02_av001/
#'
#'@examples
#'\dontrun{
#'# specify the path to GEDI Level 2A data
#'level2apath <- system.file("extdata", "GEDIexample_level02A.h5", package="rGEDI")
#'
#'# Reading GEDI level2A data
#'level2a <- readLevel2A(level2apath)
#'
#'# Get GEDI-derived Elevation and Height Metrics
#'level2AM<-getLevel2AM(level2a)
#'head(level2AM)
#'
#'#' Define your own function
#'mySetOfMetrics = function(x)
#'{
#'metrics = list(
#'    min =min(x), # Min of z
#'    max = max(x), # Max of z
#'    mean = mean(x), # Mean of z
#'    sd = sd(x)# Sd of z
#'    )
#'    return(metrics)
#'}
#'
#'#' Computing a serie of GEDI metrics
#'ZTstats<-gridStatsLevel2AM(level2AM = level2AM, func=mySetOfMetrics(ZT), res=0.0005)
#'plot(ZTstats)
#'
#'#' Computing the maximum of RH100 only
#'maxRH100<-gridStatsLevel2AM(level2AM = level2AM, func=max(rh100), res=0.0005)
#'plot(maxRH100)
#'
#'#' Computing the mean of ZG only
#'ZGmean<-gridStatsLevel2AM(level2AM = level2AM, func=mean(ZG), res=0.0005)
#'plot(ZGmean)
#'}
#'
#'@export
gridStatsLevel2AM = function(level2AM, func, res = 0.5)
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
  level2AM.dt <- level2AM[,names(level2AM) %in% c("lon_lowestmode","lat_lowestmode",vars), with=FALSE]
  level2AM.dt<-setNames(level2AM.dt,c("y","x",vars))
  layout    <- raster::raster(raster::extent(level2AM.dt), res=res)
  call      <- lazyeval::as_call(func)
  cells     <- raster::cellFromXY(layout, na.omit(level2AM.dt[,2:1]))
  metrics   <- level2AM.dt[,eval(call), by = cells]
  xy_coords <- raster::xyFromCell(layout, metrics[[1]])
  metrics[, cells := NULL]
  output.dt<-na.omit(cbind(xy_coords,metrics))
  output <- sp::SpatialPixelsDataFrame(output.dt[,1:2], output.dt[,-c(1:2)])#, proj4string = level2AM.dt@proj4string)
  if (names(metrics)[1]=="V1") {
    names(output)<-all.names(func)[2]
  } else {names(output) <- names(metrics)}
  if (length(names(metrics)) > 1 ) {output<-raster::brick(output)} else {output<-raster::raster(output)}
  rm(level2AM.dt)
  return(output)
}
