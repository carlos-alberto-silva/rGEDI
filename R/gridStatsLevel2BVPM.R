#'Compute Grids with descriptive statistics of
#'GEDI-derived Canopy Cover and Vertical Profile Metrics (Level2B)
#'
#'@description This function computes a series of user-defined descriptive statistics within
#'each grid cell for GEDI-derived Canopy Cover and Vertical Profile Metrics (Level2B)
#'
#'@usage gridStatsLevel2BVPM(level2VPM, func, res)
#'
#'@param level2BVPM A GEDI Level2AM object (output of \code{\link[rGEDI:getLevel2BVPM]{getLevel2BVPM}} function). A S4 object of class "data.table".
#'@param func the function(s) to be applied to each cell
#'@param res spatial resolution in decimal degrees for the output raster layer
#'
#'@return Returns raster layer(s) of selected GEDI Canopy Cover and Vertical Profile Metric(s)
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi02_bv001/
#'
#'@examples
#'\dontrun{
#'# specify the path to GEDI Level 2B data
#'level2bpath <- system.file("extdata", "GEDIexample_level02B.h5", package="rGEDI")
#'
#'# Reading GEDI level2B data
#'level2b <- readLevel2B(level2bpath)
#'
#'# Get GEDI-derived Canopy Cover and Vertical Profile Metrics
#'level2BVPM<-getLevel2BVPM(level2b)
#'head(level2BVPM)
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
#'#' Computing a serie of statistics of GEDI-derived canopy cover
#'cover_stats<-gridStatsLevel2BVPM(level2AM = level2AM, func=mySetOfMetrics(cover), res=0.5)
#'plot(cover_stats)
#'
#'#' Computing the max of the Total Plant Area Index only
#'pai_max<-gridStatsLevel2BVPM(level2AM = level2AM, func=max(pai), res=0.5)
#'plot(pai_max)
#'
#'#' Computing the mean of Foliage Clumping Index only
#'omega_mean<-gridStatsLevel2BVPM(level2AM = level2AM, func=mean(omega), res=0.5)
#'plot(omega_mean)
#'}
#'
gridStatsLevel2BVPM = function(level2BVPM, func, res = 0.5)
{
  requireNamespace("data.table")
  # this code has been adapted from the grid_metrics function in lidR package (Roussel et al. 2019)
  # https://github.com/Jean-Romain/lidR/blob/master/R/grid_metrics.r

  is_formula <- tryCatch(lazyeval::is_formula(func), error = function(e) FALSE)
  if (!is_formula) func <- lazyeval::f_capture(func)

  # Add data.table operator
  `:=` <- data.table::`:=`
  func<<-lazyeval::f_interp(func)
  vars<-all.names(func)[3:length(all.names(func))]
  level2a.dt <- na.omit(level2BVPM[,names(level2BVPM) %in% c("lon_lowestmode","lat_lowestmode",vars), with=FALSE])
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