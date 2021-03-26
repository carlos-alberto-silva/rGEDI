#'Compute Grids with Descriptive Statistics of
#'GEDI derived Canopy Cover and Vertical Profile Metrics (Level2B)
#'
#'@description This function computes a series of user defined descriptive statistics within
#'each grid cell for GEDI derived Canopy Cover and Vertical Profile Metrics (Level2B)
#'
#'@usage gridStatsLevel2BVPM(level2BVPM, func, res)
#'
#'@param level2BVPM A GEDI Level2AM object (output of [getLevel2BVPM()] function).
#'An S4 object of class "data.table".
#'@param func The function(s) to be applied to each cell
#'@param res Spatial resolution in decimal degrees for the output raster layer
#'
#'@return Returns a raster layer(s) of selected GEDI Canopy Cover and Vertical Profile Metric(s)
#'
#'@seealso \url{https://lpdaac.usgs.gov/products/gedi02_bv002/}
#'
#'@examples
#'# specify the path to GEDI level2B data (zip file)
#'outdir = tempdir()
#'level2B_fp_zip <- system.file("extdata",
#'                   "GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'                   package="rGEDI")
#'
#'# Unzipping GEDI level2A data
#'level2Bpath <- unzip(level2B_fp_zip,exdir = outdir)
#'
#'# Reading GEDI level2B data (h5 file)
#'level2b<-readLevel2B(level2Bpath=level2Bpath)
#'
#'# Get GEDI derived Canopy Cover and Vertical Profile Metrics
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
#'#' Computing a serie of statistics of GEDI derived canopy cover
#'cover_stats<-gridStatsLevel2BVPM(level2BVPM = level2BVPM, func=mySetOfMetrics(cover), res=0.005)
#'plot(cover_stats)
#'
#'#' Computing the max of the Total Plant Area Index only
#'pai_max<-gridStatsLevel2BVPM(level2BVPM = level2BVPM, func=max(pai), res=0.005)
#'plot(pai_max)
#'
#'#' Computing the Foliage Height Diversity Index only
#'fhd_mean<-gridStatsLevel2BVPM(level2BVPM = level2BVPM, func=mean(fhd_normal), res=0.005)
#'plot(fhd_mean)
#'
#'close(level2b)
#'@export
gridStatsLevel2BVPM = function(level2BVPM, func, res)
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
  level2b.dt <- na.omit(level2BVPM[,names(level2BVPM) %in% c("longitude_lastbin","latitude_lastbin",vars), with=FALSE])
  level2b.dt<-setNames(level2b.dt,c("y","x",vars))
  layout    <- raster::raster(raster::extent(level2b.dt), res=res)
  call      <- lazyeval::as_call(func)
  cells     <- raster::cellFromXY(layout, na.omit(level2b.dt[,2:1]))
  metrics   <- with(level2b.dt, level2b.dt[,eval(call), by = cells])
  xy_coords <- raster::xyFromCell(layout, metrics[[1]])
  metrics[, cells := NULL]
  output.dt<-na.omit(cbind(xy_coords,metrics))
  output <- sp::SpatialPixelsDataFrame(output.dt[,1:2], output.dt[,-c(1:2)])#, proj4string = level2b.dt@proj4string)
  if (names(metrics)[1]=="V1") {
    names(output)<-all.names(func)[2]
  } else {names(output) <- names(metrics)}
  if (length(names(metrics)) > 1 ) {output<-raster::brick(output)} else {output<-raster::raster(output)}
  rm(level2b.dt)
  return(output)
}
