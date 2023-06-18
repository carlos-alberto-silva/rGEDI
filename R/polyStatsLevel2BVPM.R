#' Compute descriptive statistics of GEDI Canopy Cover and Vertical Profile Metrics
#'
#' @description Computes a Series of Statistics of GEDI derived Canopy Cover and Vertical Profile metrics
#' within a given area defined or not by a polygon
#'
#' @usage polyStatsLevel2BVPM(level2BVPM, func, id=NULL)
#'
#' @param level2BVPM A GEDI Level2BVPM object (output of [getLevel2BVPM()] function).
#' An S4 object of class "data.table".
#' @param func The function to be applied for computing the defined statistics
#' @param id A vector containing the polygon id for each GEDI observation. Default is NULL
#'
#' @return Returns an S4 object of class [data.table::data.table]
#' Containing Statistics of GEDI level2BVPM defined metrics
#'
#' @seealso \url{https://lpdaac.usgs.gov/products/gedi02_bv002/}
#'
#' @examples
#' # Specifying the path to GEDI level2B data (zip file)
#' outdir <- tempdir()
#' level2B_fp_zip <- system.file("extdata",
#'   "GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'   package = "rGEDI"
#' )
#'
#' # Unzipping GEDI level2A data
#' level2Bpath <- unzip(level2B_fp_zip, exdir = outdir)
#'
#' # Reading GEDI level2B data (h5 file)
#' level2b <- readLevel2B(level2Bpath = level2Bpath)
#'
#' # Specifying the path to shapefile
#' polygon_filepath <- system.file("extdata", "stands_cerrado.shp", package = "rGEDI")
#'
#' # Reading shapefile as SpatialPolygonsDataFrame object
#' library(raster)
#' polygon_spdf <- shapefile(polygon_filepath)
#'
#' # Extracting GEDI Canopy Cover and Vertical Profile Metrics
#' level2BVPM <- getLevel2BVPM(level2b)
#' head(level2BVPM)
#'
#' # Clipping GEDI data by geometry
#' level2BVPM_clip <- clipLevel2BVPMGeometry(level2BVPM, polygon_spdf, split_by = "id")
#'
#' # Define your own function
#' mySetOfMetrics <- function(x) {
#'   metrics <- list(
#'     min = min(x), # Min of x
#'     max = max(x), # Max of x
#'     mean = mean(x), # Mean of x
#'     sd = sd(x) # Sd of x
#'   )
#'   return(metrics)
#' }
#'
#' # Computing the max of the Total Plant Area Index
#' pai_max <- polyStatsLevel2BVPM(level2BVPM_clip, func = max(pai), id = NULL)
#' pai_max
#'
#' # Computing the max of the Total Plant Area Index stratified by polygon
#' pai_max_poly <- polyStatsLevel2BVPM(level2BVPM_clip, func = max(pai), id = "poly_id")
#' head(pai_max_poly)
#'
#' # Computing the serie of statistics of canopy cover stratified by polygon
#' cover_metrics <- polyStatsLevel2BVPM(level2BVPM_clip,
#'   func = mySetOfMetrics(cover),
#'   id = level2BVPM_clip$id
#' )
#' head(cover_metrics)
#' close(level2b)
#' @export
polyStatsLevel2BVPM <- function(level2BVPM, func, id = NULL) {
  # this code has been adapted from the grid_metrics function in lidR package (Roussel et al. 2019)
  # https://github.com/Jean-Romain/lidR/blob/master/R/grid_metrics.r

  # Add data.table operator
  `:=` <- data.table::`:=`
  call <- lazy_call(func)

  if (is.null(id)) {
    metrics <- lazy_apply_dt_call(dt = level2BVPM, call = call)
    metrics <- as.data.table(metrics)
    if (ncol(metrics) < 2) {
      colnames(metrics) <- paste0(call)[1]
    }
  } else {
    metrics <- lazy_apply_dt_call(dt = level2BVPM, call = call, group.by = paste0("by = ", id))
    if (ncol(metrics) < 3) {
      colnames(metrics)[2] <- paste0(call)[1]
    }
  }

  return(metrics)
}
