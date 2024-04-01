#' Compute Grids with Descriptive Statistics of
#' GEDI derived Canopy Cover and Vertical Profile Metrics (Level2B)
#'
#' @description This function computes a series of user defined descriptive statistics within
#' each grid cell for GEDI derived Canopy Cover and Vertical Profile Metrics (Level2B)
#'
#' @param level2BVPM A GEDI Level2AM object (output of [getLevel2BVPM()] function).
#' An S4 object of class "data.table".
#' @param func The function(s) to be applied to each cell
#' @param res Spatial resolution in decimal degrees for the output stars raster layer
#'
#' @return Returns a stars raster layer(s) of selected GEDI Canopy Cover and Vertical Profile Metric(s)
#'
#' @seealso \url{https://lpdaac.usgs.gov/products/gedi02_bv002/}
#'
#' @examples
#' # specify the path to GEDI level2B data (zip file)
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
#' # Get GEDI derived Canopy Cover and Vertical Profile Metrics
#' level2BVPM <- getLevel2BVPM(level2b)
#' head(level2BVPM)
#'
#' #' Define your own function
#' mySetOfMetrics <- function(x) {
#'   metrics <- list(
#'     min = min(x), # Min of z
#'     max = max(x), # Max of z
#'     mean = mean(x), # Mean of z
#'     sd = sd(x) # Sd of z
#'   )
#'   return(metrics)
#' }
#'
#' #' Computing a serie of statistics of GEDI derived canopy cover
#' cover_stats <- gridStatsLevel2BVPM(
#'   level2BVPM = level2BVPM,
#'   func = mySetOfMetrics(cover),
#'   res = 0.005
#' )
#' plot(cover_stats)
#'
#' #' Computing the max of the Total Plant Area Index only
#' pai_max <- gridStatsLevel2BVPM(level2BVPM = level2BVPM, func = max(pai), res = 0.005)
#' plot(pai_max)
#'
#' #' Computing the Foliage Height Diversity Index only
#' fhd_mean <- gridStatsLevel2BVPM(level2BVPM = level2BVPM, func = mean(fhd_normal), res = 0.005)
#' plot(fhd_mean)
#'
#' close(level2b)
#' @export
gridStatsLevel2BVPM <- function(level2BVPM, func, res) {
  requireNamespace("data.table")
  cells <- NA
  # this code has been adapted from the grid_metrics function in lidR package (Roussel et al. 2019)
  # https://github.com/Jean-Romain/lidR/blob/master/R/grid_metrics.r

  # Add data.table operator
  `:=` <- data.table::`:=`
  `%>%` <- sf::`%>%`

  call <- lazy_call(func)


  sf <- sf::st_as_sf(level2BVPM, coords = c("longitude_bin0", "latitude_bin0"))
  bbox <- terra::ext(sf)
  layout <- terra::rast(bbox, resolution = res, vals = NA, crs = "epsg:4326")

  level2BVPM[, cells := terra::cells(layout, terra::vect(sf))[, 2]]
  metrics <- lazy_apply_dt_call(level2BVPM, call, group.by = "by = cells")
  metrics <- metrics[cells > -1]
  n_metrics <- ncol(metrics) - 1

  output <- terra::rast(
    bbox,
    resolution = res,
    vals = as.numeric(NA),
    nlyr = n_metrics,
    crs = "epsg:4326"
  )

  names(output) <- names(metrics)[-1]

  for (metric in seq_along(names(metrics)[-1])) {
    output[[metric]][metrics$cells] <- metrics[[metric + 1]]
  }

  return(output)
}
