#' Compute Grids with Descriptive Statistics of
#' GEDI derived Elevation and Height Metrics (Level2A)
#'
#' @description This function computes a series of user defined descriptive statistics within
#' each grid cell for GEDI derived Elevation and Height Metrics (Level2A)
#'
#' @usage gridStatsLevel2AM(level2AM, func, res)
#'
#' @param level2AM A GEDI Level2AM object (output of [getLevel2AM()] function).
#' An S4 object of class "data.table".
#' @param func The function(s) to be applied to each cell
#' @param res Spatial resolution in decimal degrees for the output stars raster layer
#'
#' @return Return a stars raster layer(s) of selected GEDI Elevation and Height Metric(s)
#'
#' @seealso \url{https://lpdaac.usgs.gov/products/gedi02_av002/}
#'
#' @examples
#' # specify the path to GEDI level2A data (zip file)
#' outdir <- tempdir()
#' level2A_fp_zip <- system.file("extdata",
#'   "GEDI02_A_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'   package = "rGEDI"
#' )
#'
#' # Unzipping GEDI level2A data
#' level2Apath <- unzip(level2A_fp_zip, exdir = outdir)
#'
#' # Reading GEDI level2A data (h5 file)
#' level2a <- readLevel2A(level2Apath = level2Apath)
#'
#' # Get GEDI derived Elevation and Height Metrics
#' level2AM <- getLevel2AM(level2a)
#' head(level2AM)
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
#' #' Computing a serie of GEDI metrics
#' ZTstats <- gridStatsLevel2AM(
#'   level2AM = level2AM,
#'   func = mySetOfMetrics(elev_highestreturn),
#'   res = 0.005
#' )
#' plot(ZTstats)
#'
#' #' Computing the maximum of RH100 only
#' maxRH100 <- gridStatsLevel2AM(level2AM = level2AM, func = mySetOfMetrics(rh100), res = 0.0005)
#' plot(maxRH100)
#'
#' #' Computing the mean of ZG only
#' ZGmean <- gridStatsLevel2AM(level2AM = level2AM, func = mean(elev_lowestmode), res = 0.005)
#' plot(ZGmean)
#'
#' close(level2a)
#' @importFrom stats setNames na.omit
#' @export
gridStatsLevel2AM <- function(level2AM, func, res = 0.5) {
  requireNamespace("data.table")
  cells <- NA
  # this code has been adapted from the grid_metrics function in lidR package (Roussel et al. 2019)
  # https://github.com/Jean-Romain/lidR/blob/master/R/grid_metrics.r

  # Add data.table operator
  `:=` <- data.table::`:=`
  `%>%` <- sf::`%>%`

  call <- lazy_call(func)


  sf <- sf::st_as_sf(level2AM, coords = c("lon_lowestmode", "lat_lowestmode"))
  layout <- sf::st_bbox(sf) %>%
    stars::st_as_stars(dx = res, dy = res, values = NA, crs = "epsg:4326")

  level2AM[, cells := stars::st_cells(layout, sf)]
  metrics <- lazy_apply_dt_call(level2AM, call, group.by = "by = cells")
  n_metrics <- ncol(metrics) - 1
  output <- sf::st_bbox(sf) %>%
    stars::st_as_stars(
      dx = res,
      dy = res,
      values = as.numeric(NA),
      nz = n_metrics,
      crs = "epsg:4326"
    )

  dim_data <- stars::st_dimensions(output) %>%
    setNames(c("x", "y", "bands"))
  dim_data$bands$values <- names(metrics)[-1]

  stars::st_dimensions(output) <- dim_data

  for (metric in seq_along(names(metrics)[-1])) {
    output[[1]][, , metric][metrics$cells] <- metrics[[metric]]
  }

  return(output)
}
