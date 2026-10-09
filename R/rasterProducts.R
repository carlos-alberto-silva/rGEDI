# GEDI gridded products ---------------------------------------------------

.gedi_metric_name <- function(path, product) {
  stem <- tools::file_path_sans_ext(basename(path))
  if (product == "GEDI03") {
    value <- sub("^GEDI03_", "", stem)
    value <- sub("_[0-9]{7}_[0-9]{7}_[0-9]{3}_[0-9]{2}$", "", value)
    return(value)
  }
  value <- sub("^GEDI04_B_", "", stem)
  suffix <- sub("^.*_", "", value)
  if (nzchar(suffix)) suffix else value
}

.read_gedi_raster <- function(x, product, metrics = NULL) {
  if (inherits(x, "SpatRaster")) return(x)
  if (!is.character(x) || !length(x)) stop("'x' must contain raster paths or URLs.")
  r <- terra::rast(x)
  if (terra::nlyr(r) == length(x)) names(r) <- make.unique(vapply(x, .gedi_metric_name, character(1), product = product))
  if (!is.null(metrics)) {
    missing <- setdiff(metrics, names(r))
    if (length(missing)) stop("Unknown metric(s): ", paste(missing, collapse = ", "))
    r <- r[[metrics]]
  }
  attr(r, "gedi_product") <- product
  r
}

#' Read GEDI Level 3 gridded metrics
#'
#' @param x One or more local paths, GDAL-readable URLs, or a `SpatRaster`.
#' @param metrics Optional layer names to retain.
#' @return A [`terra::SpatRaster-class`].
#' @seealso \url{https://daac.ornl.gov/GEDI/guides/GEDI_L3_LandSurface_Metrics_V3.html}
#' @export
readLevel3 <- function(x, metrics = NULL) .read_gedi_raster(x, "GEDI03", metrics)

#' Read GEDI Level 4B gridded biomass
#'
#' @inheritParams readLevel3
#' @return A [`terra::SpatRaster-class`].
#' @seealso \doi{10.3334/ORNLDAAC/2299}
#' @export
readLevel4B <- function(x, metrics = NULL) .read_gedi_raster(x, "GEDI04_B", metrics)

.gedi_rasters <- function(x) {
  if (inherits(x, "SpatRaster")) return(list(x))
  if (is.list(x) && all(vapply(x, inherits, logical(1), "SpatRaster"))) return(x)
  if (is.character(x)) return(lapply(x, terra::rast))
  stop("Expected raster paths, a SpatRaster, or a list of SpatRaster objects.")
}

.as_gedi_raster <- function(x) if (inherits(x, "SpatRaster")) x else terra::rast(x)

.mosaic_gedi <- function(x, fun = "mean", filename = "", overwrite = FALSE) {
  rasters <- .gedi_rasters(x)
  if (length(rasters) == 1L) return(rasters[[1L]])
  do.call(terra::mosaic, c(rasters, list(fun = fun, filename = filename, overwrite = overwrite)))
}

#' Mosaic GEDI Level 3 rasters
#' @param x Raster paths, a `SpatRaster`, or a list of rasters.
#' @param fun Function used for overlapping cells.
#' @param filename Optional output filename.
#' @param overwrite Logical.
#' @return A [`terra::SpatRaster-class`].
#' @export
mosaicLevel3 <- function(x, fun = "mean", filename = "", overwrite = FALSE) {
  .mosaic_gedi(x, fun, filename, overwrite)
}

#' Mosaic GEDI Level 4B rasters
#' @inheritParams mosaicLevel3
#' @return A [`terra::SpatRaster-class`].
#' @export
mosaicLevel4B <- function(x, fun = "mean", filename = "", overwrite = FALSE) {
  .mosaic_gedi(x, fun, filename, overwrite)
}

.clip_gedi_raster <- function(x, y, filename = "", overwrite = FALSE) {
  r <- .as_gedi_raster(x)
  geom <- if (is.numeric(y) && length(y) == 4L) terra::ext(y) else terra::vect(y)
  ans <- terra::crop(r, geom)
  if (!inherits(geom, "SpatExtent")) ans <- terra::mask(ans, geom)
  if (nzchar(filename)) terra::writeRaster(ans, filename, overwrite = overwrite) else ans
}

#' Clip a GEDI Level 3 raster
#' @param x A raster path or `SpatRaster`.
#' @param y An extent `c(xmin, xmax, ymin, ymax)` or polygon object.
#' @param filename Optional output filename.
#' @param overwrite Logical.
#' @return A [`terra::SpatRaster-class`].
#' @export
clipLevel3 <- function(x, y, filename = "", overwrite = FALSE) .clip_gedi_raster(x, y, filename, overwrite)

#' Clip a GEDI Level 4B raster
#' @inheritParams clipLevel3
#' @return A [`terra::SpatRaster-class`].
#' @export
clipLevel4B <- function(x, y, filename = "", overwrite = FALSE) .clip_gedi_raster(x, y, filename, overwrite)

#' Project a GEDI Level 3 raster
#' @param x A raster path or `SpatRaster`.
#' @param crs Target CRS.
#' @param method Resampling method.
#' @param ... Passed to [terra::project()].
#' @return A [`terra::SpatRaster-class`].
#' @export
projectLevel3 <- function(x, crs, method = "bilinear", ...) terra::project(.as_gedi_raster(x), crs, method = method, ...)

#' Project a GEDI Level 4B raster
#' @inheritParams projectLevel3
#' @return A [`terra::SpatRaster-class`].
#' @export
projectLevel4B <- function(x, crs, method = "bilinear", ...) terra::project(.as_gedi_raster(x), crs, method = method, ...)

#' Resample a GEDI Level 3 raster
#' @param x A raster path or `SpatRaster`.
#' @param template Target `SpatRaster`.
#' @param method Resampling method.
#' @param ... Passed to [terra::resample()].
#' @return A [`terra::SpatRaster-class`].
#' @export
resampleLevel3 <- function(x, template, method = "bilinear", ...) terra::resample(.as_gedi_raster(x), template, method = method, ...)

#' Resample a GEDI Level 4B raster
#' @inheritParams resampleLevel3
#' @return A [`terra::SpatRaster-class`].
#' @export
resampleLevel4B <- function(x, template, method = "bilinear", ...) terra::resample(.as_gedi_raster(x), template, method = method, ...)

#' Aggregate a GEDI Level 3 raster
#' @param x A raster path or `SpatRaster`.
#' @param fact Aggregation factor.
#' @param fun Aggregation function.
#' @param ... Passed to [terra::aggregate()].
#' @return A [`terra::SpatRaster-class`].
#' @export
aggregateLevel3 <- function(x, fact, fun = mean, ...) terra::aggregate(.as_gedi_raster(x), fact = fact, fun = fun, ...)

#' Aggregate a GEDI Level 4B raster
#' @inheritParams aggregateLevel3
#' @return A [`terra::SpatRaster-class`].
#' @export
aggregateLevel4B <- function(x, fact, fun = mean, ...) terra::aggregate(.as_gedi_raster(x), fact = fact, fun = fun, ...)

.extract_gedi_raster <- function(x, y, ...) terra::extract(.as_gedi_raster(x), terra::vect(y), ...)

#' Extract GEDI Level 3 raster values
#' @param x A raster path or `SpatRaster`.
#' @param y Points or polygons accepted by [terra::vect()].
#' @param ... Passed to [terra::extract()].
#' @return A data frame.
#' @export
extractLevel3 <- function(x, y, ...) .extract_gedi_raster(x, y, ...)

#' Extract GEDI Level 4B raster values
#' @inheritParams extractLevel3
#' @return A data frame.
#' @export
extractLevel4B <- function(x, y, ...) .extract_gedi_raster(x, y, ...)

.poly_stats_raster <- function(x, polygon, fun = mean, na.rm = TRUE, ...) {
  terra::extract(.as_gedi_raster(x), terra::vect(polygon), fun = fun, na.rm = na.rm, ...)
}

#' Polygon summaries of GEDI Level 3 metrics
#' @inheritParams extractLevel3
#' @param fun Summary function.
#' @param na.rm Logical.
#' @return A data frame.
#' @export
polyStatsLevel3 <- function(x, y, fun = mean, na.rm = TRUE, ...) .poly_stats_raster(x, y, fun, na.rm, ...)

#' Polygon summaries of GEDI Level 4B metrics
#' @inheritParams polyStatsLevel3
#' @return A data frame.
#' @export
polyStatsLevel4B <- function(x, y, fun = mean, na.rm = TRUE, ...) .poly_stats_raster(x, y, fun, na.rm, ...)

#' Plot GEDI Level 3 metrics
#' @param x A raster path or `SpatRaster`.
#' @param ... Passed to [terra::plot()].
#' @return Invisibly returns the raster.
#' @export
plotLevel3 <- function(x, ...) { r <- .as_gedi_raster(x); terra::plot(r, ...); invisible(r) }

#' Plot GEDI Level 4B metrics
#' @inheritParams plotLevel3
#' @return Invisibly returns the raster.
#' @export
plotLevel4B <- function(x, ...) { r <- .as_gedi_raster(x); terra::plot(r, ...); invisible(r) }

.grid_gedi_points <- function(x, lon, lat, metric, fun, res, ...) {
  if (!all(c(lon, lat, metric) %in% names(x))) stop("Coordinate or metric columns are missing.")
  pts <- terra::vect(as.data.frame(x), geom = c(lon, lat), crs = "EPSG:4326")
  template <- terra::rast(terra::ext(pts), resolution = res, crs = "EPSG:4326")
  terra::rasterize(pts, template, field = metric, fun = fun, ...)
}

.poly_gedi_points <- function(x, polygon, lon, lat, metric, fun, id = NULL, ...) {
  if (!all(c(lon, lat, metric) %in% names(x))) stop("Coordinate or metric columns are missing.")
  poly <- sf::st_transform(sf::st_as_sf(polygon), 4326)
  pts <- sf::st_as_sf(as.data.frame(x), coords = c(lon, lat), crs = 4326,
                      remove = FALSE)
  hits <- sf::st_intersects(pts, poly)
  first_hit <- vapply(hits, function(z) if (length(z)) z[[1L]] else NA_integer_,
                      integer(1))
  dt <- data.table::as.data.table(x)[!is.na(first_hit)]
  dt[[".polygon_id__"]] <- first_hit[!is.na(first_hit)]
  split_values <- split(dt[[metric]], dt[[".polygon_id__"]])
  ans <- data.table::data.table(
    .polygon_id__ = as.integer(names(split_values)),
    value = vapply(split_values, function(z) do.call(fun, c(list(z), list(...))),
                   numeric(1))
  )
  if (!is.null(id)) {
    if (!id %in% names(poly)) stop("Unknown polygon ID field: ", id)
    ans[["polygon"]] <- poly[[id]][ans[[".polygon_id__"]]]
  }
  ans[]
}
