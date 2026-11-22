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
  if (inherits(x, "gedi_gee_reference")) {
    object <- .open_gee_reference(x, product)
    if (!is.null(metrics)) object <- .gee_select(object, metrics)
    return(object)
  }
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
#' @export
readLevel3 <- function(x, metrics = NULL) .read_gedi_raster(x, "GEDI03", metrics)

#' Read GEDI Level 4B gridded biomass
#'
#' @inheritParams readLevel3
#' @return A [`terra::SpatRaster-class`] for local or GDAL-readable data, or a
#'   lazy `gedi_gee_object` for a Google Earth Engine reference.
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

.get_gedi_raster <- function(x, geometry = NULL, metrics = NULL, xy = TRUE,
                             cells = FALSE, na.rm = TRUE, global = FALSE,
                             max_cells = 5e6, ...) {
  if (inherits(x, "gedi_gee_object")) {
    if (is.null(geometry)) {
      stop("A geometry is required when retrieving pixels from Google Earth Engine.",
           call. = FALSE)
    }
    if (!is.null(metrics)) x <- .gee_select(x, metrics)
    return(extractEE(.unwrap_gee(x), geometry, ...))
  }
  raster <- .as_gedi_raster(x)
  if (!is.null(metrics)) {
    missing <- setdiff(metrics, names(raster))
    if (length(missing)) stop("Unknown metric(s): ", paste(missing, collapse = ", "), call. = FALSE)
    raster <- raster[[metrics]]
  }
  if (!is.null(geometry)) {
    geometry <- if (is.numeric(geometry) && length(geometry) == 4L) {
      terra::as.polygons(terra::ext(geometry))
    } else {
      terra::vect(geometry)
    }
    return(data.table::as.data.table(
      terra::extract(raster, geometry, ...)
    ))
  }
  requested <- terra::ncell(raster)
  if (!isTRUE(global) && requested > max_cells) {
    stop(
      "The raster has ", format(requested, big.mark = ","),
      " cells. Supply `geometry`, or set `global = TRUE` explicitly.",
      call. = FALSE
    )
  }
  data.table::as.data.table(
    terra::as.data.frame(raster, xy = xy, cells = cells, na.rm = na.rm, ...)
  )
}

#' Get GEDI Level 3 pixel measurements
#'
#' @param x A Level 3 raster returned by [readLevel3()] or a raster path.
#' @param geometry Optional points, polygons, or extent accepted by
#'   [terra::vect()]. A geometry is required for Google Earth Engine data.
#' @param metrics Optional raster layers to retrieve.
#' @param xy Include cell-center coordinates when returning all cells.
#' @param cells Include cell numbers when returning all cells.
#' @param na.rm Remove rows whose raster values are all missing.
#' @param global Explicitly allow conversion of a large complete raster.
#' @param max_cells Maximum cells converted without `global = TRUE`.
#' @param ... Passed to [terra::extract()], [terra::as.data.frame()], or
#'   [extractEE()].
#' @return A [data.table::data.table] containing pixel measurements.
#' @export
getLevel3 <- function(x, geometry = NULL, metrics = NULL, xy = TRUE,
                      cells = FALSE, na.rm = TRUE, global = FALSE,
                      max_cells = 5e6, ...) {
  .get_gedi_raster(x, geometry, metrics, xy, cells, na.rm, global, max_cells, ...)
}

#' Get GEDI Level 4B pixel measurements
#' @inheritParams getLevel3
#' @param x A Level 4B raster returned by [readLevel4B()], a raster path, or a
#'   lazy Google Earth Engine Level 4B object.
#' @return A [data.table::data.table] containing pixel-level biomass measurements.
#' @export
getLevel4B <- function(x, geometry = NULL, metrics = NULL, xy = TRUE,
                       cells = FALSE, na.rm = TRUE, global = FALSE,
                       max_cells = 5e6, ...) {
  .get_gedi_raster(x, geometry, metrics, xy, cells, na.rm, global, max_cells, ...)
}

.poly_stats_raster <- function(x, polygon, fun = mean, na.rm = TRUE, ...) {
  data.table::as.data.table(
    terra::extract(.as_gedi_raster(x), terra::vect(polygon),
                   fun = fun, na.rm = na.rm, ...)
  )
}

#' Polygon summaries for any GEDI raster
#'
#' Works with native Level 3 and Level 4B products and with rasters created
#' from Level 2A, Level 2B, or Level 4A footprint measurements.
#'
#' @param x A raster path or [terra::SpatRaster-class].
#' @param polygon Polygon geometry accepted by [terra::vect()].
#' @param fun Summary function.
#' @param na.rm Logical.
#' @param ... Passed to [terra::extract()].
#' @return A [data.table::data.table] with one row per polygon and summary
#'   columns for every raster layer.
#' @export
polyStats <- function(x, polygon, fun = mean, na.rm = TRUE, ...) {
  .poly_stats_raster(x, polygon, fun, na.rm, ...)
}

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

.grid_gedi_points <- function(x, lon, lat, metric, fun, res, ...,
                              template = NULL, input_crs = "EPSG:4326") {
  if (!all(c(lon, lat, metric) %in% names(x))) stop("Coordinate or metric columns are missing.")
  pts <- terra::vect(as.data.frame(x), geom = c(lon, lat), crs = input_crs)
  if (is.null(template)) {
    if (!terra::is.lonlat(pts)) stop("Supply a template when using projected coordinates.")
    template <- terra::rast(terra::ext(pts), resolution = res, crs = input_crs)
  } else {
    if (!inherits(template, "SpatRaster") || !nzchar(terra::crs(template))) {
      stop("The template must be a SpatRaster with a coordinate reference system.")
    }
    if (!terra::same.crs(pts, template)) pts <- terra::project(pts, terra::crs(template))
  }
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
