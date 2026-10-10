#' Clip GEDI data with a common interface
#'
#' `clip()` provides one clipping interface for open GEDI HDF5 products,
#' extracted footprint tables, and GEDI Level 3 or Level 4B rasters. It follows
#' the workflow used by ICESat2VegR while retaining rGEDI's coordinate order.
#'
#' @param x A `gedi.level1b`, `gedi.level2a`, `gedi.level2b`, or
#'   `gedi.level4a` object; a footprint `data.frame`/`data.table`; or a
#'   [`terra::SpatRaster-class`].
#' @param clip_obj A numeric extent in the order
#'   `c(xmin, xmax, ymin, ymax)`, a [`terra::SpatExtent-class`], or an `sf`,
#'   `sfc`, or [`terra::SpatVector-class`] polygon object.
#' @param ... Arguments used by the selected method. HDF5 methods accept
#'   `output` and geometry methods also accept `split_by`. Table geometry
#'   methods accept `split_by`. Raster methods accept `filename` and
#'   `overwrite`.
#'
#' @return The return type follows `x`. Open Level 1B, 2A, 2B, and 4A HDF5
#'   inputs return clipped open GEDI objects; footprint tables return a
#'   [data.table::data.table]; and Level 3/4B rasters return a
#'   [`terra::SpatRaster-class`]. Geometry clipping with `split_by` returns a
#'   named list of HDF5-backed GEDI objects, one per polygon ID.
#'
#' @details
#' Extracted tables are recognized from standard GEDI coordinates:
#' `longitude_bin0`/`latitude_bin0`, `lon_lowestmode`/`lat_lowestmode`, or
#' `longitude`/`latitude`. When the corresponding `longitude_lastbin` and
#' `latitude_lastbin` fields exist, both waveform endpoints must fall inside a
#' numeric extent.
#'
#' This function masks [graphics::clip()]. Use `graphics::clip()` explicitly
#' when changing a base graphics clipping region.
#'
#' @importClassesFrom terra SpatExtent SpatRaster
#' @examples
#' footprints <- data.frame(
#'   shot_number = 1:3,
#'   lon_lowestmode = c(-44.13, -44.12, -44.10),
#'   lat_lowestmode = c(-13.74, -13.73, -13.70)
#' )
#' clip(footprints, c(-44.14, -44.11, -13.75, -13.72))
#'
#' @export
setGeneric("clip", function(x, clip_obj, ...) standardGeneric("clip"))

.clip_bounds <- function(clip_obj) {
  if (inherits(clip_obj, "SpatExtent")) {
    bounds <- as.vector(clip_obj)
  } else {
    bounds <- clip_obj
  }
  if (!is.numeric(bounds) || length(bounds) != 4L || any(!is.finite(bounds))) {
    stop("'clip_obj' must be a finite numeric c(xmin, xmax, ymin, ymax) or SpatExtent.",
      call. = FALSE
    )
  }
  bounds <- as.numeric(bounds)
  if (bounds[1L] > bounds[2L] || bounds[3L] > bounds[4L]) {
    stop("'clip_obj' must satisfy xmin <= xmax and ymin <= ymax.", call. = FALSE)
  }
  names(bounds) <- c("xmin", "xmax", "ymin", "ymax")
  bounds
}

.gedi_table_coordinates <- function(x) {
  candidates <- list(
    c("longitude_bin0", "latitude_bin0"),
    c("lon_lowestmode", "lat_lowestmode"),
    c("longitude", "latitude")
  )
  found <- vapply(candidates, function(fields) all(fields %in% names(x)), logical(1))
  if (!any(found)) {
    stop(
      "Cannot find GEDI coordinates. Expected longitude_bin0/latitude_bin0, lon_lowestmode/lat_lowestmode, or longitude/latitude.",
      call. = FALSE
    )
  }
  candidates[[which(found)[1L]]]
}

.clip_gedi_table_extent <- function(x, clip_obj) {
  bounds <- .clip_bounds(clip_obj)
  coords <- .gedi_table_coordinates(x)
  lon <- x[[coords[1L]]]
  lat <- x[[coords[2L]]]
  keep <- lon >= bounds["xmin"] & lon <= bounds["xmax"] &
    lat >= bounds["ymin"] & lat <= bounds["ymax"]

  if (identical(coords, c("longitude_bin0", "latitude_bin0")) &&
      all(c("longitude_lastbin", "latitude_lastbin") %in% names(x))) {
    keep <- keep &
      x[["longitude_lastbin"]] >= bounds["xmin"] &
      x[["longitude_lastbin"]] <= bounds["xmax"] &
      x[["latitude_lastbin"]] >= bounds["ymin"] &
      x[["latitude_lastbin"]] <= bounds["ymax"]
  }
  keep[is.na(keep)] <- FALSE
  data.table::as.data.table(x)[keep]
}

.clip_gedi_table_geometry <- function(x, clip_obj, split_by = NULL) {
  coords <- .gedi_table_coordinates(x)
  poly <- if (inherits(clip_obj, "SpatVector")) {
    sf::st_as_sf(clip_obj)
  } else {
    tryCatch(sf::st_as_sf(clip_obj), error = function(e) {
      stop("'clip_obj' must be an sf, sfc, or SpatVector polygon object.", call. = FALSE)
    })
  }
  if (!nrow(poly)) return(data.table::as.data.table(x)[0])
  geometry_types <- unique(as.character(sf::st_geometry_type(poly)))
  if (!all(geometry_types %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop("'clip_obj' must contain polygon geometries.", call. = FALSE)
  }
  if (!is.null(split_by) && !split_by %in% names(poly)) {
    stop("'split_by' is not a polygon attribute.", call. = FALSE)
  }

  poly <- sf::st_transform(sf::st_make_valid(poly), 4326)
  points <- sf::st_as_sf(
    as.data.frame(x), coords = coords, crs = 4326, remove = FALSE
  )
  hits <- sf::st_intersects(points, poly)
  keep <- lengths(hits) > 0L
  ans <- data.table::as.data.table(x)[keep]
  if (!is.null(split_by) && nrow(ans)) {
    first_hit <- vapply(hits[keep], `[`, integer(1), 1L)
    ans[["poly_id"]] <- poly[[split_by]][first_hit]
  }
  ans[]
}

.clip_h5_extent <- function(x, clip_obj, fun, ...) {
  bounds <- .clip_bounds(clip_obj)
  fun(x, bounds["xmin"], bounds["xmax"], bounds["ymin"], bounds["ymax"], ...)
}

.clip_h5_geometry <- function(x, clip_obj, fun, ...) {
  polygon <- if (inherits(clip_obj, "sf")) clip_obj else {
    tryCatch(sf::st_as_sf(clip_obj), error = function(e) {
      stop("'clip_obj' must be an sf, sfc, or SpatVector polygon object.", call. = FALSE)
    })
  }
  dots <- list(...)
  result <- do.call(fun, c(list(x, polygon), dots))
  split_by <- if ("split_by" %in% names(dots)) dots$split_by else NULL
  if (is.null(split_by) && is.list(result) && length(result) == 1L) result[[1L]] else result
}

#' @rdname clip
#' @export
setMethod("clip", c(x = "gedi.level1b", clip_obj = "numeric"),
  function(x, clip_obj, ...) .clip_h5_extent(x, clip_obj, clipLevel1B, ...))

#' @rdname clip
#' @export
setMethod("clip", c(x = "gedi.level1b", clip_obj = "SpatExtent"),
  function(x, clip_obj, ...) .clip_h5_extent(x, clip_obj, clipLevel1B, ...))

#' @rdname clip
#' @export
setMethod("clip", c(x = "gedi.level1b", clip_obj = "ANY"),
  function(x, clip_obj, ...) .clip_h5_geometry(x, clip_obj, clipLevel1BGeometry, ...))

#' @rdname clip
#' @export
setMethod("clip", c(x = "gedi.level2a", clip_obj = "numeric"),
  function(x, clip_obj, ...) .clip_h5_extent(x, clip_obj, clipLevel2A, ...))

#' @rdname clip
#' @export
setMethod("clip", c(x = "gedi.level2a", clip_obj = "SpatExtent"),
  function(x, clip_obj, ...) .clip_h5_extent(x, clip_obj, clipLevel2A, ...))

#' @rdname clip
#' @export
setMethod("clip", c(x = "gedi.level2a", clip_obj = "ANY"),
  function(x, clip_obj, ...) .clip_h5_geometry(x, clip_obj, clipLevel2AGeometry, ...))

#' @rdname clip
#' @export
setMethod("clip", c(x = "gedi.level2b", clip_obj = "numeric"),
  function(x, clip_obj, ...) .clip_h5_extent(x, clip_obj, clipLevel2B, ...))

#' @rdname clip
#' @export
setMethod("clip", c(x = "gedi.level2b", clip_obj = "SpatExtent"),
  function(x, clip_obj, ...) .clip_h5_extent(x, clip_obj, clipLevel2B, ...))

#' @rdname clip
#' @export
setMethod("clip", c(x = "gedi.level2b", clip_obj = "ANY"),
  function(x, clip_obj, ...) .clip_h5_geometry(x, clip_obj, clipLevel2BGeometry, ...))

#' @rdname clip
#' @export
setMethod("clip", c(x = "gedi.level4a", clip_obj = "numeric"),
  function(x, clip_obj, ...) .clip_h5_extent(x, clip_obj, clipLevel4A, ...))

#' @rdname clip
#' @export
setMethod("clip", c(x = "gedi.level4a", clip_obj = "SpatExtent"),
  function(x, clip_obj, ...) .clip_h5_extent(x, clip_obj, clipLevel4A, ...))

#' @rdname clip
#' @export
setMethod("clip", c(x = "gedi.level4a", clip_obj = "ANY"),
  function(x, clip_obj, ...) .clip_h5_geometry(x, clip_obj, clipLevel4AGeometry, ...))

#' @rdname clip
#' @export
setMethod("clip", c(x = "data.frame", clip_obj = "numeric"),
  function(x, clip_obj, ...) .clip_gedi_table_extent(x, clip_obj))

#' @rdname clip
#' @export
setMethod("clip", c(x = "data.frame", clip_obj = "SpatExtent"),
  function(x, clip_obj, ...) .clip_gedi_table_extent(x, clip_obj))

#' @rdname clip
#' @export
setMethod("clip", c(x = "data.frame", clip_obj = "ANY"),
  function(x, clip_obj, ...) .clip_gedi_table_geometry(x, clip_obj, ...))

#' @rdname clip
#' @export
setMethod("clip", c(x = "SpatRaster", clip_obj = "ANY"),
  function(x, clip_obj, ...) .clip_gedi_raster(x, clip_obj, ...))
