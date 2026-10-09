# GEDI sampling ----------------------------------------------------------

.gedi_coords <- function(x, lon = NULL, lat = NULL) {
  candidates <- list(
    c("lon_lowestmode", "lat_lowestmode"), c("longitude_bin0", "latitude_bin0"),
    c("longitude", "latitude"), c("lon", "lat")
  )
  if (!is.null(lon) || !is.null(lat)) candidates <- list(c(lon, lat))
  for (pair in candidates) if (all(pair %in% names(x))) return(pair)
  stop("Could not identify longitude and latitude columns.", call. = FALSE)
}

.sampling_method <- function(name, ...) structure(list(name = name, args = list(...)), class = "gedi_sampling_method")

#' GEDI random sampling specification
#' @param size Number or fraction of observations to select.
#' @export
randomSampling <- function(size) .sampling_method("random", size = size)

#' GEDI minimum-distance sampling specification
#' @param size Number or fraction of observations to select.
#' @param radius Minimum separation in meters.
#' @param lon,lat Optional coordinate column names.
#' @export
spacedSampling <- function(size, radius, lon = NULL, lat = NULL) .sampling_method("spaced", size = size, radius = radius, lon = lon, lat = lat)

#' GEDI grid-stratified sampling specification
#' @param size Number or fraction per cell.
#' @param grid_size Grid size in decimal degrees.
#' @param lon,lat Optional coordinate column names.
#' @export
gridSampling <- function(size, grid_size, lon = NULL, lat = NULL) .sampling_method("grid", size = size, grid_size = grid_size, lon = lon, lat = lat)

#' GEDI attribute-stratified sampling specification
#' @param size Number or fraction per stratum.
#' @param variable Column used as strata, or a numeric column to bin.
#' @param breaks Optional breaks for a numeric variable.
#' @export
stratifiedSampling <- function(size, variable, breaks = NULL) .sampling_method("stratified", size = size, variable = variable, breaks = breaks)

#' GEDI polygon-stratified sampling specification
#' @param size Number or fraction per polygon.
#' @param geom Polygon object.
#' @param split_id Optional polygon ID column.
#' @param lon,lat Optional coordinate column names.
#' @export
geomSampling <- function(size, geom, split_id = NULL, lon = NULL, lat = NULL) .sampling_method("geometry", size = size, geom = geom, split_id = split_id, lon = lon, lat = lat)

#' GEDI raster-stratified sampling specification
#' @param size Number or fraction per raster class.
#' @param raster A categorical `SpatRaster`.
#' @param lon,lat Optional coordinate column names.
#' @export
rasterSampling <- function(size, raster, lon = NULL, lat = NULL) .sampling_method("raster", size = size, raster = raster, lon = lon, lat = lat)

.sample_indices <- function(index, size) {
  n <- length(index)
  if (length(size) != 1L || !is.numeric(size) || !is.finite(size) || size < 0) {
    stop("Sampling `size` must be one finite non-negative number.", call. = FALSE)
  }
  take <- if (size > 0 && size < 1) ceiling(n * size) else as.integer(size)
  take <- max(0L, min(n, take))
  if (!take) integer() else base::sample(index, take)
}

#' Sample GEDI footprint observations
#'
#' @param x A data frame or data table of GEDI footprints.
#' @param method A specification returned by one of the sampling constructors.
#' @return A [data.table::data.table].
#' @export
sampleGEDI <- function(x, method) {
  if (!inherits(x, c("data.frame", "data.table"))) stop("'x' must be a GEDI footprint table.")
  if (!inherits(method, "gedi_sampling_method")) stop("'method' must be a GEDI sampling specification.")
  dt <- data.table::as.data.table(x)
  a <- method$args
  if (method$name == "random") return(dt[.sample_indices(seq_len(nrow(dt)), a$size)])
  if (method$name == "stratified") {
    if (!a$variable %in% names(dt)) stop("Unknown stratification column: ", a$variable)
    group <- dt[[a$variable]]
    if (is.numeric(group) && !is.null(a$breaks)) group <- cut(group, a$breaks, include.lowest = TRUE)
    idx <- split(seq_len(nrow(dt)), group, drop = TRUE)
    return(dt[unlist(lapply(idx, .sample_indices, size = a$size), use.names = FALSE)])
  }
  coords <- .gedi_coords(dt, a$lon, a$lat)
  if (method$name == "grid") {
    gx <- floor(dt[[coords[1L]]] / a$grid_size)
    gy <- floor(dt[[coords[2L]]] / a$grid_size)
    idx <- split(seq_len(nrow(dt)), interaction(gx, gy, drop = TRUE))
    return(dt[unlist(lapply(idx, .sample_indices, size = a$size), use.names = FALSE)])
  }
  if (method$name == "geometry") {
    pts <- sf::st_as_sf(as.data.frame(dt), coords = coords, crs = 4326, remove = FALSE)
    poly <- sf::st_transform(sf::st_as_sf(a$geom), 4326)
    hits <- sf::st_intersects(pts, poly)
    groups <- vapply(hits, function(z) if (length(z)) z[[1L]] else NA_integer_, integer(1))
    idx <- split(which(!is.na(groups)), groups[!is.na(groups)])
    chosen <- unlist(lapply(idx, .sample_indices, size = a$size), use.names = FALSE)
    ans <- dt[chosen]
    if (!is.null(a$split_id)) ans[["sample_group"]] <- poly[[a$split_id]][groups[chosen]]
    return(ans)
  }
  if (method$name == "raster") {
    pts <- terra::vect(as.data.frame(dt), geom = coords, crs = "EPSG:4326")
    values <- terra::extract(a$raster, pts)[[2L]]
    idx <- split(which(!is.na(values)), values[!is.na(values)])
    return(dt[unlist(lapply(idx, .sample_indices, size = a$size), use.names = FALSE)])
  }
  if (method$name == "spaced") {
    if (length(a$radius) != 1L || !is.finite(a$radius) || a$radius < 0) {
      stop("Sampling `radius` must be one finite non-negative number.", call. = FALSE)
    }
    pts <- sf::st_as_sf(as.data.frame(dt), coords = coords, crs = 4326, remove = FALSE)
    pts <- sf::st_transform(pts, 6933)
    order <- base::sample(seq_len(nrow(dt)))
    selected <- integer()
    target <- if (a$size > 0 && a$size < 1) ceiling(nrow(dt) * a$size) else as.integer(a$size)
    target <- max(0L, min(nrow(dt), target))
    xy <- sf::st_coordinates(pts)
    for (i in order) {
      if (!length(selected) || all(sqrt(rowSums((xy[selected, , drop = FALSE] - matrix(xy[i, ], nrow = length(selected), ncol = 2, byrow = TRUE))^2)) >= a$radius)) selected <- c(selected, i)
      if (length(selected) >= target) break
    }
    return(dt[selected])
  }
  stop("Unknown GEDI sampling method.")
}

#' Convert GEDI footprints to a spatial vector
#' @param x GEDI footprint table.
#' @param lon,lat Optional coordinate column names.
#' @param crs Coordinate reference system.
#' @return A [`terra::SpatVector-class`].
#' @export
to_vect <- function(x, lon = NULL, lat = NULL, crs = "EPSG:4326") {
  coords <- .gedi_coords(x, lon, lat)
  terra::vect(as.data.frame(x), geom = coords, crs = crs)
}
