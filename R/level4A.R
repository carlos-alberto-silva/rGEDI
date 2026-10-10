# GEDI Level 4A -----------------------------------------------------------

.level4a_default_columns <- c(
  "beam", "shot_number", "l2_algrunflag", "l2a_quality_flag_rel3",
  "l4a_quality_flag_rel3", "degrade_include_flag",
  "elev_highestreturn_outlier_flag", "degrade_flag", "delta_time", "sensitivity",
  "solar_elevation", "surface_flag", "lat_lowestmode", "lon_lowestmode",
  "elev_lowestmode", "agbd", "agbd_se", "agbd_pi_lower",
  "agbd_pi_upper", "predict_stratum", "selected_algorithm"
)

#' Read a GEDI Level 4A granule
#'
#' Opens a GEDI04_A HDF5 granule containing footprint-level aboveground
#' biomass density estimates.
#'
#' @param level4Apath Local file path or Earthdata Cloud URL pointing to a
#'   GEDI04_A HDF5 granule.
#' @return A [`gedi.level4a-class`] object. Close it with [close()].
#' @seealso \url{https://daac.ornl.gov/GEDI/guides/GEDI_L4A_AGB_Density_V3.html}
#' @export
readLevel4A <- function(level4Apath) {
  if (inherits(level4Apath, "GEDICloudH5")) {
    h5 <- level4Apath
  } else if (.is_gedi_url(level4Apath)) {
    h5 <- .open_cloud_h5(level4Apath, product = "GEDI04_A")
  } else {
    if (!is.character(level4Apath) || length(level4Apath) != 1L || !file.exists(level4Apath)) {
      stop("'level4Apath' must be an existing file or Earthdata Cloud URL.", call. = FALSE)
    }
    h5 <- hdf5r::H5File$new(level4Apath, mode = "r")
  }
  new("gedi.level4a", h5 = h5)
}

.gedi_beams <- function(h5) {
  grep("^BEAM[0-9]{4}$", h5$ls()$name, value = TRUE)
}

.h5_read_if_present <- function(group, name) {
  if (!group$exists(name)) return(NULL)
  value <- group[[name]][]
  if (is.factor(value)) as.character(value) else value
}

#' Extract GEDI Level 4A footprint metrics
#'
#' @param level4a A [`gedi.level4a-class`] object returned by [readLevel4A()].
#' @param cols Character vector of fields to return. `NULL` returns all
#'   datasets shared by the selected beams.
#' @param quality Logical. If `TRUE`, retain observations accepted by the
#'   current Release 3 quality fields when present:
#'
#'   * `l4a_quality_flag_rel3`
#'   * `degrade_include_flag`
#'   * `elev_highestreturn_outlier_flag`
#'
#'   Legacy Release 2 quality fields remain supported.
#' @param beams Optional character vector of GEDI beam names.
#' @return A [data.table::data.table] with one row per footprint.
#' @export
getLevel4A <- function(level4a, cols = .level4a_default_columns,
                       quality = FALSE, beams = NULL) {
  if (!is(level4a, "gedi.level4a")) {
    stop("'level4a' must be returned by readLevel4A().", call. = FALSE)
  }
  h5 <- level4a@h5
  available_beams <- .gedi_beams(h5)
  if (is.null(beams)) beams <- available_beams
  unknown <- setdiff(beams, available_beams)
  if (length(unknown)) stop("Unknown beam(s): ", paste(unknown, collapse = ", "), call. = FALSE)

  out <- lapply(beams, function(beam) {
    group <- h5[[beam]]
    available <- group$ls()$name
    requested <- cols
    if (is.null(requested)) requested <- available
    requested <- unique(requested)
    values <- list()
    for (field in setdiff(requested, "beam")) {
      value <- .h5_read_if_present(group, field)
      if (!is.null(value)) values[[field]] <- value
    }
    lengths <- lengths(values)
    n <- if (length(lengths)) max(lengths) else 0L
    if (!n) return(data.table::data.table())
    values <- lapply(values, function(value) {
      if (length(value) == 1L && n > 1L) rep(value, n) else value
    })
    bad <- names(values)[lengths(values) != n]
    if (length(bad)) {
      warning("Skipping non-footprint dataset(s) in ", beam, ": ", paste(bad, collapse = ", "))
      values[bad] <- NULL
    }
    dt <- data.table::as.data.table(values)
    dt[, beam := beam]
    data.table::setcolorder(dt, c("beam", setdiff(names(dt), "beam")))
    dt
  })
  ans <- data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  missing_cols <- setdiff(if (is.null(cols)) character() else cols, names(ans))
  if (length(missing_cols)) {
    warning("Unavailable Level 4A field(s): ", paste(missing_cols, collapse = ", "))
  }
  if (isTRUE(quality) && nrow(ans)) {
    keep <- rep(TRUE, nrow(ans))
    if ("l4a_quality_flag_rel3" %in% names(ans)) {
      keep <- keep & ans$l4a_quality_flag_rel3 == 1
    } else if ("l4_quality_flag" %in% names(ans)) {
      keep <- keep & ans$l4_quality_flag == 1
    }
    if ("degrade_include_flag" %in% names(ans)) {
      keep <- keep & ans$degrade_include_flag == 1
    } else if ("degrade_flag" %in% names(ans)) {
      keep <- keep & ans$degrade_flag == 0
    }
    if ("elev_highestreturn_outlier_flag" %in% names(ans)) {
      keep <- keep & ans$elev_highestreturn_outlier_flag == 0
    }
    keep[is.na(keep)] <- FALSE
    ans <- ans[keep]
  }
  ans[]
}

#' Clip GEDI Level 4A data by an extent
#'
#' @param level4A A [`gedi.level4a-class`] object returned by [readLevel4A()]
#'   or a table returned by [getLevel4A()].
#' @param xmin,xmax,ymin,ymax Bounding coordinates in longitude/latitude.
#' @param output Optional output HDF5 path when `level4A` is a
#'   `gedi.level4a` object. A temporary file is used by default.
#' @return A clipped [`gedi.level4a-class`] object for HDF5 input, or a
#'   [data.table::data.table] for table input.
#' @export
clipLevel4A <- function(level4A, xmin, xmax, ymin, ymax, output = "") {
  bounds <- c(xmin, xmax, ymin, ymax)
  if (!is.numeric(bounds) || length(bounds) != 4L || any(!is.finite(bounds))) stop("Bounds must be finite numbers.")
  if (is(level4A, "gedi.level4a")) {
    masks <- .level4a_extent_masks(level4A, xmin, xmax, ymin, ymax)
    return(.clip_level4a_h5(level4A, masks, output))
  }
  if (!inherits(level4A, c("data.table", "data.frame"))) {
    stop("'level4A' must be a gedi.level4a object or an extracted table.")
  }
  lon <- level4A$lon_lowestmode
  lat <- level4A$lat_lowestmode
  keep <- lon >= xmin & lon <= xmax & lat >= ymin & lat <= ymax
  keep[is.na(keep)] <- FALSE
  data.table::as.data.table(level4A)[keep]
}

#' Clip GEDI Level 4A data by polygons
#'
#' @param level4A A [`gedi.level4a-class`] object returned by [readLevel4A()]
#'   or a table returned by [getLevel4A()].
#' @param polygon An `sf`, `sfc`, or `SpatVector` polygon object.
#' @param split_by Optional polygon attribute. For HDF5 input, one clipped
#'   file is returned per unique value. For table input, the value is copied
#'   to `poly_id`.
#' @param output Optional output HDF5 path when `level4A` is a
#'   `gedi.level4a` object. A temporary file is used by default.
#' @return A clipped [`gedi.level4a-class`] object (or a named list when
#'   `split_by` is used) for HDF5 input, or a [data.table::data.table] for
#'   table input.
#' @export
clipLevel4AGeometry <- function(level4A, polygon, split_by = NULL, output = "") {
  poly <- if (inherits(polygon, "SpatVector")) sf::st_as_sf(polygon) else sf::st_as_sf(polygon)
  poly <- sf::st_transform(sf::st_make_valid(poly), 4326)
  if (!is.null(split_by) && !split_by %in% names(poly)) {
    stop("'split_by' is not a polygon attribute.")
  }
  if (is(level4A, "gedi.level4a")) {
    mask_sets <- .level4a_geometry_masks(level4A, poly, split_by)
    if (is.null(split_by)) return(.clip_level4a_h5(level4A, mask_sets[[1L]], output))
    output <- checkOutput(output)
    results <- lapply(names(mask_sets), function(id) {
      path <- sub("\\.h5$", paste0("_", id, ".h5"), output)
      .clip_level4a_h5(level4A, mask_sets[[id]], path)
    })
    names(results) <- names(mask_sets)
    return(results)
  }
  if (!inherits(level4A, c("data.table", "data.frame"))) {
    stop("'level4A' must be a gedi.level4a object or an extracted table.")
  }
  pts <- sf::st_as_sf(as.data.frame(level4A),
    coords = c("lon_lowestmode", "lat_lowestmode"), crs = 4326, remove = FALSE)
  hits <- sf::st_intersects(pts, poly)
  keep <- lengths(hits) > 0L
  ans <- data.table::as.data.table(level4A)[keep]
  if (!is.null(split_by)) {
    ans[["poly_id"]] <- poly[[split_by]][vapply(hits[keep], `[`, integer(1), 1L)]
  }
  ans
}

.level4a_spatial_data <- function(level4a) {
  setNames(lapply(.gedi_beams(level4a@h5), function(beam) {
    group <- level4a@h5[[beam]]
    if (!group$exists("lon_lowestmode") || !group$exists("lat_lowestmode")) {
      stop("Level 4A beam ", beam, " does not contain lon_lowestmode/lat_lowestmode.", call. = FALSE)
    }
    data.frame(
      longitude = group[["lon_lowestmode"]][],
      latitude = group[["lat_lowestmode"]][]
    )
  }), .gedi_beams(level4a@h5))
}

.level4a_extent_masks <- function(level4a, xmin, xmax, ymin, ymax) {
  spatial <- .level4a_spatial_data(level4a)
  masks <- lapply(spatial, function(x) {
    keep <- x$longitude >= xmin & x$longitude <= xmax &
      x$latitude >= ymin & x$latitude <= ymax
    keep[is.na(keep)] <- FALSE
    which(keep)
  })
  if (!any(lengths(masks))) stop("The clipping region does not intersect the data!", call. = FALSE)
  masks
}

.level4a_geometry_masks <- function(level4a, polygon, split_by = NULL) {
  spatial <- .level4a_spatial_data(level4a)
  ids <- if (is.null(split_by)) "1" else unique(as.character(polygon[[split_by]]))
  result <- setNames(lapply(ids, function(x) setNames(vector("list", length(spatial)), names(spatial))), ids)
  for (beam in names(spatial)) {
    points <- sf::st_as_sf(spatial[[beam]], coords = c("longitude", "latitude"), crs = 4326)
    hits <- sf::st_intersects(points, polygon)
    if (is.null(split_by)) {
      result[[1L]][[beam]] <- which(lengths(hits) > 0L)
    } else {
      for (id in ids) {
        polygon_rows <- which(as.character(polygon[[split_by]]) == id)
        result[[id]][[beam]] <- which(vapply(hits, function(z) any(z %in% polygon_rows), logical(1)))
      }
    }
  }
  result <- result[vapply(result, function(x) any(lengths(x)), logical(1))]
  if (!length(result)) stop("The clipping polygon does not intersect the data!", call. = FALSE)
  result
}

.subset_level4a_dataset <- function(dataset, mask, shot_count) {
  dims <- dataset$dims
  if (!length(dims) || !any(dims == shot_count)) return(dataset[])
  index <- lapply(dims, seq_len)
  index[[which(dims == shot_count)[1L]]] <- mask
  do.call(`[`, c(list(dataset[]), index, list(drop = FALSE)))
}

.clip_level4a_h5 <- function(level4a, masks, output = "") {
  output <- checkOutput(output)
  selected_beams <- names(masks)[lengths(masks) > 0L]
  if (!length(selected_beams)) stop("The clipping region does not intersect the data!", call. = FALSE)

  source <- level4a@h5
  destination <- hdf5r::H5File$new(output, mode = "w")
  completed <- FALSE
  on.exit(if (!completed) try(destination$close_all(), silent = TRUE), add = TRUE)
  copyRootAttributes(source, destination)

  groups <- hdf5r::list.groups(source)
  keep_path <- function(path) {
    top <- strsplit(sub("^/", "", path), "/", fixed = TRUE)[[1L]][1L]
    !grepl("^BEAM[0-9]{4}$", top) || top %in% selected_beams
  }
  groups <- groups[vapply(groups, keep_path, logical(1))]
  groups <- groups[order(lengths(strsplit(groups, "/", fixed = TRUE)))]
  for (group in groups) {
    hdf5r::createGroup(destination, group)
    createAttributesWithinGroup(source, destination, group)
  }

  datasets <- hdf5r::list.datasets(source)
  datasets <- datasets[vapply(datasets, keep_path, logical(1))]
  for (path in datasets) {
    top <- strsplit(sub("^/", "", path), "/", fixed = TRUE)[[1L]][1L]
    dataset <- source[[path]]
    value <- if (top %in% selected_beams) {
      shot_count <- source[[top]][["shot_number"]]$dims[1L]
      .subset_level4a_dataset(dataset, masks[[top]], shot_count)
    } else {
      dataset[]
    }
    hdf5r::createDataSet(destination, path, value, dtype = dataset$get_type())
    createAttributesWithinGroup(source, destination, path)
  }
  destination$close_all()
  completed <- TRUE
  readLevel4A(output)
}

#' Grid GEDI Level 4A footprint metrics
#'
#' @param level4A A table returned by [getLevel4A()].
#' @param metric Metric column name.
#' @param fun Aggregation function.
#' @param res Output resolution in decimal degrees.
#' @param ... Additional arguments passed to `fun`.
#' @return A [`terra::SpatRaster-class`].
#' @export
gridStatsLevel4A <- function(level4A, metric = "agbd", fun = mean, res = 0.01, ...) {
  .grid_gedi_points(level4A, "lon_lowestmode", "lat_lowestmode", metric, fun, res, ...)
}

#' Polygon statistics for GEDI Level 4A footprints
#'
#' @param level4A A table returned by [getLevel4A()].
#' @param polygon An `sf`, `sfc`, or `SpatVector` polygon object.
#' @param metric Metric column name.
#' @param fun Aggregation function.
#' @param id Optional polygon ID field.
#' @param ... Additional arguments passed to `fun`.
#' @return A data table with one row per polygon.
#' @export
polyStatsLevel4A <- function(level4A, polygon, metric = "agbd", fun = mean,
                             id = NULL, ...) {
  .poly_gedi_points(level4A, polygon, "lon_lowestmode", "lat_lowestmode", metric, fun, id, ...)
}

#' Rasterize GEDI Level 4A footprints
#'
#' @inheritParams gridStatsLevel4A
#' @param filename Optional output GeoTIFF path.
#' @param overwrite Logical; overwrite `filename` when it exists.
#' @return A [`terra::SpatRaster-class`].
#' @export
rasterizeLevel4A <- function(level4A, metric = "agbd", fun = mean, res = 0.01,
                             filename = "", overwrite = FALSE, ...) {
  r <- gridStatsLevel4A(level4A, metric, fun, res, ...)
  if (nzchar(filename)) terra::writeRaster(r, filename, overwrite = overwrite) else r
}

#' Plot GEDI Level 4A footprint biomass
#'
#' @param level4A A table returned by [getLevel4A()].
#' @param metric Metric column to display.
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return Invisibly returns `level4A`.
#' @export
plotLevel4A <- function(level4A, metric = "agbd", ...) {
  if (!all(c("lon_lowestmode", "lat_lowestmode", metric) %in% names(level4A))) {
    stop("Required coordinate or metric columns are missing.")
  }
  metric_values <- as.numeric(level4A[[metric]])
  bins <- if (length(unique(metric_values[is.finite(metric_values)])) < 2L) {
    rep(50L, length(metric_values))
  } else {
    as.integer(cut(metric_values, 100, include.lowest = TRUE))
  }
  graphics::plot(level4A$lon_lowestmode, level4A$lat_lowestmode,
    col = grDevices::hcl.colors(100, "Viridis")[bins],
    xlab = "Longitude", ylab = "Latitude", ...)
  invisible(level4A)
}
