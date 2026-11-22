#' Accumulate the screened footprints of a cloud granule on a 1 km grid
#' @param entry One manifest row with product, version, metric and url.
#' @param from,to Optional UTC date bounds. The upper bound is exclusive.
#' @param region Optional output of [readGEDIArchiveRegion()].
#' @return A data.table of metric sums and footprint counts by EASE-Grid cell.
#' @export
accumulateGEDIGranule <- function(entry, from = NULL, to = NULL, region = NULL) {
  if (xor(is.null(from), is.null(to))) stop("Supply both UTC date bounds.")
  bounds <- if (!is.null(from)) {
    as.numeric(difftime(as.POSIXct(c(from, to), tz = "UTC"),
                        as.POSIXct("2018-01-01", tz = "UTC"), units = "secs"))
  }
  if (!is.null(bounds) && (anyNA(bounds) || bounds[1L] >= bounds[2L])) {
    stop("The UTC date bounds must be increasing.")
  }
  metric <- match.arg(entry$metric, c("rh98", "agbd"))
  cloud <- .read_gedi_cloud(entry$url, product = entry$product)
  on.exit(close(cloud))
  power <- c("BEAM0101", "BEAM0110", "BEAM1000", "BEAM1011")
  v3 <- as.integer(entry$version) == 3L
  quality <- if (metric == "rh98") {
    if (v3) "l2a_quality_flag_rel3" else "quality_flag"
  } else {
    if (v3) "l4a_quality_flag_rel3" else "l4_quality_flag"
  }
  pieces <- lapply(intersect(power, names(cloud)), function(beam) {
    group <- cloud[[beam]]
    solar <- group[["solar_elevation"]][]
    use <- which(is.finite(solar) & solar >= -90 & solar < 0 &
                   group[[quality]][] == 1 & group[["degrade_flag"]][] == 0)
    if (length(use) && !is.null(bounds)) {
      time <- group[["delta_time"]][][use]
      use <- use[which(time >= bounds[1L] & time < bounds[2L])]
    }
    if (!length(use)) return(NULL)
    lon <- group[["lon_lowestmode"]][][use]
    lat <- group[["lat_lowestmode"]][][use]
    keep <- which(data.table::between(lon, -180, 180) & data.table::between(lat, -52, 52))
    if (!is.null(region)) {
      box <- region$bounds
      keep <- keep[data.table::between(lon[keep], box[1L], box[2L]) &
                     data.table::between(lat[keep], box[3L], box[4L])]
    }
    if (!length(keep)) return(NULL)
    xy <- terra::project(cbind(lon[keep], lat[keep]), "EPSG:4326", "EPSG:6933")
    if (!is.null(region)) {
      forest <- terra::extract(region$forest, xy)[[1L]]
      covered <- which(forest == 1L)
      if (!length(covered)) return(NULL)
      points <- terra::vect(xy[covered, , drop = FALSE], type = "points", crs = "EPSG:6933")
      inside <- covered[terra::is.related(points, region$area, "intersects")]
      keep <- keep[inside]
      xy <- xy[inside, , drop = FALSE]
    }
    if (!length(keep)) return(NULL)
    use <- use[keep]
    first <- min(use) - 1L
    value <- if (metric == "rh98") {
      builtins <- reticulate::import_builtins(convert = FALSE)
      rh <- group[["rh"]]$ds
      selection <- reticulate::tuple(builtins$slice(first, max(use)), 98L)
      reticulate::py_to_r(rh$`__getitem__`(selection))[use - first]
    } else {
      group[["agbd"]][][use]
    }
    shots <- data.table::data.table(col_km = floor(xy[, 1L] / 1000),
                                    row_km = floor(xy[, 2L] / 1000), value = value)
    shots <- shots[is.finite(value) & value >= 0]
    shots[, .(total = sum(value), shots = .N), by = .(col_km, row_km)]
  })
  .archive_cell_totals(data.table::rbindlist(pieces))
}

.accumulate_gedi_granule_spec <- function(spec) {
  region <- if (is.null(spec$region_root)) spec$region else {
    readGEDIArchiveRegion(spec$region_root)
  }
  accumulateGEDIGranule(spec$entry, spec$from, spec$to, region)
}

#' Accumulate multiple cloud granules on a common 1 km grid
#'
#' Each worker opens and closes one cloud granule independently. Partial cell
#' totals are combined after all workers finish, so the result is identical to
#' summing sequential granule results.
#'
#' @param entries A data frame or data table of manifest rows accepted by
#'   [accumulateGEDIGranule()].
#' @inheritParams accumulateGEDIGranule
#' @param parallel Logical; process granules concurrently.
#' @param workers Number of workers. `NULL` uses physical cores minus one,
#'   limited by the number of granules.
#' @return A data table of combined metric sums and footprint counts by
#'   EASE-Grid cell. The `workers` attribute records actual concurrency.
#' @export
accumulateGEDIGranules <- function(entries, from = NULL, to = NULL,
                                   region = NULL, parallel = FALSE,
                                   workers = NULL) {
  entries <- data.table::as.data.table(entries)
  if (!nrow(entries)) return(.archive_cell_totals(data.table::data.table()))
  region_root <- NULL
  if (isTRUE(parallel) && !is.null(region)) {
    region_root <- region$root
    if (is.null(region_root) || !dir.exists(region_root)) {
      stop(
        "A parallel regional batch requires a region returned by readGEDIArchiveRegion().",
        call. = FALSE
      )
    }
  }
  specs <- lapply(seq_len(nrow(entries)), function(i) list(
    entry = entries[i], from = from, to = to,
    region = if (isTRUE(parallel)) NULL else region,
    region_root = region_root
  ))
  pieces <- .gedi_parallel_lapply(
    specs,
    function(spec) getFromNamespace(".accumulate_gedi_granule_spec", "rGEDI")(spec),
    parallel = parallel, workers = workers
  )
  used_workers <- attr(pieces, "workers")
  answer <- .archive_cell_totals(data.table::rbindlist(pieces))
  attr(answer, "workers") <- used_workers
  answer[]
}

#' Read the footprint boundary and forest mask for a regional archive
#' @param root Archive directory with the optional Amazon region assets.
#' @return A regional filter, or NULL for an archive without region assets.
#' @export
readGEDIArchiveRegion <- function(root) {
  area_file <- file.path(root, "region", "amazonia_core.gpkg")
  mask_file <- file.path(root, "region", "forest_2021.tif")
  if (!file.exists(area_file) && !file.exists(mask_file)) return(NULL)
  if (!file.exists(area_file) || !file.exists(mask_file)) stop("The regional archive needs a boundary and forest mask.")
  area <- terra::vect(area_file)
  list(root = normalizePath(root, winslash = "/", mustWork = TRUE),
       bounds = as.vector(terra::ext(area)),
       area = terra::makeValid(terra::project(area, "EPSG:6933")),
       forest = terra::rast(mask_file))
}

.archive_cell_totals <- function(cells) {
  if (!nrow(cells)) return(data.table::data.table(col_km = numeric(), row_km = numeric(),
                                                  total = numeric(), shots = numeric()))
  cells[, .(total = sum(total), shots = sum(shots)), by = .(col_km, row_km)]
}

.archive_raster_tile <- function(cells) {
  x0 <- 250 * floor(cells$col_km[1L] / 250)
  y0 <- 250 * floor(cells$row_km[1L] / 250)
  grid <- terra::rast(ncols = 250, nrows = 250, xmin = x0 * 1000, xmax = (x0 + 250) * 1000,
                      ymin = y0 * 1000, ymax = (y0 + 250) * 1000, crs = "EPSG:6933")
  cells[, `:=`(x = 1000 * (col_km + .5), y = 1000 * (row_km + .5))]
  totals <- rasterize(cells, "total", fun = "sum", lon = "x", lat = "y",
                          template = grid, input_crs = "EPSG:6933")
  counts <- rasterize(cells, "shots", fun = "sum", lon = "x", lat = "y",
                          template = grid, input_crs = "EPSG:6933")
  stopifnot(terra::global(counts, "sum", na.rm = TRUE)[1L, 1L] == sum(cells$shots))
  stats <- c(totals / counts, counts, totals)
  names(stats) <- c("mean", "shot_count", "sum")
  stats
}

.archive_boundary <- function(root, index, target, year) {
  out <- file.path(root, "by_year", "boundary")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  cache <- file.path(out, sprintf("%s_%d.rds", target, year))
  if (file.exists(cache)) return(readRDS(cache))
  december <- index[metric == target & month == sprintf("%d-12", year - 1L)]$file
  next_year <- sprintf("%d-01-01", year)
  nodes <- data.table::rbindlist(lapply(file.path(root, december), data.table::fread,
                                      colClasses = "character"))
  if (nrow(nodes)) nodes <- nodes[substr(time_end, 1L, 10L) >= next_year]
  later <- if (nrow(nodes)) {
    # A granule can span midnight. Use shot time for contributions after midnight.
    end <- as.character(as.Date(max(substr(nodes$time_end, 1L, 10L))) + 1L)
    region <- readGEDIArchiveRegion(root)
    data.table::rbindlist(lapply(seq_len(nrow(nodes)), function(i) {
      accumulateGEDIGranule(nodes[i], next_year, end, region)
    }))
  } else {
    data.table::data.table()
  }
  answer <- .archive_cell_totals(later)
  temp <- tempfile(tmpdir = out)
  saveRDS(answer, temp)
  stopifnot(file.rename(temp, cache))
  answer
}

.archive_read_cell_file <- function(path) {
  input <- gzfile(path, "rt")
  on.exit(close(input))
  data.table::fread(text = paste(readLines(input), collapse = "\n"))
}

.archive_write_year <- function(root, index, target, year, tier, marker_names,
                                previous = NULL) {
  output <- file.path(root, "by_year", tier, target, as.character(year))
  check <- file.path(output, "complete.csv")
  if (file.exists(check)) return(data.table::fread(check))
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  prior <- if (is.null(previous)) character() else list.files(previous, "[.]tif$", full.names = TRUE)
  if (tier == "annual") {
    folders <- list.dirs(file.path(root, "partials"), recursive = FALSE)
    folders <- folders[startsWith(basename(folders), paste0(target, "_"))]
    partials <- lapply(folders, function(folder) {
      paths <- list.files(folder, "[.]csv[.]gz$", full.names = TRUE)
      paths[sub("[.]csv[.]gz$", "", basename(paths)) %in% marker_names]
    })
    names(partials) <- basename(folders)
    incoming <- .archive_boundary(root, index, target, year)
    outgoing <- .archive_boundary(root, index, target, year + 1L)
    outgoing[, `:=`(total = -total, shots = -shots)]
    changes <- data.table::rbindlist(list(incoming, outgoing))
    if (nrow(changes)) {
      changes[, tile := paste(target, floor(col_km / 250), floor(row_km / 250), sep = "_")]
    }
    keys <- union(names(partials)[lengths(partials) > 0L], changes$tile)
  } else {
    annual_folder <- file.path(root, "by_year", "annual", target, as.character(year))
    annual <- list.files(annual_folder, "[.]tif$", full.names = TRUE)
    keys <- sub("[.]tif$", "", union(basename(annual), basename(prior)))
  }
  occupied <- shots <- 0
  for (key in keys) {
    if (tier == "annual") {
      cells <- if (nrow(changes)) changes[tile == key, .(col_km, row_km, total, shots)] else changes
      paths <- partials[[key]]
      for (chunk in split(paths, ceiling(seq_along(paths) / 32L))) {
        cells <- .archive_cell_totals(data.table::rbindlist(
          c(list(cells), lapply(chunk, .archive_read_cell_file))))
      }
      cells <- .archive_cell_totals(cells)
      if (any(cells$shots < 0)) stop("A year-boundary correction exceeded the original count.")
      cells <- cells[shots > 0]
      if (!nrow(cells)) next
      result <- .archive_raster_tile(cells)
      occupied <- occupied + nrow(cells)
      shots <- shots + sum(cells$shots)
    } else {
      files <- c(file.path(annual_folder, paste0(key, ".tif")),
                 prior[basename(prior) == paste0(key, ".tif")])
      layers <- lapply(files[file.exists(files)], terra::rast)
      sum_cells <- function(band) {
        Reduce(`+`, lapply(layers, function(r) terra::ifel(is.na(r[[band]]), 0, r[[band]])))
      }
      counts <- sum_cells("shot_count")
      totals <- sum_cells("sum")
      occupied <- occupied + terra::global(counts > 0, "sum")[1L, 1L]
      shots <- shots + terra::global(counts, "sum")[1L, 1L]
      result <- c(terra::ifel(counts > 0, totals / counts, NA),
                  terra::ifel(counts > 0, counts, NA), terra::ifel(counts > 0, totals, NA))
      names(result) <- c("mean", "shot_count", "sum")
    }
    terra::writeRaster(result, file.path(output, paste0(key, ".tif")), overwrite = TRUE,
                        datatype = "FLT8S", gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
  }
  files <- list.files(output, "[.]tif$", full.names = TRUE)
  if (length(files)) terra::vrt(files, filename = file.path(output, paste0(target, "_1km.vrt")),
                               overwrite = TRUE)
  stats <- data.table::data.table(metric = target, year = year, tier = tier,
                                  occupied_cells = occupied, shots = shots)
  temp <- tempfile(tmpdir = output)
  data.table::fwrite(stats, temp)
  stopifnot(file.rename(temp, check))
  stats
}

#' Render complete calendar years and their cumulative footprint grids
#' @param root Archive directory containing its index, manifests, markers and partials.
#' @param batch_size Batch size used during the cloud pass.
#' @return A list with the archive readiness flag and a coverage data.table.
#' @export
assembleGEDIArchiveYears <- function(root, batch_size = 16L) {
  index <- data.table::fread(file.path(root, "archive_index.csv"))
  index[, month := sub("^manifest_archive_(rh98|agbd)_([0-9]{4}-[0-9]{2})[.]csv$", "\\2", file)]
  if (any(!grepl("^[0-9]{4}-[0-9]{2}$", index$month))) stop("Unrecognized archive manifest names.")
  index[, year := as.integer(substr(month, 1L, 4L))]
  index[, markers := lapply(seq_len(.N), function(i) {
    ids <- seq_len(ceiling(granules[i] / batch_size))
    sprintf("%s_size%d_batch%06d", sub("[.]csv$", "", file[i]), batch_size, ids)
  })]
  done <- list.files(file.path(root, "done"))
  index[, ready := vapply(markers, function(m) all(m %in% done), logical(1L))]
  results <- list()
  for (target in unique(index$metric)) {
    groups <- index[metric == target, .(ready = all(ready), markers = list(unlist(markers))), by = year]
    data.table::setorder(groups, year)
    previous <- NULL
    can_accumulate <- TRUE
    for (i in seq_len(nrow(groups))) {
      yr <- groups$year[i]
      can_accumulate <- can_accumulate && groups$ready[i]
      if (!groups$ready[i]) next
      results[[length(results) + 1L]] <- .archive_write_year(
        root, index, target, yr, "annual", groups$markers[[i]])
      if (!can_accumulate) next
      results[[length(results) + 1L]] <- .archive_write_year(
        root, index, target, yr, "cumulative", character(), previous)
      previous <- file.path(root, "by_year", "cumulative", target, as.character(yr))
    }
  }
  coverage <- data.table::rbindlist(results)
  if (nrow(coverage)) {
    data.table::fwrite(coverage, file.path(root, "by_year", "annual_coverage.csv"))
  }
  list(complete = all(index$ready), coverage = coverage)
}
