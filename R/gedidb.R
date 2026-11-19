# gediDB optional backend -------------------------------------------------

.gedidb_endpoint <- "https://s3.gfz.de"
.gedidb_buckets <- c(
  `003` = "dog-ext--gedidb--gedi-l2-l4-v003",
  `002` = "dog-ext.gedidb.gedi-l2-l4-v002.0"
)
.gedidb_provider_cache <- new.env(parent = emptyenv())

.gedidb_version <- function(version = NULL) {
  if (is.null(version)) return("003")
  version <- sub("^[Vv]", "", as.character(version)[1L])
  if (version %in% c("3", "03", "003")) return("003")
  if (version %in% c("2", "02", "002", "2.0", "002.0", "2.1", "002.1")) return("002")
  stop("gediDB supports GEDI Release 3 ('003') and Release 2 ('002').", call. = FALSE)
}

.gedidb_config <- function(version = NULL) {
  version <- .gedidb_version(version)
  bucket <- unname(.gedidb_buckets[[version]])
  list(
    version = version,
    endpoint = .gedidb_endpoint,
    bucket = bucket,
    address = paste0("s3://", bucket, "/array_uri")
  )
}

#' Test availability of the public gediDB service
#'
#' Checks that the optional Python module is installed and that the selected
#' public GFZ gediDB bucket is reachable. The default check is deliberately
#' lightweight. Set `deep = TRUE` to open the TileDB array and read its
#' variable metadata; that check can take substantially longer.
#'
#' @param version GEDI database release, `"003"` or `"002"`.
#' @param deep Logical; also open the TileDB array and read its metadata.
#' @param timeout Timeout in seconds for the lightweight endpoint request.
#' @return A single logical value. A diagnostic explanation is stored in its
#'   `reason` attribute.
#' @export
gedidb_available <- function(version = "003", deep = FALSE, timeout = 5) {
  config <- .gedidb_config(version)
  backend <- getOption("rGEDI.gedidb_backend")
  if (!is.null(backend)) {
    result <- tryCatch(isTRUE(backend$available(config, deep = deep)), error = function(e) FALSE)
    attr(result, "reason") <- if (result) "available" else "mock backend unavailable"
    return(result)
  }

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    result <- FALSE
    attr(result, "reason") <- "the optional R package 'reticulate' is not installed"
    return(result)
  }
  if (!reticulate::py_module_available("gedidb")) {
    result <- FALSE
    attr(result, "reason") <- "the optional Python package 'gedidb' is not installed"
    return(result)
  }

  url <- paste0(config$endpoint, "/", config$bucket)
  response <- tryCatch(
    curl::curl_fetch_memory(url, handle = curl::new_handle(
      connecttimeout = timeout, timeout = timeout,
      followlocation = TRUE, useragent = "rGEDI/gedidb-availability-check"
    )),
    error = function(e) e
  )
  if (inherits(response, "error")) {
    result <- FALSE
    attr(result, "reason") <- conditionMessage(response)
    return(result)
  }
  # Public S3-compatible buckets commonly return 403 to a bucket-level request
  # while still allowing anonymous reads of TileDB objects. A missing bucket
  # returns 404 and must not be treated as available.
  if (!response$status_code %in% c(200L, 204L, 301L, 302L, 403L)) {
    result <- FALSE
    attr(result, "reason") <- paste("endpoint returned HTTP", response$status_code)
    return(result)
  }

  if (isTRUE(deep)) {
    check <- tryCatch({
      provider <- .gedidb_provider(config, cache = FALSE)
      on.exit(try(provider$close(), silent = TRUE), add = TRUE)
      variables <- provider$get_available_variables()
      reticulate::py_len(variables) > 0L
    }, error = function(e) e)
    if (inherits(check, "error") || !isTRUE(check)) {
      result <- FALSE
      attr(result, "reason") <- if (inherits(check, "error")) conditionMessage(check) else "TileDB metadata are unavailable"
      return(result)
    }
  }

  result <- TRUE
  attr(result, "reason") <- "available"
  result
}

.new_gedidb_reference <- function(product, bbox, daterange = NULL, version = NULL,
                                  checked = FALSE) {
  config <- .gedidb_config(version)
  structure(
    config$address,
    class = c("gedidb_reference", "character"),
    product = product,
    bbox = bbox,
    daterange = daterange,
    version = config$version,
    endpoint = config$endpoint,
    bucket = config$bucket,
    server_checked = isTRUE(checked)
  )
}

#' @export
print.gedidb_reference <- function(x, ...) {
  bbox <- attr(x, "bbox")
  dates <- attr(x, "daterange")
  cat("<gediDB remote reference>\n", sep = "")
  cat("Product: ", attr(x, "product"), "\n", sep = "")
  cat("Release: ", attr(x, "version"), "\n", sep = "")
  cat("Address: ", unclass(x)[1L], "\n", sep = "")
  cat("Bounding box: ", paste(unname(bbox), collapse = ", "), "\n", sep = "")
  if (!is.null(dates)) cat("Dates: ", paste(dates, collapse = " to "), "\n", sep = "")
  invisible(x)
}

.open_gedidb_reference <- function(x, product) {
  if (!inherits(x, "gedidb_reference")) stop("Expected a gediDB remote reference.", call. = FALSE)
  if (!identical(attr(x, "product"), product)) {
    stop("The gediDB reference contains ", attr(x, "product"), ", not ", product, ".", call. = FALSE)
  }
  suffix <- switch(product, GEDI02_A = "level2a", GEDI02_B = "level2b", GEDI04_A = "level4a")
  structure(list(reference = x), class = c(paste0("gedidb_", suffix), "gedidb_dataset"))
}

.new_gedidb_query <- function(dataset, variables, filters = list()) {
  structure(
    list(dataset = dataset, variables = unique(as.character(variables)), filters = filters),
    class = "gedidb_query"
  )
}

#' @export
print.gedidb_query <- function(x, ...) {
  reference <- x$dataset$reference
  cat("<lazy gediDB query>\n")
  cat("Product: ", attr(reference, "product"), "\n", sep = "")
  cat("Release: ", attr(reference, "version"), "\n", sep = "")
  cat("Variables: ", paste(x$variables, collapse = ", "), "\n", sep = "")
  invisible(x)
}

.gedidb_provider <- function(config, cache = TRUE) {
  key <- paste(config$endpoint, config$bucket, sep = "|")
  if (isTRUE(cache) && exists(key, envir = .gedidb_provider_cache, inherits = FALSE)) {
    return(get(key, envir = .gedidb_provider_cache, inherits = FALSE))
  }
  gdb <- reticulate::import("gedidb", convert = FALSE)
  provider <- gdb$GEDIProvider(
    storage_type = "s3",
    s3_bucket = config$bucket,
    url = config$endpoint,
    region = "eu-central-1",
    config_overrides = reticulate::dict(
      "vfs.s3.connect_timeout_ms" = "10000",
      "vfs.s3.request_timeout_ms" = "120000",
      "py.init_buffer_bytes" = as.character(16 * 1024^2)
    )
  )
  if (isTRUE(cache)) assign(key, provider, envir = .gedidb_provider_cache)
  provider
}

#' Close cached gediDB connections
#'
#' Releases the cached Python TileDB providers used to make repeated queries
#' fast. New providers are created automatically by the next query.
#'
#' @return No return value, called for the side effect of closing cached
#'   providers.
#' @export
gedidb_disconnect <- function() {
  keys <- ls(.gedidb_provider_cache, all.names = TRUE)
  for (key in keys) {
    provider <- get(key, envir = .gedidb_provider_cache, inherits = FALSE)
    try(provider$close(), silent = TRUE)
    rm(list = key, envir = .gedidb_provider_cache)
  }
  invisible(NULL)
}

.gedidb_reference_config <- function(reference) {
  list(
    version = attr(reference, "version"), endpoint = attr(reference, "endpoint"),
    bucket = attr(reference, "bucket"), address = unclass(reference)[1L]
  )
}

.gedidb_query <- function(dataset, variables, filters = list()) {
  reference <- dataset$reference
  backend <- getOption("rGEDI.gedidb_backend")
  if (!is.null(backend)) {
    result <- backend$query(reference, variables = variables, filters = filters)
  } else {
    if (!isTRUE(attr(reference, "server_checked"))) {
      available <- gedidb_available(attr(reference, "version"), deep = FALSE)
      if (!isTRUE(available)) {
        stop("gediDB is unavailable: ", attr(available, "reason"), call. = FALSE)
      }
    }
    config <- .gedidb_reference_config(reference)
    provider <- .gedidb_provider(config)
    bbox <- attr(reference, "bbox")
    dates <- attr(reference, "daterange")
    args <- list(
      variables = as.list(variables),
      lat_min = unname(bbox["ymin"]), lat_max = unname(bbox["ymax"]),
      lon_min = unname(bbox["xmin"]), lon_max = unname(bbox["xmax"])
    )
    if (!is.null(dates)) {
      args$start_time <- as.integer(as.Date(dates[1L]))
      args$end_time <- as.integer(as.Date(dates[2L]))
    }
    args <- c(args, filters)
    result <- do.call(provider$query_dataframe, args)
    if (!is.null(result) && !reticulate::py_is_null_xptr(result)) {
      columns <- reticulate::py_to_r(result$columns$tolist())
      if ("shot_number" %in% columns) {
        # pandas uint64 cannot be represented exactly by an R double. Keep the
        # identifier as text, matching its role as an ID rather than a metric.
        result[["shot_number"]] <- result[["shot_number"]]$astype("string")
      }
      result <- reticulate::py_to_r(result)
    }
  }

  if (is.null(result)) result <- data.frame()
  result <- data.table::as.data.table(result)
  renames <- c(
    latitude = "lat_lowestmode", longitude = "lon_lowestmode",
    beam_name = "beam", rh_p98 = "rh98"
  )
  present <- intersect(names(renames), names(result))
  if (length(present)) data.table::setnames(result, present, unname(renames[present]))
  profile_names <- grep("^rh_p[0-9]+$", names(result), value = TRUE)
  if (length(profile_names)) {
    data.table::setnames(result, profile_names, sub("^rh_p", "rh", profile_names))
  }
  if ("time" %in% names(result) && is.numeric(result$time)) {
    result[, time := as.Date(time, origin = "1970-01-01")]
  }
  attr(result, "gedi_source") <- "gedidb"
  attr(result, "gedi_release") <- attr(reference, "version")
  attr(result, "gedidb_address") <- unclass(reference)[1L]
  result[]
}

.gedidb_grid_tiles <- function(reference, tile_size) {
  tile_size <- as.numeric(tile_size)[1L]
  if (!is.finite(tile_size) || tile_size <= 0 || tile_size > 360) {
    stop("'tile_size' must be a positive number no greater than 360 degrees.", call. = FALSE)
  }
  bbox <- attr(reference, "bbox")
  nx <- max(1L, ceiling((bbox["xmax"] - bbox["xmin"]) / tile_size - 1e-10))
  ny <- max(1L, ceiling((bbox["ymax"] - bbox["ymin"]) / tile_size - 1e-10))
  xbreaks <- pmin(bbox["xmin"] + (0:nx) * tile_size, bbox["xmax"])
  ybreaks <- pmin(bbox["ymin"] + (0:ny) * tile_size, bbox["ymax"])
  grid <- expand.grid(ix = seq_len(length(xbreaks) - 1L),
                      iy = seq_len(length(ybreaks) - 1L))
  lapply(seq_len(nrow(grid)), function(i) {
    ix <- grid$ix[i]
    iy <- grid$iy[i]
    list(
      bbox = c(xmin = xbreaks[ix], ymin = ybreaks[iy],
               xmax = xbreaks[ix + 1L], ymax = ybreaks[iy + 1L]),
      last_x = ix == length(xbreaks) - 1L,
      last_y = iy == length(ybreaks) - 1L
    )
  })
}

.gedidb_grid_tile <- function(spec) {
  query <- spec$query
  reference <- query$dataset$reference
  attr(reference, "bbox") <- spec$tile$bbox
  dataset <- .open_gedidb_reference(reference, attr(reference, "product"))
  values <- .gedidb_query(dataset, query$variables, query$filters)
  metric <- spec$metric
  needed <- c("lon_lowestmode", "lat_lowestmode", metric)
  if (!nrow(values)) {
    return(data.table::data.table(col = numeric(), row = numeric(),
                                  total = numeric(), shots = integer()))
  }
  if (!all(needed %in% names(values))) {
    stop("gediDB did not return the requested metric and coordinates: ",
         paste(setdiff(needed, names(values)), collapse = ", "), call. = FALSE)
  }
  box <- spec$tile$bbox
  keep_x <- values$lon_lowestmode >= box["xmin"] &
    if (spec$tile$last_x) values$lon_lowestmode <= box["xmax"] else values$lon_lowestmode < box["xmax"]
  keep_y <- values$lat_lowestmode >= box["ymin"] &
    if (spec$tile$last_y) values$lat_lowestmode <= box["ymax"] else values$lat_lowestmode < box["ymax"]
  values <- values[which(keep_x & keep_y & is.finite(values[[metric]]))]
  if (!nrow(values)) {
    return(data.table::data.table(col = numeric(), row = numeric(),
                                  total = numeric(), shots = integer()))
  }
  xy <- terra::project(
    as.matrix(values[, .(lon_lowestmode, lat_lowestmode)]),
    "EPSG:4326", spec$output_crs
  )
  cells <- data.table::data.table(
    col = floor(xy[, 1L] / spec$res), row = floor(xy[, 2L] / spec$res),
    value = as.numeric(values[[metric]])
  )
  cells[, .(total = sum(value), shots = .N), by = .(col, row)]
}

.gedidb_grid_template <- function(reference, res, output_crs) {
  bbox <- attr(reference, "bbox")
  polygon <- terra::as.polygons(terra::ext(
    bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"]
  ), crs = "EPSG:4326")
  projected <- terra::project(polygon, output_crs)
  extent <- as.vector(terra::ext(projected))
  extent <- c(
    floor(extent[1L] / res) * res, ceiling(extent[2L] / res) * res,
    floor(extent[3L] / res) * res, ceiling(extent[4L] / res) * res
  )
  terra::rast(terra::ext(extent), resolution = res, crs = output_crs)
}

.gedidb_grid_query <- function(query, metric, res = 1000, fun = mean,
                               output_crs = "EPSG:6933", parallel = FALSE,
                               workers = NULL, tile_size = 5, filename = "") {
  reducer <- if (is.character(fun)) tolower(fun[1L]) else if (identical(fun, mean)) {
    "mean"
  } else if (identical(fun, sum)) {
    "sum"
  } else {
    stop("Lazy gediDB rasterization currently supports mean and sum.", call. = FALSE)
  }
  if (!reducer %in% c("mean", "sum")) {
    stop("Lazy gediDB rasterization currently supports mean and sum.", call. = FALSE)
  }
  res <- as.numeric(res)[1L]
  if (!is.finite(res) || res <= 0) stop("'res' must be a positive number.", call. = FALSE)
  metric <- as.character(metric)[1L]
  reference <- query$dataset$reference
  tiles <- .gedidb_grid_tiles(reference, tile_size)
  specs <- lapply(tiles, function(tile) list(
    query = query, tile = tile, metric = metric, res = res,
    output_crs = output_crs
  ))
  started <- proc.time()[["elapsed"]]
  if (isTRUE(parallel) && length(specs) > 1L) {
    if (is.null(workers)) {
      detected <- parallel::detectCores(logical = FALSE)
      workers <- max(1L, min(length(specs), if (is.na(detected)) 2L else detected - 1L))
    }
    workers <- max(1L, min(length(specs), as.integer(workers)[1L]))
    if (.Platform$OS.type == "windows") {
      cluster <- parallel::makePSOCKcluster(workers)
      on.exit(parallel::stopCluster(cluster), add = TRUE)
      parallel::clusterCall(cluster, function() loadNamespace("rGEDI"))
      backend <- getOption("rGEDI.gedidb_backend")
      if (!is.null(backend)) {
        parallel::clusterCall(cluster, function(value) {
          options(rGEDI.gedidb_backend = value)
          NULL
        }, backend)
      }
      pieces <- parallel::parLapply(
        cluster, specs,
        function(spec) getFromNamespace(".gedidb_grid_tile", "rGEDI")(spec)
      )
    } else {
      pieces <- parallel::mclapply(specs, .gedidb_grid_tile, mc.cores = workers)
    }
  } else {
    workers <- 1L
    pieces <- lapply(specs, .gedidb_grid_tile)
  }
  cells <- data.table::rbindlist(pieces, use.names = TRUE)
  if (nrow(cells)) {
    cells <- cells[, .(total = sum(total), shots = sum(shots)), by = .(col, row)]
    cells[, value := if (reducer == "mean") total / shots else total]
    cells[, `:=`(x = (col + 0.5) * res, y = (row + 0.5) * res)]
  }
  template <- .gedidb_grid_template(reference, res, output_crs)
  if (!nzchar(filename) && terra::ncell(template) > 1e8) {
    stop("This grid exceeds 100 million cells. Supply 'filename' for disk-backed output.", call. = FALSE)
  }
  if (!nrow(cells)) {
    result <- template
    terra::values(result) <- NA_real_
    if (nzchar(filename)) result <- terra::writeRaster(result, filename, overwrite = TRUE)
  } else {
    points <- terra::vect(as.data.frame(cells), geom = c("x", "y"), crs = output_crs)
    result <- terra::rasterize(
      points, template, field = "value", fun = "mean",
      filename = filename, overwrite = nzchar(filename)
    )
  }
  names(result) <- metric
  attr(result, "gedi_source") <- "gedidb"
  attr(result, "gedi_release") <- attr(reference, "version")
  attr(result, "gedidb_tiles") <- length(tiles)
  attr(result, "gedidb_workers") <- workers
  attr(result, "gedidb_shots") <- if (nrow(cells)) sum(cells$shots) else 0
  attr(result, "gedidb_elapsed") <- proc.time()[["elapsed"]] - started
  result
}
