# Unified Google Earth Engine backend ------------------------------------

#' Test Google Earth Engine availability
#'
#' Verifies that the optional Earth Engine Python API is installed,
#' initialized, and can access the GEDI catalogue.
#'
#' @return A single logical value with a diagnostic reason attribute.
#' @export
gee_available <- function() {
  backend <- getOption("rGEDI.gee_backend")
  if (!is.null(backend)) {
    result <- tryCatch(isTRUE(backend$available()), error = function(e) FALSE)
    attr(result, "reason") <- if (result) "available" else "mock backend unavailable"
    return(result)
  }
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    result <- FALSE
    attr(result, "reason") <- "the optional R package 'reticulate' is not installed"
    return(result)
  }
  if (!reticulate::py_module_available("ee")) {
    result <- FALSE
    attr(result, "reason") <- "the Earth Engine Python API is not installed"
    return(result)
  }
  check <- tryCatch({
    ee <- .require_ee()
    ee$data$getAsset("LARSE/GEDI/GEDI02_A_002_INDEX")
    TRUE
  }, error = function(e) e)
  if (inherits(check, "error")) {
    result <- FALSE
    attr(result, "reason") <- paste0(
      conditionMessage(check),
      "; initialize Earth Engine with ee_initialize(project = 'your-project')"
    )
    return(result)
  }
  result <- TRUE
  attr(result, "reason") <- "available"
  result
}

.gee_catalog_version <- function(product, version = NULL) {
  expected <- if (grepl("^GEDI04_A", product)) "002.1" else "002"
  if (is.null(version)) return(expected)
  normalized <- sub("^[Vv]", "", as.character(version)[1L])
  allowed <- if (grepl("^GEDI04_A", product)) {
    c("2", "02", "002", "2.1", "02.1", "002.1")
  } else {
    c("2", "02", "002")
  }
  if (!normalized %in% allowed) {
    stop(
      "Google Earth Engine currently provides ", product,
      " as GEDI Version ", expected, ", not Version ", normalized, ".",
      call. = FALSE
    )
  }
  expected
}

.new_gee_reference <- function(product, bbox, daterange = NULL, version = NULL,
                               checked = FALSE) {
  if (!product %in% names(.gedi_ee_catalog)) {
    stop(
      "Unsupported Earth Engine GEDI product: ", product, ". Choose one of: ",
      paste(names(.gedi_ee_catalog), collapse = ", "), ".",
      call. = FALSE
    )
  }
  address <- unname(.gedi_ee_catalog[[product]])
  kind <- if (product == "GEDI04_B") {
    "image"
  } else if (grepl("_MONTHLY$", product)) {
    "image_collection"
  } else {
    "features"
  }
  structure(
    address,
    class = c("gedi_gee_reference", "character"),
    product = product,
    base_product = sub("_MONTHLY$", "", product),
    bbox = bbox,
    daterange = daterange,
    version = .gee_catalog_version(product, version),
    kind = kind,
    server_checked = isTRUE(checked)
  )
}

#' @export
print.gedi_gee_reference <- function(x, ...) {
  cat("<GEDI Google Earth Engine reference>\n")
  cat("Product: ", attr(x, "product"), "\n", sep = "")
  cat("Version: ", attr(x, "version"), "\n", sep = "")
  cat("Asset: ", unclass(x)[1L], "\n", sep = "")
  cat("Type: ", attr(x, "kind"), "\n", sep = "")
  bbox <- attr(x, "bbox")
  if (!is.null(bbox)) cat("Bounding box: ", paste(unname(bbox), collapse = ", "), "\n", sep = "")
  dates <- attr(x, "daterange")
  if (!is.null(dates)) cat("Dates: ", paste(dates, collapse = " to "), "\n", sep = "")
  invisible(x)
}

.open_gee_reference <- function(x, product) {
  if (!inherits(x, "gedi_gee_reference")) {
    stop("Expected a Google Earth Engine GEDI reference.", call. = FALSE)
  }
  base_product <- attr(x, "base_product")
  if (!identical(base_product, product)) {
    stop("The Earth Engine reference contains ", base_product, ", not ", product, ".", call. = FALSE)
  }
  structure(
    list(reference = x, object = NULL, kind = attr(x, "kind")),
    class = "gedi_gee_object"
  )
}

.gee_reference_aoi <- function(reference) {
  bbox <- attr(reference, "bbox")
  if (is.null(bbox)) NULL else unname(bbox[c("xmin", "xmax", "ymin", "ymax")])
}

.gee_resolve <- function(x, quality = TRUE, max_granules = 200L) {
  if (!is.null(x$object)) return(x$object)
  backend <- getOption("rGEDI.gee_backend")
  if (!is.null(backend)) return(backend$open(x$reference, quality = quality))
  if (!isTRUE(attr(x$reference, "server_checked"))) {
    available <- gee_available()
    if (!isTRUE(available)) {
      stop("Google Earth Engine is unavailable: ", attr(available, "reason"), call. = FALSE)
    }
  }
  dates <- attr(x$reference, "daterange")
  gediEE(
    product = attr(x$reference, "product"),
    start_date = if (is.null(dates)) NULL else dates[1L],
    end_date = if (is.null(dates)) NULL else dates[2L],
    aoi = .gee_reference_aoi(x$reference),
    quality = quality,
    max_granules = max_granules
  )
}

.new_gee_object <- function(reference, object, kind = "features") {
  structure(
    list(reference = reference, object = object, kind = kind),
    class = "gedi_gee_object"
  )
}

.gee_feature_variables <- function(product, variables) {
  variables <- unique(as.character(variables))
  if (product == "GEDI04_A") {
    variables[variables == "l4a_quality_flag_rel3"] <- "l4_quality_flag"
    variables <- setdiff(variables, c(
      "l2a_quality_flag_rel3", "degrade_include_flag",
      "elev_highestreturn_outlier_flag", "delta_time"
    ))
  }
  if (product == "GEDI02_A") {
    variables[variables == "l2a_quality_flag_rel3"] <- "quality_flag"
  }
  if (product == "GEDI02_B") {
    variables[variables == "l2b_quality_flag_rel3"] <- "quality_flag"
    variables[variables == "l2_algrunflag"] <- "algorithmrun_flag"
  }
  setdiff(variables, c("lat_lowestmode", "lon_lowestmode", "time"))
}

.gee_collect <- function(x, variables, quality = TRUE, max_features = 100000L,
                         chunk_size = 5000L, max_granules = 200L) {
  if (!identical(x$kind, "features")) {
    stop("Only footprint FeatureCollections can be collected as a data.table. Use collect = FALSE for raster products.", call. = FALSE)
  }
  product <- attr(x$reference, "base_product")
  variables <- .gee_feature_variables(product, variables)
  backend <- getOption("rGEDI.gee_backend")
  if (!is.null(backend)) {
    ans <- backend$collect(
      x$reference, variables = variables, quality = quality,
      max_features = max_features, chunk_size = chunk_size
    )
  } else {
    ee <- .require_ee()
    collection <- .gee_resolve(x, quality = quality, max_granules = max_granules)
    collection <- collection$map(function(feature) {
      coordinates <- feature$geometry()$coordinates()
      feature$set(reticulate::dict(
        lon_lowestmode = coordinates$get(0L),
        lat_lowestmode = coordinates$get(1L),
        shot_number = ee$String(feature$get("shot_number"))
      ))
    })
    keep <- unique(c(
      "shot_number", "lon_lowestmode", "lat_lowestmode", variables
    ))
    collection <- collection$select(as.list(keep))
    n <- as.integer(reticulate::py_to_r(collection$size()$getInfo()))
    max_features <- as.integer(max_features)[1L]
    if (!is.finite(max_features) || max_features < 1L) {
      stop("'max_features' must be a positive integer.", call. = FALSE)
    }
    if (n > max_features) {
      stop(
        "The Earth Engine query contains ", n, " footprints, above max_features=",
        max_features, ". Use collect = FALSE for server-side processing, ",
        "or reduce the AOI/date range.", call. = FALSE
      )
    }
    if (!n) {
      ans <- data.table::data.table()
    } else {
      chunk_size <- max(1L, as.integer(chunk_size)[1L])
      offsets <- seq.int(0L, n - 1L, by = chunk_size)
      pieces <- lapply(offsets, function(offset) {
        count <- min(chunk_size, n - offset)
        subset <- ee$FeatureCollection(collection$toList(count, offset))
        .ee_to_dt(subset)
      })
      ans <- data.table::rbindlist(pieces, use.names = TRUE, fill = TRUE)
    }
  }
  ans <- data.table::as.data.table(ans)
  attr(ans, "gedi_source") <- "gee"
  attr(ans, "gedi_release") <- attr(x$reference, "version")
  attr(ans, "gee_asset") <- unclass(x$reference)[1L]
  ans[]
}

.gee_select <- function(x, variables, quality = TRUE, max_granules = 200L) {
  backend <- getOption("rGEDI.gee_backend")
  if (!is.null(backend) && !is.null(backend$select)) {
    object <- backend$select(x$reference, variables, quality = quality)
  } else {
    object <- .gee_resolve(x, quality = quality, max_granules = max_granules)
    object <- object$select(as.list(unique(variables)))
  }
  .new_gee_object(x$reference, object, x$kind)
}

.gee_clip <- function(x, clip_obj) {
  backend <- getOption("rGEDI.gee_backend")
  if (!is.null(backend)) {
    return(.new_gee_object(
      x$reference, backend$clip(x, clip_obj), x$kind
    ))
  }
  object <- .gee_resolve(x)
  geometry <- .as_ee_geom(clip_obj)
  clipped <- switch(
    x$kind,
    features = object$filterBounds(geometry),
    image_collection = object$map(function(image) image$clip(geometry)),
    image = object$clip(geometry)
  )
  .new_gee_object(x$reference, clipped, x$kind)
}

.gee_reducer <- function(fun) {
  ee <- .require_ee()
  name <- .gee_reducer_name(fun)
  switch(
    name,
    mean = ee$Reducer$mean(), sum = ee$Reducer$sum(),
    min = ee$Reducer$min(), max = ee$Reducer$max(),
    median = ee$Reducer$median(),
    stop("Unsupported GEE reducer: ", name, call. = FALSE)
  )
}

.gee_reducer_name <- function(fun) {
  if (is.character(fun)) name <- tolower(fun[1L]) else if (identical(fun, mean)) {
    name <- "mean"
  } else if (identical(fun, sum)) {
    name <- "sum"
  } else if (identical(fun, min)) {
    name <- "min"
  } else if (identical(fun, max)) {
    name <- "max"
  } else if (identical(fun, stats::median)) {
    name <- "median"
  } else {
    stop("GEE rasterization supports mean, sum, min, max, and median.", call. = FALSE)
  }
  if (!name %in% c("mean", "sum", "min", "max", "median")) {
    stop("Unsupported GEE reducer: ", name, call. = FALSE)
  }
  name
}

.gee_rasterize <- function(x, metric, res = 1000, fun = mean,
                           crs = "EPSG:3857") {
  backend <- getOption("rGEDI.gee_backend")
  if (!is.null(backend)) {
    object <- backend$rasterize(x, metric = metric, res = res, fun = fun, crs = crs)
    return(.new_gee_object(x$reference, object, "image"))
  }
  object <- .gee_resolve(x)
  if (x$kind == "features") {
    image <- object$reduceToImage(
      properties = as.list(metric), reducer = .gee_reducer(fun)
    )$rename(metric)
  } else if (x$kind == "image_collection") {
    reducer <- .gee_reducer_name(fun)
    image <- switch(
      reducer,
      mean = object$select(metric)$mean(),
      sum = object$select(metric)$sum(),
      min = object$select(metric)$min(),
      max = object$select(metric)$max(),
      median = object$select(metric)$median(),
      stop("Unsupported GEE ImageCollection reducer.", call. = FALSE)
    )
  } else {
    image <- object$select(metric)
  }
  image <- image$reproject(crs = crs, scale = as.numeric(res))
  .new_gee_object(x$reference, image, "image")
}

.unwrap_gee <- function(x) {
  if (inherits(x, "gedi_gee_object")) .gee_resolve(x) else x
}

methods::setOldClass("gedi_gee_object")

#' @rdname clip
#' @export
setMethod(
  "clip", c(x = "gedi_gee_object", clip_obj = "ANY"),
  function(x, clip_obj, ...) .gee_clip(x, clip_obj)
)
