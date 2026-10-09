# Google Earth Engine ----------------------------------------------------

.require_ee <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) stop("Install 'reticulate' for Earth Engine support.")
  if (!reticulate::py_module_available("ee")) stop("Install the Earth Engine Python API with rGEDI_configure(install = TRUE).")
  reticulate::import("ee", delay_load = FALSE, convert = FALSE)
}

#' Initialize Google Earth Engine
#' @param project Google Cloud project ID.
#' @param authenticate Authenticate when initialization fails.
#' @param quiet Logical.
#' @return Invisibly returns the Python Earth Engine module.
#' @export
ee_initialize <- function(project = Sys.getenv("EE_PROJECT", unset = ""), authenticate = TRUE, quiet = FALSE) {
  ee <- .require_ee()
  init <- function() ee$Initialize(project = if (nzchar(project)) project else NULL)
  tryCatch(init(), error = function(e) {
    if (!isTRUE(authenticate)) stop(e)
    ee$Authenticate(); init()
  })
  if (!quiet) message("Google Earth Engine initialized", if (nzchar(project)) paste0(" for ", project) else "", ".")
  invisible(ee)
}

.as_ee_geom <- function(x) {
  ee <- .require_ee()
  if (inherits(x, "python.builtin.object")) return(x)
  if (is.numeric(x) && length(x) == 4L) {
    x <- as.numeric(x)
    return(ee$Geometry$Rectangle(as.list(x[c(1L, 3L, 2L, 4L)])))
  }
  obj <- if (inherits(x, "SpatVector")) sf::st_as_sf(x) else sf::st_as_sf(x)
  path <- tempfile(fileext = ".geojson")
  on.exit(unlink(path), add = TRUE)
  suppressMessages(sf::st_write(sf::st_transform(obj, 4326), path, quiet = TRUE))
  geojson <- jsonlite::read_json(path, simplifyVector = FALSE)
  ee$FeatureCollection(reticulate::r_to_py(geojson))$geometry()
}

#' Convert a vector object to an Earth Engine feature collection
#' @param x `sf`, `sfc`, `SpatVector`, or vector file path.
#' @return An Earth Engine `FeatureCollection`.
#' @export
vect_as_ee <- function(x) {
  ee <- .require_ee()
  if (is.character(x)) x <- terra::vect(x)
  obj <- if (inherits(x, "SpatVector")) sf::st_as_sf(x) else sf::st_as_sf(x)
  path <- tempfile(fileext = ".geojson"); on.exit(unlink(path), add = TRUE)
  suppressMessages(sf::st_write(sf::st_transform(obj, 4326), path, quiet = TRUE))
  ee$FeatureCollection(reticulate::r_to_py(jsonlite::read_json(path, simplifyVector = FALSE)))
}

#' Convert an extent to an Earth Engine geometry
#' @param x A `SpatExtent`, raster/vector object, or numeric extent in
#'   `c(xmin, xmax, ymin, ymax)` order.
#' @return An Earth Engine geometry.
#' @export
ext_to_ee <- function(x) {
  if (!is.numeric(x)) x <- as.vector(terra::ext(x))
  .as_ee_geom(x)
}

#' Convert an Earth Engine rectangle to `sf`
#' @param aoi Earth Engine geometry.
#' @return An `sf` polygon.
#' @export
ee_rect_to_sf <- function(aoi) {
  info <- reticulate::py_to_r(aoi$getInfo())
  coordinates <- info$coordinates
  if (identical(info$type, "Polygon")) {
    coordinates <- lapply(coordinates, function(ring) {
      matrix(unlist(ring), ncol = 2L, byrow = TRUE)
    })
    geometry <- sf::st_polygon(coordinates)
  } else {
    stop("`aoi` must be an Earth Engine Polygon geometry.", call. = FALSE)
  }
  sf::st_as_sf(sf::st_sfc(geometry, crs = 4326))
}

.gedi_ee_catalog <- c(
  GEDI02_A = "LARSE/GEDI/GEDI02_A_002",
  GEDI02_A_MONTHLY = "LARSE/GEDI/GEDI02_A_002_MONTHLY",
  GEDI02_B = "LARSE/GEDI/GEDI02_B_002",
  GEDI02_B_MONTHLY = "LARSE/GEDI/GEDI02_B_002_MONTHLY",
  GEDI04_A = "LARSE/GEDI/GEDI04_A_002",
  GEDI04_A_MONTHLY = "LARSE/GEDI/GEDI04_A_002_MONTHLY",
  GEDI04_B = "LARSE/GEDI/GEDI04_B_002"
)

#' Search the GEDI datasets available in Earth Engine
#'
#' @param ... Words matched against the product name, catalog ID, and title.
#' @param operator Use `"and"` or `"or"` when several words are supplied.
#' @return A data table containing product, catalog ID, version, and title.
#' @export
search_datasets <- function(..., operator = c("and", "or")) {
  operator <- match.arg(tolower(operator), c("and", "or"))
  query <- unlist(list(...), use.names = FALSE)
  catalog <- data.table::data.table(
    product = names(.gedi_ee_catalog),
    id = unname(.gedi_ee_catalog),
    version = c("2", "2", "2", "2", "2.1", "2.1", "2"),
    title = c(
      "GEDI Level 2A footprint height and elevation",
      "GEDI Level 2A monthly raster",
      "GEDI Level 2B footprint canopy metrics",
      "GEDI Level 2B monthly raster",
      "GEDI Level 4A footprint biomass density",
      "GEDI Level 4A monthly raster",
      "GEDI Level 4B gridded biomass density"
    )
  )
  if (!length(query)) return(catalog[])
  text <- tolower(paste(catalog$product, catalog$id, catalog$title))
  matches <- vapply(query, function(term) grepl(tolower(term), text, fixed = TRUE),
                    logical(nrow(catalog)))
  if (is.null(dim(matches))) matches <- matrix(matches, ncol = 1L)
  keep <- if (operator == "and") rowSums(matches) == ncol(matches) else rowSums(matches) > 0L
  catalog[keep]
}

#' Return an Earth Engine catalog ID
#' @param x A GEDI product name, catalog search row, or catalog ID.
#' @return A character catalog ID.
#' @export
get_catalog_id <- function(x) {
  if (inherits(x, "data.frame")) {
    if (!nrow(x) || !"id" %in% names(x)) stop("The search result contains no catalog ID.")
    return(as.character(x$id[[1L]]))
  }
  x <- as.character(x)[1L]
  if (x %in% names(.gedi_ee_catalog)) unname(.gedi_ee_catalog[[x]]) else x
}

#' Open a GEDI product in Google Earth Engine
#' @param product GEDI Earth Engine product name.
#' @param start_date,end_date Optional dates for image collections.
#' @param aoi Optional spatial filter.
#' @param quality Apply available quality masks.
#' @param max_granules Maximum number of vector granules to merge. Vector GEDI
#'   products are stored as folders of per-granule feature collections in Earth
#'   Engine, so use dates and/or an AOI to keep the request focused.
#' @return An Earth Engine image, image collection, or feature collection.
#' @export
gediEE <- function(product = names(.gedi_ee_catalog), start_date = NULL,
                   end_date = NULL, aoi = NULL, quality = TRUE,
                   max_granules = 200L) {
  product <- match.arg(product)
  ee <- .require_ee(); id <- unname(.gedi_ee_catalog[[product]])
  if (product == "GEDI04_B") return(ee$Image(id))

  is_monthly <- grepl("_MONTHLY$", product)
  if (!is_monthly) {
    index <- ee$FeatureCollection(paste0(id, "_INDEX"))
    if (!is.null(start_date)) {
      index <- index$filter(ee$Filter$gte("time_end", as.character(start_date)))
    }
    if (!is.null(end_date)) {
      index <- index$filter(ee$Filter$lt("time_start", as.character(end_date)))
    }
    geom <- NULL
    if (!is.null(aoi)) {
      geom <- .as_ee_geom(aoi)
      index <- index$filterBounds(geom)
    }
    n <- reticulate::py_to_r(index$size()$getInfo())
    if (!n) return(ee$FeatureCollection(list()))
    max_granules <- as.integer(max_granules)[1L]
    if (is.na(max_granules) || max_granules < 1L) {
      stop("`max_granules` must be a positive integer.", call. = FALSE)
    }
    if (n > max_granules) {
      stop(sprintf(
        paste0("The query matches %s GEDI granules, above max_granules=%s. ",
               "Use a smaller AOI/date range or increase max_granules."),
        n, max_granules
      ), call. = FALSE)
    }
    ids <- reticulate::py_to_r(index$aggregate_array("table_id")$getInfo())
    collections <- lapply(ids, ee$FeatureCollection)
    collection <- Reduce(function(x, y) x$merge(y), collections)
    if (!is.null(geom)) collection <- collection$filterBounds(geom)
    if (isTRUE(quality)) {
      flag <- if (product == "GEDI04_A") "l4_quality_flag" else "quality_flag"
      collection <- collection$filter(ee$Filter$eq(flag, 1L))$
        filter(ee$Filter$eq("degrade_flag", 0L))
    }
    return(collection)
  }

  collection <- ee$ImageCollection(id)
  if (!is.null(start_date) && !is.null(end_date)) collection <- collection$filterDate(start_date, end_date)
  if (!is.null(aoi)) collection <- collection$filterBounds(.as_ee_geom(aoi))
  if (isTRUE(quality)) {
    flag <- if (grepl("04_A", product)) "l4_quality_flag" else "quality_flag"
    collection <- collection$map(function(image) image$updateMask(image$select(flag)$eq(1L)))
  }
  collection
}

.ee_to_dt <- function(collection) {
  size <- reticulate::py_to_r(collection$size()$getInfo())
  if (!size) return(data.table::data.table())
  features <- reticulate::py_to_r(collection$toList(size)$getInfo())
  data.table::rbindlist(lapply(features, `[[`, "properties"), fill = TRUE)
}

#' Extract Earth Engine image values at GEDI footprints
#' @param stack Earth Engine image or list of images.
#' @param geom GEDI table, `sf`, or `SpatVector` points.
#' @param scale Sampling scale in meters.
#' @param chunk_size Number of points per request.
#' @return A [data.table::data.table].
#' @export
extractEE <- function(stack, geom, scale = 30, chunk_size = 1000L) {
  ee <- .require_ee()
  if (inherits(geom, c("data.frame", "data.table")) &&
      !inherits(geom, c("sf", "sfc", "SpatVector"))) {
    geom <- to_vect(geom)
  }
  fc <- vect_as_ee(geom)
  image <- if (is.list(stack)) Reduce(function(a, b) a$addBands(b), stack) else stack
  n <- reticulate::py_to_r(fc$size()$getInfo())
  if (!n) return(data.table::data.table())
  out <- list(); k <- 0L
  for (first in seq.int(0L, max(0L, n - 1L), by = chunk_size)) {
    subset <- ee$FeatureCollection(fc$toList(as.integer(min(chunk_size, n - first)), as.integer(first)))
    sampled <- image$sampleRegions(collection = subset, scale = as.integer(scale), tileScale = 16L)
    k <- k + 1L; out[[k]] <- .ee_to_dt(sampled)
  }
  data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
}

#' @rdname extractEE
#' @export
extractGEDIAncillary <- extractEE

#' Build an HLS and terrain predictor stack in Earth Engine
#' @param x Area of interest.
#' @param start_date,end_date Date range.
#' @param cloud_max Maximum HLS cloud cover.
#' @param buffer_m AOI buffer in meters.
#' @return An Earth Engine image.
#' @export
ee_build_hls_s1c_terrain_stack <- function(x, start_date, end_date, cloud_max = 20, buffer_m = 30) {
  ee <- .require_ee(); aoi <- .as_ee_geom(x)$buffer(buffer_m)
  hls <- ee$ImageCollection("NASA/HLS/HLSS30/v002")$filterBounds(aoi)$filterDate(start_date, end_date)$filter(sprintf("CLOUD_COVERAGE < %s", cloud_max))
  image <- hls$median()$select(
    c("B2", "B3", "B4", "B5", "B6", "B7"),
    c("blue", "green", "red", "nir", "swir1", "swir2")
  )$clip(aoi)
  indices <- image$normalizedDifference(c("nir", "red"))$rename("ndvi")$
    addBands(image$normalizedDifference(c("green", "nir"))$rename("ndwi"))
  s1 <- ee$ImageCollection("COPERNICUS/S1_GRD")$filterBounds(aoi)$
    filterDate(start_date, end_date)$
    filter(ee$Filter$listContains("transmitterReceiverPolarisation", "VV"))$
    filter(ee$Filter$listContains("transmitterReceiverPolarisation", "VH"))$
    filter(ee$Filter$eq("instrumentMode", "IW"))$median()$
    select(c("VV", "VH"), c("vv", "vh"))$clip(aoi)
  dem <- ee$Image("NASA/NASADEM_HGT/001")$select("elevation")$clip(aoi)
  terrain <- ee$Terrain$products(dem)
  image$addBands(indices)$addBands(s1)$
    addBands(terrain$select(c("elevation", "slope", "aspect")))
}

#' Build an AlphaEarth and terrain predictor stack
#' @param geom Area of interest.
#' @param start_year,end_year Inclusive year range.
#' @param add_lonlat Add coordinate bands.
#' @return An Earth Engine image.
#' @export
ee_build_AlphaEarth_embedding_terrain_stack <- function(geom, start_year, end_year, add_lonlat = TRUE) {
  ee <- .require_ee(); aoi <- .as_ee_geom(geom)
  image <- ee$ImageCollection("GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL")$filterDate(sprintf("%04d-01-01", start_year), sprintf("%04d-12-31", end_year))$filterBounds(aoi)$median()$clip(aoi)
  dem <- ee$Image("NASA/NASADEM_HGT/001")$select("elevation")$clip(aoi)
  image <- image$addBands(ee$Terrain$products(dem)$select(c("elevation", "slope", "aspect")))
  if (isTRUE(add_lonlat)) image <- image$addBands(ee$Image$pixelLonLat()$clip(aoi))
  image
}

#' Get an Earth Engine image tile URL
#' @param image Earth Engine image.
#' @param vis Visualization parameter list.
#' @return Character tile URL.
#' @export
getTileUrl <- function(image, vis = list()) reticulate::py_to_r(image$getMapId(reticulate::r_to_py(vis))[["tile_fetcher"]]$url_format)

#' Add an Earth Engine image to a Leaflet map
#' @param map A Leaflet widget.
#' @param image Earth Engine image.
#' @param bands One or three bands.
#' @param min,max Visualization limits.
#' @param palette Optional palette.
#' @param group Layer group.
#' @param ... Passed to [leaflet::addTiles()].
#' @return A Leaflet widget.
#' @export
addEEImage <- function(map, image, bands = NULL, min = 0, max = 1, palette = NULL, group = NULL, ...) {
  if (!requireNamespace("leaflet", quietly = TRUE)) stop("Install 'leaflet'.")
  vis <- list(min = min, max = max)
  if (!is.null(bands)) vis$bands <- bands
  if (!is.null(palette)) vis$palette <- palette
  leaflet::addTiles(map, urlTemplate = getTileUrl(image, vis), group = group, ...)
}

#' Download an Earth Engine image
#' @param image Earth Engine image.
#' @param filename Destination GeoTIFF path.
#' @param region Area of interest.
#' @param scale Pixel size in meters.
#' @param crs Output CRS.
#' @param overwrite Logical.
#' @return Invisibly returns `filename`.
#' @export
map_download <- function(image, filename, region, scale = 30, crs = "EPSG:4326", overwrite = FALSE) {
  if (file.exists(filename) && !overwrite) stop("Output exists; use overwrite=TRUE.")
  params <- list(scale = scale, crs = crs, region = .as_ee_geom(region), format = "GEO_TIFF")
  url <- reticulate::py_to_r(image$getDownloadURL(reticulate::r_to_py(params)))
  curl::curl_download(url, filename, quiet = FALSE)
  invisible(filename)
}

#' Check an Earth Engine task
#' @param task Earth Engine batch task.
#' @param quiet Logical.
#' @return Task status list.
#' @export
ee_check_task_status <- function(task, quiet = TRUE) {
  status <- reticulate::py_to_r(task$status())
  if (!quiet) print(status)
  status
}

#' Train a regression forest in Earth Engine
#'
#' @param training Earth Engine `FeatureCollection` containing predictors and
#'   the response.
#' @param response Response property name, such as `"agbd"` or `"rh98"`.
#' @param predictors Predictor property names.
#' @param trees Number of trees.
#' @param seed Random seed.
#' @return A trained Earth Engine classifier configured for regression.
#' @export
build_ee_forest <- function(training, response, predictors, trees = 500L,
                            seed = 1L) {
  ee <- .require_ee()
  forest <- ee$Classifier$smileRandomForest(
    numberOfTrees = as.integer(trees), seed = as.integer(seed)
  )$setOutputMode("REGRESSION")
  forest$train(
    features = training,
    classProperty = response,
    inputProperties = reticulate::r_to_py(as.list(predictors))
  )
}

#' Apply an Earth Engine model to a predictor image
#' @param model A trained Earth Engine classifier or regressor.
#' @param stack Earth Engine predictor image.
#' @param aoi Optional area to clip.
#' @param name Output band name.
#' @return An Earth Engine image.
#' @export
map_create <- function(model, stack, aoi = NULL, name = "prediction") {
  out <- stack$classify(model, name)$toFloat()
  if (!is.null(aoi)) out <- out$clip(.as_ee_geom(aoi))
  out
}

#' Display Earth Engine layers in Leaflet
#' @param layers Named list of Earth Engine images.
#' @param vis Named list of visualization parameter lists.
#' @param aoi Optional spatial object used to fit map bounds.
#' @return A Leaflet widget.
#' @export
map_view <- function(layers, vis = list(), aoi = NULL) {
  if (!requireNamespace("leaflet", quietly = TRUE)) stop("Install 'leaflet'.")
  if (!is.list(layers)) layers <- list(GEDI = layers)
  if (is.null(names(layers))) names(layers) <- paste0("layer_", seq_along(layers))
  map <- leaflet::addProviderTiles(leaflet::leaflet(), "CartoDB.Positron")
  for (i in seq_along(layers)) {
    settings <- vis[[names(layers)[i]]] %||gee% vis[[i]] %||gee% list()
    map <- leaflet::addTiles(map, getTileUrl(layers[[i]], settings),
                             group = names(layers)[i])
  }
  if (!is.null(aoi)) {
    box <- sf::st_bbox(sf::st_transform(sf::st_as_sf(aoi), 4326))
    map <- leaflet::fitBounds(map, box[["xmin"]], box[["ymin"]],
                              box[["xmax"]], box[["ymax"]])
  }
  leaflet::addLayersControl(map, overlayGroups = names(layers))
}

`%||gee%` <- function(x, y) if (is.null(x)) y else x

#' Earth Engine terrain slope
#' @param x Earth Engine elevation image.
#' @return An Earth Engine image.
#' @export
slope <- function(x) .require_ee()$Terrain$slope(x)

#' Earth Engine terrain aspect
#' @param x Earth Engine elevation image.
#' @return An Earth Engine image.
#' @export
aspect <- function(x) .require_ee()$Terrain$aspect(x)

#' Earth Engine gray-level co-occurrence textures
#' @param x Earth Engine integer image.
#' @param size Neighborhood size.
#' @param kernel Optional Earth Engine kernel.
#' @param average Average directional metrics.
#' @return An Earth Engine image.
#' @export
glcmTexture <- function(x, size = 1L, kernel = NULL, average = TRUE) {
  x$glcmTexture(size = as.integer(size), kernel = kernel, average = average)
}
