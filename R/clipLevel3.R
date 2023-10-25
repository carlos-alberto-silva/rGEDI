#' Clip GEDI Full Waveform Geolocations by Coordinates
#'
#' @description This function clips GEDI level3 extracted geolocation ([getLevel3()])
#' data a within given bounding coordinates
#'
#' @usage clipLevel3(level3, xmin, xmax, ymin, ymax)
#'
#' @param level3 A [`rGEDI::gedi.level3-class`] resulting from [getLevel3()].
#' @param xmin Numeric. West longitude (x) coordinate of the bounding rectangle, in decimal degrees.
#' @param xmax Numeric. East longitude (x) coordinate of the bounding rectangle, in decimal degrees.
#' @param ymin Numeric. South latitude (y) coordinate of the bounding rectangle, in decimal degrees.
#' @param ymax Numeric. North latitude (y) coordinate of the bounding rectangle, in decimal degrees.
#'
#' @return Returns the clipped S4 object of class [`rGEDI::gedi.level3-class`].
#'
#' @seealso \url{https://lpdaac.usgs.gov/products/gedi01_bv002/}
#'
#' @examples
#' # Specifying the path to GEDI level3 data (zip file)
#' outdir <- tempdir()
#' level3_fp_zip <- system.file("extdata",
#'   "GEDI01_B_2019108080338_O01964_T05337_02_003_01_sub.zip",
#'   package = "rGEDI"
#' )
#'
#' # Unzipping GEDI level3 data
#' level3path <- unzip(level3_fp_zip, exdir = outdir)
#'
#' # Reading GEDI level3 data (h5 file)
#' level3 <- readLevel3(level3path = level3path)
#'
#' # Bounding rectangle coordinates
#' xmin <- -44.15036
#' xmax <- -44.10066
#' ymin <- -13.75831
#' ymax <- -13.71244
#'
#' # Clipping GEDI Full Waveform Geolocations by boundary box extent
#' level3_clip <- clipLevel3(level3, xmin, xmax, ymin, ymax)
#'
#' hasLeaflet <- require(leaflet)
#'
#' if (hasLeaflet) {
#'   leaflet() %>%
#'     addCircleMarkers(level3_clip$longitude_bin0,
#'       level3_clip$latitude_bin0,
#'       radius = 1,
#'       opacity = 1,
#'       color = "red"
#'     ) %>%
#'     addScaleBar(options = list(imperial = FALSE)) %>%
#'     addProviderTiles(providers$Esri.WorldImagery)
#' }
#'
#'
#' @import terra
#' @export
clipLevel3 <- function(level3, xmin, xmax, ymin, ymax) {
  # Create a SpatExtent object from the extents
  ext <- terra::ext(c(xmin, xmax, ymin, ymax))

  # Crop the SpatRast object by the specified extents
  clipped_rast <- terra::crop(level3@raster, ext)

  level3_out <- new("gedi.level3", list(raster = clipped_rast))
  return(level3_out)
}

setGeneric("clip", function(level3, xa, xb, xc, xd, ...) standardGeneric("clip"))

setMethod(
  "clip",
  signature("gedi.level3", "numeric", "numeric", "numeric", "numeric"),
  function(level3, xa, xb, xc, xd, ...) {
    message(ymin)
    message(ymax)
    clipLevel3(level3, xa, xb, xc, xd)
})

setMethod(
  "clip",
  signature("gedi.level3", "bbox", "missing", "missing", "missing"),
  function(level3, xa, ...) {
    clipLevel3(level3, y$xmin, y$xmax, y$ymin, y$ymax)
})

#' Clip GEDI Full Waveform Geolocations by geometry
#'
#' @description This function clips level3 extracted geolocation (level3)
#' data within a given geometry
#'
#' @param level3 A [`data.table::data.table-class`] resulting from [getLevel3()] function.
#' @param polygon_sf Polygon. An object of class [`sf`],
#' which can be loaded as an ESRI shapefile using [raster::shapefile] function in the \emph{raster} package.
#' @param split_by Polygon id. If defined, GEDI data will be clipped by each polygon using the polygon id from
#' table of attribute defined by the user.
#'
#' @return Returns a lilst of S4 object of class [`data.table::data.table-class`] containing the
#' clipped GEDI level3 extracted geolocations.
#'
#' @seealso \url{https://lpdaac.usgs.gov/products/gedi01_bv002/}
#'
#' @examples
#' # Specifying the path to GEDI level3 data (zip file)
#' outdir <- tempdir()
#' level3_fp_zip <- system.file("extdata",
#'   "GEDI01_B_2019108080338_O01964_T05337_02_003_01_sub.zip",
#'   package = "rGEDI"
#' )
#'
#' # Unzipping GEDI level3 data
#' level3path <- unzip(level3_fp_zip, exdir = outdir)
#'
#' # Reading GEDI level3 data (h5 file)
#' level3 <- readLevel3(level3path = level3path)
#'
#' # Extracting GEDI Full Waveform Geolocations
#' level3 <- getLevel3(level3)
#'
#' # Specifying the path to shapefile
#' polygon_filepath <- system.file("extdata", "stands_cerrado.shp", package = "rGEDI")
#'
#' # Reading shapefile as SpatialPolygonsDataFrame object
#' library(sf)
#' polygon_spdf <- sf::st_read(polygon_filepath)
#'
#' # Clipping GEDI Full Waveform Geolocations by Geometry
#' level3_clip <- clipLevel3Geometry(level3, polygon_spdf, split_by = "id")
#'
#' hasLeaflet <- require(leaflet)
#'
#' if (hasLeaflet) {
#'   leaflet() %>%
#'     addCircleMarkers(level3_clip$longitude_bin0,
#'       level3_clip$latitude_bin0,
#'       radius = 1,
#'       opacity = 1,
#'       color = "red"
#'     ) %>%
#'     addScaleBar(options = list(imperial = FALSE)) %>%
#'     addPolygons(
#'       data = polygon_spdf, weight = 1, col = "white",
#'       opacity = 1, fillOpacity = 0
#'     ) %>%
#'     addProviderTiles(providers$Esri.WorldImagery)
#' }
#'
#' @import terra
#' @export
clipLevel3Geometry <- function(level3, polygon_sf, split_by = "id") {
  # Split the SpatialPolygonsDataFrame into a list of polygons by ID
  polygons_list <- split(polygon_sf, polygon_sf[[split_by]])

  # Initialize an empty list to store cropped raster objects
  cropped_rasts <- list()

  # Iterate through each polygon, crop the raster, and store the results
  for (polygon_id in names(polygons_list)) {
    polygon <- polygons_list[[polygon_id]]
    extent <- terra::ext(polygon) # Calculate extent for the current polygon
    cropped <- terra::crop(level3@raster, extent) # Crop the raster
    cropped_rasts[[polygon_id]] <- cropped # Store the cropped raster in the list
  }
  
  return(cropped_rasts)
}
