#' Clip GEDI Full Waveform Geolocations by Coordinates
#'
#' @description This function clips GEDI level3 extracted geolocation ([readLevel3()])
#' data a within given bounding coordinates
#'
#' @usage clipLevel3(level3, xmin, xmax, ymin, ymax)
#'
#' @param level3 A [`rGEDI::gedi.level3-class`] resulting from [readLevel3()].
#' @param xmin Numeric. West longitude (x) coordinate of the bounding rectangle, in decimal degrees.
#' @param xmax Numeric. East longitude (x) coordinate of the bounding rectangle, in decimal degrees.
#' @param ymin Numeric. South latitude (y) coordinate of the bounding rectangle, in decimal degrees.
#' @param ymax Numeric. North latitude (y) coordinate of the bounding rectangle, in decimal degrees.
#'
#' @return Returns the clipped S4 object of class [`rGEDI::gedi.level3-class`].
#'
#' @seealso \url{https://daac.ornl.gov/GEDI/guides/GEDI_L3_LandSurface_Metrics_V2.html}
#'
#' @examples
#' # Specifying the path to GEDI level3 data (zip file)
#' outdir <- tempdir()
#'
#' level3path <- "C:/Users/caiohamamura/Downloads/GEDI03_rh100_mean_2019108_2020287_002_02.tif"
#' level3 <- readLevel3(level3path = level3path)
#' xmin <- -5888627
#' xmax <- -5878615
#' ymin <- -990000
#' ymax <- -976400
#'
#' # Clipping GEDI Full Waveform Geolocations by boundary box extent
#' level3_crop <- crop(level3, xmin, xmax, ymin, ymax)
#'
#' if (require(ggplot2)) {
#'   ggplot() +
#'     geom_spatraster(data = level3_crop)
#' }
#'
#' @export
clipLevel3 <- function(level3, xmin, xmax, ymin, ymax) {
  rast <- level3@raster
  # Create a SpatExtent object from the extents
  ext <- terra::ext(c(xmin, xmax, ymin, ymax))

  # Crop the SpatRast object by the specified extents
  clipped_rast <- terra::crop(rast, ext)

  level3_out <- new("gedi.level3", raster = clipped_rast)
  return(level3_out)
}

#' Plot function for gedi.level3 data
#'
#' @description Plots gedi.level3 raster data
#'
#' @param x A [`gedi.level3-class`] object
#' @param y Not used, inherited by generic method
#' @param ... Additional parameters for plotting
#'
#' @examples
#' level3path <- "C:/Users/caiohamamura/Downloads/GEDI03_rh100_mean_2019108_2020287_002_02.tif"
#' level3 <- readLevel3(level3path = level3path)
#' xmin <- -5888627
#' xmax <- -5878615
#' ymin <- -990000
#' ymax <- -976400
#'
#' # Clipping GEDI Full Waveform Geolocations by boundary box extent
#' level3_crop <- crop(level3, xmin, xmax, ymin, ymax)
#'
#' plot(level3_crop)
#'
#' @return Nothing
#' @export
setMethod(
  f = "plot",
  signature("gedi.level3", y = "missing"),
  definition = function(x, ...) {
    plot(x@raster, ...)
  }
)

#' Generic function for geom_spatraster
#' @param mapping Set of aesthetic mappings created by [`ggplot2::aes()`] or
#' [`ggplot2::aes()`]. See Aesthetics specially in the use of fill aesthetic.
#' @param data A `SpatRaster` object
#' @param ... Other arguments passed on to [`ggplot2::layer()`]. These are often
#' aesthetics, used to set an aesthetic to a fixed value, like colour = "red" or
#' size = 3. They may also be parameters to the paired geom/stat.
#'
#' @return A ggplot2 layer
#'
#' @examples
#' level3path <- "C:/Users/caiohamamura/Downloads/GEDI03_rh100_mean_2019108_2020287_002_02.tif"
#' level3 <- readLevel3(level3path = level3path)
#' xmin <- -5888627
#' xmax <- -5878615
#' ymin <- -990000
#' ymax <- -976400
#'
#' # Clipping GEDI Full Waveform Geolocations by boundary box extent
#' level3_crop <- crop(level3, xmin, xmax, ymin, ymax)
#'
#' if (require(ggplot2)) {
#'   ggplot() +
#'     geom_spatraster(data = level3_crop)
#' }
#' @seealso [`tidyterra::geom_spatraster`]
#' @importFrom tidyterra geom_spatraster
setGeneric(
  "geom_spatraster",
  function(mapping = tidyterra::aes(), data, ...) standardGeneric("geom_spatraster")
)

#' Maps geom_spatraster function from tidyterra
#'
#' @param mapping Set of aesthetic mappings created by [`ggplot2::aes()`] or
#' [`ggplot2::aes()`]. See Aesthetics specially in the use of fill aesthetic.
#' @param data A `SpatRaster` object
#' @param ... Other arguments passed on to [`ggplot2::layer()`]. These are often
#' aesthetics, used to set an aesthetic to a fixed value, like colour = "red" or
#' size = 3. They may also be parameters to the paired geom/stat.
#'
#' @return A ggplot2 layer
#' @seealso [`tidyterra::geom_spatraster`]
#' @importFrom tidyterra geom_spatraster
#' @export
setMethod(
  f = "geom_spatraster",
  signature(mapping = "ANY", data = "gedi.level3"),
  definition = function(mapping = tidyterra::aes(), data, ...) {
    tidyterra::geom_spatraster(mapping = mapping, data = data, ...)
  }
)

#' ggplot2 function for gedi.level3 data
#'
#' @description Creates the ggplot2 layer for gedi.level3
#'
#' @param mapping Set of aesthetic mappings created by [`ggplot2::aes()`] or
#' [`ggplot2::aes()`]. See Aesthetics specially in the use of fill aesthetic.
#' @param data A [`gedi.level3-class`] object
#' @param ... Other arguments passed on to [`ggplot2::layer()`]. These are often
#' aesthetics, used to set an aesthetic to a fixed value, like colour = "red" or
#' size = 3. They may also be parameters to the paired geom/stat.
#'
#' @return A ggplot2 layer
#' @export
setMethod(
  f = "geom_spatraster",
  signature(mapping = "ANY", data = "gedi.level3"),
  definition = function(mapping = tidyterra::aes(), data, ...) {
    tidyterra::geom_spatraster(mapping = mapping, data = data@raster, ...)
  }
)

#' @export
crop <- terra::crop

#' Crop gedi.level3 using by extents
#'
#' @param x A [`gedi.level3-class`] object
#' @param y numeric. The lowest longitude xmin
#' @param xb numeric. The maximum longitude xmax
#' @param xc numeric. The mininum latitude ymin
#' @param xd numeric. The maximum latitude ymax
#'
#' @return The cropped gedi.level3
#'
#' @export
setMethod(
  "crop",
  signature("gedi.level3", "numeric"),
  function(x, y, xb, xc, xd) {
    clipLevel3(x, y, xb, xc, xd)
  }
)

#' Crop gedi.level3 using by extents
#'
#' @param x A [`gedi.level3-class`] object
#' @param y SpatVect object from loaded with [`terra::vect`].
#' @param xb Not used.
#' @param xc Not used.
#' @param xd Not used.
#'
#' @return The cropped gedi.level3
#'
#' @examples
#' level3path <- "C:/Users/caiohamamura/Downloads/GEDI03_rh100_mean_2019108_2020287_002_02.tif"
#' level3 <- readLevel3(level3path = level3path)
#'
#'
#' # Specifying the path to shapefile
#' polygon_filepath <- system.file("extdata", "stands_cerrado.shp", package = "rGEDI")
#'
#' # Reading shapefile as SpatialPolygonsDataFrame object
#' library(terra)
#' polygon <- terra::vect(polygon_filepath)
#'
#' # Clipping GEDI Full Waveform Geolocations by Geometry
#' level3_clip <- crop(level3, polygon, "id")
#'
#' @export
setMethod(
  "crop",
  signature("gedi.level3", "SpatVector"),
  function(x, y, xb = NULL, xc = NULL, xd = NULL) {
    clipLevel3Geometry(x, y, xb)
  }
)




#' Clip GEDI Full Waveform Geolocations by geometry
#'
#' @description This function clips level3 extracted geolocation (level3)
#' data within a given geometry
#'
#' @param level3 A [`gedi.level3-class`] resulting from [`readLevel3()`] function.
#' @param polygon Polygon. An object of class `SpatVect`,
#' which can be loaded as an ESRI shapefile using [terra::vect] function in the \emph{terra} package.
#' @param split_by Polygon id. If defined, GEDI data will be clipped by each polygon using the polygon id from
#' table of attribute defined by the user.
#'
#' @return Returns a list of S4 object of class [`data.table::data.table-class`] containing the
#' clipped GEDI level3 extracted geolocations.
#'
#' @seealso \url{https://daac.ornl.gov/GEDI/guides/GEDI_L3_LandSurface_Metrics_V2.html}
#'
#' @examples
#' # Specifying the path to GEDI level3 data (zip file)
#' outdir <- tempdir()
#'
#' level3path <- "C:/Users/caiohamamura/Downloads/GEDI03_rh100_mean_2019108_2020287_002_02.tif"
#' level3 <- readLevel3(level3path = level3path)
#'
#'
#' # Specifying the path to shapefile
#' polygon_filepath <- system.file("extdata", "stands_cerrado.shp", package = "rGEDI")
#'
#' # Reading shapefile as SpatialPolygonsDataFrame object
#' library(terra)
#' polygon <- terra::vect(polygon_filepath)
#'
#' # Clipping GEDI Full Waveform Geolocations by Geometry
#' level3_clip <- clipLevel3Geometry(level3, polygon, split_by = "id")
#'
#' # Calculate global means
#' lapply(level3_clip, terra::global, mean)
#'
#' # Plot all in a grid
#' n <- length(level3_clip)
#' sqrt_n <- ceiling(sqrt(n))
#' oldpar <- par()
#' par(mfrow = c(sqrt_n, sqrt_n))
#'
#' for (l in level3_clip) {
#'   plot(l)
#' }
#' par(oldpar)
#'
#' @export
clipLevel3Geometry <- function(level3, polygon, split_by = "id") {
  # Initialize an empty list to store cropped raster objects
  cropped_rasts <- list()

  # Polygon
  polygon <- terra::project(polygon, terra::crs(level3@raster))

  # Iterate through each polygon, crop the raster, and store the results
  for (pol_id in seq_along(polygon)) {
    pol <- polygon[pol_id]
    cropped <- terra::crop(level3@raster, pol) # Crop the raster
    cropped_rasts[[pol[[split_by]][[1]]]] <- cropped # Store the cropped raster in the list
  }

  return(cropped_rasts)
}
