#' Clip GEDI Elevation and Height Metrics by Coordinates
#'
#' @description This function clips GEDI Level2A extracted Elevation and Height Metrics (Level2AM)
#' within a given bounding coordinates
#'
#' @usage clipLevel2AM(level2AM, xmin, xmax, ymin, ymax)
#'
#' @param level2AM A GEDI Level2A object (output of [readLevel2A()] function).
#' An S4 object of class "gedi.level2a".
#' @param xmin Numeric. West longitude (x) coordinate of bounding rectangle, in decimal degrees.
#' @param xmax Numeric. East longitude (x) coordinate of bounding rectangle, in decimal degrees.
#' @param ymin Numeric. South latitude (y) coordinate of bounding rectangle, in decimal degrees.
#' @param ymax Numeric. North latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'
#' @return Returns an S4 object of class [data.table::data.table]
#' containing the clipped elevation and relative heights metrics.
#'
#' @seealso \url{https://lpdaac.usgs.gov/products/gedi02_av002/}
#'
#' @examples
#' # Specifying the path to GEDI level2A data (zip file)
#' outdir <- tempdir()
#' level2A_fp_zip <- system.file("extdata",
#'   "GEDI02_A_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'   package = "rGEDI"
#' )
#'
#' # Unzipping GEDI level2A data
#' level2Apath <- unzip(level2A_fp_zip, exdir = outdir)
#'
#' # Reading GEDI level2A data (h5 file)
#' level2a <- readLevel2A(level2Apath = level2Apath)
#'
#' # Extracting GEDI Elevation and Height Metrics
#' level2AM <- getLevel2AM(level2a)
#'
#' # Bounding rectangle coordinates
#' xmin <- -44.15036
#' xmax <- -44.10066
#' ymin <- -13.75831
#' ymax <- -13.71244
#'
#' # Clipping GEDI data by boundary box extent
#' level2AM_clip <- clipLevel2AM(level2AM, xmin, xmax, ymin, ymax)
#'
#' close(level2a)
#' @import hdf5r stats
#' @export
clipLevel2AM <- function(level2AM, xmin, xmax, ymin, ymax) {
  # xmin ymin xmax ymax
  mask <-
    level2AM$lon_lowestmode >= xmin &
      level2AM$lon_lowestmode <= xmax &
      level2AM$lat_lowestmode >= ymin &
      level2AM$lat_lowestmode <= ymax &
      level2AM$lon_lowestmode >= xmin &
      level2AM$lon_lowestmode <= xmax &
      level2AM$lat_lowestmode >= ymin &
      level2AM$lat_lowestmode <= ymax
  mask[!stats::complete.cases(mask)] <- FALSE
  mask <- (seq_along(level2AM$lat_lowestmode))[mask]
  newFile <- level2AM[mask, ]
  # newFile<- new("gedi.level1b.dt", dt = level1bdt[mask,])
  return(newFile)
}

#' Clip GEDI Elevation and Height Metrics by Coordinates
#'
#' @description This function clips GEDI Level2A extracted Elevation and Height Metrics (Level2AM)
#' within a given bounding coordinates
#'
#' @param level2AM A GEDI Level2A object (output of [readLevel2A()] function).
#' An S4 object of class "data.table".
#' @param polygon Polygon. An object of class [`sf::sf`],
#' which can be loaded as an ESRI shapefile using [sf::st_read] function in the
#' \emph{sf} package.
#' @param split_by Polygon id. If defined, GEDI data will be clipped by each polygon using
#' the polygon id from table of attribute defined by the user
#'
#' @return Returns an S4 object of class [data.table::data.table]
#' containing the clipped elevation and relative heights metrics.
#'
#' @examples
#' # Specifying the path to GEDI level2A data (zip file)
#' outdir <- tempdir()
#' level2A_fp_zip <- system.file("extdata",
#'   "GEDI02_A_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'   package = "rGEDI"
#' )
#'
#' # Unzipping GEDI level2A data
#' level2Apath <- unzip(level2A_fp_zip, exdir = outdir)
#'
#' # Reading GEDI level2A data (h5 file)
#' level2a <- readLevel2A(level2Apath = level2Apath)
#'
#' # Extracting GEDI Elevation and Height Metrics
#' level2AM <- getLevel2AM(level2a)
#'
#' # Specifying the path to shapefile
#' polygon_filepath <- system.file("extdata", "stands_cerrado.shp", package = "rGEDI")
#'
#' # Reading shapefile as sf object
#' library(sf)
#' polygon <- sf::st_read(polygon_filepath)
#'
#' # Clipping GEDI data by Geometry
#' level2AM_clip <- clipLevel2AMGeometry(level2AM, polygon, split_by = "id")
#'
#' hasLeaflet <- require(leaflet)
#'
#' if (hasLeaflet) {
#'   leaflet() %>%
#'     addCircleMarkers(level2AM_clip$lat_lowestmode,
#'       level2AM_clip$lon_lowestmode,
#'       radius = 1,
#'       opacity = 1,
#'       color = "red"
#'     ) %>%
#'     addScaleBar(options = list(imperial = FALSE)) %>%
#'     addPolygons(
#'       data = polygon, weight = 1, col = "white",
#'       opacity = 1, fillOpacity = 0
#'     ) %>%
#'     addProviderTiles(providers$Esri.WorldImagery)
#' }
#'
#' close(level2a)
#' @export
clipLevel2AMGeometry <- function(level2AM, polygon, split_by = "id") {
  exshp <- sf::st_bbox(polygon)
  level2adt <- clipLevel2AM(
    level2AM,
    xmin = exshp$xmin,
    xmax = exshp$xmax,
    ymin = exshp$ymin,
    ymax = exshp$ymax
  )
  if (nrow(level2adt) == 0) {
    print("The polygon does not overlap the GEDI data")
  } else {
    .I <- data.table::.I
    level2adt[, id := .I]
    points <- sf::st_as_sf(
      level2adt,
      coords = c("lon_lowestmode", "lat_lowestmode"),
      crs = sf::st_crs(polygon)
    )

    names(points) <- gsub("^(?!geometry)", "x_\\1", names(points), perl = TRUE)
    pts <- sf::st_intersection(sf::st_make_valid(points), sf::st_make_valid(polygon))

    if (!is.null(split_by)) {
      if (any(names(polygon) == split_by)) {
        mask <- as.integer(pts$x_id)
        newFile <- level2adt[mask, ]
        newFile$poly_id <- pts[[split_by]]
      } else {
        stop(paste("The", split_by, "is not included in the attribute table.
                       Please check the names in the attribute table"))
      }
    } else {
      mask <- as.integer(rownames(pts))
      newFile <- level2adt[mask, ]
    }
    # newFile<- new("gedi.level1b.dt", dt = level2adt2@dt[mask,])
    return(newFile)
  }
}
