#' Clip GEDI Level 4A by Coordinates
#'
#' @description This function clips GEDI Level 4A data within given bounding coordinates
#'
#'
#' @usage clipLevel4A(level4A, xmin, xmax, ymin, ymax)
#'
#'
#' @param level4A A GEDI Level4A object (output of [readLevel4A()] function).
#' An S4 object of class "data.table".
#' @param xmin Numeric. West longitude (x) coordinate of the bounding rectangle, in decimal degrees.
#' @param xmax Numeric. East longitude (x) coordinate of the bounding rectangle, in decimal degrees.
#' @param ymin Numeric. South latitude (y) coordinate of the bounding rectangle, in decimal degrees.
#' @param ymax Numeric. North latitude (y) coordinate of the bounding rectangle, in decimal degrees.
#'
#' @return Returns an S4 object of class [data.table::data.table]
#' containing the AGBD and standard error estimates.
#'
#' @seealso \url{https://daac.ornl.gov/GEDI/guides/GEDI_L4A_AGB_Density.html}
#'
#' @examples
#' # Specifying the path to GEDI level4A data (zip file)
#' outdir <- tempdir()
#' level4A_fp_zip <- system.file("extdata",
#'   "GEDI04_A_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'   package = "rGEDI"
#' )
#'
#' # Unzipping GEDI leve42A data
#' level4Apath <- unzip(level4A_fp_zip, exdir = outdir)
#'
#' # Reading GEDI level2B data (h5 file)
#' level4a <- readLevel4A(level4Apath = level4Apath)
#'
#' # Extracting canopy cover and vertical profile metrics
#' level4AVPM <- getLevel4AVPM(level4a)
#'
#' # Bounding rectangle coordinates
#' xmin <- -44.15036
#' xmax <- -44.10066
#' ymin <- -13.75831
#' ymax <- -13.71244
#'
#' # Clipping level4a by extent boundary box
#' level4a_clip <- clipLevel4a(level4a, xmin, xmax, ymin, ymax)
#'
#' hasLeaflet <- require(leaflet)
#'
#' if (hasLeaflet) {
#'   leaflet() %>%
#'     addCircleMarkers(level4a_clip$lon_lowestmode,
#'       level4a_clip$lat_lowestmode,
#'       radius = 1,
#'       opacity = 1,
#'       color = "red"
#'     ) %>%
#'     addScaleBar(options = list(imperial = FALSE)) %>%
#'     addProviderTiles(providers$Esri.WorldImagery)
#' }
#' close(level4a)
#' @export
clipLevel4A <- function(level4A, xmin, xmax, ymin, ymax) {
  # xmin ymin xmax ymax
  mask <-
    level4A$lon_lowestmode >= xmin &
      level4A$lon_lowestmode <= xmax &
      level4A$lat_lowestmode >= ymin &
      level4A$lat_lowestmode <= ymax


  mask[!stats::complete.cases(mask)] <- FALSE
  mask <- seq_len(length(level4A$lon_lowestmode))[mask]
  newFile <- level4A[mask, ]
  return(newFile)
}

#' Clip GEDI Level 4A by geometry
#'
#' @description This function clips GEDI Level 4A data within a given geometry
#'
#' @usage clipLevel4AGeometry(level4A, polygon_spdf, split_by)
#'
#'
#' @param level4A A GEDI Level4A object (output of [readLevel4A()] function).
#' An S4 object of class "gedi.level4a".
#' @param polygon_spdf Polygon. An object of class [`sp::SpatialPolygonsDataFrame-class`],
#' which can be loaded as an ESRI shapefile using [raster::shapefile] function in the \emph{raster} package.
#' @param split_by Polygon id. If defined, GEDI data will be clipped by each polygon using the attribute
#' specified by `split_by` from the attribute table.
#'
#' @return Returns an S4 object of class [data.table::data.table]
#' containing the AGBD and standard error estiates.
#'
#' @seealso \url{https://daac.ornl.gov/GEDI/guides/GEDI_L4A_AGB_Density.html}
#'
#' @examples
#' # Specifying the path to GEDI level4A data (zip file)
#' outdir <- tempdir()
#' level4A_fp_zip <- system.file("extdata",
#'   "GEDI04_A_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'   package = "rGEDI"
#' )
#'
#' # Unzipping GEDI leve42A data
#' level4Apath <- unzip(level4A_fp_zip, exdir = outdir)
#'
#' # Reading GEDI level2B data (h5 file)
#' level4a <- readLevel4A(level4Apath = level4Apath)
#'
#' # Extracting canopy cover and vertical profile metrics
#' level4AVPM <- getLevel4AVPM(level4a)
#'
#' # Specifying the path to shapefile
#' polygon_filepath <- system.file("extdata", "stands_cerrado.shp", package = "rGEDI")
#'
#' # Reading shapefile as SpatialPolygonsDataFrame object
#' library(raster)
#' polygon_spdf <- shapefile(polygon_filepath)
#'
#' # Clipping level4A by geometry
#' level4a_clip_geometry <- clipLevel4AGeometry(level4A, polygon_spdf, split_by = "id")
#'
#' hasLeaflet <- require(leaflet)
#'
#' if (hasLeaflet) {
#'   leaflet() %>%
#'     addCircleMarkers(level4a_clip$lon_lowestmode,
#'       level4a_clip$lat_lowestmode,
#'       radius = 1,
#'       opacity = 1,
#'       color = "red"
#'     ) %>%
#'     addScaleBar(options = list(imperial = FALSE)) %>%
#'     addProviderTiles(providers$Esri.WorldImagery)
#' }
#' close(level4a)
#' @export
clipLevel4AGeometry <- function(level4A, polygon_spdf, split_by = NULL) {
  exshp <- raster::extent(polygon_spdf)
  level4adt <- clipLevel4A(level4A, xmin = exshp[1], xmax = exshp[2], ymin = exshp[3], ymax = exshp[4])

  if (nrow(level4adt) == 0) {
    print("The polygon does not overlap the GEDI data")
  } else {
    points <- sp::SpatialPointsDataFrame(
      coords = matrix(c(level4adt$lon_lowestmode, level4adt$lat_lowestmode), ncol = 2),
      data = data.frame(id = seq_len(length(level4adt$lat_lowestmode))), proj4string = polygon_spdf@proj4string
    )
    pts <- raster::intersect(points, polygon_spdf)
    colnames(pts@data) <- c("rowids", names(polygon_spdf))

    if (!is.null(split_by)) {
      if (any(names(polygon_spdf) == split_by)) {
        mask <- as.integer(pts@data$rowids)
        newFile <- level4adt[mask, ]
        newFile$poly_id <- pts@data[, split_by]
      } else {
        stop(paste("The", split_by, "is not included in the attribute table.
                       Please check the names in the attribute table"))
      }
    } else {
      mask <- as.integer(pts@data$rowids)
      newFile <- level4adt[mask, ]
    }
    return(newFile)
  }
}