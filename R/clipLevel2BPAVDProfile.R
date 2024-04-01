#' Clip GEDI Plant Area Volume Density profile by Coordinates
#'
#' @description This function clips GEDI level2B derived
#' Plant Area Volume Density profile within a given bounding coordinates
#'
#' @usage clipLevel2BPAVDProfile(level2BPAVDProfile, xmin, xmax, ymin, ymax)
#'
#' @param level2BPAVDProfile A GEDI Level2B object (output of [getLevel2BPAVDProfile()] function).
#' An S4 object of class "data.table".
#' @param xmin Numeric. West longitude (x) coordinate of the bounding rectangle, in decimal degrees.
#' @param xmax Numeric. East longitude (x) coordinate of the bounding rectangle, in decimal degrees.
#' @param ymin Numeric. South latitude (y) coordinate of the bounding rectangle, in decimal degrees.
#' @param ymax Numeric. North latitude (y) coordinate of the bounding rectangle, in decimal degrees.
#'
#' @return Returns an S4 object of class [data.table::data.table]
#' containing the Plant Area Volume Density profile data.
#'
#' @seealso \url{https://lpdaac.usgs.gov/products/gedi02_bv002/}
#'
#' @examples
#' # specify the path to GEDI level2B data (zip file)
#' outdir <- tempdir()
#' level2B_fp_zip <- system.file("extdata",
#'   "GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'   package = "rGEDI"
#' )
#'
#' # Unzipping GEDI level2A data
#' level2Bpath <- unzip(level2B_fp_zip, exdir = outdir)
#'
#' # Reading GEDI level2B data (h5 file)
#' level2b <- readLevel2B(level2Bpath = level2Bpath)
#'
#' # Extracting GEDI Plant Area Volume Density profile
#' level2BPAVDProfile <- getLevel2BPAVDProfile(level2b)
#'
#' # Bounding rectangle coordinates
#' xmin <- -44.15036
#' xmax <- -44.10066
#' ymin <- -13.75831
#' ymax <- -13.71244
#'
#' # Clipping GEDI Plant Area Volume Density profile by boundary box extent
#' level2BPAVDProfile_clip <- clipLevel2BPAVDProfile(level2BPAVDProfile, xmin, xmax, ymin, ymax)
#'
#' close(level2b)
#' @export
clipLevel2BPAVDProfile <- function(level2BPAVDProfile, xmin, xmax, ymin, ymax) {
  # xmin ymin xmax ymax
  mask <-
    level2BPAVDProfile$lon_lowestmode >= xmin &
      level2BPAVDProfile$lon_lowestmode <= xmax &
      level2BPAVDProfile$lat_lowestmode >= ymin &
      level2BPAVDProfile$lat_lowestmode <= ymax

  mask[!stats::complete.cases(mask)] <- FALSE
  mask <- (seq_along(level2BPAVDProfile$lon_lowestmode))[mask]
  newFile <- level2BPAVDProfile[mask, ]
  # newFile<- new("gedi.level1b.dt", dt = level2bdt[mask,])
  if (nrow(newFile) == 0) {
    print("The polygon does not overlap the GEDI data")
  } else {
    return(newFile)
  }
}

#' Clip GEDI Plant Area Volume Density profile by geometry
#'
#' @description This function clips GEDI level2B derived
#' Plant Area Index profile within a given geometry
#'
#' @param level2BPAVDProfile A GEDI Level2B object (output of [getLevel2BPAIProfile()] function).
#' An S4 object of class "gedi.level2b".
#' @param polygon Polygon. An object of class [`sf::sf`],
#' which can be loaded as an ESRI shapefile using [sf::st_read] function in the \emph{sf} package.
#' @param split_by Polygon id. If defined, GEDI data will be clipped by each polygon using the attribute specified by `split_by` from the attribute table.
#'
#' @return Returns an S4 object of class [data.table::data.table]
#' containing the Plant Area Volume Density profile data.
#'
#' @seealso \url{https://lpdaac.usgs.gov/products/gedi02_bv002/}
#'
#' @examples
#' # Specifying the path to GEDI level2B data (zip file)
#' outdir <- tempdir()
#' level2B_fp_zip <- system.file("extdata",
#'   "GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'   package = "rGEDI"
#' )
#'
#' # Unzipping GEDI level2A data
#' level2Bpath <- unzip(level2B_fp_zip, exdir = outdir)
#'
#' # Reading GEDI level2B data (h5 file)
#' level2b <- readLevel2B(level2Bpath = level2Bpath)
#'
#' # Extracting GEDI Plant Area Volume Density profile
#' level2BPAVDProfile <- getLevel2BPAVDProfile(level2b)
#'
#' # Specifying the path to shapefile
#' polygon_filepath <- system.file("extdata", "stands_cerrado.shp", package = "rGEDI")
#'
#' # Reading shapefile as sf object
#' library(sf)
#' polygon <- sf::st_read(polygon_filepath)
#'
#' # Clipping GEDI Plant Area Volume Density profile by geometry
#' level2BPAVDProfile_clip <- clipLevel2BPAVDProfileGeometry(
#'   level2BPAVDProfile,
#'   polygon,
#'   split_by = "id"
#' )
#'
#' close(level2b)
#' @export
clipLevel2BPAVDProfileGeometry <- function(level2BPAVDProfile, polygon, split_by = NULL) {
  exshp <- sf::st_bbox(polygon)
  level2bdt <- clipLevel2BPAIProfile(
    level2BPAVDProfile,
    xmin = exshp$xmin,
    xmax = exshp$xmax,
    ymin = exshp$ymin,
    ymax = exshp$ymax
  )

  if (nrow(level2bdt) == 0) {
    print("The polygon does not overlap the GEDI data")
  } else {
    .I <- data.table::.I
    level2bdt[, id := .I]

    points <- sf::st_as_sf(
      level2bdt,
      coords = c("lon_lowestmode", "lat_lowestmode"),
      crs = sf::st_crs(polygon)
    )

    names(points) <- gsub("^(?!geometry)", "x_\\1", names(points), perl = TRUE)
    pts <- sf::st_intersection(sf::st_make_valid(points), sf::st_make_valid(polygon))

    if (!is.null(split_by)) {
      if (any(names(polygon) == split_by)) {
        mask <- as.integer(pts$x_id)
        newFile <- level2bdt[mask, ]
        newFile$poly_id <- pts[[split_by]]
      } else {
        stop(paste("The", split_by, "is not included in the attribute table.
                       Please check the names in the attribute table"))
      }
    } else {
      mask <- as.integer(rownames(pts))
      newFile <- level2bdt[mask, ]
    }
    # newFile<- new("gedi.level1b.dt", dt = level2bdt2@dt[mask,])
    return(newFile)
  }
}
