#' Clip GEDI Level2A data by Coordinates
#'
#' @description This function clips GEDI Level2A data within a given bounding coordinates
#'
#' @usage clipLevel2A(level2a, xmin, xmax, ymin, ymax, output)
#'
#' @param level2a A GEDI Level2A object (output of [readLevel2A()] function).
#' An S4 object of class "gedi.level2a".
#' @param xmin Numeric. West longitude (x) coordinate of the bounding rectangle, in decimal degrees.
#' @param xmax Numeric. East longitude (x) coordinate of the bounding rectangle, in decimal degrees.
#' @param ymin Numeric. South latitude (y) coordinate of the bounding rectangle, in decimal degrees.
#' @param ymax Numeric. North latitude (y) coordinate of the bounding rectangle, in decimal degrees.
#' @param output Optional character path where to save the new hdf5file. The default stores a
#' temporary file only.
#'
#' @return Returns a list of S4 objects of class "gedi.level2a" containing clipped
#' GEDI Level2A data.
#'
#' @seealso \url{https://lpdaac.usgs.gov/products/gedi02_av002/}
#'
#' @examples
#' \donttest{
#' outdir <- tempdir()
#'
#' # Specifying the path to GEDI level2A data (zip file)
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
#' # Bounding rectangle coordinates
#' xmin <- -44.13
#' xmax <- -44.12
#' ymin <- -13.74
#' ymax <- -13.73
#'
#' print(level2a)
#'
#' # Specifying output file and path
#' output <- file.path(outdir, "GEDI02_A_2019108080338_O01964_T05337_02_001_01_clip.h5")
#'
#' # Clipping GEDI Level2A data by boundary box extent
#' level2a_clip <- clipLevel2A(level2a, xmin, xmax, ymin, ymax, output)
#'
#' close(level2a)
#' close(level2a_clip)
#' }
#' @export
clipLevel2A <- function(level2a, xmin, xmax, ymin, ymax, output = "") {
  output <- checkOutput(output)
  checkClipExtentInputs(level2a, "gedi.level2a", xmin, xmax, ymin, ymax)

  # Get all spatial data as a list of dataframes with spatial information
  spData <- getSpatialData2A(level2a)

  masks <- clipSpDataByExtentLevel2A(spData, xmin, xmax, ymin, ymax)
  print(masks)

  newFile <- clipByMask2A(
    level2a,
    masks,
    output
  )
  output <- newFile@h5$filename
  close(newFile)
  result <- readLevel2A(output)

  return(result)
}

#' Clip GEDI Level2A data by geometry
#'
#' @description This function clips GEDI Level2A data within a given geometry
#'
#' @param level2a A GEDI Level2A object (output of [readLevel2A()] function).
#' An S4 object of class "gedi.level2a".
#' @param polygon Polygon. An object of class [`sf::sf`],
#' which can be loaded as an ESRI shapefile using [sf::st_read()] function in the \emph{sf} package.
#' @param output optional character path where to save the new h5file. Default "" (temporary file).
#' @param split_by Polygon id. If defined, GEDI data will be clipped by each polygon using the
#' attribute specified by `split_by` from the attribute table.
#'
#' @return Returns a list of S4 object of class "gedi.level2a" containing clipped GEDI Level2A data.
#'
#' @seealso \url{https://lpdaac.usgs.gov/products/gedi02_av002/}
#'
#' @examples
#' \donttest{
#' outdir <- tempdir()
#'
#' # Specifying the path to GEDI level2A data (zip file)
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
#' # Specifying the path to shapefile
#' polygon_filepath <- system.file("extdata", "stands_cerrado.shp", package = "rGEDI")
#'
#' # Reading shapefile as sf object
#' library(sf)
#' polygon <- sf::st_read(polygon_filepath)
#'
#' # Specifying output file and path
#' output <- file.path(outdir, "GEDI02_A_2019108080338_O01964_T05337_02_001_01_clip")
#'
#' # Clipping GEDI Level2A data by geometry
#' level2a_clip <- clipLevel2AGeometry(level2a,
#'   polygon = polygon,
#'   output = output,
#'   split_by = "id"
#' )
#' close(level2a)
#' lapply(level2a_clip, close)
#' }
#' @export
clipLevel2AGeometry <- function(level2a, polygon, output = "", split_by = NULL) {
  output <- checkOutput(output)
  checkClipGeoInputs(level2a, "gedi.level2a", polygon, split_by)

  spData <- getSpatialData2A(level2a)

  bbox <- sf::st_bbox(polygon)
  xmin <- bbox$xmin
  xmax <- bbox$xmax
  ymin <- bbox$ymin
  ymax <- bbox$ymax

  masks <- clipSpDataByExtentLevel2A(spData, xmin, xmax, ymin, ymax)

  polygon_masks <- getPolygonMaskLevel2A(spData, masks, polygon, split_by)

  results <- clipByMasks(level2a, polygon_masks, output, split_by, clipByMask2A)

  return(results)
}


# Helper function to return spatial data within a dataframe
getSpatialData2A <- function(level2a) {
  level2a.h5 <- level2a@h5
  groups_id <- grep("BEAM\\d{4}$", gsub(
    "/", "",
    hdf5r::list.groups(level2a.h5, recursive = FALSE)
  ), value = TRUE)

  beams_spdf <- list()

  for (i in gsub("/", "", groups_id)) {
    algs_coordinates <- list()

    n_algs <- level2a.h5[[i]][["ancillary/l2a_alg_count"]][1]
    beams_spdf[[i]] <- algs_coordinates

    for (j in 1:n_algs) {
      beams_spdf[[i]][[j]] <- data.frame(
        latitude_highest = level2a.h5[[paste0(i, "/geolocation/lat_highestreturn_a", j)]][],
        latitude_lowest = level2a.h5[[paste0(i, "/geolocation/lat_lowestreturn_a", j)]][],
        longitude_highest = level2a.h5[[paste0(i, "/geolocation/lon_highestreturn_a", j)]][],
        longitude_lowest = level2a.h5[[paste0(i, "/geolocation/lon_lowestreturn_a", j)]][]
      )
    }

    beams_spdf[[i]][["main"]] <- data.frame(
      latitude_highest = level2a.h5[[i]][["lat_highestreturn"]][],
      latitude_lowest = level2a.h5[[i]][["lat_lowestmode"]][],
      longitude_highest = level2a.h5[[i]][["lon_highestreturn"]][],
      longitude_lowest = level2a.h5[[i]][["lon_lowestmode"]][]
    )
  }

  return(beams_spdf)
}

clipByMask2A <- function(level2a, masks, output = "") {
  newFile <- hdf5r::H5File$new(output, mode = "w")
  if (length(masks) == 0) {
    message("No intersection found!")
    newFile$close_all()
    newFile <- hdf5r::H5File$new(output, mode = "r")
    result <- new("gedi.level2a", h5 = newFile)
    return(result)
  }

  for (attr in hdf5r::list.attributes(level2a@h5)) {
    hdf5r::h5attr(newFile, attr) <- hdf5r::h5attr(level2a@h5, attr)
  }


  all_groups <- hdf5r::list.groups(level2a@h5)


  # Check if the beam has any intersecting area
  beams_with_value <- lapply(lapply(masks, function(x) sapply(x, length)), sum) > 0
  beams_with_value <- names(which(beams_with_value))
  beams_with_value <- c(beams_with_value, "METADATA")
  which_groups <- gsub("([^/]*).*", "\\1", all_groups) %in% beams_with_value
  groups_with_value <- all_groups[which_groups]

  # Setup progress bar
  all_datasets <- hdf5r::list.datasets(level2a@h5)
  which_datasets <- gsub("([^/]*).*", "\\1", all_datasets) %in% beams_with_value
  datasets_with_value <- all_datasets[which_datasets]

  total <- length(datasets_with_value)
  pb <- utils::txtProgressBar(min = 0, max = total, style = 3)
  progress <- 0

  for (group in groups_with_value) {
    beam_id <- strsplit(group, "/")[[1]][1]

    hdf5r::createGroup(newFile, group)
    createAttributesWithinGroup(level2a@h5, newFile, group)

    # Datasets to loop
    datasets <- hdf5r::list.datasets(level2a@h5[[group]], recursive = FALSE, full.names = TRUE)

    # Create list of algorithm ids for the datasets to choose right mask
    alg_ids <- as.list(as.integer(
      sapply(regmatches(datasets, regexec("_a(\\d+)", datasets)), function(x) x[2])
    ))
    names(alg_ids) <- datasets

    for (dt in datasets) {
      # Get right mask
      alg_id <- alg_ids[[dt]]

      if (is.na(alg_id)) {
        mask <- masks[[beam_id]][["main"]]
        beam_shot_n <- level2a@h5[[beam_id]][["shot_number"]]$dims
      } else {
        mask <- masks[[beam_id]][[alg_id]]
        beam_shot_n <- level2a@h5[[
          sprintf("%s/geolocation/elev_highestreturn_a%d", beam_id, alg_id)
        ]]$dims
      }
      mask_size <- length(mask)

      if (mask_size == 0) {
        # Update progress
        progress <- progress + 1
        utils::setTxtProgressBar(pb, progress)
        next
      }

      applyMaskAndCreateDataset(level2a@h5, newFile, dt, mask, beam_shot_n)

      # Update progress
      progress <- progress + 1
      utils::setTxtProgressBar(pb, progress)
    }
  }

  hdf5r::h5flush(newFile)
  newFile$close_all()
  newFile <- hdf5r::H5File$new(output, mode = "r")
  result <- new("gedi.level2a", h5 = newFile)
  close(pb)
  return(result)
}
