#' Clip GEDI Level1B data by Coordinates
#'
#' @description This function clips GEDI Level1B data (geolocated waveforms) within a given
#' bounding coordinates
#'
#' @usage clipLevel1B(level1b, xmin, xmax, ymin, ymax, output)
#'
#' @param level1b A [`gedi.level1b-class`] object (output of [readLevel1B()] function).
#' An S4 object of class [`gedi.level1b-class`].
#' @param xmin Numeric. West longitude (x) coordinate of the bounding rectangle, in decimal degrees.
#' @param xmax Numeric. East longitude (x) coordinate of the bounding rectangle, in decimal degrees.
#' @param ymin Numeric. South latitude (y) coordinate of the bounding rectangle, in decimal degrees.
#' @param ymax Numeric. North latitude (y) coordinate of the bounding rectangle, in decimal degrees.
#' @param output Optional character path where to save the new hdf5file. The default stores a
#' temporary file only.
#'
#' @return Returns a list of S4 objects of class [`gedi.level1b-class`] containing
#' clipped GEDI Level1B data.
#'
#' @seealso \url{https://lpdaac.usgs.gov/products/gedi01_bv002/}
#'
#' @examples
#' \donttest{
#' # Specifying the path to GEDI level1B data (zip file)
#' outdir <- tempdir()
#'
#' level1B_fp_zip <- system.file("extdata",
#'   "GEDI01_B_2019108080338_O01964_T05337_02_003_01_sub.zip",
#'   package = "rGEDI"
#' )
#'
#' # Unzipping GEDI level1B data
#' level1Bpath <- unzip(level1B_fp_zip, exdir = outdir)
#'
#' # Reading GEDI level1B data (h5 file)
#' level1b <- readLevel1B(level1Bpath = level1Bpath)
#'
#' # Bounding rectangle coordinates
#' xmin <- -44.13
#' xmax <- -44.12
#' ymin <- -13.74
#' ymax <- -13.73
#'
#' # Specifying output file and path
#' output <- file.path(outdir, "GEDI01_B_2019108080338_O01964_T05337_02_003_01_clip")
#'
#' # Clipping GEDI Level1B data by extent boundary box
#' level1b_clip <- clipLevel1B(level1b, xmin, xmax, ymin, ymax, output)
#'
#' close(level1b)
#' close(level1b_clip)
#' }
#' @import hdf5r
#' @export
clipLevel1B <- function(level1b, xmin, xmax, ymin, ymax, output = "") {
  output <- checkOutput(output)
  checkClipExtentInputs(level1b, "gedi.level1b", xmin, xmax, ymin, ymax)

  # Get all spatial data as a list of dataframes with spatial information
  spData <- getSpatialData1B(level1b)

  masks <- clipSpDataByExtentLevelB(spData, xmin, xmax, ymin, ymax)

  newFile <- clipByMask1B(
    level1b,
    masks,
    output
  )
  output <- newFile@h5$filename
  close(newFile)
  result <- readLevel1B(output)

  return(result)
}

#' Clip GEDI Level1B data by geometry
#'
#' @description This function clips GEDI Level1B (geolocated waveforms) data within a
#' given bounding geometry
#'
#'
#' @param level1b A [`gedi.level1b-class`] object (output of [readLevel1B()] function).
#' An S4 object of class "gedi.level1b".
#' @param polygon SpatVect. An object opened with `terra::vect`,
#' @param output Optional character path where to save the new
#' [`hdf5r::H5File`][hdf5r::H5File-class]. The default stores a temporary file only.
#' @param split_by Polygon id. If defined, GEDI data will be clipped by each polygon using
#' the attribute specified by `split_by` from the attribute table.
#'
#' @return Returns a list of S4 object of class [`gedi.level1b-class`] containing clipped
#' GEDI Level1B data.
#'
#' @examples
#' \donttest{
#' outdir <- tempdir()
#'
#' # Specifying the path to GEDI level1B data (zip file)
#' level1B_fp_zip <- system.file("extdata",
#'   "GEDI01_B_2019108080338_O01964_T05337_02_003_01_sub.zip",
#'   package = "rGEDI"
#' )
#'
#' # Unzipping GEDI level1B data
#' level1Bpath <- unzip(level1B_fp_zip, exdir = outdir)
#'
#' # Reading GEDI level1B data (h5 file)
#' level1b <- readLevel1B(level1Bpath = level1Bpath)
#'
#' # Specifying the path to shapefile
#' polygon_filepath <- system.file("extdata", "stands_cerrado.shp", package = "rGEDI")
#'
#' # Reading shapefile as SpatVect object
#' library(terra)
#' polygon <- terra::vect(polygon_filepath)
#'
#' # Spepecifing output file and path
#' output <- file.path(outdir, "GEDI01_B_2019108080338_O01964_T05337_02_003_01_clip")
#'
#' # clipping GEDI Level1B data by extent boundary box
#' level1b_clip <- clipLevel1BGeometry(level1b,
#'   polygon = polygon,
#'   output = output,
#'   split_by = "id"
#' )
#'
#' close(level1b)
#' lapply(level1b_clip, close)
#' }
#' @import hdf5r
#' @export
clipLevel1BGeometry <- function(level1b, polygon, output = "", split_by = NULL) {
  output <- checkOutput(output)
  checkClipGeoInputs(level1b, "gedi.level1b", polygon, split_by)

  spData <- getSpatialData1B(level1b)

  bbox <- terra::ext(polygon)
  xmin <- bbox$xmin
  xmax <- bbox$xmax
  ymin <- bbox$ymin
  ymax <- bbox$ymax

  masks <- clipSpDataByExtentLevelB(spData, xmin, xmax, ymin, ymax)
  polygon_masks <- getPolygonMaskLevelB(spData, masks, polygon, split_by)
  results <- clipByMasks(level1b, polygon_masks, output, split_by, clipByMask1B)

  return(results)
}


# Helper function to return spatial data within a dataframe
getSpatialData1B <- function(level1b) {
  level1b.h5 <- level1b@h5
  groups_id <- grep("BEAM\\d{4}$", gsub(
    "/", "",
    hdf5r::list.groups(level1b.h5, recursive = FALSE)
  ), value = TRUE)

  beams_spdf <- list()

  for (i in gsub("/", "", groups_id)) {
    beams_spdf[[i]] <- data.frame(
      latitude_bin0 = level1b.h5[[paste0(i, "/geolocation/latitude_bin0")]][],
      latitude_lastbin = level1b.h5[[paste0(i, "/geolocation/latitude_lastbin")]][],
      longitude_bin0 = level1b.h5[[paste0(i, "/geolocation/longitude_bin0")]][],
      longitude_lastbin = level1b.h5[[paste0(i, "/geolocation/longitude_lastbin")]][]
    )
  }

  return(beams_spdf)
}

handleNoIntersection <- function(h5File, output, type) {
  message("No intersection found!")
  h5File$close_all()
  h5File <- hdf5r::H5File$new(output, mode = "r")
  result <- new(type, h5 = h5File)
  return(result)
}

copyRootAttributes <- function(inH5, outH5) {
  for (attr in hdf5r::list.attributes(inH5)) {
    hdf5r::h5attr(outH5, attr) <- hdf5r::h5attr(inH5, attr)
  }
}

getGroupsWithValue <- function(h5, masks) {
  all_groups <- hdf5r::list.groups(h5)
  beams_with_value <- sapply(masks, length) > 0
  beams_with_value <- names(which(beams_with_value))
  beams_with_value <- c(beams_with_value, "METADATA")
  which_groups <- gsub("([^/]*).*", "\\1", all_groups) %in% beams_with_value

  return(all_groups[which_groups])
}

countDatasetsWithData <- function(h5, groups_with_value) {
  all_datasets <- hdf5r::list.datasets(h5)
  which_datasets <- gsub("([^/]*).*", "\\1", all_datasets) %in% groups_with_value
  datasets_with_value <- all_datasets[which_datasets]

  return(length(datasets_with_value))
}

clipByMask1B <- function(level1b, masks, output = "") {
  newFile <- hdf5r::H5File$new(output, mode = "w")
  if (length(masks) == 0) {
    return(handleNoIntersection(newFile, output, "gedi.level1b"))
  }

  copyRootAttributes(level1b@h5, newFile)
  groups_with_value <- getGroupsWithValue(level1b@h5, masks)

  # Setup progress bar
  n_datasets_with_data <- countDatasetsWithData(level1b@h5, groups_with_value)
  pb <- utils::txtProgressBar(min = 0, max = n_datasets_with_data, style = 3)
  progress <- 0

  for (group in groups_with_value) {
    beam_id <- strsplit(group, "/")[[1]][1]
    mask <- masks[[beam_id]]
    mask_size <- length(mask)

    hdf5r::createGroup(newFile, group)
    createAttributesWithinGroup(level1b@h5, newFile, group)

    if (beam_id != "METADATA") {
      beam_shot_n <- level1b@h5[[beam_id]][["shot_number"]]$dims
      total_waveforms <- list()
      total_waveforms[["rx"]] <- sum(level1b@h5[[beam_id]][["rx_sample_count"]][])
      total_waveforms[["tx"]] <- sum(level1b@h5[[beam_id]][["tx_sample_count"]][])
    }


    for (dt in hdf5r::list.datasets(level1b@h5[[group]], recursive = FALSE, full.names = TRUE)) {
      if (grepl("[rt]x_sample_start_index$", dt)) {
        progress <- progress + 1
        utils::setTxtProgressBar(pb, progress)
        next
      }
      h5_dt <- level1b@h5[[dt]]
      dt_dim <- h5_dt$dims
      dtype <- h5_dt$get_type()
      if (is.na(all(h5_dt$chunk_dims))) {
        chunkdims <- NULL
      } else {
        chunkdims <- h5_dt$chunk_dims
      }


      if (length(dt_dim) == 1) {
        if (dt_dim == 1) {
          hdf5r::createDataSet(newFile, dt, h5_dt[], dtype = dtype, chunk_dim = chunkdims)
        } else if (dt_dim == beam_shot_n) {
          hdf5r::createDataSet(newFile, dt, h5_dt[mask], dtype = dtype, chunk_dim = chunkdims)
        } else if (dt_dim %in% total_waveforms || dt_dim %% beam_shot_n == 0) {
          prefix <- ifelse(substr(basename(dt), 1, 2) == "rx", "rx", "tx")
          sampleCount <- sprintf("%s_sample_count", prefix)
          sampleStartIndex <- sprintf("%s_sample_start_index", prefix)
          countPath <- file.path(beam_id, sampleCount, fsep = "/")
          startIdxPath <- file.path(beam_id, sampleStartIndex, fsep = "/")

          counts <- level1b@h5[[countPath]][mask]
          startIdx <- level1b@h5[[startIdxPath]][mask]
          v.seq <- Vectorize(
            seq.default,
            vectorize.args = c("from", "length.out"),
            SIMPLIFY = TRUE
          )
          mask_list_waveform <- v.seq(startIdx, length.out = counts)
          mask_waveform <- Reduce(c, mask_list_waveform)
          newStartIdx <- c(1, cumsum(counts) + 1)[-mask_size - 1]
          newFile$create_dataset(
            name       = startIdxPath,
            robj       = newStartIdx,
            dims       = length(newStartIdx),
            chunk_dims = level1b@h5[[startIdxPath]]$chunk_dims,
            dtype      = level1b@h5[[startIdxPath]]$get_type()
          )

          total_size <- length(mask_waveform)
          chunk_part <- 1
          dt_res <- hdf5r::createDataSet(
            newFile,
            dt,
            dtype = dtype,
            chunk_dim = chunkdims,
            dims = total_size
          )
          while (chunk_part < total_size) {
            end <- chunk_part + chunkdims - 1
            if (end > total_size) {
              end <- total_size
            }
            get_part <- mask_waveform[(chunk_part):(end)]
            dt_res[chunk_part:end] <- h5_dt[get_part]
            chunk_part <- end + 1
          }
        }
      } else if (length(dt_dim) == 2 && dt_dim[1] == beam_shot_n) {
        if (length(mask) == 1) {
          chunkdims <- chunkdims[[2]]
        }
        hdf5r::createDataSet(
          newFile,
          dt,
          level1b@h5[[dt]][mask, ],
          dtype = dtype,
          chunk_dim = chunkdims
        )
      } else {
        stop(
          paste0(
            "Don't know how to handle dataset: ",
            dt,
            "\nContact the maintainer of the package!"
          )
        )
      }

      # Update progress
      progress <- progress + 1
      utils::setTxtProgressBar(pb, progress)
    }
  }

  newFile$close_all()
  newFile <- hdf5r::H5File$new(output, mode = "r")
  result <- new("gedi.level1b", h5 = newFile)
  close(pb)
  return(result)
}
