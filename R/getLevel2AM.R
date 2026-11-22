#'Get GEDI Elevation and Height Metrics (GEDI Level2A)
#'
#'@description This function extracts Elevation and Relative Height (RH) metrics from GEDI Level2A data.
#'
#'@param level2a A GEDI Level2A object (output of [readLevel2A()] function).
#'An S4 object of class "gedi.level2a".
#'@param beams Optional beam names. `NULL` reads every beam.
#'@param include_rh Logical. Include relative-height columns. Set to `FALSE`
#'for a faster metadata-only read.
#'@param rh Optional RH percentiles from 0 through 100. `NULL` retains all RH
#'metrics for HDF5 data and requests the full RH profile from gediDB.
#'@param quality Logical; for gediDB, retain Release 3 quality observations.
#'@param collect Logical; collect a remote footprint query into R. Set to
#'`FALSE` to retain a lazy GEE or gediDB object for rasterization.
#'@param max_features Maximum GEE footprints allowed when collecting.
#'@param chunk_size Number of GEE features transferred per request.
#'@param max_granules Maximum vector granules merged by GEE.
#'
#'@return A [data.table::data.table] containing elevation and relative-height
#'metrics, or a lazy remote query when `collect = FALSE`.
#'
#'@seealso \url{https://www.earthdata.nasa.gov/data/catalog/lpcloud-gedi02-a-002}
#'
#'@details Characteristics. Flag indicating likely invalid waveform (1=valid, 0=invalid).
#'\itemize{
#'\item \emph{beam} Beam identifier
#'\item \emph{shot_number} Shot number
#'\item \emph{degrade_flag} Flag indicating degraded state of pointing and/or positioning information
#'\item \emph{quality_flag} Flag simplifying selection of most useful data
#'\item \emph{delta_time} Transmit time of the shot since Jan 1 00:00 2018
#'\item \emph{sensitivity} Maxmimum canopy cover that can be penetrated
#'\item \emph{solar_elevation} Solar elevation
#'\item \emph{lat_lowestmode} Latitude of center of lowest mode
#'\item \emph{lon_lowestmode} Longitude of center of lowest mode
#'\item \emph{elev_highestreturn} Elevation of highest detected return relative to reference ellipsoid Meters
#'\item \emph{elev_lowestmode} Elevation of center of lowest mode relative to reference ellipsoid
#'\item \emph{rh} Relative height metrics at 1% interval
#'}
#'
#'@examples
#'
#'# Specifying the path to GEDI level2A data (zip file)
#'outdir = tempdir()
#'level2A_fp_zip <- system.file("extdata",
#'                   "GEDI02_A_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'                   package="rGEDI")
#'
#'# Unzipping GEDI level2A data
#'level2Apath <- unzip(level2A_fp_zip,exdir = outdir)
#'
#'# Reading GEDI level2A data (h5 file)
#'level2a<-readLevel2A(level2Apath=level2Apath)
#'
#'# Extracting GEDI Elevation and Height Metrics
#'level2AM<-getLevel2AM(level2a)
#'head(level2AM)
#'
#'close(level2a)
#' @param parallel Logical; process a multi-granule collection concurrently.
#' @param workers Number of parallel workers, or `NULL` for automatic selection.
#'@export
getLevel2AM<-function(level2a, beams = NULL, include_rh = TRUE, rh = NULL,
                      quality = FALSE, collect = TRUE,
                      max_features = 100000L, chunk_size = 5000L,
                      max_granules = 200L, parallel = FALSE, workers = NULL){
  if (inherits(level2a, "gedi_granule_collection")) {
    return(.extract_gedi_collection(
      level2a,
      list(beams = beams, include_rh = include_rh, rh = rh, quality = quality,
           collect = collect, max_features = max_features,
           chunk_size = chunk_size, max_granules = max_granules),
      parallel, workers
    ))
  }
  if (inherits(level2a, "gedi_gee_object")) {
    if (!is.null(beams)) warning("'beams' is ignored for a Google Earth Engine query; filter the returned table instead.")
    variables <- c(
      "beam", "shot_number", "degrade_flag", "quality_flag",
      "sensitivity", "solar_elevation", "elev_lowestmode"
    )
    if (isTRUE(include_rh)) {
      selected_rh <- if (is.null(rh)) 0:100 else unique(as.integer(rh))
      if (anyNA(selected_rh) || any(selected_rh < 0L | selected_rh > 100L)) {
        stop("'rh' must contain percentiles from 0 through 100.", call. = FALSE)
      }
      variables <- c(variables, paste0("rh", selected_rh))
    }
    if (!isTRUE(collect)) {
      return(.gee_select(level2a, variables, quality = quality, max_granules = max_granules))
    }
    return(.gee_collect(
      level2a, variables, quality = quality, max_features = max_features,
      chunk_size = chunk_size, max_granules = max_granules
    ))
  }
  if (inherits(level2a, "gedidb_level2a")) {
    if (!is.null(beams)) warning("'beams' is ignored for a gediDB query; filter the returned table instead.")
    variables <- c("beam_name", "shot_number", "degrade_flag",
                   "l2a_quality_flag_rel3", "sensitivity", "solar_elevation",
                   "elev_lowestmode")
    if (isTRUE(include_rh)) {
      if (is.null(rh)) {
        variables <- c(variables, "rh")
      } else {
        rh <- unique(as.integer(rh))
        if (anyNA(rh) || any(rh < 0L | rh > 100L)) stop("'rh' must contain percentiles from 0 through 100.", call. = FALSE)
        variables <- c(variables, paste0("rh:", rh))
      }
    }
    filters <- if (isTRUE(quality)) list(l2a_quality_flag_rel3 = "== 1") else list()
    if (!isTRUE(collect)) return(.new_gedidb_query(level2a, variables, filters))
    return(.gedidb_query(level2a, variables, filters))
  }
  level2a<-level2a@h5
  groups_id<-grep("BEAM\\d{4}$",gsub("/","",
                                     .gedi_list_groups(level2a, recursive = FALSE)), value = T)
  if (!is.null(beams)) groups_id <- intersect(groups_id, beams)
  if (!length(groups_id)) return(data.table::data.table())
  rh.dt<-data.table::data.table()
  pb <- utils::txtProgressBar(min = 0, max = length(groups_id), style = 3)
  i.s=0

  for ( i in groups_id){
    i.s<-i.s+1
    utils::setTxtProgressBar(pb, i.s)
    level2a_i<-level2a[[i]]

    if (any(.gedi_list_datasets(level2a_i)=="shot_number")){

    available <- .gedi_list_datasets(level2a_i)
    quality_name <- intersect(c("l2a_quality_flag_rel3", "quality_flag",
                                "l2a_quality_flag_rel2"), available)[1L]
    values <- list(
      beam=rep(i,length(level2a_i[["shot_number"]][])),
      shot_number=level2a_i[["shot_number"]][],
      degrade_flag=level2a_i[["degrade_flag"]][],
      delta_time=level2a_i[["delta_time"]][],
      sensitivity=level2a_i[["sensitivity"]][],
      solar_elevation=level2a_i[["solar_elevation"]][],
      lat_lowestmode=level2a_i[["lat_lowestmode"]][],
      lon_lowestmode=level2a_i[["lon_lowestmode"]][],
      elev_highestreturn=level2a_i[["elev_highestreturn"]][],
      elev_lowestmode=level2a_i[["elev_lowestmode"]][])
    if (!is.na(quality_name)) values[[quality_name]] <- level2a_i[[quality_name]][]
    rhs <- data.table::as.data.table(values)
    if (isTRUE(include_rh)) {
      if(length(level2a_i[["rh"]]$dims)==2) {
        rh_values=t(level2a_i[["rh"]][,])
      } else {
        rh_values=t(level2a_i[["rh"]][])
      }
      rhs <- cbind(rhs, rh_values)
    }
    rh.dt<-rbind(rh.dt,rhs)
    }
  }

  if (isTRUE(include_rh) && ncol(rh.dt) >= 101L) {
    tail_names <- tail(names(rh.dt), 101L)
    data.table::setnames(rh.dt, tail_names, paste0("rh",seq(0,100)))
  }
  if (isTRUE(include_rh) && !is.null(rh)) {
    rh <- unique(as.integer(rh))
    if (anyNA(rh) || any(rh < 0L | rh > 100L)) stop("'rh' must contain percentiles from 0 through 100.", call. = FALSE)
    keep_rh <- paste0("rh", rh)
    rh.dt <- rh.dt[, c(setdiff(names(rh.dt), grep("^rh[0-9]+$", names(rh.dt), value = TRUE)), keep_rh), with = FALSE]
  }
  close(pb)
  return(rh.dt)
}


