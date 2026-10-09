#'Get GEDI Elevation and Height Metrics (GEDI Level2A)
#'
#'@description This function extracts Elevation and Relative Height (RH) metrics from GEDI Level2A data.
#'
#'@param level2a A GEDI Level2A object (output of [readLevel2A()] function).
#'An S4 object of class "gedi.level2a".
#'@param beams Optional beam names. `NULL` reads every beam.
#'@param include_rh Logical. Include the 101 relative-height columns. Set to
#'`FALSE` for a faster metadata-only cloud read.
#'
#'@return Returns an S4 object of class [data.table::data.table]
#'containing the elevation and relative heights metrics.
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
#'@export
getLevel2AM<-function(level2a, beams = NULL, include_rh = TRUE){
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
        rh=t(level2a_i[["rh"]][,])
      } else {
        rh=t(level2a_i[["rh"]][])
      }
      rhs <- cbind(rhs, rh)
    }
    rh.dt<-rbind(rh.dt,rhs)
    }
  }

  if (isTRUE(include_rh) && ncol(rh.dt) >= 101L) {
    tail_names <- tail(names(rh.dt), 101L)
    data.table::setnames(rh.dt, tail_names, paste0("rh",seq(0,100)))
  }
  close(pb)
  return(rh.dt)
}


