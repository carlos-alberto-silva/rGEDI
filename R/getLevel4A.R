`:=` <- data.table::`:=`
l4a.var.map <- list()
l4a.var.map[["beam"]] <- "beam"
l4a.var.map[["shot_number"]] <- "shot_number"
l4a.var.map[["algorithm_run_flag"]] <- "algorithm_run_flag"
l4a.var.map[["l4_quality_flag"]] <- "l4_quality_flag"
l4a.var.map[["delta_time"]] <- "delta_time"
l4a.var.map[["sensitivity"]] <- "sensitivity"
l4a.var.map[["solar_elevation"]] <- "solar_elevation"
l4a.var.map[["lat_lowestmode"]] <- "lat_lowestmode"
l4a.var.map[["lon_lowestmode"]] <- "lon_lowestmode"
l4a.var.map[["elev_lowestmode"]] <- "elev_lowestmode"
l4a.var.map[["agbd"]] <- "agbd"
l4a.var.map[["agbd_pi_upper"]] <- "agbd_pi_upper"
l4a.var.map[["agbd_pi_lower"]] <- "agbd_pi_lower"
l4a.var.map[["agbd_se"]] <- "agbd_se"


#' Get GEDI-derived Aboveground biomass density (GEDI Level4A)
#'
#' @description This function extracts the aboveground biomass density (AGBD; in Mg/ha) and
#' standard error estimates within each sampled geolocated laser footprint.
#'
#' @usage getLevel4A(level4a)
#'
#' @param level4a A GEDI Level4A object (output of [readLevel4A()] function).
#' An S4 object of class "gedi.level4a".
#'
#' @return Returns an S4 object of class [data.table::data.table]
#' containing the aboveground biomass density (AGBD; in Mg/ha) and
#' standard error estimates.
#'
#' @seealso \url{https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1907}
#'
#' @details Characteristics. Flag indicating likely invalid waveform (1=valid, 0=invalid).
#' \itemize{
#' \item \emph{beam} Beam identifier
#' \item \emph{shot_number} Shot number
#' \item \emph{algorithm_run_flag} The L4A algorithm is run if this flag is set to 1.
#' This flag selects data that have sufficient waveform fidelity for AGBD estimation.
#' \item \emph{l4_quality_flag} Flag simplifying selection of most useful biomass predictions
#' \item \emph{delta_time} Transmit time of the shot since Jan 1 00:00 2018
#' \item \emph{sensitivity} Maxmimum canopy cover that can be penetrated
#' \item \emph{solar_elevation} Solar elevation
#' \item \emph{lat_lowestmode} Latitude of center of lowest mode
#' \item \emph{lon_lowestmode} Longitude of center of lowest mode
#' \item \emph{elev_lowestmode} Elevation of center of lowest mode relative to reference ellipsoid
#' \item \emph{agbd} Aboveground biomass density (Mg/ha)
#' \item \emph{agbd_pi_upper} Lower prediction interval (see alpha attribute for the level)
#' \item \emph{agbd_pi_lower} Upper prediction interval (see alpha attribute for the level)
#' \item \emph{agbd_se} Aboveground biomass density (Mg/ha) prediction standard error
#' }
#'
#'
#' @examples
#'
#' # Specifying the path to GEDI level4A data (zip file)
#' outdir <- tempdir()
#' level4A_fp_zip <- system.file("extdata",
#'   "GEDI04_A_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'   package = "rGEDI"
#' )
#'
#' # Unzipping GEDI level4A data
#' level4Apath <- unzip(level4A_fp_zip, exdir = outdir)
#'
#' # Reading GEDI level4A data (h5 file)
#' level4a <- readLevel2A(level4Apath = level4Apath)
#'
#' # Extracting GEDI AGBD and standard error estimates.
#' level4AM <- getLevel4A(level4a)
#' head(level4AM)
#'
#' close(level4a)
#' @export
getLevel4A <- function(level4a, cols = c(
                         "beam",
                         "shot_number",
                         "algorithm_run_flag",
                         "l4_quality_flag",
                         "delta_time",
                         "sensitivity",
                         "solar_elevation",
                         "lat_lowestmode",
                         "lon_lowestmode",
                         "elev_lowestmode",
                         "agbd",
                         "agbd_pi_upper",
                         "agbd_pi_lower",
                         "agbd_se"
                       )) {
  level4a <- level4a@h5
  groups_id <- grep("BEAM\\d{4}$", gsub(
    "/", "",
    hdf5r::list.groups(level4a, recursive = F)
  ), value = T)
  m.dt <- data.table::data.table()
  pb <- utils::txtProgressBar(min = 0, max = length(groups_id), style = 3)
  i.s <- 0

  # i = groups_id[1]
  for (i in groups_id) {
    i.s <- i.s + 1
    utils::setTxtProgressBar(pb, i.s)
    level4a_i <- level4a[[i]]
    m <- data.table::data.table()
    # col = cols[1]
    for (col in cols) {
      h5.address <- l4a.var.map[[col]]
      if (is.null(h5.address)) {
        if (i.s == 1) warning(sprintf("The column '%s' is not available in the GEDI Level 4A product!", col))
        next
      }
      m[, eval(col) := level4a_i[[h5.address]][]]
    }
    m.dt <- rbind(m.dt, m)
  }
  close(pb)
  return(m.dt)
}