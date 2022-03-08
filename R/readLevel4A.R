#' Read GEDI Level4A data
#'
#' @description This function reads GEDI Level 4A (L4A) products: Predictions of the
#' aboveground biomass density (AGBD; in Mg/ha) and estimates of the prediction standard
#' error within each sampled geolocated laser footprint.
#'
#'
#' @usage readLevel4A(level4Apath)
#'
#' @param level2Apath File path pointing to GEDI level2A data. Data in HDF5 Hierarchical Data Format (.h5).
#'
#' @return Returns an S4 object of class [`gedi.level4a-class`] containing GEDI level4A data.
#'
#' @seealso \url{https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1907}
#'
#' @examples
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
#' level2a <- readLevel4A(level4Apath = level4Apath)
#'
#' close(level4a)
#' @import hdf5r
#' @export
readLevel4A <- function(level4Apath) {
  level4a_h5 <- hdf5r::H5File$new(level4Apath, mode = "r")
  level4a <- new("gedi.level4a", h5 = level4a_h5)
  return(level4a)
}