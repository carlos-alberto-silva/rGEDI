#'Read GEDI Level1B data (Geolocated Waveforms)
#'
#'@description This function reads GEDI level1B products: geolocated Waveforms
#'
#'@usage readLevel1B(level1Bpath)
#'
#'@param level1Bpath Local file path or Earthdata Cloud URL pointing to a
#'GEDI Level 1B HDF5 granule.
#'
#'@return Returns an S4 object of class [`gedi.level1b-class`] containing GEDI level1B data.
#'
#'@seealso [`hdf5r::H5File-class`] in the \emph{hdf5r} package and
#'\url{https://www.earthdata.nasa.gov/data/catalog/lpcloud-gedi01-b-002}
#'
#'@examples
#'# Specifying the path to GEDI level1B data (zip file)
#'outdir = tempdir()
#'level1B_fp_zip <- system.file("extdata",
#'                   "GEDI01_B_2019108080338_O01964_T05337_02_003_01_sub.zip",
#'                   package="rGEDI")
#'
#'# Unzipping GEDI level1B data
#'level1Bpath <- unzip(level1B_fp_zip,exdir = outdir)
#'
#'# Reading GEDI level1B data (h5 file)
#'level1b<-readLevel1B(level1Bpath=level1Bpath)
#'
#'close(level1b)
#'@import hdf5r
#'@export
readLevel1B <-function(level1Bpath) {
  level1b_h5 <- if (inherits(level1Bpath, "GEDICloudH5")) level1Bpath else if (.is_gedi_url(level1Bpath)) {
    .open_cloud_h5(level1Bpath, product = "GEDI01_B")
  } else {
    hdf5r::H5File$new(level1Bpath, mode = 'r')
  }
  level1b<- new("gedi.level1b", h5 = level1b_h5)
  return(level1b)
}
