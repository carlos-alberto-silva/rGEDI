#'Read GEDI Level2A data (Basic Full Waveform derived Metrics)
#'
#'@description This function reads GEDI level2A products: ground elevation, canopy top height, and relative heights (RH).
#'
#'
#'@usage readLevel2A(level2Apath)
#'
#'@param level2Apath One or more local paths or Earthdata Cloud URLs pointing
#'to GEDI Level 2A HDF5 granules, or a gediDB reference returned by
#'[gedifinder()] with `source = "gedidb"`.
#'
#'@return An S4 object of class [`gedi.level2a-class`] for HDF5 input, or a
#'lazy `gedidb_level2a` object for a gediDB reference. Multiple paths return a
#'lazy collection consumed by [getLevel2AM()].
#'
#'@seealso \url{https://www.earthdata.nasa.gov/data/catalog/lpcloud-gedi02-a-002}
#'
#'@examples
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
#'close(level2a)
#'@import hdf5r
#'@export
readLevel2A <-function(level2Apath) {
  if (is.character(level2Apath) && length(level2Apath) > 1L) {
    return(.new_gedi_granule_collection(level2Apath, "GEDI02_A"))
  }
  if (inherits(level2Apath, "gedi_gee_reference")) {
    return(.open_gee_reference(level2Apath, "GEDI02_A"))
  }
  if (inherits(level2Apath, "gedidb_reference")) {
    return(.open_gedidb_reference(level2Apath, "GEDI02_A"))
  }
  level2a_h5 <- if (inherits(level2Apath, "GEDICloudH5")) level2Apath else if (.is_gedi_url(level2Apath)) {
    .open_cloud_h5(level2Apath, product = "GEDI02_A")
  } else {
    hdf5r::H5File$new(level2Apath, mode = 'r')
  }
  level2a<- new("gedi.level2a", h5 = level2a_h5)
  return(level2a)
}
