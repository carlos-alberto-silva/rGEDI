#'Read GEDI Level2B data (Biophysical Variables)
#'
#'@description This function reads GEDI level2B products: canopy cover, Plant Area Index (PAI), Plant Area Volume Density (PAVD), and Foliage Height Diversity (FHD).
#'
#'@usage readLevel2B(level2Bpath)
#'
#'@param level2Bpath One or more local paths or Earthdata Cloud URLs pointing
#'to GEDI Level 2B HDF5 granules, or a gediDB reference returned by [gedifinder()].
#'
#'@return A [`gedi.level2b-class`] object for one HDF5 granule, a lazy remote
#'object, or a multi-granule collection consumed by [getLevel2BVPM()].
#'
#'
#'@examples
#'# Specifying the path to GEDI level2B data (zip file)
#'outdir = tempdir()
#'level2B_fp_zip <- system.file("extdata",
#'                   "GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'                   package="rGEDI")
#'
#'# Unzipping GEDI level2A data
#'level2Bpath <- unzip(level2B_fp_zip,exdir = outdir)
#'
#'# Reading GEDI level2B data (h5 file)
#'level2b<-readLevel2B(level2Bpath=level2Bpath)
#'
#'close(level2b)
#'@import hdf5r
#'@export
readLevel2B <-function(level2Bpath) {
  if (is.character(level2Bpath) && length(level2Bpath) > 1L) {
    return(.new_gedi_granule_collection(level2Bpath, "GEDI02_B"))
  }
  if (inherits(level2Bpath, "gedi_gee_reference")) {
    return(.open_gee_reference(level2Bpath, "GEDI02_B"))
  }
  if (inherits(level2Bpath, "gedidb_reference")) {
    return(.open_gedidb_reference(level2Bpath, "GEDI02_B"))
  }
  level2b_h5 <- if (inherits(level2Bpath, "GEDICloudH5")) level2Bpath else if (.is_gedi_url(level2Bpath)) {
    .open_cloud_h5(level2Bpath, product = "GEDI02_B")
  } else {
    hdf5r::H5File$new(level2Bpath, mode = 'r')
  }
  level2b<-new("gedi.level2b", h5 = level2b_h5)
  return(level2b)
}
