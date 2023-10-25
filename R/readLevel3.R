#'Read GEDI Level3 data (Rasterized metrics)
#'
#'@description This function reads GEDI level3 products: rasterized metrics
#'
#'@usage readLevel3(level3path)
#'
#'@param level3path File path pointing to GEDI level3 data. Data in HDF5 Hierarchical Data Format (.h5).
#'
#'@return Returns an S4 object of class [`gedi.level3-class`] containing GEDI level3 data.
#'
#'@seealso [`hdf5r::H5File-class`] in the \emph{hdf5r} package and
#'\url{https://lpdaac.usgs.gov/products/gedi01_bv002/}
#'
#'@examples
#'# Specifying the path to GEDI level3 data (zip file)
#'outdir = tempdir()
#'level3_fp_zip <- system.file("extdata",
#'                   "GEDI01_B_2019108080338_O01964_T05337_02_003_01_sub.zip",
#'                   package="rGEDI")
#'
#'# Unzipping GEDI level3 data
#'level3path <- unzip(level3_fp_zip,exdir = outdir)
#'
#'# Reading GEDI level3 data (h5 file)
#'level3<-readLevel3(level3path=level3path)
#'
#'close(level3)
#'@import hdf5r
#'@export
readLevel3 <-function(level3path) {
  level3_rast <- terra::rast(level3path)
  level3<- new("gedi.level3", raster = level3_rast)
  return(level3)
}
