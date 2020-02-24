#'Read Geolocated Waveforms (GEDI Level1B)
#'
#'@description This function reads GEDI level1B products: geolocated Waveforms
#'
#'@usage readLevel1B(level1Bpath)
#'
#'@param level1Bpath file path pointing to GEDI level1B data. Data in HDF5 Hierarchical Data Format (.h5).
#'
#'@return S4 object of class "gedi.level1b".
#'@seealso \code{\link[hdf5r]{hdf5rfile}} in the \emph{hdf5r} package and
#'https://lpdaac.usgs.gov/products/gedi01_bv001/
#'
#'@examples
#'# specify the path and data file to be read
#'level1bpath <- system.file("extdata", "GEDIexample_level01B.h5", package="rGEDI")
#'
#'# read the file
#'gedilevel1b<-readLevel1B(level1bpath)
#'
#'
#'@import hdf5r
#'@export
readLevel1B <-function(level1Bpath) {
  level1b_h5 <- hdf5r::H5File$new(level1Bpath, mode = 'r')
  level1b<- new("gedi.level1b", h5 = level1b_h5)
  return(level1b)
}
