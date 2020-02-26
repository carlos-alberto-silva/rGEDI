#'Read GEDI Basic Waveform Metrics (Level2A data)
#'
#'@description This function reads GEDI level2A products: ground elevation, canopy top height, and relative heights (RH).
#'
#'
#'@usage readLevel2A(level2Apath)
#'
#'@param level2Apath file path pointing to GEDI level2A data. Data in HDF5 Hierarchical Data Format (.h5).
#'
#'@return S4 object of class "gedi.level1a".
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi02_av001/
#'
#'#'@examples
#'# specify the path and data file to be read
#'level2apath <- system.file("extdata", "GEDIexample_level02A.h5", package="rGEDI")
#'
#'# read the file
#'gedilevel2a<-readLevel2A(level2apath)
#'
#'@import hdf5r
#'@export
readLevel2A <-function(level2Apath) {
  level2a_h5 <- hdf5r::H5File$new(level2Apath, mode = 'r')
  level2a<- new("gedi.level2a", h5 = level2a_h5)
  return(level2a)
}