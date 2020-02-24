#'Read GEDI Biophysical Metrics (Level2B data)
#'
#'@description This function reads GEDI level2B products: canopy cover, Plant Area Index (PAI), Plant Area Volume Density (PAVD), and Foliage Height Diversity (FHD).
#'
#'#'
#'@usage readLevel2B(level2Bpath)
#'
#'@param level2Bpath file path pointing to GEDI level2B data. Data in HDF5 Hierarchical Data Format (.h5).
#'
#'@return S4 object of class "gedi.level1b".
#'
#'
#'#'@examples
#'# specify the path and data file to be read
#'level2bpath <- system.file("extdata", "GEDIexample_level02B.h5", package="rGEDI")
#'
#'# read the file
#'gedilevel2b<-readLevel2B(level2bpath)
#'
#'@import hdf5r
#'@export
readLevel2B <-function(level2Bpath) {
  level2b_h5 <- hdf5r::H5File$new(level2Bpath, mode = 'r')
  level2b<-new("gedi.level2b", h5 = level2b_h5)
  return(level2b)
}
