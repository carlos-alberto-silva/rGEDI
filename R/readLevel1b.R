#'Read GEDI Level1B data
#'
#'@description This function reads GEDI level1B data
#'
#'@param level1Bpath file path pointing to GEDI level1B data (H5 format)
#'@return level1BFile; S4 object of class "gedi.level1b";
#'
#'@import hdf5r
#'@export
readLevel1B <-function(level1Bpath) {
  level1B_h5 <- hdf5r::H5File$new(level1Bpath, mode = 'r')
  level2<- new("gedi.level1b", h5 = level1B_h5)
  return(level1B)
}
