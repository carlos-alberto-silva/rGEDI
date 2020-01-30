#'Read GEDI Level2B data
#'
#'@description This function reads GEDI level2B data
#'
#'@param level2bpath file path pointing to GEDI level2B data (H5 format)
#'@return level2bFile; S4 object of class "gedi.level2b";
#'
#'@import hdf5r
#'@export
#list.datasets(level1b., recursive = T))
readLevel2B <-function(level2Bpath) {
  level2B_h5 <- hdf5r::H5File$new(level2Bpath, mode = 'r')
  level2B<- new("gedi.level2b", h5 = level2B_h5)
  return(level2B)
}
