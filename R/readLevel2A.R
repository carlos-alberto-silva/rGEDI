#'Read GEDI Level2A data
#'
#'@description This function reads GEDI level2A data
#'
#'@param level2bpath file path pointing to GEDI level2B data (H5 format)
#'@return level2AFile; S4 object of class "gedi.level2a";
#'
#'@import hdf5r
#'@export
#list.datasets(level1b., recursive = T))
readLevel2A <-function(level2Apath) {
  level2A_h5 <- hdf5r::H5File$new(level2Apath, mode = 'r')
  level2A<- new("gedi.level2a", h5 = level2A_h5)
  return(level2A)
}
