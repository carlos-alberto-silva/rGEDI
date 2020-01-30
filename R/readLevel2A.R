#'Read GEDI Level1A data
#'
#'@description This function reads LVIS level1b data
#'
#'@param level1bpath file path pointing to LVIS level1b data (H5 format)
#'@return H5File; S4 object of class H5File;
#'@author Carlos Alberto Silva. This function calls \emph{h5file} function from h5 package (Author: Mario Annau)
#'@seealso \code{\link[h5]{h5file}} in the \emph{h5} package.
#'
#'@importFrom data.table data.table
#'@import hdf5r
#'@export
#list.datasets(level1b., recursive = T))
readLevel2A <-function(level2Apath) {
  level2A_h5 <- hdf5r::H5File$new(level2Apath, mode = 'r')
  level2A<- new("gedi.level2a", h5 = level2A_h5)
  return(level2A)
}
