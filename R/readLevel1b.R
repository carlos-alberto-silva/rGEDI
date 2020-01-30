#'Read LVIS Level1b data
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
readLevel1B <-function(level1bpath) {
  level1b_h5 <- hdf5r::H5File$new(level1bpath, mode = 'r')
  level1b<- new("gedi.level1b", h5 = level1b_h5)
  return(level1b)
}
