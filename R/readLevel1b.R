#'Read LVIS Level1b data
#'
#'@description This function reads LVIS level1b data
#'
#'@usage readLevel1b(level1bpath)
#'
#'@param level1bpath file path pointing to LVIS level1b data (H5 format)
#'@return H5File; S4 object of class H5File;
#'@author Carlos Alberto Silva. This function calls \emph{h5file} function from h5 package (Author: Mario Annau)
#'@seealso \code{\link[h5]{h5file}} in the \emph{h5} package.
#'@examples
#'
#'# LVIS level1b file path
#'level1b_filepath_zip <- system.file("extdata", "LVIS_Mondah_level1b.zip", package="rLVIS")
#'unzip(level1b_filepath_zip, exdir = tempdir())
#'level1b_filepath <- file.path(tempdir(), "LVIS_Mondah_level1b.h5")
#'
#'# Reading LVIS level1b file
#'level1b<-readLevel1b(level1bpath=level1b_filepath)
#'
#'
#'@importFrom hdf5r H5File
#'@importFrom data.table data.table
#'@export

#list.datasets(level1b., recursive = T))

readLevel1b<-function(level1bpath) {
  level1b_h5<- hdf5r::H5File$new(level1bpath, mode = 'r')
  level1b<- new("gedi.level1b", h5 = level1b_h5)
  return(level1b_h5)
}
class(level1b_h5)
