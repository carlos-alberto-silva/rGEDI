#'Read GEDI Level2B data (Biophysical Variables)
#'
#'@description This function reads GEDI level2B products: canopy cover, Plant Area Index (PAI), Plant Area Volume Density (PAVD), and Foliage Height Diversity (FHD).
#'
#'@usage readLevel2B(level2Bpath)
#'
#'@param level2Bpath File path pointing to GEDI level2B data. Data in HDF5 Hierarchical Data Format (.h5).
#'
#'@return Returns an S4 object of class "gedi.level2b" containing GEDI level2B data.
#'
#'
#'@examples
#'# Specifying the path to GEDI level2B data (zip file)
#'level2B_fp_zip <- system.file("extdata",
#'                   "GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'                   package="rGEDI")
#'
#'# Unzipping GEDI level2A data
#'level2Bpath <- unzip(level2B_fp_zip,exdir = dirname(level2B_fp_zip))
#'
#'# Reading GEDI level2B data (h5 file)
#'level2b<-readLevel2B(level2Bpath=level2Bpath)
#'
#'close(level2b)
#'@import hdf5r
#'@export
readLevel2B <-function(level2Bpath) {
  level2b_h5 <- hdf5r::H5File$new(level2Bpath, mode = 'r')
  level2b<-new("gedi.level2b", h5 = level2b_h5)
  return(level2b)
}
