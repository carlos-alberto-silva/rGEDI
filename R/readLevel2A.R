#'Read GEDI Level2A data (Basic Full Waveform derived Metrics)
#'
#'@description This function reads GEDI level2A products: ground elevation, canopy top height, and relative heights (RH).
#'
#'
#'@usage readLevel2A(level2Apath)
#'
#'@param level2Apath File path pointing to GEDI level2A data. Data in HDF5 Hierarchical Data Format (.h5).
#'
#'@return Returns an S4 object of class "gedi.level2a" containing GEDI level2A data.
#'
#'@seealso \href{https://lpdaac.usgs.gov/products/gedi02_av002/}{https://lpdaac.usgs.gov/products/gedi02_av002/}
#'
#'@examples
#'# Specifying the path to GEDI level2A data (zip file)
#'outdir = tempdir()
#'level2A_fp_zip <- system.file("extdata",
#'                   "GEDI02_A_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'                   package="rGEDI")
#'
#'# Unzipping GEDI level2A data
#'level2Apath <- unzip(level2A_fp_zip,exdir = outdir)
#'
#'# Reading GEDI level2A data (h5 file)
#'level2a<-readLevel2A(level2Apath=level2Apath)
#'
#'close(level2a)
#'@import hdf5r
#'@export
readLevel2A <-function(level2Apath) {
  level2a_h5 <- hdf5r::H5File$new(level2Apath, mode = 'r')
  level2a<- new("gedi.level2a", h5 = level2a_h5)
  return(level2a)
}
