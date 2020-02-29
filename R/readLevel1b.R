#'Read Geolocated Waveforms (GEDI Level1B)
#'
#'@description This function reads GEDI level1B products: geolocated Waveforms
#'
#'@usage readLevel1B(level1Bpath)
#'
#'@param level1Bpath file path pointing to GEDI level1B data. Data in HDF5 Hierarchical Data Format (.h5).
#'
#'@return S4 object of class "gedi.level1b".
#'
#'@seealso \code{\link[hdf5r]{hdf5rfile}} in the \emph{hdf5r} package and
#'https://lpdaac.usgs.gov/products/gedi01_bv001/
#'
#'@examples
#'\dontrun{
#'# specify the path to download GEDI example dataset
#'outdir<-getwd()
#'
#'# downloading GEDI example dataset (zip file)
#'download.file("https://github.com/carlos-alberto-silva/rGEDI/releases/download/examples/examples.zip",destfile=paste0(outdir,"/examples.zip"))
#'
#'# unzip the file
#'unzip(paste0(outdir,"\\examples.zip"))
#'
#'# specify the path to GEDI level1B data
#'level1bpath = paste0(outdir,"\\GEDI01_B_2019108080338_O01964_T05337_02_003_01_sub.h5"))
#'
#'# Reading GEDI level1B file
#'level1b<-readLevel1b(gedilevel1b)
#'
#'}
#'@import hdf5r
#'@export
readLevel1B <-function(level1Bpath) {
  level1b_h5 <- hdf5r::H5File$new(level1Bpath, mode = 'r')
  level1b<- new("gedi.level1b", h5 = level1b_h5)
  return(level1b)
}
