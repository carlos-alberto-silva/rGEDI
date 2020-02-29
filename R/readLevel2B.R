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
#'# specify the path to GEDI level2B data
#'level2bpath = paste0(outdir,"\\GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub.h5"))
#'
#'# Reading GEDI level1B file
#'level2b<-readLevel2b(gedilevel2b)
#'
#'}
#'@import hdf5r
#'@export
readLevel2B <-function(level2Bpath) {
  level2b_h5 <- hdf5r::H5File$new(level2Bpath, mode = 'r')
  level2b<-new("gedi.level2b", h5 = level2b_h5)
  return(level2b)
}
