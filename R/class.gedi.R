#setClass("gedi.level1b", representation(h5="H5File",level1b.spdf='SpatialPointsDataFrame'))
#' @importFrom hdf5r H5File
setRefClass("H5File")

#' Class for GEDI level1B
#'
#' @slot h5 Object of class H5File from hdf5r package
#'
#' @import methods
#' @export
gedi.level1b <- setClass(
  Class="gedi.level1b",
  slots = list(h5 = "H5File")
)

#' Class for GEDI level2A
#'
#' @slot h5 Object of class H5File from hdf5r package
#'
#' @import methods
#' @export
gedi.level2a <- setClass(
  Class="gedi.level2a",
  slots = list(h5 = "H5File")
)

#' Class for GEDI level2B
#'
#' @slot h5 Object of class H5File from hdf5r package
#'
#' @import methods
#' @export
gedi.level2b <- setClass(
  Class="gedi.level2b",
  slots = list(h5 = "H5File")
)

setMethod("plot", signature("gedi.level1b", y = "missing"), function(x,shot_number,relative=TRUE,polygon=FALSE,...) {
  PlotWaveform(x,shot_number=shot_number,relative=relative,polygon=polygon,return=FALSE,plotWave=TRUE,...)
})

