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

#' Class for GEDI level1B Full Waveform
#'
#' @slot dt Object of class data.table from data.table package
#'
#' @import methods
#' @export
gedi.level2b <- setClass(
  Class="gedi.fullwaveform",
  slots = list(dt = "data.table")
)

setMethod("plot", signature("gedi.fullwaveform", y = "missing"), function(x,relative=FALSE,polygon=FALSE,...) {

  x0<-as.data.frame(x@dt)
  x<-x0[,1]
  z<-x0[,2]

  if (relative==TRUE){
    x=c(x-min(x))/(max(x)-min(x))*100
  } else{
    x=x
  }

  if (polygon==TRUE){

      xstart<-x[which(z==min(z, na.rm=T))]
      xend<-x[which(z==max(z, na.rm=T))]

      xl<-c(min(x),min(x),xstart,rev(x),xend,min(x))
      yl<-c(max(z, na.rm=T),min(z, na.rm=T),min(z, na.rm=T),rev(z),max(z, na.rm=T),max(z, na.rm=T))

      plot(xl,yl,...)
      suppressWarnings(polygon(xl,yl,...))
    } else {
      plot(x=x,y=z)
    }

})

