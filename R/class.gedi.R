#setClass("gedi.level1b", representation(h5="H5File",level1b.spdf='SpatialPointsDataFrame'))
#' @importFrom hdf5r H5File
setRefClass("H5File")
requireNamespace("data.table")

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


#'Plot GEDI full-waveform
#'
#'Plots a single GEDI full-waveform (level1b)
#'
#'@param x An object of class "gedi.fullwaveform". (output of \code{\link[rGEDI:getLevel1BWF]{getLevel1BWF}} function)
#'@param relative if TRUE, the Wavform Amplitude will be showed in percentage (%)
#'@param polygon if TRUE, polygon will be added to the plot
#'
#'@examples
#'#'# specify the path to GEDI Level 1B data
#'level1bpath <- system.file("extdata", "GEDIexample_level01B.h5", package="rGEDI")
#'
#'# Reading GEDI level1B data
#'level1b <- readLevel1B(level1bpath)
#'
#'# extract the desired information into a dataframe
#'wf <- getLevel1BWF(level1b, shot_number="19850022900500000")
#'
#'# Plot Full-waveform
#'par(mfrow = c(1,2), cex.axis = 1.5)
#'plot(wf, relative=FALSE, polygon=TRUE, type="l", lwd=2, col="forestgreen",
#'xlab="", ylab="Elevation (m)")
#'
#'plot(wf, relative=TRUE, polygon=TRUE, type="l", lwd=2, col="forestgreen",
#'xlab="Waveform Amplitude (%)", ylab="Elevation (m)")
#' @export
#' @method plot gedi.fullwaveform
setGeneric("plot", function(x, y, ...)
  standardGeneric("plot"))

#' @export
#' @rdname plot
setMethod("plot", signature("gedi.fullwaveform", y = "missing"), function(x,relative=FALSE,polygon=FALSE,...) {

  if (!class(x)=="gedi.fullwaveform"){

    print("Invalid input file. It should be an object of class 'gedi.fullwaveform' ")
  } else {



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

      suppressWarnings({plot(xl,yl,...)})
      suppressWarnings({polygon(xl,yl,...)})
    } else {
      suppressWarnings({plot(x=x,y=z,...)})
    }
  }
})

