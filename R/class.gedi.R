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

#' Class for GEDI Full-Waveform Simulation
#'
#' @slot h5 Object of class H5File from hdf5r package
#'
#' @import methods
#' @export
gedi.level1bSim <- setClass(
  Class="gedi.level1bSim",
  slots = list(h5 = "H5File")
)

#' Class for GEDI level1B Full Waveform
#'
#' @slot dt Object of class data.table from data.table package
#'
#' @import methods
#' @export
gedi.fullwaveform <- setClass(
  Class="gedi.fullwaveform",
  slots = list(dt = "data.table")
)


#'Plot GEDI* object
#'
#'@param x An object of class "gedi.fullwaveform". (output of \code{\link[rGEDI:getLevel1BWF]{getLevel1BWF}} function)
#'@param y not used (inherited from R base)
#'@param ... will be passed to the main plot
#'
#'@param relative if TRUE, the Wavform Amplitude will be showed in percentage (\%)
#'@param polygon if TRUE, polygon will be added to the plot
#'
#'@param method methods used for simulating the GEDI full-waveform ("RXWAVEINT","RXWAVEINT" or "RXWAVEINT"). Default is "RXWAVECOUNT".
#'
#' @export
#' @method plot gedi.fullwaveform
setGeneric("plot", function(x, y, ...)
  standardGeneric("plot"))

#'@description for gedi.fullwaveform: will plot the full waveform\cr\cr
#'@examples
#'\donttest{
#'# specify the path to download GEDI example dataset
#'outdir<-getwd()
#'
#'# downloading GEDI example dataset (zip file)
#'download.file(
#'              paste0(
#'                     "https://github.com/carlos-alberto-silva/rGEDI/",
#'                     "releases/download/examples/examples.zip"
#'              ),
#'              destfile=file.path(outdir, "examples.zip"))
#'
#'# unzip the file
#'unzip(file.path(outdir, "examples.zip"))
#'
#'# specify the path to GEDI level1B data
#'level1bpath = file.path(outdir, "GEDI01'_B_2019108080338_O01964_T05337_02_003_01_sub.h5")
#'
#'# Reading GEDI level1B file
#'level1b<-readLevel1b(gedilevel1b)
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
#'}
#'
#'
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


#'@description for gedi.level1bSim: will plot the simulated waveform
#'
#'@examples
# specify the path to ALS data
#'lasfile_amazon <- system.file("extdata", "Amazon.las", package="rGEDI")
#'lasfile_cerrado <- system.file("extdata", "Cerrado.las", package="rGEDI")
#'
#'# Reading and plot ALS file
#'libsAvailable = require(lidR) && require(plot3D)
#'if (libsAvailable) {
#'las_amazon<-readLAS(lasfile_amazon)
#'las_cerrado<-readLAS(lasfile_cerrado)
#'
#'# Extracting plot center geolocations
#'xcenter_amazon = mean(las_amazon@bbox[1,])
#'ycenter_amazon = mean(las_amazon@bbox[2,])
#'xcenter_cerrado = mean(las_cerrado@bbox[1,])
#'ycenter_cerrado = mean(las_cerrado@bbox[2,])
#'
#'# Simulating GEDI full-waveform
#'wf_amazon<-gediWFSimulator(
#'                           input=lasfile_amazon,
#'                           output=paste0(
#'                                         getwd(),
#'                                         "//gediWF_amazon_simulation.h5"
#'                                         ),
#'                           coords = c(xcenter_amazon, ycenter_amazon))
#' wf_cerrado<-gediWFSimulator(
#'                             input=lasfile_cerrado,
#'                             output=paste0(
#'                                           getwd(),
#'                                           "//gediWF_cerrado_simulation.h5"
#'                                           ),
#'                             coords = c(xcenter_cerrado, ycenter_cerrado))
#'# Plot Full-waveform
#'par(mfrow=c(2,2), mar=c(4,4,0,0), oma=c(0,0,1,1),cex.axis = 1.2)
#'scatter3D(
#'          las_amazon@data$X,
#'          las_amazon@data$Y,
#'          las_amazon@data$Z,
#'          pch = 16, colkey = FALSE, main="",
#'          cex = 0.5, bty = "u", col.panel ="gray90",
#'          phi = 30, alpha=1, theta=45, col.grid = "gray50",
#'          xlab="UTM Easting (m)", ylab="UTM Northing (m)", zlab="Elevation (m)"
#'          )
#'
#'plot(wf_amazon, relative=TRUE, polygon=TRUE, type="l", lwd=2, col="forestgreen",
#'     xlab="", ylab="Elevation (m)", ylim=c(90,140))
#'grid()
#'scatter3D(
#'          las_cerrado@data$X,las_cerrado@data$Y,las_cerrado@data$Z,
#'          pch = 16,colkey = FALSE, main="",
#'          cex = 0.5,bty = "u",col.panel ="gray90",
#'          phi = 30,alpha=1,theta=45,col.grid = "gray50",
#'          xlab="UTM Easting (m)", ylab="UTM Northing (m)", zlab="Elevation (m)"
#'          )
#'
#'plot(wf_cerrado, relative=TRUE, polygon=TRUE, type="l", lwd=2, col="green",
#'     xlab="Waveform Amplitude (%)", ylab="Elevation (m)", ylim=c(815,835))
#'grid()
#'}
#' @rdname plot
setMethod("plot", signature("gedi.level1bSim", y = "missing"), function(x,relative=FALSE,polygon=FALSE,method="RXWAVEINT",...) {

  if (!class(x)=="gedi.level1bSim"){

    print("Invalid input file. It should be an object of class 'gedi.fullwaveform' ")
  } else {

    wfh5<-x@h5
    z = seq(wfh5[["Z0"]][],wfh5[["ZN"]][], length.out = wfh5[["NBINS"]][])

    if  (method=="RXWAVEINT") { x0<-wfh5[["RXWAVEINT"]][,]}
    if  (method=="RXWAVECOUNT") { x0<-wfh5[["RXWAVECOUNT"]][,]}
    if  (method=="RXWAVEFRAC") { x0<-wfh5[["RXWAVEFRAC"]][,]}

    x<-x0

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

h5closeall = function(con, ...) {
  try(con@h5$close_all(), silent=TRUE)
}


#'Close hdf5 connections from gedi* objects
#'
#'@param con An object of class gedi*
#'@param ... Inherited from base
#'
#' @export
#' @rdname close
#' @method close gedi.level1b
setGeneric("close", function(con, ...)
  standardGeneric("close"))

#'@rdname close
setMethod("close", signature = c("gedi.level1b"), h5closeall)
#'@rdname close
setMethod("close", signature = c("gedi.level2a"), h5closeall)
#'@rdname close
setMethod("close", signature = c("gedi.level2b"), h5closeall)
#'@rdname close
setMethod("close", signature = c("gedi.level1bSim"), h5closeall)
