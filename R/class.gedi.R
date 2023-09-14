# setClass("gedi.level1b", representation(h5="H5File",level1b.spdf='SpatialPointsDataFrame'))
#' @importFrom hdf5r H5File
setRefClass("H5File")
requireNamespace("data.table")

#' Class for GEDI level1B
#'
#' @slot h5 Object of class [`H5File`][hdf5r::H5File-class] from `hdf5r` package containing the
#' GEDI level1B products: geolocated Waveforms
#'
#' @seealso [`H5File`][hdf5r::H5File-class] in the `hdf5r` package and
#' \url{https://lpdaac.usgs.gov/products/gedi01_bv002/}
#'
#' @import methods
#' @export
gedi.level1b <- setClass(
  Class = "gedi.level1b",
  slots = list(h5 = "H5File")
)

#' Class for GEDI level2A
#'
#' @slot h5 Object of class H5File from `hdf5r` package containing the
#' GEDI level2A products: ground elevation, canopy top height, and relative heights (RH).
#'
#' @seealso [`H5File`][hdf5r::H5File-class] in the `hdf5r` package and
#' \url{https://lpdaac.usgs.gov/products/gedi02_av002/}
#'
#' @import methods
#' @export
gedi.level2a <- setClass(
  Class = "gedi.level2a",
  slots = list(h5 = "H5File")
)

#' Class for GEDI level2B
#'
#' @slot h5 Object of class [`H5File`][hdf5r::H5File-class] from `hdf5r` package containing the
#' GEDI level2B products: canopy cover, Plant Area Index (PAI), Plant Area Volume Density (PAVD),
#' and Foliage Height Diversity (FHD).
#'
#' @seealso [`H5File`][hdf5r::H5File-class] in the `hdf5r` package and
#' \url{https://lpdaac.usgs.gov/products/gedi02_bv002/}
#'
#' @import methods
#' @export
gedi.level2b <- setClass(
  Class = "gedi.level2b",
  slots = list(h5 = "H5File")
)


#' Class for GEDI level1B Full Waveform
#'
#' @slot dt Object of class data.table from \emph{data.table} package containing
#' the extracted GEDI full-waveform elevation and amplitude.
#'
#' @import methods
#' @export
gedi.fullwaveform <- setClass(
  Class = "gedi.fullwaveform",
  slots = list(dt = "data.table")
)


#' Plot GEDI* object
#'
#' @param x An object of class [`gedi.fullwaveform-class`] (output of [getLevel1BWF()] function)
#' @param y not used (inherited from R base)
#'
#' @param relative if TRUE, the Waveform Amplitude will be showed in percentage (%)
#' @param polygon if TRUE, the polygon will be added to the plot
#'
#' @param method methods used for simulating the GEDI full-waveform ("RXWAVEINT", "RXWAVECOUNT" or "RXWAVEFRAC"). Default is "RXWAVECOUNT".
#' @param ... will be passed to the main plot
#' @return No return value
#'
#' @export
#' @method plot gedi.fullwaveform
setGeneric("plot", function(x, y, ...) {
  standardGeneric("plot")
})

#' @description For [`gedi.fullwaveform-class`]: will plot the full waveform
#' @examples
#' # Specifying the path to GEDI level1B data (zip file)
#' outdir <- tempdir()
#' level1B_fp_zip <- system.file("extdata",
#'   "GEDI01_B_2019108080338_O01964_T05337_02_003_01_sub.zip",
#'   package = "rGEDI"
#' )
#'
#' # Unzipping GEDI level1B data
#' level1Bpath <- unzip(level1B_fp_zip, exdir = outdir)
#'
#' # Reading GEDI level1B data (h5 file)
#' level1b <- readLevel1B(level1Bpath = level1Bpath)
#'
#' # Extracting GEDI Full-Waveform
#' wf <- getLevel1BWF(level1b, shot_number = "19640521100108408")
#'
#' # Plotting GEDI Full-waveform
#' oldpar <- par()
#' par(mfrow = c(1, 2), cex.axis = 1.5)
#' plot(wf,
#'   relative = FALSE, polygon = TRUE, type = "l", lwd = 2, col = "forestgreen",
#'   xlab = "", ylab = "Elevation (m)"
#' )
#'
#' plot(wf,
#'   relative = TRUE, polygon = TRUE, type = "l", lwd = 2, col = "forestgreen",
#'   xlab = "Waveform Amplitude (%)", ylab = "Elevation (m)"
#' )
#'
#' par(oldpar)
#' close(level1b)
#' @rdname plot
setMethod("plot", signature("gedi.fullwaveform", y = "missing"), function(x, relative = FALSE, polygon = FALSE, ...) {
  if (!is(x, "gedi.fullwaveform")) {
    print("Invalid input file. It should be an object of class 'gedi.fullwaveform' ")
  } else {
    x0 <- as.data.frame(x@dt)
    x <- x0[, 1]
    z <- x0[, 2]

    if (relative == TRUE) {
      x <- c(x - min(x)) / (max(x) - min(x)) * 100
    } else {
      x <- x
    }

    if (polygon == TRUE) {
      xstart <- x[which(z == min(z, na.rm = T))]
      xend <- x[which(z == max(z, na.rm = T))]

      xl <- c(min(x), min(x), xstart, rev(x), xend, min(x))
      yl <- c(max(z, na.rm = T), min(z, na.rm = T), min(z, na.rm = T), rev(z), max(z, na.rm = T), max(z, na.rm = T))

      suppressWarnings({
        plot(xl, yl, ...)
      })
      suppressWarnings({
        polygon(xl, yl, ...)
      })
    } else {
      suppressWarnings({
        plot(x = x, y = z, ...)
      })
    }
  }
})


h5closeall <- function(con, ...) {
  try(con@h5$close_all(), silent = TRUE)
}


#' Close hdf5 connections from gedi* objects
#'
#' @description
#' Closing files will avoid locking HDF5 GEDI files.
#'
#' @param con An object of class `gedi.level*`
#' @param ... Inherited from base
#'
#' @export
#' @rdname close
#' @method close gedi.level1b
setGeneric("close", function(con, ...) {
  standardGeneric("close")
})

#' Handles the [`rGEDI::gedi.level1b-class`].
#' @rdname close
setMethod("close", signature = c("gedi.level1b"), h5closeall)
#' Handles the [`rGEDI::gedi.level2a-class`].
#' @rdname close
setMethod("close", signature = c("gedi.level2a"), h5closeall)
#' Handles the [`rGEDI::gedi.level2b-class`].
#' @rdname close
setMethod("close", signature = c("gedi.level2b"), h5closeall)
