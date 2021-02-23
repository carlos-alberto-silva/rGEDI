Rcpp::loadModule("gdal_module", TRUE)

.onUnload <- function (libpath) {
  library.dynam.unload("rGEDI", libpath)
  invisible()
}

.onLoad <- function(libname, pkgname) {
  .RGEDI_CACHE <- new.env(FALSE, parent=globalenv())
  assign("old.PROJ_LIB", Sys.getenv("PROJ_LIB"), envir=.RGEDI_CACHE)
}

.onAttach <- function(libname, pkgname) {
  #Setup copied from rgdal package
  Sys.setenv("PROJ_LIB"=system.file("proj", package = pkgname)[1])
}

#' @export
GDALDataType = list(
  GDT_Byte = 1,
  GDT_UInt16 = 2,
  GDT_Int16 = 3,
  GDT_UInt32 = 4,
  GDT_Int32 = 5,
  GDT_Float32 = 6,
  GDT_Float64 = 7
)

#' @export CPP_GDALDataset CPP_GDALRasterBand create_dataset
Rcpp::loadModule("gdal_module", TRUE)


#' @export
GDALDataset <- R6::R6Class("GDALDataset",
                    private=list(
                     ds = NULL,
                     datatype = NULL
                    ),
                    public=list(
                    initialize = function(name, nbands, datatype, proj4string, ul_lat, ul_lon, lr_lat, lr_lon, res, nodata=0, co=NULL) {
                      private$ds = create_dataset(name, nbands, datatype, proj4string, lr_lat, ul_lat, ul_lon, lr_lon, res, nodata, co)
                      private$datatype = datatype
                    },
                    close = function() {
                      private$ds$Close()
                    },
                    GetRasterBand = function(x) {
                      GDALRasterBand$new(private$ds$GetRasterBand(x), private$datatype)
                    },
                    GetRasterXSize = function() private$ds$GetRasterXSize(),
                    GetRasterYSize = function() private$ds$GetRasterYSize()
                  )
)

#' @export
'[[.GDALDataset' = function(x, slice) {
  return (x$GetRasterBand(slice))
}

#' @export
GDALRasterBand <- R6::R6Class("GDALRasterBand",
                    private=list(
                     band = NULL,
                     datatype = NULL
                    ),
                    public=list(
                    initialize = function(band, datatype) {
                      private$band = band
                      private$datatype = datatype
                    },
                    close = function() {
                      private$band = NULL
                    },
                    ReadBlock = function(iXBlock, iYBlock) {
                      funcs = list(
                        private$band$ReadBlock1,
                        private$band$ReadBlock2,
                        private$band$ReadBlock3,
                        private$band$ReadBlock4,
                        private$band$ReadBlock5,
                        private$band$ReadBlock6,
                        private$band$ReadBlock7
                      )
                      funcs[[private$datatype]](iXBlock, iYBlock)
                    },
                    WriteBlock = function(iXBlock, iYBlock, buffer) {
                      funcs = list(
                        private$band$WriteBlock1,
                        private$band$WriteBlock2,
                        private$band$WriteBlock3,
                        private$band$WriteBlock4,
                        private$band$WriteBlock5,
                        private$band$WriteBlock6,
                        private$band$WriteBlock7
                      )
                      funcs[[private$datatype]](iXBlock, iYBlock, buffer)
                    }
                  )
)
