.RGEDI_CACHE <- new.env(FALSE, parent=globalenv())

.onUnload <- function (libpath) {
  Sys.setenv("PROJ_LIB"=get("old.PROJ_LIB", envir=.RGEDI_CACHE))
  library.dynam.unload("rGEDI", libpath)
  invisible()
}

.onLoad <- function(libname, pkgname) {
  assign("old.PROJ_LIB", Sys.getenv("PROJ_LIB"), envir=.RGEDI_CACHE)
}

.onAttach <- function(libname, pkgname) {
  #Setup copied from rgdal package
  Sys.setenv("PROJ_LIB"=system.file("proj", package = pkgname)[1])
}

#' List of datatypes supported by the GDALDataset R6 class
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

defaultProjection = 'GEOGCS["WGS 84",
    DATUM["WGS_1984",
        SPHEROID["WGS 84",6378137,298.257223563,
            AUTHORITY["EPSG","7030"]],
        AUTHORITY["EPSG","6326"]],
    PRIMEM["Greenwich",0,
        AUTHORITY["EPSG","8901"]],
    UNIT["degree",0.01745329251994328,
        AUTHORITY["EPSG","9122"]],
    AUTHORITY["EPSG","4326"]]'

Rcpp::loadModule("gdal_module", TRUE)

#' R6 Class GDALDataset wrapping
#' 
#' @description 
#' Wrapping class for GDALDataset C++ API exporting GetRasterBand, GetRasterXSize, GetRasterYSize
#' @export
GDALDataset <- R6::R6Class("GDALDataset",
                    private=list(
                     ds = NULL,
                     datatype = NULL
                    ),
                    public=list(
                    #' @description
                    #' Create a new raster file based on specified data. It will output a *.tif file.
                    #' @param raster_path Character. The output path for the raster data
                    #' @param ul_lat Numeric. The upper left latitude.
                    #' @param ul_lon Numeric. The upper left longitude.
                    #' @param lr_lat Numeric. The lower right latitude.
                    #' @param lr_lon Numeric. The lower right longitude.
                    #' @param res Numeric. The resolution of the output raster
                    #' @param datatype GDALDataType. The GDALDataType to use for the raster, use (GDALDataType$) to find the options. Default GDALDataType$GDT_Float64
                    #' @param nbands Integer. Number of bands. Default 1.
                    #' @param projstring The projection string, either proj or WKT is accepted. 
                    #' @param nodata Numeric. The no data value for the raster. Default 0.
                    #' @param co CharacterVector. A CharacterVector of creation options for GDAL. Default NULL
                    #' @return An object from GDALDataset R6 class.
                    initialize = function(raster_path, ul_lat, ul_lon, lr_lat, lr_lon, res, datatype = GDALDataType$GDT_Float64, nbands = 1, projstring = defaultProjection, nodata=0, co=NULL) {
                      private$ds = create_dataset(raster_path, nbands, datatype, projstring, lr_lat, ul_lat, ul_lon, lr_lon, res, nodata, co)
                      private$datatype = datatype
                    },
                    #' @description
                    #' Function to retrieve the GDALRasterBand R6 Object.
                    #' @param x Integer. The band index, starting from 1 to number of bands.
                    #' @return An object of GDALRasterBand R6 class.
                    GetRasterBand = function(x) {
                      GDALRasterBand$new(private$ds$GetRasterBand(x), private$datatype)
                    },
                    #' @description
                    #' Get the width for the raster
                    #' @return An integer indicating the raster width
                    GetRasterXSize = function() private$ds$GetRasterXSize(),
                    #' @description
                    #' Get the height for the raster
                    #' @return An integer indicating the raster height
                    GetRasterYSize = function() private$ds$GetRasterYSize()
                  )
)

#' GDALDataset [[]] accessor
#' @description
#' This function gives access to the GDALRasterBand using [[i]], where i is the band index to return.
#' @param slice Integer. The index for the band to access.
#' @return An object of GDALRasterBand R6 class.
#' @export
'[[.GDALDataset' = function(x, slice) {
  return (x$GetRasterBand(slice))
}

#' R6 Class GDALRasterBand wrapping
#' \loadmathjax 
#' 
#' @description 
#' Wrapping class for GDALRasterBand C++ API exporting ReadBlock, WriteBlock for better IO speed.
#' @export
GDALRasterBand <- R6::R6Class("GDALRasterBand",
                    private=list(
                     band = NULL,
                     datatype = NULL
                    ),
                    public=list(
                    #' @description
                    #' Creates a new GDALRasterBand
                    #' @param band The C++ pointer to the GDALRasterBandR object.
                    #' @param datatype The GDALDataType for this band
                    #' @return An object of GDALRasterBand R6 class
                    #' 
                    #' @note This constructor must not be called at all, this is automatically called from GDALDataset$GetRasterBand function.
                    initialize = function(band, datatype) {
                      private$band = band
                      private$datatype = datatype
                    },
                    #' @description
                    #' Efficiently reads a raster block
                    #' 
                    #' @param iXBlock Integer. The i-th column block to access. The `iXBlock` will be offset \mjeqn{ BLOCKXSIZE \times iXBlock }{BLOCKXSIZE x iXBlock} from the origin.
                    #' @param iYBlock Integer. The i-th row block to access. The `iYBlock` will be offset \mjeqn{ BLOCKYSIZE \times iYBlock }{BLOCKYSIZE x iYBlock} from the origin.
                    #' 
                    #' @return RawVector for GDALDataType$GDT_Byte, IntegerVector for int types and NumericVector for floating point types.
                    #' @details 
                    #' The returned Vector will be single dimensional with the length \mjeqn{ BLOCKXSIZE \times BLOCKYSIZE }{BLOCKXSIZE x BLOCKYSIZE}. If you use matrix(, ncol=BLOCKXSIZE) the matrix returned will actually be transposed. You should either transpose it or you can calculate the indices using \mjeqn{ y \cdot xsize + x }{y*xsize + x}
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
                    #' @description
                    #' Efficiently writes a raster block
                    #' 
                    #' @param iXBlock Integer. The i-th column block to write. The `iXBlock` will be offset \mjeqn{ BLOCKXSIZE \times iXBlock }{BLOCKXSIZE x iXBlock} from the origin.
                    #' @param iYBlock Integer. The i-th row block to write. The `iYBlock` will be offset \mjeqn{ BLOCKYSIZE \times iYBlock }{BLOCKYSIZE x iYBlock} from the origin.
                    #' @param buffer RawVector/IntegerVector/NumericVector depending on the GDALDataType. This should be a 1D vector with size equal to raster \mjeqn{ BLOCKXSIZE \times BLOCKYSIZE }{BLOCKXSIZE x BLOCKYSIZE}.
                    #'
                    #' @return Nothing
                    #' @details
                    #' The returned Vector will be single dimensional with the length \mjeqn{ BLOCKXSIZE \times BLOCKYSIZE }{BLOCKXSIZE x BLOCKYSIZE}. If you use matrix(, ncol=BLOCKXSIZE) the matrix returned will actually be transposed. You should either transpose it or you can calculate the indices using \mjeqn{ y \cdot xsize + x }{ y*xsize + x }.
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
