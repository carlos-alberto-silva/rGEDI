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


#' Creates a new GDALDataset
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
#' 
#' @examples 
#' # Parameters
#' raster_path = file.path(tempdir(), "output.tif")
#' ul_lat = -15
#' ul_lon = -45
#' lr_lat = -25
#' lr_lon = -35
#' res = c(0.01, -0.01)
#' datatype = GDALDataType$GDT_Int32
#' nbands = 1
#' projstring = "EPSG:4326"
#' nodata = -1
#' co = c("TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512", "COMPRESSION=LZW")
#'
#' # Create a new raster dataset
#' ds = createDataset(
#' raster_path = raster_path, 
#' nbands = nbands, 
#' datatype = datatype, 
#' projstring = projstring, 
#' lr_lat = lr_lat, 
#' ul_lat = ul_lat, 
#' ul_lon = ul_lon, 
#' lr_lon = lr_lon, 
#' res = res, 
#' nodata = nodata, 
#' co = co
#' )
#'
#' # Get the GDALRasterBand for ds
#' band = ds[[1]]
#' 
#' # Set some dummy values
#' band[[0,0]] = 1:(512*512)
#'
#' # Calculate the square - 10
#' formulaCalculate(
#'     formula = "x*2 - 10", 
#'     data = list(x = band), 
#'     updateBand = band
#' )
#'
#' ds$Close()
#' @export
createDataset = function(raster_path, nbands, datatype, projstring, lr_lat, ul_lat, ul_lon, lr_lon, res, nodata, co = NULL) {
  ds = create_dataset(raster_path, nbands, datatype, projstring, lr_lat, ul_lat, ul_lon, lr_lon, res, nodata, co)
  GDALDataset$new(ds)
}

#' R6 Class GDALDataset wrapping
#'
#' @description
#' Wrapping class for GDALDataset C++ API exporting GetRasterBand, GetRasterXSize, GetRasterYSize
#' @export
GDALDataset <- R6::R6Class("GDALDataset",
                           private = list(
                             ds = NULL,
                             datatype = NULL
                           ),
                           public = list(
#' @description
#' Create a new raster file based on specified data. It will output a *.tif file.
#' @param ds GDALDatasetR pointer. Should not be used
#' @param datatype GDALDataType. The GDALDataType to use for the raster, use (GDALDataType$) to find the options. Default GDALDataType$GDT_Float64
#' @return An object from GDALDataset R6 class.
                             initialize = function(ds) {
                              private$ds = ds
                             },
#' @description
#' Function to retrieve the GDALRasterBand R6 Object.
#' @param x Integer. The band index, starting from 1 to number of bands.
#' @return An object of GDALRasterBand R6 class.
                             GetRasterBand = function(x) {
                              GDALRasterBand$new(private$ds$GetRasterBand(x))
                             },
#' @description
#' Get the width for the raster
#' @return An integer indicating the raster width
                             GetRasterXSize = function() private$ds$GetRasterXSize(),
#' @description
#' Get the height for the raster
#' @return An integer indicating the raster height
                             GetRasterYSize = function() private$ds$GetRasterYSize(),

#' @description
#' Closes the GDALDataset
#' @return An integer indicating the raster width
                             Close = function() private$ds$Close()
                           )
)


#' Open GDAL raster
#' @description
#' Function to open GDAL Dataset
#'
#' @param filename Character. The path to a GDAL dataset.
#' @param readonly Logical. Flag to open a read only GDALDataset with GA_ReadOnly or GA_Update. Default TRUE.
#' @return An R6 object of GDALDataset class.
#' @export
GDALOpen = function(filename, readonly = TRUE) {
  ds = RGDALOpen(filename, readonly)
  GDALDataset$new(ds)
}

#' GDALDataset [[]] accessor
#' @description
#' This function gives access to the GDALRasterBand using [[i]], where i is the band index to return.
#' @param x GDALDatset. Automatically obtained from GDALDataset[[]] call.
#' @param slice Integer. The index for the band to access.
#' @return An object of GDALRasterBand R6 class.
#' @export
'[[.GDALDataset' = function(x, slice) {
  return(x$GetRasterBand(slice))
}

#' R6 Class GDALRasterBand wrapping
#'
#' @description
#' Wrapping class for GDALRasterBand C++ API exporting GetBlockXSize, GetBlockYSize, ReadBlock, WriteBlock for better IO speed.
#' @export
GDALRasterBand <- R6::R6Class("GDALRasterBand",
                              private = list(
                                band = NULL,
                                datatype = NULL
                              ),
                              public = list(
#' @description
#' Creates a new GDALRasterBand
#' @param band The C++ pointer to the GDALRasterBandR object.
#' @param datatype The GDALDataType for this band
#' @return An object of GDALRasterBand R6 class
#'
#' @note This constructor must not be called at all, this is automatically called from GDALDataset$GetRasterBand function.
                                initialize = function(band) {
                                  private$band = band
                                },
#' @description
#' Efficiently reads a raster block
#'
#' @param iXBlock Integer. The i-th column block to access. The `iXBlock` will be offset \eqn{ BLOCKXSIZE \times iXBlock }{BLOCKXSIZE x iXBlock} from the origin.
#' @param iYBlock Integer. The i-th row block to access. The `iYBlock` will be offset \eqn{ BLOCKYSIZE \times iYBlock }{BLOCKYSIZE x iYBlock} from the origin.
#'
#' @return RawVector for GDALDataType$GDT_Byte, IntegerVector for int types and NumericVector for floating point types.
#' @details
#' The returned Vector will be single dimensional with the length \eqn{ BLOCKXSIZE \times BLOCKYSIZE }{BLOCKXSIZE x BLOCKYSIZE}. If you use matrix(, ncol=BLOCKXSIZE) the matrix returned will actually be transposed. You should either transpose it or you can calculate the indices using \eqn{ y \cdot xsize + x }{y*xsize + x}
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
                                  result = funcs[[self$GetRasterDataType()]](iXBlock, iYBlock)
                                  nodata = self$GetNoDataValue()
                                  result[result == nodata] = NA
                                  return(result)
                                },
#' @description
#' Efficiently writes a raster block
#'
#' @param iXBlock Integer. The i-th column block to write. The `iXBlock` will be offset \eqn{ BLOCKXSIZE \times iXBlock }{BLOCKXSIZE x iXBlock} from the origin.
#' @param iYBlock Integer. The i-th row block to write. The `iYBlock` will be offset \eqn{ BLOCKYSIZE \times iYBlock }{BLOCKYSIZE x iYBlock} from the origin.
#' @param buffer RawVector/IntegerVector/NumericVector depending on the GDALDataType. This should be a 1D vector with size equal to raster \eqn{ BLOCKXSIZE \times BLOCKYSIZE }{BLOCKXSIZE x BLOCKYSIZE}.
#'
#' @return Nothing
#' @details
#' The returned Vector will be single dimensional with the length \eqn{ BLOCKXSIZE \times BLOCKYSIZE }{BLOCKXSIZE x BLOCKYSIZE}. If you use matrix(, ncol=BLOCKXSIZE) the matrix returned will actually be transposed. You should either transpose it or you can calculate the indices using \eqn{ y \cdot xsize + x }{ y*xsize + x }.
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
                                  nodata = self$GetNoDataValue()
                                  buffer[is.na(buffer)] = nodata
                                  funcs[[self$GetRasterDataType()]](iXBlock, iYBlock, buffer)
                                },
#' @description
#' Get the block width
#' @return An integer indicating block width
                                GetBlockXSize = function() {
                                  private$band$GetBlockXSize()
                                },
#' @description
#' Get the block height
#' @return An integer indicating block height
                                GetBlockYSize = function() {
                                  private$band$GetBlockYSize()
                                },

#' @description
#' Get the band width
#' @return An integer indicating band width
                                GetXSize = function() {
                                  private$band$GetXSize()
                                },

#' @description
#' Get the band height
#' @return An integer indicating band height
                                GetYSize = function() {
                                  private$band$GetYSize()
                                },

#' @description
#' Get band no data value
#' @return Numeric indicating no data value
                                GetNoDataValue = function() {
                                  private$band$GetNoDataValue()
                                },

#' @description
#' Get band datatype
#' @return Numeric indicating the datatype
                                GetRasterDataType = function() {
                                  private$band$GetRasterDataType()
                                }
                              )
)


#' GDALRasterBand [[]]= setter
#' @description
#' This function gives access to the GDALRasterBand using [[i]], where i is the band index to return.
#' @param x GDALRasterBand. Automatically obtained from GDALDataset[[]] call.
#' @param blockX Integer. The x index for block to access.
#' @param blockY Integer. The y index for block to access.
#' @param value Integer. The value buffer to write
#'
#' @return Nothing, this is a setter
#' @export
'[[<-.GDALRasterBand' = function(x, blockX, blockY, value) {
  x$WriteBlock(blockX, blockY, value)
  return(x)
}

#' GDALRasterBand [[]] getter
#' @description
#' This function gives access to the GDALRasterBand using [[i]], where i is the band index to return.
#' @param x GDALRasterBand. Automatically obtained from GDALDataset[[]] call.
#' @param blockX Integer. The x index for block to access.
#' @param blockY Integer. The y index for block to access.
#'
#' @return Nothing, this is a setter
#' @export
'[[.GDALRasterBand' = function(x, blockX, blockY) {
  x$ReadBlock(blockX, blockY)
}

#' Calculate raster values based on a formula
#'
#' @param formula Character. A textual formula to apply to the RasterBands from `data`
#' @param data List. A named list with the used variables in the textual formula
#' @param updateBand GDALRasterBand. The GDALRasterBand which will be updated with the calculated values.
#'
#' @return Nothing, it just updates the band of interest.
#'
#' @examples 
#' # Parameters
#' raster_path = file.path(tempdir(), "output.tif")
#' ul_lat = -15
#' ul_lon = -45
#' lr_lat = -25
#' lr_lon = -35
#' res = c(0.01, -0.01)
#' datatype = GDALDataType$GDT_Int32
#' nbands = 1
#' projstring = "EPSG:4326"
#' nodata = -1
#' co = c("TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512", "COMPRESSION=LZW")
#'
#' # Create a new raster dataset
#' ds = createDataset(
#' raster_path = raster_path, 
#' nbands = nbands, 
#' datatype = datatype, 
#' projstring = projstring, 
#' lr_lat = lr_lat, 
#' ul_lat = ul_lat, 
#' ul_lon = ul_lon, 
#' lr_lon = lr_lon, 
#' res = res, 
#' nodata = nodata, 
#' co = co
#' )
#'
#' # Get the GDALRasterBand for ds
#' band = ds[[1]]
#' 
#' # Set some dummy values
#' band[[0,0]] = 1:(512*512)
#'
#' # Calculate the square - 10
#' formulaCalculate(
#'     formula = "x*2 - 10", 
#'     data = list(x = band), 
#'     updateBand = band
#' )
#'
#' ds$Close()
#' 
#' @export
formulaCalculate =  function(formula, data, updateBand) {
  first = data[[1]]
  blocksize1 = c(first$GetBlockXSize(), first$GetBlockYSize())
  bandsize1 = c(first$GetXSize(), first$GetYSize())
  iters = floor(bandsize1 / blocksize1)
  nodata = first$GetNoDataValue()
  form = as.formula(paste0("~I(",formula,")"))
  totalBlocks = (iters[1] + 1) * (iters[2] + 1)
  blockCounter = 0
  termNames = names(data)
  thisTerms = termNames[sapply(Vectorize(grep, vectorize.args = c("pattern"))(termNames, form), length) > 0]
  thisData = data[thisTerms]

  for (xblock in 0:iters[1]) {
      for (yblock in 0:iters[2]) {
        blockCounter = blockCounter + 1
        message(
          sprintf(
            "\rCalculating block %d/%d (%.2f%%)           ", 
            blockCounter, 
            totalBlocks, 
            100*blockCounter/totalBlocks
          ), appendLF = FALSE
        )
        vals = as.data.table(lapply(thisData, function(x) x[[xblock, yblock]]))
        mask = !is.na(vals[[1]])
        if (!any(mask)) next
        vals[[1]][mask] = model.frame(form, vals[mask], na.action = na.pass)[[1]]
        vals[[1]][!mask] = updateBand$GetNoDataValue()
        updateBand[[xblock, yblock]] = vals[[1]]
      }
  }
  message()
}

