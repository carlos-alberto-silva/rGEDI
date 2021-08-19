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
    InitializeGDAL()
}

#' @export CPP_GDALDataset CPP_GDALRasterBand create_dataset
Rcpp::loadModule("gdal_module", TRUE)
