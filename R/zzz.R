Rcpp::loadModule("gdal_module", TRUE)

.onUnload <- function (libpath) {
  library.dynam.unload("rGEDI", libpath)
  invisible()
}

.onLoad <- function(libname, pkgname) {
}

.onAttach <- function(libname, pkgname) {
  InitializeGDAL()
}

