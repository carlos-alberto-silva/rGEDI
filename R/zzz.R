Rcpp::loadModule("gdal_module", TRUE)
.RGEDI_CACHE <- new.env(FALSE, parent=globalenv())

.onUnload <- function (libpath) {
  library.dynam.unload("rGEDI", libpath)
  Sys.setenv("PROJ_LIB"=get("old.PROJ_LIB", envir=.RGEDI_CACHE))
  invisible()
}

.onLoad <- function(libname, pkgname) {
  assign("old.PROJ_LIB", Sys.getenv("PROJ_LIB"), envir=.RGEDI_CACHE)
}

.onAttach <- function(libname, pkgname) {
  #Setup copied from rgdal package
  Sys.setenv("PROJ_LIB"=system.file("proj", package = pkgname)[1])
}

