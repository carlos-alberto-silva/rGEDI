unloadLibrary = function() {
  if (isLibraryLoaded()) {
    library.dynam.unload("rGEDI", find.package("rGEDI"))
  }
  unloadNamespace("rGEDI")
  require("rGEDI", quietly = T)
  invisible()
}


isLibraryLoaded = function() {
  is.loaded("C_gediSimulator", "rGEDI")
}
