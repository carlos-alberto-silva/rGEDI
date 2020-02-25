unloadLibrary = function() {
  if (isLibraryLoaded()) {
    unloadNamespace("rGEDI")
  }
  require("rGEDI", quietly = T)
  invisible()
}


isLibraryLoaded = function() {
  is.loaded("C_gediSimulator", "rGEDI")
}
