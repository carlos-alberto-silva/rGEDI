unloadLibrary = function() {
  if (isLibraryLoaded()) {
    library.dynam.unload("rGEDI", find.package("rGEDI"))
  }
  invisible()
}

loadLibrary = function() {
  library.dynam("rGEDI", "rGEDI", paste0(find.package("rGEDI"), "/.."))
  invisible()
}

isLibraryLoaded = function() {
  is.loaded("C_gediSimulator", "rGEDI")
}