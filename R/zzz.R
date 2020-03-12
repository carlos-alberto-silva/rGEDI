.onUnload <- function (libpath) {
  objs = c("gedi.level2a", "gedi.level1b", "gedi.level2b", "gedi.level1bSim")
  classes = sapply(ls(.GlobalEnv), function(x) class(get(x)))
  classes = classes[classes %in% objs]
  objNames = names(classes)
  for (o in objNames) {
    try(close(get(o, envir = .GlobalEnv)), silent=T)
  }
  library.dynam.unload("rGEDI", libpath)
  invisible()
}

.onLoad <- function(libname, pkgname){
  # This will reopen closed files which are in the environment,
  # necessary for cleaning up between calls of gedisimulator.
  # 'gedisimulator' is a foreign C library, which have some
  # memory issues. Unloading and reloading solves the problem,
  # but unloading should handle hdf5 opened files to avoid
  # IO locking, so that when reloading, the objects needs to
  # be opened again.
  objs = c("gedi.level2a", "gedi.level1b", "gedi.level2b", "gedi.level1bSim")
  classes = sapply(ls(.GlobalEnv), function(x) class(get(x)))
  classes = classes[classes %in% objs]
  objNames = names(classes)
  for (i in 1:length(objNames)) {
    if (length(classes) > 0) {
      assign(
        objNames[i], 
        new(
          classes[i][[1]],
          h5=hdf5r::h5file(get(objNames[i])@h5$filename)), 
        envir=.GlobalEnv)
    }
  }
  
  invisible()
}