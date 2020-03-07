.onUnload <- function (libpath) {
  objs = c("gedi.level1a", "gedi.level1b", "gedi.level2a", "gedi.level1bSim")
  classes = sapply(ls(.GlobalEnv), function(x) class(get(x)))
  classes = classes[classes %in% objs]
  objNames = names(classes)
  for (o in objNames) {
    try(close(get(o)), silent=T)
  }
  rm(list=as.vector(objNames))
  library.dynam.unload("rGEDI", libpath)
  invisible()
}
