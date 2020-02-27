stopifnotMessage = function(...) {
  ok = TRUE
  errors = list()
  listargs = list(...)
  for (i in 1:length(listargs)) {
    if (listargs[i] == FALSE) {
      errors[[""]] = names(listargs)[i]
      ok = FALSE
    }
  }
  if (ok == FALSE) {
    stop(paste0("\n\nError validating the arguments:\n    ", paste(errors, collapse="\n    ")))
  }
}

checkNumericLength = function(x, len) {
  return (is.null(x) || (length(x) == len && is.numeric(x)))
}

checkNumeric = function(x) {
  return (checkNumericLength(x, 1))
}

checkLogical = function(x) {
  return (is.null(x) || (length(x) == 1 && is.logical(x)))
}

checkInteger = function(x) {
  x_int = as.integer(x)
  return (is.null(x) || (length(x_int) == 1 && is.integer(x_int) && !is.na(x_int)))
}

checkCharacter = function(x) {
  return (is.null(x) || (length(x) == 1 && is.character(x)))
}

checkFilepath = function(x, newFile=TRUE, optional = TRUE) {
  exists = TRUE
  if (is.null(x)) {
    if (optional)
      return (TRUE)
    else
      return (FALSE)
  }
  if (!is.character(x) || length(x) != 1)
    return (FALSE)


  if (!newFile)
    return (file.exists(x))

  return (TRUE)
}

checkParentDir = function(x, optional=FALSE) {
  if (optional && is.null(x)) {
    return (TRUE)
  }
  dirName = fs::path_dir(x)
  return (fs::dir_exists(dirName)[[1]])
}

inputOrInList = function(input) {
  inList=NULL
  if (length(input) > 1) {
    inList = tempfile(fileext=".txt")
    fileHandle = file(inList, "w")
    writeLines(input, fileHandle)
    close(fileHandle)
    inFile = fileHandle
    return (list(NULL, inList))
  }
  return (list(input, NULL))
}

cleanInList = function(x) {
  if (!is.null(x[[2]]) && file.exists(x[[2]])) {
    file.remove(x[[2]])
  }
}
