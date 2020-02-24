stopifnotMessage = function(msg, ...) {
  tryCatch(
    stopifnot(...),
    error=function(e) {
      stop(paste0("\n\nInput arguments are invalid!\n",
                  e$message))
    }
  )
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
  return (is.null(x) || (length(x)) == 1 && is.character(x))
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

  return TRUE
}
