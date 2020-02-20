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
