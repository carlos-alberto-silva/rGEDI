# Earthdata Cloud --------------------------------------------------------

#' @importFrom R6 R6Class
NULL

.earthdata_netrc <- function() {
  candidates <- unique(c(
    Sys.getenv("NETRC", unset = ""),
    file.path(Sys.getenv("USERPROFILE", unset = ""), ".netrc"),
    file.path(Sys.getenv("HOME", unset = ""), ".netrc"),
    path.expand("~/.netrc")
  ))
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (length(candidates)) candidates[[1L]] else ""
}

.earthdata_env_from_netrc <- function(path) {
  if (!nzchar(path) || !file.exists(path)) return(FALSE)
  tokens <- scan(path, what = character(), quiet = TRUE, comment.char = "")
  machines <- which(tolower(tokens) == "machine")
  for (start in machines) {
    if (start == length(tokens) ||
        tolower(tokens[start + 1L]) != "urs.earthdata.nasa.gov") next
    finish <- c(machines[machines > start] - 1L, length(tokens))[[1L]]
    block <- tokens[start:finish]
    login_at <- match("login", tolower(block))
    password_at <- match("password", tolower(block))
    if (!is.na(login_at) && login_at < length(block) &&
        !is.na(password_at) && password_at < length(block)) {
      Sys.setenv(EARTHDATA_USERNAME = block[login_at + 1L],
                 EARTHDATA_PASSWORD = block[password_at + 1L])
      return(TRUE)
    }
  }
  FALSE
}

.earthdata_sync_python_env <- function() {
  values <- Sys.getenv(c("EARTHDATA_USERNAME", "EARTHDATA_PASSWORD",
                         "EARTHDATA_TOKEN"), unset = "")
  values <- as.list(values[nzchar(values)])
  if (length(values) && requireNamespace("reticulate", quietly = TRUE)) {
    os <- reticulate::import("os", convert = FALSE)
    os$environ$update(reticulate::r_to_py(values))
  }
  invisible(NULL)
}

.is_gedi_url <- function(x) {
  is.character(x) && length(x) == 1L &&
    grepl("^(https?|s3)://", x, ignore.case = TRUE)
}

.infer_gedi_product <- function(x, product = NULL) {
  if (!is.null(product) && nzchar(product)) return(toupper(product))
  source <- if (inherits(x, "GEDICloudH5")) x$source else as.character(x)[1L]
  name <- toupper(basename(sub("[?#].*$", "", source)))
  hit <- regmatches(name, regexpr("GEDI0[124]_[AB]", name))
  if (!length(hit) || !nzchar(hit)) {
    stop("Cannot infer a supported GEDI product; supply `product`.", call. = FALSE)
  }
  hit
}

# Cloud wrappers deliberately mirror the small part of the hdf5r interface
# used by rGEDI's extractors. Data stay remote until a dataset is indexed.
GEDICloudH5 <- R6Class("GEDICloudH5", public = list(
  h5 = NULL,
  handles = NULL,
  source = NULL,
  product = NULL,
  initialize = function(h5, handles = NULL, source = NULL, product = NULL) {
    self$h5 <- h5
    self$handles <- handles
    self$source <- source
    self$product <- product
  },
  ls = function(recursive = FALSE) {
    if (isTRUE(recursive)) return(self$dt_datasets(TRUE))
    keys <- reticulate::py_to_r(
      reticulate::import_builtins(convert = FALSE)$list(self$h5$keys())
    )
    data.frame(name = as.character(keys), stringsAsFactors = FALSE)
  },
  ls_groups = function(recursive = FALSE) {
    out <- character()
    if (isTRUE(recursive)) {
      self$h5$visititems(function(name, obj) {
        if (inherits(obj, "h5py._hl.group.Group")) out <<- c(out, as.character(name))
        NULL
      })
    } else {
      items <- reticulate::py_to_r(
        reticulate::import_builtins(convert = FALSE)$list(self$h5$items())
      )
      for (item in items) {
        if (inherits(item[[2L]], "h5py._hl.group.Group")) out <- c(out, as.character(item[[1L]]))
      }
    }
    out
  },
  dt_datasets = function(recursive = FALSE) {
    out <- list()
    add <- function(name, obj) {
      if (inherits(obj, "h5py._hl.dataset.Dataset")) {
        shape <- as.numeric(reticulate::py_to_r(obj$shape))
        out[[length(out) + 1L]] <<- data.frame(
          name = as.character(name),
          dataset.dims = paste(rev(shape), collapse = " x "),
          dataset.rank = length(shape), stringsAsFactors = FALSE
        )
      }
      NULL
    }
    if (isTRUE(recursive)) {
      self$h5$visititems(add)
    } else {
      items <- reticulate::py_to_r(
        reticulate::import_builtins(convert = FALSE)$list(self$h5$items())
      )
      for (item in items) add(item[[1L]], item[[2L]])
    }
    data.table::rbindlist(out, fill = TRUE)
  },
  exists = function(path) {
    tryCatch({ self$h5[[path]]; TRUE }, error = function(e) FALSE)
  },
  close_all = function() {
    try(self$h5$close(), silent = TRUE)
    self$h5 <- NULL
    self$handles <- NULL
    invisible(NULL)
  }
))

GEDICloudDataset <- R6Class("GEDICloudDataset", public = list(
  ds = NULL,
  dims = NULL,
  chunk_dims = NULL,
  initialize = function(ds) {
    self$ds <- ds
    self$dims <- rev(as.numeric(reticulate::py_to_r(ds$shape)))
    chunks <- reticulate::py_to_r(ds$chunks)
    self$chunk_dims <- if (is.null(chunks)) NULL else rev(as.numeric(chunks))
  }
))

#' @export
`[[.GEDICloudH5` <- function(x, i, ...) {
  obj <- x$h5[[i]]
  if (inherits(obj, "h5py._hl.dataset.Dataset")) GEDICloudDataset$new(obj) else
    GEDICloudH5$new(obj, handles = x$handles, source = x$source, product = x$product)
}

#' @export
names.GEDICloudH5 <- function(x) x$ls()$name

.cloud_dataset_value <- function(x) {
  value <- x$ds$`__getitem__`(reticulate::tuple())
  dtype <- as.character(x$ds$dtype)
  integer64 <- grepl("^(u?int64|<u8|<i8|>u8|>i8)", dtype)
  if (integer64) value <- value$astype("str")
  value <- reticulate::py_to_r(value)
  if (integer64) value <- bit64::as.integer64(value)
  if (!is.null(dim(value)) && length(dim(value)) > 1L) {
    value <- aperm(value, rev(seq_along(dim(value))))
  }
  value
}

#' @export
`[.GEDICloudDataset` <- function(x, ...) {
  value <- .cloud_dataset_value(x)
  call <- match.call()
  call[[1L]] <- quote(`[`)
  call[[2L]] <- value
  eval(call, parent.frame())
}

#' @export
length.GEDICloudDataset <- function(x) prod(x$dims)

.gedi_list_groups <- function(h5, recursive = FALSE) {
  if (inherits(h5, "GEDICloudH5")) h5$ls_groups(recursive) else
    hdf5r::list.groups(h5, recursive = recursive)
}

.gedi_list_datasets <- function(h5, recursive = FALSE) {
  if (inherits(h5, "GEDICloudH5")) h5$dt_datasets(recursive)$name else
    hdf5r::list.datasets(h5, recursive = recursive)
}

.open_cloud_h5 <- function(x, product = NULL) {
  if (inherits(x, "GEDICloudH5")) return(x)
  if (!requireNamespace("reticulate", quietly = TRUE) ||
      !reticulate::py_module_available("earthaccess") ||
      !reticulate::py_module_available("h5py")) {
    stop("Cloud HDF5 access requires Python earthaccess and h5py; run rGEDI_configure(install = TRUE).")
  }
  earthdata_login()
  ea <- reticulate::import("earthaccess", convert = FALSE)
  h5py <- reticulate::import("h5py", convert = FALSE)
  handles <- ea$open(reticulate::r_to_py(list(x)))
  h5 <- h5py$File(handles[[0L]], "r")
  GEDICloudH5$new(h5, handles = handles, source = as.character(x)[1L], product = product)
}

#' Configure optional rGEDI cloud and Earth Engine support
#'
#' @param python Optional Python executable or environment name.
#' @param install Logical; install missing Python modules.
#' @param earth_engine_project Optional Google Earth Engine project ID.
#' @param quiet Logical.
#' @return Invisibly returns a capability table.
#' @export
rGEDI_configure <- function(python = NULL, install = FALSE,
                            earth_engine_project = Sys.getenv("EE_PROJECT", unset = ""),
                            quiet = FALSE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Install the 'reticulate' package to configure cloud capabilities.")
  }
  if (!is.null(python)) {
    if (file.exists(python)) reticulate::use_python(python, required = TRUE) else reticulate::use_condaenv(python, required = TRUE)
  }
  modules <- c("earthaccess", "h5py", "ee")
  available <- vapply(modules, reticulate::py_module_available, logical(1))
  if (isTRUE(install) && any(!available)) {
    packages <- c(earthaccess = "earthaccess", h5py = "h5py",
                  ee = "earthengine-api")
    reticulate::py_install(unname(packages[!available]), pip = TRUE)
    available <- vapply(modules, reticulate::py_module_available, logical(1))
  }
  if (nzchar(earth_engine_project)) Sys.setenv(EE_PROJECT = earth_engine_project)
  ans <- data.frame(module = modules, available = unname(available))
  if (!quiet) print(ans, row.names = FALSE)
  invisible(ans)
}

#' Configure temporary NASA Earthdata S3 credentials
#'
#' Requests short-lived AWS credentials from Earthaccess and exports them for
#' GDAL/terra. Direct S3 access to GEDI is available from AWS `us-west-2`.
#'
#' @param daac NASA DAAC short name. GEDI products are hosted by ORNL DAAC.
#' @return Invisibly returns the credential list.
#' @export
earthdata_s3_credentials <- function(daac = "ORNL") {
  netrc <- .earthdata_netrc()
  if (nzchar(netrc)) {
    Sys.setenv(NETRC = normalizePath(netrc, winslash = "/"))
    .earthdata_env_from_netrc(netrc)
  }
  if (!requireNamespace("reticulate", quietly = TRUE) ||
      !reticulate::py_module_available("earthaccess")) {
    stop("Install Python 'earthaccess' with rGEDI_configure(install = TRUE).")
  }
  .earthdata_sync_python_env()
  ea <- reticulate::import("earthaccess", convert = FALSE)
  strategy <- if (nzchar(Sys.getenv("EARTHDATA_TOKEN")) ||
                  nzchar(Sys.getenv("EARTHDATA_USERNAME"))) {
    "environment"
  } else if (nzchar(Sys.getenv("NETRC"))) {
    "netrc"
  } else {
    "interactive"
  }
  auth <- ea$login(strategy = strategy, persist = FALSE)
  raw_credentials <- if (toupper(daac) %in% c("ORNL", "ORNLDAAC")) {
    auth$get_s3_credentials(
      endpoint = "https://data.ornldaac.earthdata.nasa.gov/s3credentials"
    )
  } else {
    auth$get_s3_credentials(daac = daac)
  }
  credentials <- reticulate::py_to_r(raw_credentials)
  if (!length(credentials)) stop("Earthdata returned no temporary S3 credentials.")
  value <- function(...) {
    keys <- c(...)
    hit <- keys[keys %in% names(credentials)]
    if (length(hit)) as.character(credentials[[hit[[1L]]]]) else ""
  }
  Sys.setenv(
    AWS_ACCESS_KEY_ID = value("accessKeyId", "AccessKeyId", "aws_access_key_id"),
    AWS_SECRET_ACCESS_KEY = value("secretAccessKey", "SecretAccessKey", "aws_secret_access_key"),
    AWS_SESSION_TOKEN = value("sessionToken", "SessionToken", "aws_session_token"),
    AWS_REGION = "us-west-2",
    AWS_DEFAULT_REGION = "us-west-2",
    GDAL_DISABLE_READDIR_ON_OPEN = "EMPTY_DIR"
  )
  if (requireNamespace("terra", quietly = TRUE)) {
    try(terra::setGDALconfig("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR"),
        silent = TRUE)
  }
  invisible(credentials)
}

#' Authenticate with NASA Earthdata
#'
#' Uses an existing bearer token, environment credentials, or `.netrc` file.
#'
#' @param persist Logical; allow Earthaccess to persist an interactive login.
#' @param netrc Path to an Earthdata `.netrc` file.
#' @return Invisibly returns `TRUE` after successful authentication.
#' @export
earthdata_login <- function(persist = TRUE, netrc = Sys.getenv("NETRC", unset = "")) {
  if (!nzchar(netrc)) netrc <- .earthdata_netrc()
  if (nzchar(netrc)) {
    Sys.setenv(NETRC = normalizePath(netrc, winslash = "/", mustWork = TRUE))
    .earthdata_env_from_netrc(netrc)
  }
  if (!requireNamespace("reticulate", quietly = TRUE) || !reticulate::py_module_available("earthaccess")) {
    if (nzchar(Sys.getenv("EARTHDATA_TOKEN")) ||
        (nzchar(Sys.getenv("EARTHDATA_USERNAME")) && nzchar(Sys.getenv("EARTHDATA_PASSWORD"))) ||
        (nzchar(Sys.getenv("NETRC")) && file.exists(Sys.getenv("NETRC")))) return(invisible(TRUE))
    stop("Install Python 'earthaccess' with rGEDI_configure(install = TRUE), or configure EARTHDATA_TOKEN/NETRC.")
  }
  .earthdata_sync_python_env()
  ea <- reticulate::import("earthaccess", convert = FALSE)
  strategy <- if (nzchar(Sys.getenv("EARTHDATA_TOKEN")) || nzchar(Sys.getenv("EARTHDATA_USERNAME"))) "environment" else if (nzchar(Sys.getenv("NETRC"))) "netrc" else "interactive"
  auth <- ea$login(strategy = strategy, persist = persist)
  ok <- isTRUE(reticulate::py_to_r(auth$authenticated))
  if (!ok) stop("NASA Earthdata authentication failed.")
  invisible(TRUE)
}

#' Open GEDI data from local storage or Earthdata Cloud
#'
#' Raster products are opened with GDAL through [terra::rast()]. HDF5 cloud
#' granules are opened lazily with Python `earthaccess` and `h5py`.
#'
#' @param x Local path, HTTPS URL, S3 URI, OPeNDAP URL, or an Earthaccess granule.
#' @param product Optional product identifier.
#' @param stream Logical; stream raster data rather than download it.
#' @param cache_dir Download directory for non-streamed data.
#' @return A `SpatRaster`, `gedi.cloud_h5`, or local GEDI object.
#' @export
openGEDI <- function(x, product = NULL, stream = TRUE, cache_dir = tempdir()) {
  if (inherits(x, c("SpatRaster", "gedi.cloud_h5", "gedi.level1b", "gedi.level2a", "gedi.level2b", "gedi.level4a"))) return(x)
  is_raster <- is.character(x) && grepl("\\.(tif|tiff)(\\?|$)", x, ignore.case = TRUE)
  if (is_raster) {
    source <- x
    if (isTRUE(stream) && grepl("^https?://", source)) source <- paste0("/vsicurl/", source)
    if (isTRUE(stream) && grepl("^s3://", source)) source <- sub("^s3://", "/vsis3/", source)
    return(terra::rast(source))
  }
  if (is.character(x) && file.exists(x)) return(.read_local_gedi(x, product))
  product <- .infer_gedi_product(x, product)
  switch(product,
    GEDI01_B = readLevel1B(x), GEDI02_A = readLevel2A(x),
    GEDI02_B = readLevel2B(x), GEDI04_A = readLevel4A(x),
    readGEDICloud(x, product = product, cache_dir = cache_dir)
  )
}

.read_local_gedi <- function(path, product = NULL) {
  if (is.null(product)) {
    name <- toupper(basename(path))
    product <- sub("^(GEDI0[124]_[AB]).*$", "\\1", name)
  }
  switch(toupper(product),
    GEDI01_B = readLevel1B(path), GEDI02_A = readLevel2A(path),
    GEDI02_B = readLevel2B(path), GEDI04_A = readLevel4A(path),
    stop("Cannot infer a supported GEDI HDF5 product from the filename."))
}

#' Open a GEDI HDF5 granule through Earthaccess
#'
#' @param x URL, S3 URI, or Earthaccess granule result.
#' @param product Optional product identifier.
#' @param cache_dir Used when cloud modules are unavailable and a URL is downloaded.
#' @return An object of class `gedi.cloud_h5`.
#' @export
readGEDICloud <- function(x, product = NULL, cache_dir = tempdir()) {
  if (!requireNamespace("reticulate", quietly = TRUE) ||
      !reticulate::py_module_available("earthaccess") ||
      !reticulate::py_module_available("h5py")) {
    if (!is.character(x)) stop("Cloud HDF5 access requires Python earthaccess and h5py.")
    gediDownload(x, cache_dir)
    return(.read_local_gedi(file.path(cache_dir, basename(x)), product))
  }
  ans <- .open_cloud_h5(x, product = product)
  class(ans) <- unique(c("gedi.cloud_h5", class(ans)))
  ans
}

#' @export
print.gedi.cloud_h5 <- function(x, ...) {
  keys <- x$ls()$name
  cat("GEDI cloud HDF5", if (!is.null(x$product)) paste0(" (", x$product, ")"), "\n", sep = "")
  cat(paste(keys, collapse = "\n"), "\n")
  invisible(x)
}

#' Close a cloud HDF5 connection
#' @param con A `gedi.cloud_h5` object returned by [readGEDICloud()].
#' @param ... Reserved for compatibility with [close()].
#' @return Invisibly returns `NULL`.
#' @method close gedi.cloud_h5
#' @export
close.gedi.cloud_h5 <- function(con, ...) {
  con$close_all()
  invisible(NULL)
}
