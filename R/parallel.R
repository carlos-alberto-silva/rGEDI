# Shared parallel helpers ------------------------------------------------

.gedi_worker_count <- function(workers = NULL, tasks = Inf) {
  tasks <- as.integer(tasks)[1L]
  if (is.na(tasks) || tasks < 1L) return(0L)
  if (is.null(workers)) {
    detected <- parallel::detectCores(logical = FALSE)
    workers <- if (is.na(detected)) 2L else max(1L, detected - 1L)
  }
  workers <- as.integer(workers)[1L]
  if (is.na(workers) || workers < 1L) {
    stop("'workers' must be a positive integer.", call. = FALSE)
  }
  min(workers, tasks)
}

.gedi_parallel_lapply <- function(x, fun, parallel = FALSE, workers = NULL,
                                  worker_options = list()) {
  if (!length(x)) return(list())
  workers <- .gedi_worker_count(workers, length(x))
  if (!isTRUE(parallel) || workers < 2L) {
    answer <- lapply(x, fun)
    attr(answer, "workers") <- 1L
    return(answer)
  }
  if (.Platform$OS.type == "windows") {
    cluster <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterCall(cluster, function() loadNamespace("rGEDI"))
    if (length(worker_options)) {
      parallel::clusterCall(cluster, function(value) {
        options(value)
        NULL
      }, worker_options)
    }
    answer <- parallel::parLapply(cluster, x, fun)
  } else {
    answer <- parallel::mclapply(x, fun, mc.cores = workers, mc.preschedule = TRUE)
  }
  attr(answer, "workers") <- workers
  answer
}

.extract_gedi_granule <- function(spec) {
  product <- spec$product
  object <- switch(
    product,
    GEDI01_B = readLevel1B(spec$path),
    GEDI02_A = readLevel2A(spec$path),
    GEDI02_B = readLevel2B(spec$path),
    GEDI04_A = readLevel4A(spec$path)
  )
  on.exit(close(object), add = TRUE)
  values <- switch(
    product,
    GEDI01_B = do.call(getLevel1BGeo, c(list(object), spec$args)),
    GEDI02_A = do.call(getLevel2AM, c(list(object), spec$args)),
    GEDI02_B = do.call(getLevel2BVPM, c(list(object), spec$args)),
    GEDI04_A = do.call(getLevel4A, c(list(object), spec$args))
  )
  values <- data.table::as.data.table(values)
  values[, `:=`(granule = spec$name, granule_index = spec$index)]
  values[]
}

.new_gedi_granule_collection <- function(sources, product) {
  if (!is.character(sources) || length(sources) < 2L || anyNA(sources)) {
    stop("A granule collection requires at least two paths or URLs.", call. = FALSE)
  }
  structure(
    list(sources = unclass(sources), product = product),
    class = c("gedi_granule_collection", "list")
  )
}

.extract_gedi_collection <- function(x, args, parallel = FALSE, workers = NULL) {
  specs <- lapply(seq_along(x$sources), function(i) list(
    path = x$sources[[i]],
    name = basename(sub("[?#].*$", "", x$sources[[i]])),
    index = i, product = x$product, args = args
  ))
  pieces <- .gedi_parallel_lapply(
    specs,
    function(spec) getFromNamespace(".extract_gedi_granule", "rGEDI")(spec),
    parallel = parallel, workers = workers
  )
  used_workers <- attr(pieces, "workers")
  answer <- data.table::rbindlist(pieces, use.names = TRUE, fill = TRUE)
  data.table::setorder(answer, granule_index)
  attr(answer, "workers") <- used_workers
  answer[]
}

#' @export
print.gedi_granule_collection <- function(x, ...) {
  cat("<", x$product, " granule collection>\n", sep = "")
  cat("Granules: ", length(x$sources), "\n", sep = "")
  invisible(x)
}

#' @method close gedi_granule_collection
#' @export
close.gedi_granule_collection <- function(con, ...) invisible(NULL)
