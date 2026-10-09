# Modeling and prediction ------------------------------------------------

#' Model fit metrics
#' @param observed,predicted Numeric vectors.
#' @param plotstat Draw an observed-versus-predicted plot.
#' @param unit Unit label.
#' @param ... Passed to [graphics::plot()].
#' @return A data frame of accuracy statistics.
#' @export
fit_metrics <- function(observed, predicted, plotstat = FALSE, unit = "", ...) {
  ok <- stats::complete.cases(observed, predicted)
  observed <- as.numeric(observed[ok]); predicted <- as.numeric(predicted[ok])
  if (length(observed) < 2L) stop("At least two complete pairs are required.")
  error <- predicted - observed
  denom <- mean(observed)
  fit <- stats::lm(observed ~ predicted)
  value <- c(
    rmse = sqrt(mean(error^2)), mae = mean(abs(error)), bias = mean(error),
    r = stats::cor(observed, predicted), adj_r2 = summary(fit)$adj.r.squared)
  relative <- if (isTRUE(all.equal(denom, 0))) rep(NA_real_, 3L) else 100 * value[c("rmse", "mae", "bias")] / denom
  ans <- data.frame(stat = c("rmse", "rmseR", "mae", "maeR", "bias", "biasR", "r", "adj_r2"),
    value = c(value["rmse"], relative[1L], value["mae"], relative[2L], value["bias"], relative[3L], value["r"], value["adj_r2"]),
    unit = c(unit, "%", unit, "%", unit, "%", "", ""), row.names = NULL)
  if (isTRUE(plotstat)) {
    graphics::plot(observed, predicted, ...)
    graphics::abline(0, 1, col = "red", lwd = 2)
    graphics::abline(fit, col = "black", lwd = 2)
  }
  ans
}

.make_folds <- function(n, k) {
  k <- as.integer(k)
  if (!is.finite(k) || k < 2L || k > n) stop("`k` must be between 2 and the number of observations.")
  split(base::sample.int(n), rep(seq_len(k), length.out = n))
}

#' Fit a GEDI statistical or random-forest model
#'
#' @param x Predictor data frame.
#' @param y Response vector or response column name in `x`.
#' @param method `"randomForest"` or `"lm"`.
#' @param test Validation method: `"none"`, `"split"`, `"kfold"`,
#'   `"loocv"`, or `"bootstrap"`.
#' @param k Number of folds.
#' @param test_size Fraction held out for split validation.
#' @param iterations Bootstrap iterations.
#' @param seed Optional random seed.
#' @param ... Model arguments.
#' @return An object of class `gedi_model`.
#' @export
fit_model <- function(x, y, method = c("randomForest", "lm"),
                      test = c("none", "split", "kfold", "loocv", "bootstrap"),
                      k = 5L, test_size = 0.3, iterations = 100L,
                      seed = NULL, ...) {
  method <- match.arg(method); test <- match.arg(test)
  if (!is.null(seed)) set.seed(seed)
  x <- as.data.frame(x)
  if (is.character(y) && length(y) == 1L) { response <- x[[y]]; x[[y]] <- NULL } else response <- as.numeric(y)
  ok <- stats::complete.cases(x, response); x <- x[ok, , drop = FALSE]; response <- response[ok]
  if (nrow(x) < 3L) stop("At least three complete observations are required.")
  trainer <- function(xx, yy) {
    if (method == "lm") return(stats::lm(.response ~ ., data = cbind(.response = yy, xx), ...))
    if (!requireNamespace("randomForest", quietly = TRUE)) stop("Install 'randomForest' for this method.")
    randomForest::randomForest(x = xx, y = yy, ...)
  }
  predictor <- function(model, xx) as.numeric(stats::predict(model, newdata = xx))
  full <- trainer(x, response)
  fitted <- predictor(full, x)
  validation <- rep(NA_real_, nrow(x))
  groups <- list()
  if (test == "split") {
    ntest <- max(1L, round(nrow(x) * test_size)); groups <- list(base::sample.int(nrow(x), ntest))
  } else if (test == "kfold") groups <- .make_folds(nrow(x), k)
  else if (test == "loocv") groups <- as.list(seq_len(nrow(x)))
  if (length(groups)) for (idx in groups) {
    train <- setdiff(seq_len(nrow(x)), idx)
    validation[idx] <- predictor(trainer(x[train, , drop = FALSE], response[train]), x[idx, , drop = FALSE])
  }
  if (test == "bootstrap") {
    values <- vector("list", nrow(x))
    for (i in seq_len(iterations)) {
      train <- base::sample.int(nrow(x), replace = TRUE)
      hold <- setdiff(seq_len(nrow(x)), unique(train))
      if (length(hold)) {
        pred <- predictor(trainer(x[train, , drop = FALSE], response[train]), x[hold, , drop = FALSE])
        for (j in seq_along(hold)) values[[hold[j]]] <- c(values[[hold[j]]], pred[j])
      }
    }
    validation <- vapply(values, function(z) if (length(z)) mean(z) else NA_real_, numeric(1))
  }
  structure(list(model = full, method = method, predictors = names(x), response = response,
    fitted = fitted, validation = validation, stats_train = fit_metrics(response, fitted),
    stats_test = if (sum(is.finite(validation)) >= 2L) fit_metrics(response[is.finite(validation)], validation[is.finite(validation)]) else NULL,
    test = test), class = "gedi_model")
}

#' @export
predict.gedi_model <- function(object, newdata, ...) stats::predict(object$model, newdata = as.data.frame(newdata)[, object$predictors, drop = FALSE], ...)

#' Select predictor variables for a GEDI model
#' @param x Predictor data frame.
#' @param y Response vector.
#' @param method Selection method.
#' @param threshold Minimum absolute correlation for correlation selection.
#' @param ... Passed to `randomForest` when requested.
#' @return An object of class `gedi_var_selection`.
#' @export
varSel <- function(x, y, method = c("correlation", "randomForest"), threshold = 0.1, ...) {
  method <- match.arg(method); x <- as.data.frame(x)
  if (method == "correlation") {
    scores <- vapply(x, function(z) abs(stats::cor(z, y, use = "complete.obs")), numeric(1))
  } else {
    if (!requireNamespace("randomForest", quietly = TRUE)) stop("Install 'randomForest'.")
    fit <- randomForest::randomForest(x = x, y = y, importance = TRUE, ...)
    imp <- randomForest::importance(fit)
    scores <- imp[, ncol(imp)]
    scores <- scores / max(scores, na.rm = TRUE)
  }
  structure(list(selected = names(scores)[is.finite(scores) & scores >= threshold], scores = sort(scores, decreasing = TRUE), method = method), class = "gedi_var_selection")
}

#' Predict GEDI footprint values
#' @param model A `gedi_model` or fitted R model.
#' @param data Predictor table.
#' @param name Output column name.
#' @return A data table containing input data and predictions.
#' @export
predictGEDI <- function(model, data, name = "prediction") {
  pred <- if (inherits(model, "gedi_model")) predict(model, data) else stats::predict(model, newdata = as.data.frame(data))
  ans <- data.table::as.data.table(data)
  ans[, (name) := as.numeric(pred)]
  ans
}

#' Predict values and write them to HDF5
#' @param model A `gedi_model` or fitted R model.
#' @param data Predictor table.
#' @param output Output HDF5 filename.
#' @param dataset Dataset name.
#' @return Invisibly returns `output`.
#' @export
predictGEDIH5 <- function(model, data, output, dataset = "prediction") {
  values <- predictGEDI(model, data, dataset)[[dataset]]
  h5 <- hdf5r::H5File$new(output, mode = "w"); on.exit(h5$close_all(), add = TRUE)
  h5[[dataset]] <- values
  invisible(output)
}

#' Rasterize GEDI footprint data
#' @param x GEDI footprint table.
#' @param metric Metric column.
#' @param res Resolution in decimal degrees.
#' @param fun Aggregation function.
#' @param lon,lat Optional coordinate columns.
#' @param ... Passed to [terra::rasterize()].
#' @return A [`terra::SpatRaster-class`].
#' @export
rasterizeGEDI <- function(x, metric, res = 0.01, fun = mean, lon = NULL, lat = NULL, ...) {
  coords <- .gedi_coords(x, lon, lat)
  .grid_gedi_points(x, coords[1L], coords[2L], metric, fun, res, ...)
}
