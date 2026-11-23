# An independent reader for the public R decision-tree text format.
predict_tree_text <- function(text, newdata) {
  tokens <- strsplit(trimws(strsplit(text, "\n", fixed = TRUE)[[1L]]), "[[:space:]]+")
  ids <- vapply(tokens, function(x) as.numeric(sub(")", "", x[1L], fixed = TRUE)), numeric(1))
  nodes <- stats::setNames(tokens, ids)
  vapply(seq_len(nrow(newdata)), function(i) {
    id <- 1
    repeat {
      node <- nodes[[as.character(id)]]
      if (tail(node, 1L) == "*") return(as.numeric(node[length(node) - 1L]))
      left <- nodes[[as.character(2 * id)]]
      stopifnot(left[3L] == "<=")
      id <- if (newdata[[left[2L]]][i] <= as.numeric(left[4L])) 2 * id else 2 * id + 1
    }
  }, numeric(1))
}

test_that("serialized trees reproduce the fitted numeric forest including split boundaries", {
  skip_if_not_installed("randomForest")
  set.seed(918)
  x <- data.frame(A00 = runif(120, -10, 10), elevation = runif(120, 100, 300))
  y <- x$A00^2 + x$elevation / 10
  rf <- randomForest::randomForest(x, y, ntree = 7, maxnodes = 8)
  tree <- randomForest::getTree(rf, 1, labelVar = FALSE)
  variable <- tree[1L, "split var"]
  split <- tree[tree[, "split var"] == variable, "split point"]
  boundary <- data.frame(A00 = rep(0, 3 * length(split)), elevation = 180)
  boundary[[variable]] <- rep(split, each = 3) + rep(c(-1e-12, 0, 1e-12), length(split))
  probe <- rbind(x, data.frame(A00 = c(-100, 0, 100), elevation = c(0, 200, 1000)), boundary)
  strings <- rGEDI:::.rf_tree_strings(rf)
  decoded <- rowMeans(vapply(strings, predict_tree_text, numeric(nrow(probe)), newdata = probe))
  expect_equal(decoded, as.numeric(predict(rf, probe)), tolerance = 1e-12)
  expect_length(strings, rf$ntree)
  expect_identical(strings, rGEDI:::.rf_tree_strings(structure(list(model = rf), class = "gedi_model")))
})

test_that("constant regression trees keep their leaf predictions", {
  skip_if_not_installed("randomForest")
  x <- data.frame(A00 = seq_len(20))
  rf <- suppressWarnings(randomForest::randomForest(x, rep(7.25, 20), ntree = 3))
  strings <- rGEDI:::.rf_tree_strings(rf)
  expect_equal(predict_tree_text(strings[1], x), rep(7.25, 20))
})

test_that("unsupported forests fail before Python or Earth Engine initialization", {
  skip_if_not_installed("randomForest")
  local_mocked_bindings(.require_ee = function() stop("Unexpected EE initialization"), .package = "rGEDI")
  expect_error(build_ee_forest(lm(mpg ~ wt, mtcars)), "regression model")
  rf <- randomForest::randomForest(Species ~ ., iris, ntree = 2)
  expect_error(build_ee_forest(rf), "regression model")
  x <- data.frame(group = factor(rep(letters[1:2], 10)), value = seq_len(20))
  rf <- randomForest::randomForest(x, seq_len(20), ntree = 2)
  expect_error(build_ee_forest(rf), "numeric predictors")
  expect_error(map_create(list(), NULL), "fitted random forest")
})

test_that("map_create converts local models internally and keeps EE models usable", {
  skip_if_not_installed("randomForest")
  set.seed(32)
  model <- fit_model(data.frame(A00 = seq_len(20)), seq_len(20)^2,
    ntree = 3, maxnodes = 4)
  seen <- new.env(parent = emptyenv())
  ee_model <- structure(list(id = "converted"), class = "python.builtin.object")
  local_mocked_bindings(.require_ee = function() list(Classifier = list(
    decisionTreeEnsemble = function(strings) {
      seen$strings <- strings
      list(setOutputMode = function(mode) {seen$mode <- mode; ee_model})
    }
  )), .package = "rGEDI")
  image <- list(toFloat = function() "mapped")
  stack <- list(classify = function(model, name) {seen$model <- model; seen$name <- name; image})
  expect_identical(map_create(model, stack, name = "AGBD"), "mapped")
  expect_identical(seen$mode, "REGRESSION")
  expect_identical(seen$model, ee_model)
  expect_identical(seen$name, "AGBD")
  expect_length(seen$strings, 3)
  expect_identical(map_create(ee_model, stack), "mapped")
})
