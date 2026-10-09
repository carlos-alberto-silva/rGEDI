test_that("GEDI sampling methods select valid footprints", {
  set.seed(10)
  x <- data.frame(longitude = runif(100, -1, 1), latitude = runif(100, -1, 1),
                  stratum = rep(letters[1:4], each = 25), value = rnorm(100))
  expect_equal(nrow(sampleGEDI(x, randomSampling(10))), 10)
  expect_equal(nrow(sampleGEDI(x, stratifiedSampling(2, "stratum"))), 8)
  expect_gt(nrow(sampleGEDI(x, gridSampling(1, 1))), 0)
  expect_lte(nrow(sampleGEDI(x, spacedSampling(10, 1000))), 10)

  polygon <- sf::st_as_sf(data.frame(zone = "all", wkt =
    "POLYGON((-2 -2,2 -2,2 2,-2 2,-2 -2))"), wkt = "wkt", crs = 4326)
  expect_equal(nrow(sampleGEDI(x, geomSampling(5, polygon, "zone"))), 5)
  r <- terra::rast(nrows = 2, ncols = 2, xmin = -2, xmax = 2,
                   ymin = -2, ymax = 2, crs = "EPSG:4326")
  terra::values(r) <- c(1, 1, 2, 2)
  expect_lte(nrow(sampleGEDI(x, rasterSampling(3, r))), 6)
  expect_s4_class(to_vect(x), "SpatVector")
})

test_that("GEDI modeling, prediction, and rasterization work", {
  set.seed(1)
  x <- data.frame(a = rnorm(40), b = runif(40))
  y <- 3 + 2 * x$a - x$b + rnorm(40, sd = .1)
  model <- fit_model(x, y, method = "lm", test = "kfold", k = 4, seed = 2)
  expect_s3_class(model, "gedi_model")
  expect_true(all(c("rmse", "adj_r2") %in% fit_metrics(y, predict(model, x))$stat))
  expect_true("a" %in% varSel(x, y, method = "correlation", threshold = .2)$selected)
  pred <- predictGEDI(model, x)
  expect_equal(nrow(pred), 40)
  h5file <- tempfile(fileext = ".h5")
  expect_equal(predictGEDIH5(model, x, h5file), h5file)
  h5 <- hdf5r::H5File$new(h5file, mode = "r")
  expect_equal(length(h5[["prediction"]][]), 40)
  h5$close_all()

  points <- data.frame(longitude = runif(40), latitude = runif(40), value = y)
  expect_s4_class(rasterizeGEDI(points, "value", res = .25), "SpatRaster")
})

test_that("random-forest GEDI modeling works when installed", {
  skip_if_not_installed("randomForest")
  set.seed(9)
  x <- data.frame(a = rnorm(35), b = runif(35))
  y <- x$a * 2 + x$b + rnorm(35, sd = .2)
  model <- fit_model(x, y, method = "randomForest", test = "split",
                     test_size = .2, seed = 4, ntree = 20)
  expect_s3_class(model, "gedi_model")
  expect_length(predict(model, x), 35)
  expect_true(length(varSel(x, y, method = "randomForest", threshold = 0,
                            ntree = 20)$selected) > 0)
})

test_that("portable waveform simulation and metrics work", {
  set.seed(3)
  cloud <- data.frame(
    X = rnorm(500, sd = 3), Y = rnorm(500, sd = 3),
    Z = c(rnorm(250, 1, .3), rnorm(250, 18, 2))
  )
  f <- tempfile(fileext = ".h5")
  sim <- gediWFSimulator(cloud, output = f, coords = c(0, 0), seed = 1)
  on.exit(close(sim), add = TRUE)
  expect_s4_class(sim, "gedi.level1b")
  expect_s4_class(getLevel1BWF(sim, 0), "gedi.fullwaveform")
  metrics <- gediWFMetrics(sim)
  expect_equal(nrow(metrics), 1)
  expect_true(all(c("cover", "rh100", "waveEnergy") %in% names(metrics)))

  txt <- tempfile(fileext = ".txt")
  ascii <- gediWFSimulator(cloud, output = txt, coords = c(0, 0), ascii = TRUE)
  expect_true(file.exists(txt))
  expect_s3_class(ascii, "data.table")
})

test_that("orbit animation writes self-contained HTML", {
  orbit <- data.frame(longitude = seq(-60, -50, length.out = 12),
                      latitude = seq(-5, 5, length.out = 12),
                      delta_time = 1:12, beam = rep(c("A", "B"), each = 6))
  f <- tempfile(fileext = ".html")
  expect_equal(plot_gedi_orbit_animation(orbit, output_file = f, launch = FALSE),
               normalizePath(f, winslash = "/"))
  html <- paste(readLines(f, warn = FALSE), collapse = "\n")
  expect_match(html, "Interactive GEDI ground-track playback", fixed = TRUE)
  expect_false(grepl("<script src=", html, fixed = TRUE))
})
