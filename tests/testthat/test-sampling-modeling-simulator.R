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
  selected <- varSel(x, y, method = "rfe", threshold = 0,
                     seed = 4, ntree = 20)
  expect_s3_class(selected, "gedi_var_selection")
  expect_true(length(selected$selvars) > 0)
  expect_equal(nrow(selected$test), ncol(x))
  png(tempfile(fileext = ".png")); plot(selected); dev.off()
})

test_that("portable waveform simulation and metrics work", {
  set.seed(3)
  cloud <- data.frame(
    X = rnorm(500, sd = 3), Y = rnorm(500, sd = 3),
    Z = c(rnorm(250, 1, .3), rnorm(250, 18, 2)),
    Classification = c(rep(2L, 250), rep(5L, 250))
  )
  f <- tempfile(fileext = ".h5")
  sim <- gediWFSimulator(cloud, output = f, coords = c(0, 0), seed = 1)
  on.exit(close(sim), add = TRUE)
  expect_s4_class(sim, "gedi.level1b")
  expect_s4_class(getLevel1BWF(sim, 0), "gedi.fullwaveform")
  metrics <- gediWFMetrics(sim)
  expect_equal(nrow(metrics), 1)
  expect_true(all(c("cover", "rh100", "waveEnergy") %in% names(metrics)))
  expect_equal(metrics$ground_method, "classified_ground")
  expect_equal(metrics$waveEnergy, 1, tolerance = 0.01)
  expect_equal(metrics$gHeight, 1, tolerance = 0.5)
  expect_true(metrics$cover > 0 && metrics$cover < 1)

  trimmed_file <- tempfile(fileext = ".h5")
  trimmed <- gediWFSimulator(cloud, output = trimmed_file, coords = c(0, 0),
                             maxBins = 64, res = 0.15, seed = 1)
  on.exit(close(trimmed), add = TRUE)
  trimmed_wave <- getLevel1BWF(trimmed, 0)@dt
  expect_equal(nrow(trimmed_wave), 64)
  expect_equal(abs(diff(trimmed_wave$elevation)), rep(0.15, 63),
               tolerance = 1e-10)

  noisy_file <- tempfile(fileext = ".h5")
  noisy <- gediWFSimulator(cloud, output = noisy_file, coords = c(0, 0),
                           noise = 0.03, seed = 9)
  on.exit(close(noisy), add = TRUE)
  expect_false(isTRUE(all.equal(
    getLevel1BWF(sim, 0)@dt$rxwaveform,
    getLevel1BWF(noisy, 0)@dt$rxwaveform
  )))

  txt <- tempfile(fileext = ".txt")
  ascii <- gediWFSimulator(cloud, output = txt, coords = c(0, 0), ascii = TRUE)
  expect_true(file.exists(txt))
  expect_s3_class(ascii, "data.table")
  expect_error(gediWFSimulator(cloud, coords = c(0, 0), density_res = 0),
               "density_res")
  expect_error(gediWFSimulator(cloud, coords = c(0, 0), intensity_threshold = 2),
               "intensity_threshold")
})

test_that("orbit animation writes ICESat2-style interactive HTML", {
  orbit <- data.frame(longitude = seq(-60, -50, length.out = 12),
                      latitude = seq(-5, 5, length.out = 12),
                      delta_time = 1:12, beam = rep(c("A", "B"), each = 6))
  f <- tempfile(fileext = ".html")
  expect_equal(plot_gedi_orbit_animation(orbit, output_file = f, launch = FALSE),
               normalizePath(f, winslash = "/"))
  html <- paste(readLines(f, warn = FALSE), collapse = "\n")
  expect_match(html, "ISS + GEDI", fixed = TRUE)
  expect_match(html, "0xff1744", fixed = TRUE)
  expect_match(html, "Earth rotation speed", fixed = TRUE)
  expect_match(html, "OrbitControls", fixed = TRUE)
  expect_match(html, "data:image/png;base64", fixed = TRUE)
  expect_match(html, "data:image/jpeg;base64", fixed = TRUE)
  expect_match(html, '"reference":true', fixed = TRUE)
  expect_match(html, "cdn.jsdelivr.net/npm/three", fixed = TRUE)
  expect_error(
    plot_gedi_orbit_animation(orbit, output_file = f, track_speed = 0),
    "between 1 and 15"
  )
  expect_error(
    plot_gedi_orbit_animation(orbit, output_file = f,
                              earth_rotation_speed = 21),
    "between 0 and 20"
  )
})

test_that("GEDI tracks are standardized and thinned", {
  x <- data.frame(lon_lowestmode = 1:10, lat_lowestmode = 11:20,
                  delta_time = 10:1, beam = "BEAM0000")
  track <- getGEDITrack(x, every = 2)
  expect_s3_class(track, "data.table")
  expect_named(track, c("longitude", "latitude", "sequence", "track",
                        "delta_time", "beam"))
  expect_equal(nrow(track), 5)
  expect_true(all(diff(track$sequence) >= 0))

  split <- getGEDITrack(data.frame(
    longitude = c(0, .001, .002, 1, 1.001, 1.002),
    latitude = rep(0, 6), delta_time = 1:6, beam = "BEAM0000"
  ))
  expect_equal(length(unique(split$track)), 2)

  unsplit <- getGEDITrack(data.frame(
    longitude = c(0, .001, .002, 1, 1.001, 1.002),
    latitude = rep(0, 6), delta_time = 1:6, beam = "BEAM0000"
  ), segment_gaps = FALSE)
  expect_equal(unique(unsplit$track), "BEAM0000")

  dateline <- getGEDITrack(data.frame(
    longitude = c(179.8, 179.9, -179.9, -179.8),
    latitude = c(0, .01, .02, .03), delta_time = 1:4,
    beam = "BEAM0000"
  ))
  expect_equal(unique(dateline$track), "BEAM0000")
})
