test_that("clip subsets standard GEDI footprint tables by extent", {
  x <- data.table::data.table(
    shot_number = 1:4,
    lon_lowestmode = c(-44.14, -44.12, -44.10, NA_real_),
    lat_lowestmode = c(-13.74, -13.73, -13.70, -13.73)
  )

  ans <- clip(x, c(-44.13, -44.11, -13.75, -13.72))
  expect_s3_class(ans, "data.table")
  expect_equal(ans$shot_number, 2L)

  ans_ext <- clip(x, terra::ext(-44.13, -44.11, -13.75, -13.72))
  expect_equal(ans_ext$shot_number, 2L)
})

test_that("clip requires both waveform endpoints inside an extent", {
  x <- data.table::data.table(
    shot_number = 1:3,
    longitude_bin0 = c(0.2, 0.2, 0.2),
    latitude_bin0 = c(0.2, 0.2, 0.2),
    longitude_lastbin = c(0.3, 1.2, 0.3),
    latitude_lastbin = c(0.3, 0.3, NA_real_)
  )
  ans <- clip(x, c(0, 1, 0, 1))
  expect_equal(ans$shot_number, 1L)
})

test_that("clip subsets footprint tables with polygon geometry", {
  x <- data.frame(
    shot_number = 1:3,
    lon_lowestmode = c(0.25, 1.25, 3),
    lat_lowestmode = c(0.25, 0.25, 3)
  )
  polygons <- sf::st_sf(
    stand = c("west", "east"),
    geometry = sf::st_sfc(
      sf::st_polygon(list(matrix(c(0, 0, 1, 0, 1, 1, 0, 1, 0, 0), ncol = 2, byrow = TRUE))),
      sf::st_polygon(list(matrix(c(1, 0, 2, 0, 2, 1, 1, 1, 1, 0), ncol = 2, byrow = TRUE))),
      crs = 4326
    )
  )

  ans <- clip(x, polygons, split_by = "stand")
  expect_s3_class(ans, "data.table")
  expect_equal(ans$shot_number, 1:2)
  expect_equal(ans$poly_id, c("west", "east"))
})

test_that("clip dispatches for GEDI rasters", {
  r <- terra::rast(ncols = 10, nrows = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  terra::values(r) <- seq_len(terra::ncell(r))
  ans <- clip(r, c(2, 6, 3, 8))
  expect_s4_class(ans, "SpatRaster")
  expect_equal(unname(as.vector(terra::ext(ans))), c(2, 6, 3, 8))
})

test_that("clip registers methods for open GEDI products", {
  expect_s4_class(methods::selectMethod("clip", c("gedi.level1b", "numeric")), "MethodDefinition")
  expect_s4_class(methods::selectMethod("clip", c("gedi.level2a", "SpatExtent")), "MethodDefinition")
  expect_s4_class(methods::selectMethod("clip", c("gedi.level2b", "ANY")), "MethodDefinition")
  expect_s4_class(methods::selectMethod("clip", c("gedi.level4a", "numeric")), "MethodDefinition")
})

test_that("clip validates extents and coordinate fields", {
  expect_error(clip(data.frame(longitude = 0, latitude = 0), c(1, 0, 0, 1)),
    "xmin <= xmax")
  expect_error(clip(data.frame(x = 0, y = 0), c(0, 1, 0, 1)),
    "Cannot find GEDI coordinates")
})
