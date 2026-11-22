mock_gedidb_backend <- function(available = TRUE) {
  list(
    available = function(config, deep = FALSE) available,
    query = function(reference, variables, filters) {
      n <- 12L
      ans <- data.frame(
        latitude = seq(-13.799, -13.788, length.out = n),
        longitude = seq(-44.199, -44.188, length.out = n),
        time = rep(as.integer(as.Date("2023-06-01")), n),
        shot_number = sprintf("23748010010010%03d", seq_len(n)),
        stringsAsFactors = FALSE
      )
      if (any(grepl("^rh", variables))) ans$rh_p98 <- seq(10, 32, length.out = n)
      if ("agbd" %in% variables) ans$agbd <- seq(20, 130, length.out = n)
      if ("agbd_se" %in% variables) ans$agbd_se <- rep(4, n)
      if ("l2a_quality_flag_rel3" %in% variables) ans$l2a_quality_flag_rel3 <- 1L
      if ("l4a_quality_flag_rel3" %in% variables) ans$l4a_quality_flag_rel3 <- 1L
      if ("beam_name" %in% variables) ans$beam_name <- "BEAM0101"
      ans
    }
  )
}

test_that("gedifinder creates lazy gediDB references and fails early", {
  old <- options(rGEDI.gedidb_backend = mock_gedidb_backend())
  on.exit(options(old), add = TRUE)

  ref <- gedifinder(
    "GEDI02_A", ul_lat = -13.7, ul_lon = -44.2,
    lr_lat = -13.8, lr_lon = -44.1,
    daterange = c("2023-01-01", "2023-12-31"), source = "gedidb"
  )
  expect_s3_class(ref, "gedidb_reference")
  expect_match(unclass(ref), "gedi-l2-l4-v003")
  expect_s3_class(readLevel2A(ref), "gedidb_level2a")

  options(rGEDI.gedidb_backend = mock_gedidb_backend(FALSE))
  expect_error(gedifinder(
    "GEDI02_A", ul_lat = -13.7, ul_lon = -44.2,
    lr_lat = -13.8, lr_lon = -44.1, source = "gedidb"
  ), "unavailable")
})

test_that("gediDB RH98 and L4A data rasterize on a 1 km grid", {
  old <- options(rGEDI.gedidb_backend = mock_gedidb_backend())
  on.exit(options(old), add = TRUE)

  make_ref <- function(product) gedifinder(
    product, ul_lat = -13.7, ul_lon = -44.2,
    lr_lat = -13.8, lr_lon = -44.1,
    daterange = c("2023-01-01", "2023-12-31"), source = "gedidb"
  )

  l2a <- getLevel2AM(readLevel2A(make_ref("GEDI02_A")), rh = 98, quality = TRUE)
  expect_s3_class(l2a, "data.table")
  expect_true(all(c("lon_lowestmode", "lat_lowestmode", "rh98") %in% names(l2a)))
  expect_identical(attr(l2a, "gedi_source"), "gedidb")

  l4a <- getLevel4A(
    readLevel4A(make_ref("GEDI04_A")),
    cols = c("shot_number", "agbd", "agbd_se", "l4a_quality_flag_rel3"),
    quality = TRUE
  )
  expect_true(all(c("lon_lowestmode", "lat_lowestmode", "agbd") %in% names(l4a)))

  all_points <- rbind(
    l2a[, .(lon_lowestmode, lat_lowestmode, value = rh98)],
    l4a[, .(lon_lowestmode, lat_lowestmode, value = agbd)]
  )
  points <- terra::vect(
    as.data.frame(all_points),
    geom = c("lon_lowestmode", "lat_lowestmode"), crs = "EPSG:4326"
  )
  projected <- terra::project(points, "EPSG:6933")
  bounds <- as.vector(terra::ext(projected))
  bounds <- c(
    floor(bounds[1L] / 1000) * 1000, ceiling(bounds[2L] / 1000) * 1000,
    floor(bounds[3L] / 1000) * 1000, ceiling(bounds[4L] / 1000) * 1000
  )
  template <- terra::rast(
    terra::ext(bounds), resolution = 1000, crs = "EPSG:6933"
  )

  rh98_grid <- rasterize(l2a, "rh98", template = template)
  agbd_grid <- rasterize(l4a, "agbd", template = template)
  expect_s4_class(rh98_grid, "SpatRaster")
  expect_s4_class(agbd_grid, "SpatRaster")
  expect_equal(terra::res(rh98_grid), c(1000, 1000))
  expect_equal(terra::res(agbd_grid), c(1000, 1000))
  expect_true(any(is.finite(terra::values(rh98_grid))))
  expect_true(any(is.finite(terra::values(agbd_grid))))
})

test_that("lazy gediDB queries rasterize by spatial tiles", {
  old <- options(rGEDI.gedidb_backend = mock_gedidb_backend())
  on.exit(options(old), add = TRUE)

  ref <- gedifinder(
    "GEDI02_A", ul_lat = -13.7, ul_lon = -44.2,
    lr_lat = -13.8, lr_lon = -44.1,
    daterange = c("2023-01-01", "2023-12-31"), source = "gedidb"
  )
  lazy <- getLevel2AM(readLevel2A(ref), rh = 98, quality = TRUE, collect = FALSE)
  expect_s3_class(lazy, "gedidb_query")

  grid <- rasterize(
    lazy, "rh98", res = 1000, output_crs = "EPSG:6933",
    tile_size = 0.05, parallel = FALSE
  )
  expect_s4_class(grid, "SpatRaster")
  expect_equal(terra::res(grid), c(1000, 1000))
  expect_true(any(is.finite(terra::values(grid))))
  expect_identical(attr(grid, "gedidb_tiles"), 4L)
  expect_identical(attr(grid, "gedidb_workers"), 1L)
  expect_equal(attr(grid, "gedidb_shots"), 12)

  # Counts must conserve observations across retrieval tiles, enabling
  # correctly weighted means when spatial exports are assembled.
  for (count_fun in list("length", base::length)) {
    counts <- rasterize(lazy, "rh98", fun = count_fun, res = 1000,
                        output_crs = "EPSG:6933", tile_size = 0.05)
    expect_equal(sum(terra::values(counts), na.rm = TRUE), 12)
  }
})
