test_that("Level 3 and Level 4B raster workflows operate end to end", {
  r <- terra::rast(nrows = 8, ncols = 8, xmin = -2, xmax = 2,
                   ymin = -2, ymax = 2, crs = "EPSG:4326")
  terra::values(r) <- seq_len(terra::ncell(r))
  names(r) <- "metric"
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(r, f)

  l3 <- readLevel3(f)
  l4b <- readLevel4B(f)
  expect_s4_class(l3, "SpatRaster")
  expect_s4_class(l4b, "SpatRaster")
  expect_equal(terra::ncell(clipLevel3(l3, c(-1, 1, -1, 1))), 16)
  expect_equal(terra::ncell(clipLevel4B(l4b, c(-1, 1, -1, 1))), 16)
  expect_equal(terra::ncell(aggregateLevel3(l3, 2)), 16)
  expect_equal(terra::ncell(aggregateLevel4B(l4b, 2)), 16)

  template <- terra::rast(nrows = 4, ncols = 4, xmin = -2, xmax = 2,
                          ymin = -2, ymax = 2, crs = "EPSG:4326")
  expect_equal(terra::ncell(resampleLevel3(l3, template)), 16)
  expect_equal(terra::ncell(resampleLevel4B(l4b, template)), 16)
  expect_s4_class(projectLevel3(l3, "EPSG:3857"), "SpatRaster")
  expect_s4_class(projectLevel4B(l4b, "EPSG:3857"), "SpatRaster")
  expect_s4_class(mosaicLevel3(list(l3, l3)), "SpatRaster")
  expect_s4_class(mosaicLevel4B(list(l4b, l4b)), "SpatRaster")

  p <- sf::st_as_sf(data.frame(id = 1, x = 0, y = 0),
                    coords = c("x", "y"), crs = 4326)
  poly <- sf::st_as_sf(data.frame(id = "center", wkt =
    "POLYGON((-1 -1,1 -1,1 1,-1 1,-1 -1))"), wkt = "wkt", crs = 4326)
  expect_equal(nrow(extractLevel3(l3, p)), 1)
  expect_equal(nrow(extractLevel4B(l4b, p)), 1)
  expect_equal(nrow(polyStatsLevel3(l3, poly)), 1)
  expect_equal(nrow(polyStatsLevel4B(l4b, poly)), 1)

  png(tempfile(fileext = ".png")); expect_s4_class(plotLevel3(l3), "SpatRaster"); dev.off()
  png(tempfile(fileext = ".png")); expect_s4_class(plotLevel4B(l4b), "SpatRaster"); dev.off()
  expect_s4_class(openGEDI(f), "SpatRaster")
})

test_that("Level 4A HDF5 workflows operate end to end", {
  f <- tempfile(pattern = "GEDI04_A_", fileext = ".h5")
  h5 <- hdf5r::H5File$new(f, mode = "w")
  b <- h5$create_group("BEAM0101")
  values <- list(
    shot_number = 1:6, delta_time = 11:16, l4_quality_flag = c(1, 1, 0, 1, 1, 1),
    degrade_flag = c(0, 0, 0, 1, 0, 0), lat_lowestmode = seq(0, .05, .01),
    lon_lowestmode = seq(0, .05, .01), agbd = 10:15, agbd_se = rep(1, 6)
  )
  for (nm in names(values)) b[[nm]] <- values[[nm]]
  h5$close_all()

  x <- readLevel4A(f)
  on.exit(close(x), add = TRUE)
  all <- getLevel4A(x, cols = c("beam", names(values)))
  good <- getLevel4A(x, cols = c("beam", names(values)), quality = TRUE)
  expect_equal(nrow(all), 6)
  expect_equal(nrow(good), 4)
  expect_equal(nrow(clipLevel4A(all, 0, .03, 0, .03)), 4)

  poly <- sf::st_as_sf(data.frame(zone = "a", wkt =
    "POLYGON((-0.01 -0.01,0.031 -0.01,0.031 0.031,-0.01 0.031,-0.01 -0.01))"),
    wkt = "wkt", crs = 4326)
  expect_equal(nrow(clipLevel4AGeometry(all, poly, "zone")), 4)
  expect_equal(nrow(polyStatsLevel4A(all, poly, id = "zone")), 1)
  expect_s4_class(gridStatsLevel4A(all, res = .01), "SpatRaster")
  expect_s4_class(rasterizeLevel4A(all, res = .01), "SpatRaster")
  png(tempfile(fileext = ".png")); expect_equal(plotLevel4A(all), all); dev.off()
})

test_that("Level 4A Release 3 quality fields are applied", {
  f <- tempfile(pattern = "GEDI04_A_", fileext = ".h5")
  h5 <- hdf5r::H5File$new(f, mode = "w")
  b <- h5$create_group("BEAM0000")
  values <- list(
    shot_number = 1:5,
    lat_lowestmode = seq(0, .04, .01),
    lon_lowestmode = seq(0, .04, .01),
    agbd = 11:15,
    l4a_quality_flag_rel3 = c(1, 1, 0, 1, 1),
    degrade_include_flag = c(1, 0, 1, 1, 1),
    elev_highestreturn_outlier_flag = c(0, 0, 0, 1, 0)
  )
  for (nm in names(values)) b[[nm]] <- values[[nm]]
  h5$close_all()

  x <- readLevel4A(f)
  on.exit(close(x), add = TRUE)
  good <- getLevel4A(x, cols = c("beam", names(values)), quality = TRUE)
  expect_equal(good$shot_number, c(1, 5))
})

