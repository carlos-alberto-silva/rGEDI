make_gee_backend <- function(available = TRUE) {
  list(
    available = function() available,
    open = function(reference, quality = TRUE) {
      structure(list(asset = unclass(reference)[1L], quality = quality), class = "mock_ee")
    },
    collect = function(reference, variables, quality, max_features, chunk_size) {
      values <- list(
        shot_number = c("1", "2"),
        lon_lowestmode = c(-44.5, -44.4),
        lat_lowestmode = c(-13.5, -13.4)
      )
      for (variable in setdiff(variables, names(values))) {
        values[[variable]] <- if (variable == "beam") c("BEAM0000", "BEAM0001") else c(10, 20)
      }
      data.frame(values, check.names = FALSE)
    },
    select = function(reference, variables, quality) {
      structure(list(asset = unclass(reference)[1L], variables = variables), class = "mock_ee")
    },
    clip = function(x, clip_obj) {
      structure(list(input = x$object, clip = clip_obj), class = "mock_ee_clip")
    },
    rasterize = function(x, metric, res, fun, crs) {
      structure(list(metric = metric, res = res, crs = crs), class = "mock_ee_image")
    }
  )
}

test_that("gedifinder creates a lazy GEE reference", {
  old <- options(rGEDI.gee_backend = make_gee_backend())
  on.exit(options(old), add = TRUE)

  ref <- gedifinder(
    "GEDI02_A", 10, 20, 9, 21,
    daterange = c("2023-01-01", "2023-01-31"), source = "gee"
  )
  expect_s3_class(ref, "gedi_gee_reference")
  expect_identical(unclass(ref)[1L], "LARSE/GEDI/GEDI02_A_002")
  expect_identical(attr(ref, "version"), "002")

  tab <- gedifinder("GEDI04_A", 10, 20, 9, 21, source = "gee", return = "table")
  expect_s3_class(tab, "data.table")
  expect_identical(tab$source, "gee")
  expect_identical(tab$version, "002.1")
})

test_that("GEE readers collect footprints or keep them lazy", {
  old <- options(rGEDI.gee_backend = make_gee_backend())
  on.exit(options(old), add = TRUE)

  l2a <- readLevel2A(gedifinder("GEDI02_A", 10, 20, 9, 21, source = "gee"))
  expect_s3_class(l2a, "gedi_gee_object")
  values <- getLevel2AM(l2a, rh = 98, quality = TRUE)
  expect_s3_class(values, "data.table")
  expect_true(all(c("shot_number", "rh98") %in% names(values)))
  expect_identical(attr(values, "gedi_source"), "gee")

  lazy <- getLevel2AM(l2a, rh = 98, collect = FALSE)
  expect_s3_class(lazy, "gedi_gee_object")
  expect_identical(lazy$object$variables, c(
    "beam", "shot_number", "degrade_flag", "quality_flag",
    "sensitivity", "solar_elevation", "elev_lowestmode", "rh98"
  ))

  l4a <- readLevel4A(gedifinder("GEDI04_A", 10, 20, 9, 21, source = "gee"))
  biomass <- getLevel4A(l4a, cols = c("shot_number", "agbd", "agbd_se"))
  expect_true(all(c("shot_number", "agbd", "agbd_se") %in% names(biomass)))
})

test_that("Level 4B reader accepts its GEE raster reference", {
  old <- options(rGEDI.gee_backend = make_gee_backend())
  on.exit(options(old), add = TRUE)

  l4b <- readLevel4B(gedifinder("GEDI04_B", 10, 20, 9, 21, source = "gee"))
  expect_s3_class(l4b, "gedi_gee_object")
})

test_that("generic clip and rasterization preserve lazy GEE processing", {
  old <- options(rGEDI.gee_backend = make_gee_backend())
  on.exit(options(old), add = TRUE)

  l2a <- readLevel2A(gedifinder("GEDI02_A", 10, 20, 9, 21, source = "gee"))
  clipped <- clip(l2a, c(20.1, 20.9, 9.1, 9.9))
  expect_s3_class(clipped, "gedi_gee_object")
  expect_s3_class(clipped$object, "mock_ee_clip")

  image <- rasterize(getLevel2AM(l2a, rh = 98, collect = FALSE), "rh98")
  expect_s3_class(image, "gedi_gee_object")
  expect_identical(image$kind, "image")
  expect_identical(image$object$res, 1000)
  expect_identical(image$object$crs, "EPSG:3857")
})

test_that("GEE validation fails early and explains version limits", {
  old <- options(rGEDI.gee_backend = make_gee_backend(FALSE))
  on.exit(options(old), add = TRUE)

  expect_error(
    gedifinder("GEDI02_A", 10, 20, 9, 21, source = "gee"),
    "unavailable"
  )
  expect_error(
    gedifinder("GEDI02_A", 10, 20, 9, 21, source = "gee", version = "003", check = FALSE),
    "Version 002"
  )
})
