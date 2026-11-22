test_that("parallel helpers cap workers and preserve order", {
  expect_identical(rGEDI:::.gedi_worker_count(8, 3), 3L)
  expect_error(rGEDI:::.gedi_worker_count(0, 3), "positive integer")
  result <- rGEDI:::.gedi_parallel_lapply(
    1:4, function(x) x^2, parallel = FALSE, workers = 4
  )
  expect_identical(unlist(result), c(1, 4, 9, 16))
  expect_identical(attr(result, "workers"), 1L)
})

test_that("multi-granule extraction closes files and preserves input order", {
  archive <- system.file(
    "extdata", "GEDI02_A_2019108080338_O01964_T05337_02_001_01_sub.zip",
    package = "rGEDI"
  )
  path <- unzip(archive, exdir = tempdir())
  collection <- readLevel2A(c(path, path))
  expect_s3_class(collection, "gedi_granule_collection")
  result <- getLevel2AM(collection, beams = "BEAM0001", rh = 98,
                        parallel = FALSE)
  expect_s3_class(result, "data.table")
  expect_identical(unique(result$granule_index), 1:2)
  expect_true(all(c("rh98", "granule", "granule_index") %in% names(result)))
  expect_identical(attr(result, "workers"), 1L)
  halves <- split(result$rh98, result$granule_index)
  expect_equal(halves[[1L]], halves[[2L]])
})

test_that("simulator accepts the shared parallel controls", {
  cloud <- data.frame(
    X = c(-1, 0, 1, -1, 0, 1), Y = c(-1, 0, 1, 1, -1, 0),
    Z = c(0, 1, 2, 8, 10, 12),
    Classification = c(2L, 2L, 2L, 5L, 5L, 5L)
  )
  result <- gediWFSimulator(
    cloud, output = tempfile(fileext = ".txt"), ascii = TRUE,
    coords = c(0, 0, 0.5, 0.5), parallel = FALSE, workers = 2
  )
  expect_s3_class(result, "data.table")
  expect_identical(attr(result, "workers"), 1L)
  expect_identical(unique(result$waveID), 0:1)
})

test_that("archive batches combine independently accumulated cells", {
  testthat::local_mocked_bindings(
    accumulateGEDIGranule = function(entry, from = NULL, to = NULL, region = NULL) {
      data.table::data.table(
        col_km = c(1, 2), row_km = c(3, 4),
        total = c(entry$value, 2 * entry$value), shots = c(1, 2)
      )
    },
    .package = "rGEDI"
  )
  entries <- data.frame(value = c(10, 20))
  result <- accumulateGEDIGranules(entries, parallel = FALSE)
  expect_equal(result$total, c(30, 60))
  expect_equal(result$shots, c(2, 4))
  expect_identical(attr(result, "workers"), 1L)
})
