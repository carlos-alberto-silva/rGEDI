test_that("GEDI product and Earth Engine catalogs expose modern products", {
  expect_equal(rGEDI:::.gedi_short_name("GEDI01_B", NULL)$version, "003")
  expect_equal(rGEDI:::.gedi_short_name("GEDI03", NULL)$version, "3")
  expect_equal(rGEDI:::.gedi_short_name("GEDI04_A", NULL)$version, "3")
  expect_equal(rGEDI:::.gedi_short_name("GEDI04_B", NULL)$version, "2.1")
  catalog <- search_datasets("GEDI", "biomass", operator = "and")
  expect_true(all(c("GEDI04_A", "GEDI04_B") %in% catalog$product))
  expect_equal(get_catalog_id("GEDI04_A"), "LARSE/GEDI/GEDI04_A_002")
})

test_that("cloud helpers validate unsupported local files and optional modules", {
  f <- tempfile(pattern = "unknown_", fileext = ".h5")
  writeBin(as.raw(1:4), f)
  expect_error(openGEDI(f), "Cannot infer")
  expect_error(earthdata_login(netrc = tempfile()), "path|configure|Install")
})

test_that("cloud product inference and URL subsets retain product semantics", {
  urls <- structure(
    c("https://example.test/GEDI04_A_example_V003.h5",
      "https://example.test/GEDI04_A_second_V003.h5"),
    class = c("gedi.granules_cloud", "character")
  )
  expect_equal(rGEDI:::.infer_gedi_product(urls[1]), "GEDI04_A")
  expect_s3_class(urls[1], "gedi.granule_cloud")
  expect_s3_class(urls[1:2], "gedi.granules_cloud")
})
