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
