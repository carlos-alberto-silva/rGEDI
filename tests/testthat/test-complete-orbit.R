test_that("gedifinder retrieves every granule for an orbit", {
  seen_url <- NULL
  local_mocked_bindings(
    .cmr_get = function(url) {
      seen_url <<- url
      entry <- function(part) list(
        collection_concept_id = "C000000-LPCLOUD",
        producer_granule_id = paste0(
          "GEDI02_A_2019108080339_O01964_", part,
          "_T05337_02_003_01_V002"
        ),
        time_start = "2019-04-18T08:03:39Z",
        time_end = "2019-04-18T09:25:12Z",
        links = list(list(href = paste0("https://example.test/orbit_", part, ".h5")))
      )
      list(feed = list(entry = lapply(c("04", "02", "01", "03"), entry)))
    },
    .package = "rGEDI"
  )

  result <- gedifinder(
    "GEDI02_A", version = "002", orbit = "O01964", return = "table"
  )
  expect_equal(result$granule_part, 1:4)
  expect_true(all(result$orbit == "01964"))
  expect_match(seen_url, "orbit_number=1964", fixed = TRUE)
  expect_false(grepl("bounding_box", seen_url, fixed = TRUE))
})

test_that("gedifinder validates spatial and orbit searches", {
  expect_error(gedifinder("GEDI02_A"), "bounding coordinates")
  expect_error(gedifinder("GEDI02_A", orbit = "invalid"), "positive orbit")
  expect_error(gedifinder("GEDI02_A", orbit = c(1, 2)), "one orbit")
})

test_that("gediDownload returns existing downloaded paths", {
  outdir <- tempfile("rgedi-download-")
  dir.create(outdir)
  target <- file.path(outdir, "example.h5")
  writeBin(as.raw(1:3), target)
  netrc <- tempfile("rgedi-netrc-")
  writeLines(c("machine urs.earthdata.nasa.gov", "login test", "password test"), netrc)

  result <- suppressMessages(gediDownload(
    "https://example.test/example.h5?token=not-used",
    outdir = outdir, netrc = netrc
  ))
  expect_equal(result, normalizePath(target, winslash = "/"))
  expect_error(gediDownload(character(), outdir = outdir, netrc = netrc),
               "at least one")
})
