<!-- NEWS.md is maintained by https://cynkra.github.io/fledge, do not edit -->

# rGEDI 0.6.0 (2026-10-09)

* Add Level 3, Level 4A, and Level 4B readers and spatial processing tools
* Add Earthdata Cloud and Google Earth Engine integration
* Stream GEDI01_B, GEDI02_A, GEDI02_B, and GEDI04_A HDF5 granules through
  Earthaccess with the same typed readers used for local files
* Update Level 2A and Level 4A extraction for the current Release 3 quality
  fields and add beam/column controls for efficient cloud reads
* Resolve Earth Engine vector products through their granule indexes and add a
  complete real-data README workflow with reproducible figures and animation
* Add spatial sampling, modeling, prediction, and mapping helpers
* Add an interactive GEDI orbit animation
* Add static and GIF orbit products, recursive feature elimination diagnostics,
  Earth Engine Drive export, and a fully ordered end-to-end README workflow
* Restore waveform simulation and waveform metrics with a portable R/HDF5
  implementation aligned with the focused Gaussian-footprint workflow in
  Steven Hancock's simulator while avoiding its non-portable native library
  stack


# rGEDI 0.5.8 (2026-10-09)

* Update the designated maintainer and contact email
* Mark the long-running grid statistics example for optional example checks


# rGEDI 0.5.7 (2026-10-09)

* Fixes #68: direct call to S3 method from bit64
* Add testthat automatic testing
* Update url links from README.md and add tests
* Fix grid statistics functions after sf stopped exporting the pipe operator
* Update the GEDI data products URL after its permanent redirect
* Keep optional visualizations out of non-interactive example checks
* Make URL tests tolerate transient network and server failures


# rGEDI 0.5.6 (2025-09-23)

* Update url links from README.md and add tests


# rGEDI 0.5.5 (2025-09-23)

* Removed unused script tools
* Fixed redirected urls
* Added testthat for testing urls
* Fix minor problems with checks:


# rGEDI 0.5.4 (2025-09-23)

* Fixed redirected urls
* Added testthat for testing urls
* Fix minor problems with checks:


# rGEDI 0.5.3 (2025-09-23)

* Fix minor problems with checks:


# rGEDI 0.5.2 (2025-09-22)

* Update documentation with tests


# rGEDI 0.5.1.9000 (2025-09-22)

- Same as previous version.


# rGEDI 0.5.1 (2025-09-22)

* Update polyStatsLevel2AM.R: fix example
* Update gedifinder.R to new NASA EarthCloud IDs
* Fix incompatibilities with sp and raster
* Use sf instead of raster::shapefile
* Replaced every usage of raster with stars
* ClipLevel1B use only sf and terra
* Split gedisimulator apart
* Update gedifinder to find level3 and level4 products
* Use user defined functions from parent env instead global


# rGEDI 0.5.0 (2023-10-31)

* Replaced every usage of raster with stars
* Split gedisimulator apart
* Update gedifinder to find level3 and level4 products
* Use user defined functions from parent env instead global


# rGEDI 0.4.0 (2023-06-15)

* gediMetrics: now working again, had a bug, allow beam filter
* gediSimulator: allow ascii output
* Update Hancocks's gedisimulator
* Makevars.ucrt: remove GDAL from requirements


# rGEDI 0.3.1 (2022-08-31)

* Use Makevars.ucrt for R 4.2
* Remove @bbox from lidR: deprecated slot
* Add Makevars.ucrt for the new Rtools42
* Use markdown in functions documentation


# rGEDI 0.3.0 (2021-08-20)

* Add inst/proj and inst/gdal to .gitignore
* Use \href for external links
* update links to v002
* Use \doi macro for referencing


# rGEDI 0.2.1 (2021-07-02)

* gedifinder: use GEDI version="002" in examples


# rGEDI 0.2.0 (2021-07-02)

* Update gedifinder to use CMR for fetching v2
- Same as previous version.
