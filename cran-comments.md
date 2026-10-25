## Resubmission

This is a resubmission. In response to the CRAN review:

* Corrected the README GitHub Actions link to reference the repository's
  actual `.github/workflows/r.yml` workflow file.
* Replaced README links written as relative file paths with absolute public
  GitHub URLs for the example-building scripts and generated reference manual.

* Expanded the `Description` field to explain the package functionality,
  implemented methods, supported GEDI products, and intended workflows.
* Added explicit value documentation to all usage-bearing help files,
  including `close.Rd`, and documented both the output class and meaning.
* Replaced every `\dontrun{}` example with an executable example or
  `\donttest{}` for operations that take longer than five seconds or download
  data. The download example only runs when NASA credentials are available.
* Removed unconditional informational `print()` calls. Functions now return
  usable objects silently, use suppressible messages when requested, or use
  console output only in registered print and plotting methods.

## R CMD check results

0 errors | 0 warnings | 1 note

The incoming-feasibility note records that rGEDI was archived on 2021-11-05
and identifies version 0.6.0 as a new submission. There are no code,
documentation, example, test, URL, or timing issues in the check result.

* This release adds support for GEDI Level 3, Level 4A, and Level 4B products,
  optional Earthdata Cloud and Google Earth Engine workflows, sampling and
  modeling helpers, an orbit animation, and a portable waveform simulator.
* Optional online integrations fail with an informative message when their
  suggested dependencies or credentials are unavailable.
* The package was tested with R 4.6.1 on Windows 11.

## Additional checks

All standard examples and testthat tests pass. Live NASA CMR discovery,
Earthdata authentication, temporary ORNL S3 credentials, lazy remote Level 1B,
Level 2A, Level 2B, and Level 4A HDF5 access, and the complete Google Earth
Engine sampling, modeling, mapping, and GeoTIFF download workflow were verified.
