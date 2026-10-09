## R CMD check results

0 errors | 0 warnings | 0 notes

* This release adds support for GEDI Level 3, Level 4A, and Level 4B products,
  optional Earthdata Cloud and Google Earth Engine workflows, sampling and
  modeling helpers, an orbit animation, and a portable waveform simulator.
* Optional online integrations fail with an informative message when their
  suggested dependencies or credentials are unavailable.
* The package was tested with R 4.6.1 on Windows 11.

## Additional checks

All standard examples and testthat tests pass. Live NASA CMR discovery,
Earthdata authentication, temporary ORNL S3 credentials, and lazy remote
Level 4A HDF5 access were also verified.
