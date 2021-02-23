# Code inspired by rgdal
# Download GSL 2.4; hdf5; libgeotiff;
VERSION <- commandArgs(TRUE)
HDF_VERSION = VERSION[1]
GDAL_VERSION = VERSION[2]
GDAL_MAJOR = substr(GDAL_VERSION, 1, 1)
R_ARCH = Sys.getenv("R_ARCH")


# Download gdal
if (Sys.getenv("GDAL_HOME") == "") {
  testfile = sprintf("../windows/gdal%1$s-%2$s/include/gdal-%2$s/gdal.h", GDAL_MAJOR, GDAL_VERSION)
  download.link = sprintf("https://github.com/rwinlib/gdal%s/archive/v%s.zip", GDAL_MAJOR, GDAL_VERSION)
  if (Sys.getenv("LIB_GSL") == "") {
    if (!file.exists(testfile)) {
      print("Downloading and installing GDAL...")
      download.file(download.link, "lib.zip", quiet = FALSE)
      dir.create("../windows", showWarnings = FALSE)
      unzip("lib.zip", exdir = "../windows")
      unlink("lib.zip")
    }
  }
}
#
if (Sys.getenv("LIB_GSL") == "") {
  if (!file.exists("../windows/gsl-2.4/include/gsl/gsl_blas.h")) {
    print("Downloading and unpacking GSL...")
    download.file("https://github.com/rwinlib/gsl/archive/v2.4.zip", "lib.zip", quiet = TRUE)
    dir.create("../windows", showWarnings = FALSE)
    unzip("lib.zip", exdir = "../windows")
    unlink("lib.zip")
  }
}

if (!file.exists(sprintf("../windows/hdf5-%s/include/hdf5.h", HDF_VERSION))) {
  print("Downloading and unpacking HDF5...")
  download.file(sprintf("https://github.com/rwinlib/hdf5/archive/v%s.zip", HDF_VERSION),
                "lib.zip", quiet = TRUE)
  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
}

# if(!file.exists("../windows/libgeotiff-1.4.3/geo_config.h")){
#   print("Downloading and installing libgeotiff...")
#   download.file("https://github.com/OSGeo/libgeotiff/releases/download/1.4.3/libgeotiff-1.4.3.zip", "lib.zip", quiet = FALSE)
#   dir.create("../windows", showWarnings = FALSE)
#   unzip("lib.zip", exdir = "../windows")
#   unlink("lib.zip")
#   file.rename("../windows/libgeotiff-1.4.3/geo_config.h.vc", "../windows/libgeotiff-1.4.3/geo_config.h")
# }
#
#
# if(!file.exists("../windows/libtiff-4.0.9/include/tiffio.h")){
#   print("Downloading and installing libtiff...")
#   download.file("https://github.com/rwinlib/libtiff/archive/v4.0.9.zip", "lib.zip", quiet = FALSE)
#   dir.create("../windows", showWarnings = FALSE)
#   unzip("lib.zip", exdir = "../windows")
#   unlink("lib.zip")
# }
