# Download GSL 2.4; hdf5; libgeotiff;

# Download gdal
if(!file.exists("../windows/gdal2-2.2.3/include/gdal/gdal.h")){
  print("Downloading and installing GDAL...")
  download.file("https://github.com/rwinlib/gdal2/archive/v2.2.3.zip", "lib.zip", quiet = FALSE)
  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
}

if(!file.exists("../windows/gsl-2.4/include/gsl/gsl_blas.h")){
  print("Downloading and installing GSL...")
  download.file("https://github.com/rwinlib/gsl/archive/v2.4.zip", "lib.zip", quiet = FALSE)
  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
}

if(!file.exists("../windows/mingw64-libhdf5-dev-1.8.20/include/hdf5.h")){
  print("Downloading and installing HDF5...")
  download.file("https://github.com/caiohamamura/mingw64-libhdf5-dev/archive/v1.8.20.zip", "lib.zip", quiet = FALSE)
  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
}

if(!file.exists("../windows/libgeotiff-1.4.3/geo_config.h")){
  print("Downloading and installing libgeotiff...")
  download.file("https://github.com/OSGeo/libgeotiff/releases/download/1.4.3/libgeotiff-1.4.3.zip", "lib.zip", quiet = FALSE)
  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
  file.rename("../windows/libgeotiff-1.4.3/geo_config.h.vc", "../windows/libgeotiff-1.4.3/geo_config.h")
}


if(!file.exists("../windows/libtiff-4.0.9/include/tiffio.h")){
  print("Downloading and installing libtiff...")
  download.file("https://github.com/rwinlib/libtiff/archive/v4.0.9.zip", "lib.zip", quiet = FALSE)
  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
}
