
AC_INIT([rgdal], 1.6-7)

# Find pkg-config
PKG_PROG_PKG_CONFIG

# Check for GDAL
PKG_CHECK_MODULES([GDAL], [gdal], [
    PKG_LIBS="$GDAL_LIBS"
    PKG_CFLAGS="$GDAL_CFLAGS"
], [
    # Fallback if pkg-config fails
    AC_CHECK_PROG([GDAL_CONFIG], [gdal-config], [yes], [no])
    if test "$GDAL_CONFIG" = "yes"; then
        PKG_LIBS="`gdal-config --libs`"
        PKG_CFLAGS="`gdal-config --cflags`"
    else
        AC_MSG_ERROR([gdal-config not found in path
        Please install libgdal package:
            Debian: sudo apt install libgdal-dev
            Fedora: sudo dnf install gdal-devel
            Arch:   sudo pacman -S gdal
            Alpine: sudo apk add gdal-dev
            Mac:    brew install gdal])
    fi
])

# Check for geotiff
PKG_CHECK_MODULES([GEOTIFF], [libgeotiff], [
    PKG_LIBS="$PKG_LIBS $GEOTIFF_LIBS"
    PKG_CFLAGS="$PKG_CFLAGS $GEOTIFF_CFLAGS"
], [
    # Fallback if pkg-config fails
    AC_CHECK_PROG([GEOTIFF_CONFIG], [geotiff-config], [yes], [no])
    if test "$GEOTIFF_CONFIG" = "yes"; then
        PKG_LIBS="$PKG_LIBS `geotiff-config --libs`"
        PKG_CFLAGS="$PKG_CFLAGS `geotiff-config --cflags`"
    else
        AC_MSG_ERROR([geotiff-config not found in path
        Please install libgeotiff package:
            Debian: sudo apt install libgeotiff-dev
            Fedora: sudo dnf install libgeotiff-devel
            Arch:   sudo pacman -S libgeotiff
            Alpine: sudo apk add libgeotiff-dev
            Mac:    brew install libgeotiff])
    fi
])



# Check for HDF5
PKG_CHECK_MODULES([HDF5], [hdf5], [
    PKG_LIBS="$PKG_LIBS $HDF5_LIBS"
    PKG_CFLAGS="$PKG_CFLAGS $HDF5_CFLAGS"
], [
    # Fallback if pkg-config fails
    AC_CHECK_PROG([HDF5_CONFIG], [h5cc], [yes], [no])
    if test "$HDF5_CONFIG" = "yes"; then
        PKG_LIBS="$PKG_LIBS `h5cc -show | sed -e 's/.*\(.-L.*\)$/\1/'`"
        PKG_CFLAGS="$PKG_CFLAGS `h5cc -show | sed -e 's/\(.-I.[^ ]*\).*/\1/'`"
    else
        AC_MSG_ERROR([h5cc not found in path
        Please install hdf5 package:
            Debian: sudo apt install libhdf5-dev
            Fedora: sudo dnf install hdf5-devel
            Arch:   sudo pacman -S hdf5
            Alpine: sudo apk add hdf5-dev
            Mac:    brew install hdf5])
    fi
])


# Check for libtiff
PKG_CHECK_MODULES([TIFF], [libtiff-4], [
    PKG_LIBS="$PKG_LIBS $TIFF_LIBS"
    PKG_CFLAGS="$PKG_CFLAGS $TIFF_CFLAGS"
], [
    # Fallback if pkg-config fails
    AC_CHECK_PROG([TIFF_CONFIG], [tiff-config], [yes], [no])
    if test "$TIFF_CONFIG" = "yes"; then
        PKG_LIBS="$PKG_LIBS `tiff-config --libs`"
        PKG_CFLAGS="$PKG_CFLAGS `tiff-config --cflags`"
    else
        AC_MSG_ERROR([tiff-config not found in path
        Please install libtiff package:
            Debian: sudo apt install libtiff-dev
            Fedora: sudo dnf install libtiff-devel
            Arch:   sudo pacman -S libtiff
            Alpine: sudo apk add tiff-dev
            Mac:    brew install libtiff])
    fi
])


# Check for GSL
PKG_CHECK_MODULES([GSL], [gsl], [
    PKG_LIBS="$PKG_LIBS $GSL_LIBS"
    PKG_CFLAGS="$PKG_CFLAGS $GSL_CFLAGS"
], [
    # Fallback if pkg-config fails
    AC_CHECK_PROG([GSL_CONFIG], [gsl-config], [yes], [no])
    if test "$GSL_CONFIG" = "yes"; then
        PKG_LIBS="$PKG_LIBS `gsl-config --libs`"
        PKG_CFLAGS="$PKG_CFLAGS `gsl-config --cflags`"
    else
        AC_MSG_ERROR([gsl-config not found in path
        Please install libgsl package:
            Debian: sudo apt install libgsl-dev
            Fedora: sudo dnf install gsl-devel
            Arch:   sudo pacman -S gsl
            Alpine: sudo apk add gsl-dev
            Mac:    brew install gsl])
    fi
])

AC_SUBST([PKG_LIBS])
AC_SUBST([PKG_CFLAGS])
AC_CONFIG_FILES([
    src/Makefile
])
AC_OUTPUT