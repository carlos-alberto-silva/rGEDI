## Process this file with autoconf to produce a configure script.
m4_include([inst/m4/gsl.m4])
m4_include([inst/m4/m4_ax_lib_gdal.m4])


# The version set here will propagate to other files from here
AC_INIT([gedisimulator], 0.0.1)

# Checks for common programs using default macros
AC_PROG_CC


##################
## GSL
##################
AX_PATH_GSL(0.1,,
	    AC_MSG_ERROR(could not find GSL:
			 Debian: sudo apt install libgsl-dev
			 RPM-based: sudo yum install gsl-devel
			 MacOS: brew install gsl))


##################
## HDF5
##################
AX_LIB_HDF5()

if test "$HAVE_HDF5" = "no" ; then
	AC_MSG_ERROR(could not find hdf5:
                         Debian: sudo apt install libhdf5-dev
                         RPM-based: sudo yum install hdf5-devel
                         MacOS: brew install hdf5)
fi

##################
## GEOTIFF
##################

GEOTIFF_INCLUDE=

dnl Now seach for headers
if test -r /usr/include/geotiff.h ; then
	GEOTIFF_INCLUDE=
dnl Debian
elif test -r /usr/include/geotiff/geotiff.h ; then
	GEOTIFF_INCLUDE="-I/usr/include/geotiff"
dnl Fedore and OpenSuse
elif test -r /usr/include/libgeotiff/geotiff.h ; then
	GEOTIFF_INCLUDE="-I/usr/include/libgeotiff"
else
	AC_CHECK_HEADERS([geotiff.h])
	if test "$ac_cv_header_geotiff_h" = "no" ; then
		AC_MSG_ERROR([cannot find geotiff.h])
	fi
fi

GEOTIFF_LIBS="-lgeotiff $LIBS"


AC_SUBST(GEOTIFF_INCLUDE,$GEOTIFF_INCLUDE)
AC_SUBST(GEOTIFF_LIBS, $GEOTIFF_LIBS)


##################
## SZIP
##################
AC_CHECK_LIB( sz, SZ_BufftoBuffCompress, 
	     HAVE_LIBSZ="yes",
	     AC_MSG_ERROR(could not find szip:
                         Debian: sudo apt install libsz2
                         RPM-based: sudo yum install szip
                         MacOS: brew install szip))




##################
## LIBZ
##################
AX_CHECK_ZLIB(,
	      AC_MSG_ERROR(could not find libz:
                         Debian: sudo apt install zlib1g-dev
                         RPM-based: sudo yum install zlib-devel
                         MacOS: brew install zlib))


##################
## END
##################


AC_OUTPUT(src/Makevars)
