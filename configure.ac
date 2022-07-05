## Process this file with autoconf to produce a configure script.
m4_include([inst/m4/gsl.m4])
m4_include([inst/m4/m4_ax_lib_gdal.m4])
m4_include([inst/m4/zlib.m4])


# The version set here will propagate to other files from here
AC_INIT([gedisimulator],[0.0.1])

# Checks for common programs using default macros
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`
AC_PROG_CC


##################
## GSL
##################
AX_PATH_GSL(0.1,,
	    AC_MSG_ERROR(could not find GSL:
			 Debian: sudo apt install libgsl-dev
			 RPM-based: sudo dnf install gsl-devel
			 MacOS: brew install gsl))


##################
## HDF5
##################
AX_LIB_HDF5()

if test "$with_hdf5" != "yes"; then
    AC_MSG_ERROR(could not find hdf5:
                         Debian: sudo apt install libhdf5-dev
                         RPM-based: sudo dnf install hdf5-devel
                         MacOS: brew install hdf5)
fi

##################
## SZIP
##################
AC_CHECK_LIB( sz, SZ_BufftoBuffCompress,
	     HAVE_LIBSZ="yes",
	     AC_MSG_ERROR(could not find szip:
                         Debian: sudo apt install libsz2
                         RPM-based: sudo dnf install szip
                         MacOS: brew install szip))




##################
## LIBZ
##################
AX_CHECK_ZLIB(,
	      AC_MSG_ERROR(could not find zlib:
                         Debian: sudo apt install zlib1g-dev
                         RPM-based: sudo dnf install zlib-devel
                         MacOS: brew install zlib))


##################
## END
##################


AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

