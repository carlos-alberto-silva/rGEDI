# Makefile for GEDI simulator tools
LIBS = -lm -lgdal -lgeotiff 
EXE=$(THIS)
ARCH=$(shell uname -m)
LINKER=
COPY=cp
THIS=gediRat
OUTPUT=/usr/bin/$(THIS)

ifeq ("$(OS)","Windows_NT")
	LIBS = -lgdal -lgeotiff -lgdal -lsqlite3 -lspatialite -lproj -lgeos_c -lgeos  \
		-ljson-c -lnetcdf -lmariadbclient -lpq -lintl -lwebp -lcurl -lssh2 -lssl -lcrypto \
		-lkea -lhdf5_cpp -lhdf5_hl -lhdf5 -lexpat -lfreexl -lcfitsio \
		-lmfhdf -ldf -lxdr \
		-lopenjp2 -ljasper -lpng16 -ljpeg -ltiff -lgeotiff -lgif -lxml2 -llzma -lszip -lz \
		-lodbc32 -lodbccp32 -liconv -lpsapi -lws2_32 -lcrypt32 -lwldap32 -lsecur32 -lgdi32 \
		-lgsl -lgslcblas
	CFLAGS += -D_USE_MATH_DEFINES -D_WIN32 -DWIN32
	EXE=$(THIS).exe
	ARCH=$(PROCESSOR_ARCHITECTURE)
	LINKER=$(where link)
	COPY=copy /Y
	OUTPUT = $(CONDA_PREFIX)\Library\bin\$(THIS)$(BIN)
endif


LIBS += -lgsl -lgslcblas -ltiff -lhdf5 -L${GSL_ROOT} -L${HDF5_LIB}/lib \
	-L"C:\Users\caioh\lsrc\rLVIS\windows\gdal2-2.2.3\lib\x64" \
		-L"C:\Users\caioh\lsrc\rLVIS\windows\gsl-2.4\lib\x64" 
INCLS = -I/usr/local/include -I$(HANCOCKTOOLS_ROOT) -I$(CMPFIT_ROOT) -I${LIBCLIDAR_ROOT} -I. -I/usr/include/libgeotiff -I/usr/include/gdal \
		-IC:/Users/caioh/Documents/R/win-library/3.5/Rcpp/include \
		-IC:/Users/caioh/lsrc/rLVIS/windows/mingw64-libhdf5-dev-1.8.20/include \
		-I. \
		-I"C:\Users\caioh\lsrc\rLVIS\windows\gdal2-2.2.3\include\gdal" \
		-I"C:\Users\caioh\lsrc\rLVIS\windows\gsl-2.4\include" \
		-I"C:\Users\caioh\lsrc\rLVIS\windows\libgeotiff-1.4.3" \
		-I"C:\Users\caioh\lsrc\rLVIS\windows\libtiff-4.0.9\include" \
		-I"C:\Users\caioh\lsrc\rLVIS\windows\libgeotiff-1.4.3\libxtiff" \
		-I"C:\Program Files\R\R-3.5.3\include" \
		-I"C:\Program Files\R\R-3.5.3\include\R_ext" \
		-fPIC

CFLAGS += $(INCLS)
LIBFILES = $(HANCOCKTOOLS_ROOT)/msgHandling.o \
		   $(LIBCLIDAR_ROOT)/libLasProcess.o \
		   $(LIBCLIDAR_ROOT)/libLasRead.o \
		   $(LIBCLIDAR_ROOT)/tiffWrite.o \
		   $(LIBCLIDAR_ROOT)/gaussFit.o \
		   $(LIBCLIDAR_ROOT)/libLidVoxel.o \
		   $(LIBCLIDAR_ROOT)/libTLSread.o \
		   $(LIBCLIDAR_ROOT)/libLidarHDF.o \
		   gedisimulator/gediIO.o \
		   $(LIBCLIDAR_ROOT)/libOctree.o \
		   gedisimulator/gediNoise.o \
		   $(CMPFIT_ROOT)/mpfit.o \
		   $(HANCOCKTOOLS_ROOT)/tools.o \
		   gedisimulator/gediRat.o \
		   gediRat.o
CPPFLAGS = $(CFLAGS)

ARCH=$(shell uname -m)
obj = $(LIBFILES)
	  

CC=C:/msys64/mingw64/bin/gcc.exe
CXX=C:/msys64/mingw64/bin/g++.exe

debug: CFLAGS += -Wall -g -D_DEBUG -DDEBUG -DDLLEXPORT -D__USE_MINGW_ANSI_STDIO
debug: $(EXE)

release: CFLAGS += -O3
debug: $(EXE)

$(EXE): $(THIS).o $(obj)
	$(CXX) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
	rm -f *% *~ *.o $(obj)

install: $(THIS)$(BIN)
	$(COPY) "$(THIS)$(BIN)" "$(OUTPUT)"

.PHONY: clean install