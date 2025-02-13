#!/bin/bash -f


######################
# Scripts to install #
# GEDI software and  #
# set up environemnt #
# This will modify   #
# your .bashrc file  #
######################


HOMDIR="$HOME"


# set up environment variables
export ARCH=`uname -m`
export PATH=$PATH:./:$HOMDIR/bin/$ARCH:$HOMDIR/bin/csh
export GEDIRAT_ROOT=$HOMDIR/src/gedisimulator
export CMPFIT_ROOT=$HOMDIR/src/cmpfit-1.2
export GSL_ROOT=/usr/local/lib
export LIBCLIDAR_ROOT=$HOMDIR/src/libclidar
export HANCOCKTOOLS_ROOT=$HOMDIR/src/tools
export HDF5_LIB=/apps/hdf5/1.8.15/patch1

envFile="$HOMDIR/.bashrc"
echo "export ARCH=`uname -m`" >> $envFile
echo "export PATH=$PATH:./:$HOMDIR/bin/$ARCH:$HOMDIR/bin/csh" >> $envFile
echo "export GEDIRAT_ROOT=$HOMDIR/src/gedisimulator" >> $envFile
echo "export CMPFIT_ROOT=$HOMDIR/src/cmpfit-1.2" >> $envFile
echo "export GSL_ROOT=/usr/local/lib" >> $envFile
echo "export LIBCLIDAR_ROOT=$HOMDIR/src/libclidar" >> $envFile
echo "export HANCOCKTOOLS_ROOT=$HOMDIR/src/tools" >> $envFile
echo "export HDF5_LIB=/apps/hdf5/1.8.15/patch1" >> $envFile


# set up directory structure
if [ ! -e $HOMDIR/src ];then
  mkdir $HOMDIR/src
fi
if [ ! -e $HOMDIR/bin ];then
  mkdir $HOMDIR/bin
fi
if [ ! -e $HOMDIR/bin/$ARCH ];then
  mkdir $HOMDIR/bin/$ARCH
fi
if [ ! -e $HOMDIR/bin/csh ];then
  mkdir $HOMDIR/bin/csh
fi

pushd $HOMDIR/src
wget https://www.physics.wisc.edu/~craigm/idl/down/cmpfit-1.2.tar.gz
tar -xvf cmpfit-1.2.tar.gz
popd

pushd $HOMDIR/src
git clone https://bitbucket.org/StevenHancock/libclidar
git clone https://bitbucket.org/StevenHancock/tools
git clone https://bitbucket.org/StevenHancock/gedisimulator


programList="gediRat gediMetric mapLidar collocateWaves lasPoints fitTXpulse"
cd $GEDIRAT_ROOT/
make clean

for program in $programList;do
  make THIS=$program
  make THIS=$program install
done

programList="gediRatList.csh listGediWaves.csh overlapLasFiles.csh filtForR.csh"
for program in $cshList;do
  cp $program $HOME/bin/csh/
done

#cp *.csh $HOMDIR/src/csh/
#cp *.bash $HOMDIR/src/csh/

popd

