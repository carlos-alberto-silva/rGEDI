#!/bin/tcsh -f

#########################
# Runs GEDI simulations #
# over a list of coords #
#########################


set temp="/tmp/boundFiles.$$.dat"

# default settings
set pBuff=0.5
set outRoot="teast"
set LVIS=" "
set pSigma=" "
set pFWHM=" "
set fSigma=" "
set ground=" "
set sideLobe=" "
set lobeAng=" "
set topHat=" "
set noNorm=" "
set checkCove=" "
set maxScanAng=" "
set polyGround=" "
set pFile=" "
set res=" "
set hdf=" "
set l1b=" "
set aEPSG=" "
@ useHDF=0
@ maxPer=40000
set ending="wave"
set wavefront=" "
set octree=" "
set octLevels=" "
set nOctPix=" "
set countOnly=" "
set pulseAfter=" "
set decimate=" "
set seed=" "
set pcl=" "


# read options
while ($#argv>0)
  switch("$argv[1]")

  case -inList
    set inList="$argv[2]"
  shift argv;shift argv
  breaksw

  case -outRoot
    set outRoot="$argv[2]"
  shift argv;shift argv
  breaksw

  case -coordList
    set coordList="$argv[2]"
    @ nCoords=`wc -l < $coordList`
  shift argv;shift argv
  breaksw

  case -coord
    set coords="$argv[2] $argv[3]"
    @ nCoords=1
  shift argv;shift argv;shift argv
  breaksw

  case -pBuff
    set pBuff="$argv[2]"
  shift argv;shift argv
  breaksw

  case -LVIS
    set LVIS="-LVIS"
  shift argv
  breaksw

  case -pSigma
    set pSigma="-pSigma $argv[2]"
  shift argv;shift argv
  breaksw

  case -pFWHM
    set pFWHM="-pFWHM $argv[2]"
  shift argv;shift argv
  breaksw

  case -fSigma
    set fSigma="-fSigma $argv[2]"
  shift argv;shift argv
  breaksw

  case -ground
    set ground="-ground"
  shift argv
  breaksw

  case -sideLobe
    set sideLobe="-sideLobe"
  shift argv
  breaksw

  case -lobeAng
    set lobeAng="-lobeAng $argv[2]"
  shift argv;shift argv
  breaksw

  case -topHat 
     set topHat="-topHat"
  shift argv
  breaksw

  case -noNorm
    set noNorm="-noNorm"
  shift argv
  breaksw

  case -checkCover
    set checkCover="-checkCover"
  shift argv
  breaksw

  case -maxScanAng
    set maxScanAng="-maxScanAng $argv[2]"
  shift argv;shift argv
  breaksw

  case -pFile
    set pFile="-readPulse $argv[2]"
  shift argv;shift argv
  breaksw

  case -res
    set res="-res $argv[2]"
  shift argv;shift argv
  breaksw

  case -polyGround
    set polyGround="-polyGround"
  shift argv
  breaksw

  case -hdf
    set hdf="-hdf"
    @ useHDF=1
  shift argv
  breaksw

  case -l1b
    set l1b="-l1b"
    @ useHDF=1
  shift argv
  breaksw

  case -aEPSG
    set aEPSG="-aEPSG $argv[2]"
  shift argv;shift argv
  breaksw

  case -maxPer
    @ maxPer=$argv[2]
  shift argv;shift argv
  breaksw

  case -wavefront
    set wavefront="-wavefront $argv[2]"
  shift argv;shift argv
  breaksw

  case -noOctree
    set octree="-noOctree"
  shift argv
  breaksw

  case -octLevels
    set octLevels="-octLevels $argv[2]"
  shift argv;shift argv
  breaksw

  case -nOctPix
    set nOctPix="-nOctPix $argv[2]"
  shift argv;shift argv
  breaksw

  case -countOnly
    set countOnly="-countOnly"
  shift argv
  breaksw

  case -pulseAfter
    set pulseAfter="-pulseAfter"
  shift argv
  breaksw

  case -decimate
    set decimate="-decimate $argv[2]"
  shift argv;shift argv
  breaksw

  case -seed
    set seed="-seed $argv[2]"
  shift argv;shift argv
  breaksw

  case -pcl
    set pcl="-pcl"
  shift argv
  breaksw

  case -help
    echo " "
    echo "-inList name;      name of list with las file names"
    echo "-outRoot root;     output filename root"
    echo "-coordList name;   file containing list of coordinates and waveIDs (x y waveID)"
    echo "-coord x y;        coordinate of a single footprint to simulate"
    echo "-pBuff size;       RAM buffer size"
    echo "-LVIS;             use LVIS pulse size"
    echo "-pSigma sigma;     pulse width, sigma in m"
    echo "-pFWHM fwhm;       pulse FWHM in ns"
    echo "-pFile name;       read the pulse from an ASCII file"
    echo "-pulseAfter;       apply pulse shape after binning for speed"
    echo "-wavefront file;   read wavefront shape for non-Gaussian footprints"
    echo "-fSigma sigma;     footprint width, sigma in m"
    echo "-res res;          output range resolution, in metres"
    echo "-ground;           output ground and canopy waveforms"
    echo "-countOnly;        only simulate using the count method"
    echo "-sideLobe;         use the old side-lobes"
    echo "-lobeAng ang;      side lobe major axis azimuth, degrees"
    echo "-topHat;           use top hat rather than Gaussian footprint"
    echo "-noNorm;           do not normalise for footprint density"
    echo "-checkCover;       check that at least 2/3 of footprint is covered by ALS"
    echo "-maxScanAng ang;   maximimum scan angle to use, degrees"
    echo "-polyGround;       find the ground by fitting polynomial"
    echo "-hdf;              output in HDF5"
    echo "-l1b;              output in L1B HDF5 format"
    echo "-aEPSG epsg;       EPSG code of ALS data if reprojecting for L1B HDF5 files"
    echo "-maxPer n;         maximum number of runs per processor"
    echo "-decimate x;       decimate ALS beams by a factor"
    echo "-seed n;           random number seed"
    echo "-pcl;              will be PCL, so don't pad pulse"
    echo "# octree"
    echo "-noOctree;         do not use an octree"
    echo "-octLevels n;      number of octree levels to use"
    echo "-nOctPix n;        number of octree pixels along a side for the top level"
    echo " "
    exit

  default:
    echo "Unknown argument $argv[1]"
    echo "Who knows"
    echo "Type gediRatList.csh -help"
    exit;

  endsw
end

set grabDir="butabe"
if( ! -e $grabDir )mkdir $grabDir


# split into sub files
set tempRoot="/tmp/gediRatList.$$"
@ nReps=`echo "$nCoords $maxPer"|gawk '{print int($1/$2+1)}'`

# if in single footprint or group of footprint mode
if( $nCoords > 1 )then
  gawk -v maxPer=$maxPer -f $GEDIRAT_ROOT/awk/splitCoords.awk root="$tempRoot" < $coordList
else
  echo "$coords" > $tempRoot.0.coords
endif

@ j=0
while( $j <= $nReps )
  set input="$tempRoot.$j.coords"
  if( ! -e $input )then
    @ j++
    continue
  endif

  set output="$outRoot.$j"
  if( $useHDF )set output="$output.h5"
  set grab="$grabDir/$output.grab"

  if( ! -e $grab )then
    touch $grab
    overlapLasFiles.csh -input $inList -coordList $input -rad 100 -output $temp
    gediRat -inList $temp -output $output -listCoords $input $input -pBuff $pBuff $LVIS $pSigma $pFWHM $fSigma $ground $sideLobe $lobeAng $topHat $noNorm $checkCove $maxScanAng $pFile $res $polyGround $hdf $l1b $aEPSG $wavefront $octree $octLevels $nOctPix $countOnly $pulseAfter $decimate $seed $pcl

  endif

  # delete zero sized files
  if( -e $temp )rm $temp
  if( -e $input )rm $input
  @ j++
end

echo "Ping"

