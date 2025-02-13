#!/bin/tcsh -f

# directory containing awk scripts needed
set bin="$GEDIRAT_ROOT/extractTracks"

# defaults
@ findLat=1
@ epsg=32732
@ seed=0
set orbitDir="/gpfs/data1/vclgp/data/gedi/ancillary/orbits/test"
set orbitAng=51
set cloudFrac=0.5
set latRes=5
set alongTrack=60
set res=1000
@ readBounds=1
@ readALS=0
@ readMetric=0
@ definedBounds=0
set minSep=30
set output="teast.dat"
@ maxLines=2000    # reduces swap space needed
set gridRes=2000
@ leafOn=0         # ignore phenology


# temporary workspace files
set tempTrack="/tmp/gediTrackSpace.$$.dat"
set workSpace="/tmp/trackWorkSpace.$$.dat"


# read command line
while ($#argv>0)
  switch("$argv[1]")

  case -metricFile
    set metricFile="$argv[2]"
    @ readMetric=1
  shift argv;shift argv
  breaksw

  case -alsFile
    set alsFile="$argv[2]"
    @ readMetric=0
    @ readALS=1
    @ readBounds=1
  shift argv;shift argv
  breaksw

  case -lat
    set meanLat=$argv[2]
    @ findLat=0
  shift argv;shift argv
  breaksw

  case -epsg
    @ epsg=$argv[2]
  shift argv;shift argv
  breaksw

  case -orbitDir
    set orbitDir="$argv[2]"
  shift argv;shift argv
  breaksw

  case -seed
    @ seed=$argv[2]
  shift argv;shift argv
  breaksw

  case -cloud
    set cloudFrac="$argv[2]"
  shift argv;shift argv
  breaksw

  case -output
    set output="$argv[2]"
  shift argv;shift argv
  breaksw

  case -minSep
    set minSep="$argv[2]"
  shift argv;shift argv
  breaksw

  case -gridRes
    set gridRes="$argv[2]"
  shift argv;shift argv
  breaksw

  case -leafOn
    @ leafOn=1
  shift argv
  breaksw

  case -bound
    set bounds=`echo "$argv[2] $argv[3] $argv[4] $argv[5]"|gawk '{for(i=1;i<=NF;i++)print $i}'`
    @ definedBounds=1
    @ readBounds=0
  shift argv;shift argv;shift argv;shift argv;shift argv
  breaksw

  case -help
    echo " "
    echo "-output name;      output filename"
    echo "-metricFile name;  input metric file, if sampling from metrics"
    echo "-alsFile name;     input ALS bound file, if setting coordinates"
    echo "-lat y;            latitude in degrees EPSG:4326"
    echo "-epsg n;           EPSG code of input data"
    echo "-minSep sep;       minimum horizontal separation to accept"
    echo "-orbitDir dir;     directory with GEDI track density files"
    echo "-cloud frac;       cloud fraction"
    echo "-seed n;           random number seed"
    echo "-gridRes r;        search grid resolution, metres"
    echo "-leafOn;           only include leaf on tracks"
    echo "-bound minX maxX minY maxY;   define bounds"
    echo " "
    exit

  default:
    echo "Unknown argument $argv[1]"
    echo "Type listGediWaves.csh -help"
    exit;

  endsw
end




# do we need to find the mean latitude
if( $definedBounds )then
  set meanLon=`echo "$bounds[1] $bounds[2]"|gawk '{print ($1+$2)/2}'`
  set meanLat=`echo "$bounds[3] $bounds[4]"|gawk '{print ($1+$2)/2}'`
  set meanLat=`echo $meanLon $meanLat|gdaltransform -s_srs EPSG:$epsg -t_srs EPSG:4326|gawk '{print $2}'`
else if( $findLat && $readMetric )then   # from metric fle
  set meanLat=`gawk -f $bin/meanCoord.awk < $metricFile|gdaltransform -s_srs EPSG:$epsg -t_srs EPSG:4326|gawk '{print $2}'`
else if( $findLat )then            # from ALS bounds file
  set meanLat=`gawk 'BEGIN{x=y=0;n=0}($0&&($1!="#")){x+=$2+$5;y+=$3+$6;n+=2}END{print x/n,y/n}' < $alsFile|gdaltransform -s_srs EPSG:$epsg -t_srs EPSG:4326|gawk '{print $2}'`
endif


# convert along and across track to angles
set resAng=`echo $res $meanLat|gawk '{pi=4*atan2(1,1);lat=$2*pi/180;print ($1/(cos(lat)*6371000))*180/pi}'`

# GEDI track angle at this point
set trackAngle=`echo $meanLat $orbitAng|gawk 'BEGIN{pi=4*atan2(1,1);scale=pi/180}{lat=$1*scale;a=$2*scale;print a*cos(lat*pi/(2*a))}'`

# idensitfy relevant track file
set minLat=`echo $meanLat|gawk -v res=$latRes '{lat=(int($1/res+0.5+100)-100)*res;printf("%f",lat)}'`
set trackRoot=`echo $minLat|gawk '{if($1>=0)printf("_lat%dN",$1);else printf("_lat%dS",-1*$1)}'`
set trackFile=`ls $orbitDir/*$trackRoot.txt`

# find bounds if needed
if( $readBounds && $readMetric )then # read from metric file
  set bounds=`gawk -f $bin/metricBounds.awk` < $metricFile
else if( $readBounds && $readALS )then # read from ALS file
  set bounds=`gawk -f $bin/alsBounds.awk` < $alsFile
endif  # otherwise they have already been defined

# translate bounds to same system as tracks
set min=`echo "$bounds[1] $bounds[3]"|gdaltransform -t_srs EPSG:4326 -s_srs EPSG:$epsg|gawk '{for(i=1;i<=NF;i++)print $i}'`
set max=`echo "$bounds[2] $bounds[4]"|gdaltransform -t_srs EPSG:4326 -s_srs EPSG:$epsg|gawk '{for(i=1;i<=NF;i++)print $i}'`
set minX=$min[1]
set maxX=$max[1]
set minY=$min[2]
set maxY=$max[2]

# determine track start locations
gawk -f $bin/trackStart.awk -v seed=$seed -v cloudFrac=$cloudFrac -v minX=$minX -v maxX=$maxX -v minY=$minY -v maxY=$maxY -v res=$resAng -v ang=$trackAngle -v usePhen=$leafOn < $trackFile > $workSpace
gawk '{print $1,$2}' < $workSpace|gdaltransform -s_srs EPSG:4326 -t_srs EPSG:$epsg|gawk '{print $1,$2}' > $workSpace.1
gawk '{print $3,$4}' < $workSpace> $workSpace.2
paste  $workSpace.1 $workSpace.2 > $workSpace
if( -e $workSpace.1 )rm $workSpace.1
if( -e $workSpace.2 )rm $workSpace.2
# place  footprints
gawk -f $bin/placeGediTracks.awk -v minX=$bounds[1] -v maxX=$bounds[2] -v minY=$bounds[3] -v maxY=$bounds[4] -v alongTrack=$alongTrack -v res=$resAng -v ang=$trackAngle -v usePhen=$leafOn < $workSpace > $tempTrack
if( -e $workSpace )rm $workSpace

# alter minsep if in degrees
if( $epsg == 4326 )then
  # reproject minSep
  set deg1=`echo 0 0|gdaltransform -t_srs EPSG:4326 -s_srs EPSG:32632|gawk '{print $1}'`
  set deg2=`echo $minSep 0|gdaltransform -t_srs EPSG:4326 -s_srs EPSG:32632|gawk '{print $1}'`
  set minSep=`echo $deg2 $deg1|gawk '{print $1-$2}'`
  # reproject grid resolution
  set deg1=`echo 0 0|gdaltransform -t_srs EPSG:4326 -s_srs EPSG:32632|gawk '{print $1}'`
  set deg2=`echo $gridRes 0|gdaltransform -t_srs EPSG:4326 -s_srs EPSG:32632|gawk '{print $1}'`
  set gridRes=`echo $deg2 $deg1|gawk '{print $1-$2}'`
endif

# output results
if( $readMetric )then    # extract tracks from metric file
  $bin/chooseMetricPrints -minSep $minSep -metric $metricFile -tracks $tempTrack -output $output -gridRes $gridRes
else if( $readALS || $definedBounds )then  # generate list of footprints to simulate
  gawk '{printf("%.10f %.10f %s\n",$1,$2,$3)}' < $tempTrack > $output
else
  echo "Then why did you run this?"
  exit(1)
endif


# tidy up workspace
if( -e $tempTrack )rm $tempTrack
if( -e $workSpace )rm $workSpace

echo "Written to $output"

