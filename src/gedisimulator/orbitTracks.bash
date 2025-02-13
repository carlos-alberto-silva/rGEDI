#!/bin/bash -f

# directory containing awk scripts needed
bin="$GEDIRAT_ROOT/extractTracks"
alongTrack=60
readALS=0
seed=0
cloudFrac=0
bounds=`echo "-180 180 -51 51"|gawk '{for(i=1;i<=NF;i++)print $i}'`
powList="1 2 3 4"
skipList="9 10"
output="teast.coords"
mission_length=395
orbit_sim_dir="/gpfs/data1/vclgp/htang/GEDI/ISS/BaselineSIM"
gridRes=1000
usePhen="--usePhen 1"
useWeak=" "   # do not use weak day by default
pointErr=0

# read command line
while [[ $# -gt 0 ]]
do
key="$1"
  case $key in
    -alsFile)
      alsFile="$2"
      readALS=1
      shift # past argument
      shift # past value
      ;;
    -missLen)
      mission_length="$2"
      shift # past argument
      shift # past value
      ;;
    -noPhen)
      usePhen="--usePhen 0"
      shift # past argument
      ;;
    -useWeak)
      useWeak="--useWeak 1"
      shift # past argument
      ;;
    -epsg)
      epsg="$2"
      shift # past argument
      shift # past value
      ;;
    -seed)
      seed="$2"
      shift # past argument
      shift # past value
      ;;
    -cloud)
      cloudFrac="$2"
      shift # past argument
      shift # past value
      ;;
    -output)
      output="$2"
      shift # past argument
      shift # past value
      ;;
    -gridRes)
      gridRes="$2"
      shift # past argument
      shift # past value
      ;;
    -bound)
      bounds=`echo "$2 $3 $4 $5"|gawk '{for(i=1;i<=NF;i++)print $i}'`
      readALS=0
      shift;shift;shift;shift;shift
      ;;
    -orbit_sim_dir)
      orbit_sim_dir="$2"
      shift # past argument
      shift # past value
      ;;
    -pointErr)
      pointErr="$2"
      shift # past argument
      shift # past value
      ;;
    -help)
      echo " "
      echo "-output name;        output filename"
      echo "-orbit_sim_dir name; directory containing GEDI tracks"
      echo "-alsFile name;       input ALS bound file, if setting coordinates"
      echo "-epsg n;             EPSG code of input data"
      echo "-cloud frac;         cloud fraction"
      echo "-seed n;             random number seed"
      echo "-missLen len;        mission length in days, accounting for allocation"
      echo "-bound minX maxX minY maxY;   define bounds"
      echo "-noPhen;             use leaf-off data"
      echo "-useWeak;            use daytime weak beam"
      echo "-pointErr sig;       pointing error, 1 sigma in metres"
      echo "-gridRes res;        grid resolution, in EPSG units"
      echo " "
      exit
      ;;
    *)
      echo $"Unrecognised option: $key"
      exit 1
  esac
done


# determine bounds
if [ $readALS == 1 ]; then
  echo "Reading ALS $readALS"
  bounds=`gawk -f $bin/alsBounds.awk < $alsFile`
fi

# read data
temp="/tmp/tempCoords.$$.dat"
python $bin/orbitTracks.py --bounds $bounds --oEPSG $epsg --output $temp --seed $seed --mission_length $mission_length --skip_track $skipList --power_beams $powList --cloud $cloudFrac --orbit_sim_dir $orbit_sim_dir $usePhen $useWeak --pointErr $pointErr


# collate into patches
if [ $res > 0 ];then
  python $bin/gridPrints.py --input $temp --output $output --res $gridRes
  rm $temp
else
  mv $temp $output
fi

