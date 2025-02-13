#!/bin/tcsh -f

set inDir="./"
set output="teast.list"

# read options
while ($#argv>0)
  switch("$argv[1]")

  case -inDir
    set inDir="$argv[2]"
  shift argv;shift argv
  breaksw

  case -output
    set output="$argv[2]"
  shift argv;shift argv
  breaksw

  case -help
    echo " "
    echo "-inDir name;      directory containing ALS data"
    echo "-output name;     output list filename"
    echo " "
    exit

  default:
    echo "Unknown argument $argv[1]"
    echo "Who knows"
    echo "Type listALS.csh -help"
    exit;

  endsw
end

pushd $inDir/
set inDir=`pwd`
ls -l|sed -e s%"*"%""%g|gawk -F. '{if(($NF=="las")||($NF=="LAS"))print $0}'|gawk '{printf("%s/%s\n",dir,$NF)}' dir="$inDir" > $output
popd

echo "Written to $output"

