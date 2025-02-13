#!/bin/tcsh -f

@ nCores=1

# read command line
while ($#argv>0)
  switch("$argv[1]")

  case -progRoot
    set progRoot="$argv[2]"
  shift argv;shift argv
  breaksw

  case -nCores
    @ nCores=$argv[2]
  shift argv;shift argv
  breaksw

  case -help
    echo " "
    echo "-progRoot name;  progress filename root"
    echo "-nCores n;       number of cores to use"
    echo " "
    exit

  default:
    echo "Unknown argument $argv[1]"
    echo "A mistake has been made, by you."
    echo "Type chopPlots.csh -help"
    exit;

  endsw
end

@ notDone=1
while( $notDone )
  @ notDone=0
  @ i=0
  while( $i < $nCores )
    set progFile="$progRoot.$i"
    if( -e $progFile )@ notDone=1
    @ i++
  end
  sleep 5
end

