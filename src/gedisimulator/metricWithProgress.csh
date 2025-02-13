#!/bin/tcsh -f

########################
# Runs gediMetric with #
# progress check       #
########################

# defaults
set outRoot="teast"
set writeFit=" "
set ground=" "
set useInt=" "
set useFrac=" "
set rhRes=5
set bayesGround=" "
set gTol=" "
set noRHgauss=" "
set dcBias=" "
set nSig=" "
set seed=" "
set hNoise=" "
set linkNoise=" "
set trueSig=" "
set renoise=" "
set newPsig=" "
set oldPsig=" "
set missGround=" "
set minGap=" "
set maxDN=" "
set bitRate=" "
set meanN=" "
set thresh=" "
set sWidth=" "
set psWidth=" "
set gWidth=" "
set minWidth=" "
set varNoise=" "
set varScale=" "
set statsLen=" "
set medNoise=" "
set noiseTrack=" "
set rhoG=" "
set rhoC=" "
set offset=" "
set dontTrustGround=" "

# read command line
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

  case -writeFit
    set writeFit="-writeFit"
  shift argv
  breaksw

  case -ground
    set ground="-ground"
  shift argv
  breaksw

  case -useInt
    set useInt="-useInt"
  shift argv
  breaksw

  case -useFrac
    set useFrac="-useFrac"
  shift argv
  breaksw

  case -rhRes
    set rhRes="$argv[2]"
  shift argv;shift argv
  breaksw

  case -bayesGround
    set bayesGround="-bayesGround"
  shift argv
  breaksw

  case -gTol
    set gTol="-gTol $argv[2]"
  shift argv;shift argv
  breaksw

  case -noRHgauss
    set noRHgauss="-noRHgauss"
  shift argv
  breaksw

  case -dcBias
    set dcBias="-dcBias $argv[2]"
  shift argv;shift argv
  breaksw

  case -nSig
    set nSig="-nSig $argv[2]"
  shift argv;shift argv
  breaksw

  case -seed
    set seed="-seed $argv[2]"
  shift argv;shift argv
  breaksw

  case -hNoise
    set hNoise="-hNoise $argv[2]"
  shift argv;shift argv
  breaksw

  case -linkNoise
    set linkNoise="-linkNoise $argv[2] $argv[3]"
  shift argv;shift argv;shift argv
  breaksw

  case -trueSig
    set trueSig="-trueSig $argv[2]"
  shift argv;shift argv
  breaksw

  case -renoise
    set renoise="-renoise"
  shift argv
  breaksw

  case -oldPsig
    set oldPsig="-oldPsig $argv[2]"
  shift argv;shift argv
  breaksw

  case -newPsig
    set newPsig="-newPsig $argv[2]"
  shift argv;shift argv
  breaksw

  case -missGround
    set missGround="-missGround"
  shift argv
  breaksw

  case -minGap
    set minGap="-minGap $argv[2]"
  shift argv;shift argv
  breaksw

  case -maxDN
    set maxDN="-maxDN $argv[2]"
  shift argv;shift argv
  breaksw

  case -bitRate
    set bitRate="-bitRate $argv[2]"
  shift argv;shift argv
  breaksw

  case -meanN
    set meanN="-meanN $argv[2]"
  shift argv;shift argv
  breaksw

  case -thresh
    set thresh="-thresh $argv[2]"
  shift argv;shift argv
  breaksw

  case -sWidth
    set sWidth="-sWidth $argv[2]"
  shift argv;shift argv
  breaksw

  case -psWidth
    set psWidth="-psWidth $argv[2]"
  shift argv;shift argv
  breaksw

  case -gWidth
    set gWidth="-gWidth $argv[2]"
  shift argv;shift argv
  breaksw

  case -minWidth
    set minWidth="-minWidth $argv[2]"
  shift argv;shift argv
  breaksw

  case -varNoise
    set varNoise="-varNoise"
  shift argv
  breaksw

  case -varScale
    set varScale="-varScale $argv[2]"
  shift argv;shift argv
  breaksw

  case -statsLen
    set statsLen="-statsLen $argv[2]"
  shift argv;shift argv
  breaksw

  case -medNoise
    set medNoise="-medNoise"
  shift argv
  breaksw

  case -noiseTrack
    set noiseTrack="-noiseTrack"
  shift argv
  breaksw

  case -rhoG
    set rhoG="-rhoG $argv[2]"
  shift argv;shift argv
  breaksw

  case -rhoC
    set rhoC="-rhoC $argv[2]"
  shift argv;shift argv
  breaksw

  case -offset
    set offset="-offset $argv[2]"
  shift argv;shift argv
  breaksw

  case -progRoot
    set progRoot="$argv[2]"
  shift argv;shift argv
  breaksw

  case -dontTrustGround
    set dontTrustGround="-dontTrustGround"
  shift argv
  breaksw

  default:
    echo "Unknown argument $argv[1]"
    echo "Who knows"
    echo "Type metricWithProgress.csh -help"
    exit;

  endsw
end

# progress check
touch $progRoot

# extract metrics
gediMetric -inList $inList -outRoot $outRoot $writeFit $ground $useInt $useFrac -rhRes $rhRes $bayesGround $gTol $noRHgauss $dcBias $nSig $seed $hNoise $linkNoise $trueSig $renoise $newPsig $oldPsig $missGround $minGap $maxDN $bitRate $meanN $thresh $sWidth $psWidth $gWidth $minWidth $varNoise $varScale $statsLen $medNoise $noiseTrack $rhoG $rhoC $offset $dontTrustGround

# mark finishing
rm $progRoot

