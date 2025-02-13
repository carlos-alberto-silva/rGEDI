#!/bin/tcsh -f

#############################
# Pararellises gediMetric.c #
#############################


set bin="$GEDIRAT_ROOT"
set progRoot="/tmp/gediMetricProgress.$$.dat"

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
set minGsig=" "
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

@ nCores=1



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

  case -minGsig
    set minGsig="-minGsig $argv[2]"
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

  case -nCores
    @ nCores=$argv[2]
  shift argv;shift argv
  breaksw

  case -dontTrustGround
    set dontTrustGround="-dontTrustGround"
  shift argv
  breaksw

  case -help
    echo " "
    echo "-inList list;      input filename list"
    echo "-outRoot name;     output filename root"
    echo "-writeFit;         write out fitted waveforms"
    echo "-ground;           waveform files contain ground estimates"
    echo "-useInt;           use waveform from intensity instead of count"
    echo "-useFrac;          use waveform from fraction instead of count"
    echo "-rhRes res;        RH metric resolution in %"
    echo "-bayesGround;      turn robust ground decision on"
    echo "-gTol tol;         can't remember"
    echo "-noRHgauss;        don't so Gaussian ground finding"
    echo "-nCores n;         number of cores to split job over"
    echo "-dontTrustGround;  do not trust polynomial ground"
    echo " "
    echo "# Adding noise:"
    echo "-dcBias n;       mean noise level"
    echo "-nSig sig;       noise sigma"
    echo "-seed n;         random number seed"
    echo "-hNoise n;       hard threshold noise as a fraction of integral"
    echo "-linkNoise linkM cov;     apply Gaussian noise based on link margin at a cover"
    echo "-trueSig sig;    true sigma of background noise"
    echo "-renoise;        remove noise feom truth"
    echo "-oldPsig sig;    sigma of existing pulse when renoising"
    echo "-newPsig sig;    sigma of new pulse when renoising"
    echo "-missGround;     assume ground is missed to assess RH metrics"
    echo "-minGap gap;     delete signal beneath min detectable gap fraction"
    echo "-maxDN max;      maximum DN"
    echo "-bitRate n;      DN bit rate"
    echo " "
    echo "# Denoising:"
    echo "-meanN n;        mean noise level"
    echo "-thresh n;       noise threshold"
    echo "-sWidth sig;     smoothing width"
    echo "-psWidth sigma;  pre-smoothing width"
    echo "-gWidth sig;     Gaussian paremter selection width"
    echo "-minGsig sig;    minimum Gaussian width to fit"
    echo "-minWidth n;     minimum feature width in bins"
    echo "-varNoise;       variable noise threshold"
    echo "-varScale x;     variable noise threshold scale"
    echo "-statsLen len;   length to calculate noise stats over"
    echo "-medNoise;       use median stats rather than mean"
    echo "-noiseTrack;     use noise tracking"
    echo "-rhoG rho;       ground reflectance"
    echo "-rhoC rho;       canopy reflectance"
    echo "-offset y;       waveform DN offset"
    echo " "
    exit

  default:
    echo "Unknown argument $argv[1]"
    echo "Who knows"
    echo "Type gediMetric.csh -help"
    exit;

  endsw
end



# split input into multiple files
set tempList="/tmp/tempList.$$.dat"
set tempOut="/tmp/tempOut.$$.dat"
@ nFiles=`wc -l` < $inList
@ nPer=`echo $nCores $nFiles|gawk '{print int($2/$1+1)}'` < $inList
@ i=1
@ j=0
while( $i <= $nFiles )
  set thisList="$tempList.$j"
  gawk -v i=$i -v nPer=$nPer '{if((NR>=i)&&(NR<(i+nPer)))print $0}' < $inList > $thisList
  @ i+=$nPer
  @ j++
end


# run over all
@ i=0
while( $i < $nCores )
  set thisList="$tempList.$i"
  set thisOut="$tempOut.$i"

  $bin/metricWithProgress.csh -inList $thisList -outRoot $thisOut $writeFit $ground $useInt $useFrac -rhRes $rhRes $bayesGround $gTol $noRHgauss $dcBias $nSig $seed $hNoise $linkNoise $trueSig $renoise $newPsig $oldPsig $missGround $minGap $maxDN $bitRate $meanN $thresh $sWidth $psWidth $gWidth $minWidth $varNoise $varScale $statsLen $medNoise $noiseTrack $rhoG $rhoC $offset $minGsig -progRoot $progRoot.$i $dontTrustGround & 
  @ i++
end

sleep 2

# check progress
$bin/checkProgress.csh -progRoot $progRoot -nCores $nCores
sleep 2

# combine multiple files
set output="$outRoot.metric.txt"
@ i=0
while( $i < $nCores )
  set thisList="$tempList.$i"
  set thisOut="$tempOut.$i.metric.txt"
  set thisGauss="$tempOut.$i.gauss.txt"

  if( $i == 0 )then
    cat $thisOut > $output
  else
    gawk '(NR>1){print $0}' < $thisOut >> $output
  endif

  if( -e $thisList )rm $thisList
  if( -e $thisGauss )rm $thisGauss
  if( -e $thisOut )rm $thisOut
  @ i++
end

echo "Written to $output"

