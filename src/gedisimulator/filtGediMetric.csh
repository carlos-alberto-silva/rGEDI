#!/bin/tcsh -f

set list="/tmp/filtList.$$.dat"
ls *.metric.txt|gawk '{for(i=1;i<=NF;i++)print $i}' > $list

@ nFiles=`wc -l` < $list
@ i=1
while( $i <= $nFiles )
  set file=`gawk -v i=$i '{if(i==NR)print $1}'` < $list
  set filt="$file:r.filt"

  gawk '{if(($2>0)&&($2<10000)&&($6>0)&&($6<10000))print $0}' < $file > $filt
  echo "Written to $filt"

  @ i++
end


if( -e $list )rm $list

