BEGIN{
  readMetric=0;
  thresh=minSep*minSep;
  numb=0;
}

($0){
  if($1=="###")readMetric=1;
  else if(readMetric==1){
    if($1!="#"){
      for(i=0;i<numb;i++){
        sepSq=(($lonInd-x[i])^2)*($latInd-y[i])^2;
        if((sepSq<minSepSq[i])&&(sepSq<=thresh)){
          line[i]=$0;
          minSepSq[i]=sepSq;
          found[i]=1;
        }
      }
    }else{ # read and write header
      print $0;
      for(i=1;i<=NF;i++){
        if($i=="lon,")lonInd=$(i-1);
        else if($i=="lat,")latInd=$(i-1);
      }
    }
  }else{
    x[numb]=$1;
    y[numb]=$2;
    minSepSq[numb]=10000000000.0;
    found[numb]=0
    numb++;
  }
}

END{
  # write out data
  for(i=0;i<numb;i++){
    if(found[i]==1)print line[i];
  }
}

