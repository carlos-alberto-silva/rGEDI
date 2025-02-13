BEGIN{
  meanLat=meanLon=0.0;
  numb=0;
}

($0){
  if($1!="#"){  # read data
    meanLon+=$lonInd;
    meanLat+=$latInd;
    numb++;
  }else{  # read header
    for(i=1;i<=NF;i++){
      if($i=="lon,")lonInd=$(i-1);
      else if($i=="lat,")latInd=$(i-1);
    }
  }
}

END{
  if(numb>0){
    meanLon/=numb;
    meanLat/=numb;
  }
  printf("%.8f %.8f\n",meanLon,meanLat);
}

