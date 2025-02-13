BEGIN{
  minX=minY=10000000000.0;
  maxX=maxY=-10000000000.0;
}

($0){
  if($1!="#"){
    if(NF>=latInd){  # to avoid truncated files
      x=$lonInd;
      y=$latInd;
      if(x<minX)minX=x;
      if(x>maxX)maxX=x;
      if(y<minY)minY=y;
      if(y>maxY)maxY=y;
    }
  }else{  # read header
    for(i=1;i<=NF;i++){
      if($i=="lon,")lonInd=$(i-1);
      else if($i=="lat,")latInd=$(i-1);
    }
  }
}

END{
  printf("%f\n%f\n%f\n%f\n",minX,maxX,minY,maxY);
}

