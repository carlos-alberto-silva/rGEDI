BEGIN{
  readCoords=1;
  nCoords=0;
}

($0&&($1!="#")){
  if($1=="###")readCoords=0;
  else if(readCoords==1){   # read coordinates
    minX[nCoords]=$1-rad;
    maxX[nCoords]=$1+rad;
    minY[nCoords]=$2-rad;
    maxY[nCoords]=$2+rad;
    if(NF>=3)waveID[nCoords]=$3;
    nCoords++;
  }else{   # read LAS files
    if(NF>1){
      tMinX=$2;
      tMinY=$3;
      tMaxX=$5;
      tMaxY=$6;
      file=$1;
      for(i=0;i<nCoords;i++){
        if((tMinX<=maxX[i])&&(tMaxX>=minX[i])&&(tMinY<=maxY[i])&&(tMaxY>=minY[i])){
          print file,tMinX,tMinY;
          break;
        }
      }
    }
  }
}

