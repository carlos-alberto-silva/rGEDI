
($0&&($1!="#")){
  file=$1;
  if(NF>1){
    tMinX=$2;
    tMinY=$3;
    tMaxX=$5;
    tMaxY=$6;
    if((tMinX<=maxX)&&(tMaxX>=minX)&&(tMinY<=maxY)&&(tMaxY>=minY))print file;
  }else print file;
}

