BEGIN{
  minX=minY=10000000000.0;
  maxX=maxY=-10000000000.0;
}

($0&&($1!="#")){
  x=$2;
  y=$3;;
  if(x<minX)minX=x;
  if(y<minY)minY=y;
  x=$5;
  y=$6;;
  if(x>maxX)maxX=x;
  if(y>maxY)maxY=y;
}

END{
  printf("%f\n%f\n%f\n%f\n",minX,maxX,minY,maxY);
}

