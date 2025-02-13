BEGIN{
  minX=minY=10000000000.0;
  maxX=maxY=-10000000000.0;
}

($0&&($1!="#")){
  x=$1;
  y=$2;

  if((x-rad)<minX)minX=x-rad;
  if((x+rad)>maxX)maxX=x+rad;
  if((y-rad)<minY)minY=y-rad;
  if((y+rad)>maxY)maxY=y+rad;
}

END{
  printf("%f\n%f\n%f\n%f\n",minX,minY,maxX,maxY);
}

