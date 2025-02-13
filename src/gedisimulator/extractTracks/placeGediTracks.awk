BEGIN{
  numb=0;
  tol=0.000001;
}

{
  xS[numb]=$1;
  yS[numb]=$2;
  zen[numb]=$3;
  waveID[numb]=$4;
  numb++;
}

END{
  # write them out
  for(i=0;i<numb;i++){
    y=yS[i];
    d=0.0;
    j=0;
    while(y<=(maxY+tol)){
      x=xS[i]+d*sin(zen[i]);
      y=yS[i]+d*cos(zen[i]);

      if((x>=minX)&&(x<=maxX)&&(y>=minY)&&(y<=maxY)){
        printf("%.2f %.2f %s.%d\n",x,y,waveID[i],j);
      }
      d+=alongTrack;
      j++;
    }
  }
}

