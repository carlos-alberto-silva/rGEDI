BEGIN{
  numb=-1;
  srand(seed);
  tol=0.000001;
}

{
  n=split($1,bits,"_");
  if(bits[1]=="C"){  # It is a new entry
    if(n==2){  # new cell
      numb++;
      nIn[numb]=-1;
    }else if(n==4){  # new track within cell
      if((usePhen==0)||(leafon[numb,nIn[numb]]>0.01))nIn[numb]++;
    }
  }else{  # read a track
    if($1=="latenter")y0[numb,nIn[numb]]=$3;
    else if($1=="lonenter")x0[numb,nIn[numb]]=$3;
    else if($1=="latexit")y1[numb,nIn[numb]]=$3;
    else if($1=="lonexit")x1[numb,nIn[numb]]=$3;
    else if($1=="LEAFON")leafon[numb,nIn[numb]]=$3;
    else if($1=="beamid")beamid[numb,nIn[numb]]=$3;
    else if($1=="lon")pixX[numb]=$3;     
    else if($1=="lat")pixY[numb]=$3;
    else if($1=="sunel")sunel[numb,nIn[numb]]=$3;
    else if($1=="tenter"){
      split($3,bits,"");
      year[numb,nIn[numb]]=sprintf("%d%d",bits[1],bits[2]);
      doy[numb,nIn[numb]]=sprintf("%d%d%d",bits[3],bits[4],bits[5]);
    }
  }
}

END{
  # add the last one on
  numb++;
  for(i=0;i<numb;i++)nIn[i]++;

  # remove gaps in data by compressing x.
  for(i=0;i<numb;i++){
    xStart=-180+i*res;
    for(j=0;j<nIn[i];j++){
      x0[i,j]+=xStart-pixX[i];
      x1[i,j]+=xStart-pixX[i];
    }
  }

  # width of target area
  buffX=(maxY-minY)*cos(ang)/sin(ang);
  buffBins=int(buffX/res+1);
  nX=int((maxX-minX)/res+1)+2*buffBins;

  # pick a random start point
  start=int(rand()*numb)-buffBins;
  if(start<0)start=0;
  else if((start+nX)>=numb)start=numb-nX;  # make sure there's enough to the right
  ending=start+nX;

  # determine bounds
  minLon=minLat=100000.0;
  maxLon=maxLat=-1000000.0;
  for(i=start;i<=ending;i++){
    for(j=0;j<nIn[i];j++){
      if(x0[i,j]<minLon)minLon=x0[i,j];
      if(x0[i,j]>maxLon)maxLon=x0[i,j];
      if(x1[i,j]<minLon)minLon=x1[i,j];
      if(x1[i,j]>maxLon)maxLon=x1[i,j];
      if(y0[i,j]<minLat)minLat=y0[i,j];
      if(y0[i,j]>maxLat)maxLat=y0[i,j];
      if(y1[i,j]<minLat)minLat=y1[i,j];
      if(y1[i,j]>maxLat)maxLat=y1[i,j];
    }
  }

  # choose footprint coords
  nUse=0;
  for(i=start;i<=ending;i++){

    for(j=0;j<nIn[i];j++){
      # sort orders
      if(y0[i,j]<y1[i,j]){
        xBot=x0[i,j];
        yBot=y0[i,j];
        xTop=x1[i,j];
        yTop=y1[i,j];
      }else{
        xBot=x1[i,j];
        yBot=y1[i,j];
        xTop=x0[i,j];
        yTop=y0[i,j];
      }

      if((yBot-minLat)<tol){  # accept bottom only
        xS[nUse]=xBot;
        yS[nUse]=yBot;
        zen[nUse]=atan2(xTop-xBot,yTop-yBot);
        beamType[nUse]=beamid[i,j];
        sunZen[nUse]=sunel[i,j];
        ye[nUse]=year[i,j];
        d[nUse]=doy[i,j];
        nUse++;
      }else if(i==start){   # accept left edge
        #if((xBot-minLon)<tol){
        xS[nUse]=xBot;
        yS[nUse]=yBot;
        zen[nUse]=atan2(xTop-xBot,yTop-yBot);
        beamType[nUse]=beamid[i,j];
        ye[nUse]=year[i,j];
        d[nUse]=doy[i,j];
        nUse++;
        #}
      }else if(i==ending){  # accept right edge
        if((maxLon-xBot)<tol){
          xS[nUse]=xBot;
          yS[nUse]=yBot;
          zen[nUse]=atan2(xTop-xBot,yTop-yBot);
          beamType[nUse]=beamid[i,j];
          ye[nUse]=year[i,j];
          d[nUse]=doy[i,j];
          nUse++;
        }
      }
    }
  }

  # write them out
  for(i=0;i<nUse;i++){
    # do we accept from cloud?
    if(rand()>=cloudFrac){
      y=minY;
      x=xS[i]+minX-(minLon+buffX);
      if(sunZen[i]>=0.0)tim="day"
      else              tim="night"
      waveID=sprintf("%s.%d.%d.%s.zen.%d.y.%d.doy.%d",beamType[i],i,j,tim,int(sunZen[i]),ye[i],d[i]);
      printf("%.10f %.10f %f %s\n",x,y,zen[i],waveID);
    }# cloud test
  }
}

