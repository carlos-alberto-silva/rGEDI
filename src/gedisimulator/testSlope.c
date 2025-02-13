#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "tools.h"
#include "tools.c"

/*#######################################*/
/*# Copyright 2015-2016, Steven Hancock #*/
/*# The program is distributed under    #*/
/*# the terms of the GNU General Public #*/
/*# License.    svenhancock@gmail.com   #*/
/*#######################################*/


/*########################################################################*/
/*# This file is part of the NASA GEDI simulator, gediRat.               #*/
/*#                                                                      #*/
/*# gediRat is free software: you can redistribute it and/or modify      #*/
/*# it under the terms of the GNU General Public License as published by #*/
/*# the Free Software Foundation, either version 3 of the License, or    #*/
/*#  (at your option) any later version.                                 #*/
/*#                                                                      #*/
/*# gediRat is distributed in the hope that it will be useful,           #*/
/*# but WITHOUT ANY WARRANTY; without even the implied warranty of       #*/
/*#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #*/
/*#   GNU General Public License for more details.                       #*/
/*#                                                                      #*/
/*#    You should have received a copy of the GNU General Public License #*/
/*#    along with gediRat.  If not, see <http://www.gnu.org/licenses/>.  #*/
/*########################################################################*/

#define TOL 0.00000001


/*#############################*/
/*control structure*/

typedef struct{
  char outNamen[1000];   /*output filename*/
  char writeWave;       /*write waveform switch*/
  float res;            /*instrument resolution*/
  float A;              /*ground amplitude*/
  float pSig;           /*pulse length (sigma)*/
  float fSig;           /*footprint width (sigma)*/
  float maxSlope;       /*maximum slope, degrees*/
  float step;           /*slope step size, derees*/
}control;


/*#############################*/
/*main*/

int main(int argc,char **argv)
{
  int nBins;
  float slope=0,fSlope=0;
  float *wave=NULL;
  float *makeWave(float,float,float,float,int *,float);
  float stdev=0;
  float waveStdev(float *,int,float);
  float findSlope(float,float,float);
  float findMax(float *,int);
  float maxAmp=0;
  control *dimage=NULL;
  control *readCommands(int,char **);
  void writeWave(float,float *,int,control *);
  FILE *opoo=NULL;


  /*read command line*/
  dimage=readCommands(argc,argv);

  /*open output and write header*/
  if((opoo=fopen(dimage->outNamen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",dimage->outNamen);
    exit(1);
  }
  fprintf(opoo,"# 1 slope, 2 stdev, 3 fSlope, 4 pSig, 5 fSig, 6 maxAmp\n");


  /*loop over slopes*/
  slope=0.0;
  while(slope<=dimage->maxSlope){
    /*make the waveform*/
    wave=makeWave(slope,dimage->pSig,dimage->fSig,dimage->A,&nBins,dimage->res);

    if(dimage->writeWave)writeWave(slope,wave,nBins,dimage);

    /*calculate stdev*/
    stdev=waveStdev(wave,nBins,dimage->res);

    /*determine maximum*/
    maxAmp=findMax(wave,nBins);
    TIDY(wave);

    /*calculate slope*/
    fSlope=findSlope(stdev,dimage->pSig,dimage->fSig);

    fprintf(opoo,"%f %f %f %f %f %f\n",slope*180.0/M_PI,stdev,fSlope,dimage->pSig,dimage->fSig,maxAmp);
    slope+=dimage->step;
  }

  fprintf(stdout,"Written to %s\n",dimage->outNamen);

  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  TIDY(dimage);
  return(0);
}/*main*/


/*#############################*/
/*write waveform*/

void writeWave(float slope, float *wave,int nBins,control *dimage)
{
  int i=0;
  char namen[200];
  FILE *opoo=NULL;


  sprintf(namen,"%s.%g.wave",dimage->outNamen,slope*180.0/M_PI);
  if((opoo=fopen(namen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",namen);
    exit(1);
  }

  for(i=0;i<nBins;i++){
    fprintf(opoo,"%f %f\n",(float)i*dimage->res,wave[i]);
  }

  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Written to %s\n",namen);
  return;
}/*writeWave*/


/*#############################*/
/*determine maximum amplitude*/

float findMax(float *wave,int nBins)
{
  int i=0;
  float maxAmp=0;

  maxAmp=-109.0;
  for(i=0;i<nBins;i++){
    if(wave[i]>maxAmp)maxAmp=wave[i];
  }


  return(maxAmp);
}/*findMax*/


/*#############################*/
/*slope from stdev*/

float findSlope(float stdev,float pSig,float fSig)
{
  float slope=0;

  slope=atan2(sqrt(stdev*stdev-pSig*pSig),fSig)*180.0/M_PI;

  return(slope);
}/*findSlope*/


/*#############################*/
/*determine standard deviation*/

float waveStdev(float *wave,int nBins,float res)
{
  int i=0;
  float x=0,tot=0;
  float mean=0,stdev=0;

  mean=tot=0.0;
  for(i=0;i<nBins;i++){
    x=(float)i*res;
    mean+=wave[i]*x;
    tot+=wave[i];
  }
  mean/=tot;

  stdev=0.0;
  for(i=0;i<nBins;i++){
    x=(float)i*res;
    stdev+=(x-mean)*(x-mean)*wave[i];
  }
  stdev=sqrt(stdev/tot);

  return(stdev);
}/*waveStdev*/


/*#############################*/
/*make waveform*/

float *makeWave(float slope,float pSig,float fSig,float A,int *nBins,float res)
{
  int bin=0;
  float *wave=NULL,aRes=0;
  float footRad=0;
  float x=0,y=0,z=0;
  float tot=0,*smoothed=NULL;;
  float *smooth(float,int,float *,float);
  double gaussian(double,double,double);

  *nBins=20000;
  wave=falloc((uint64_t)(*nBins),"wave",0);

  /*determine footprint radius*/
  footRad=0.0;
  aRes=0.05;
  do{
    y=A*(float)gaussian((double)footRad,(double)fSig,0.0);
    footRad+=aRes;
  }while(y>=TOL);

  /*calculate floor return*/
  x=-1.0*footRad;
  do{
    z=x*sin(slope)/cos(slope);
    bin=(int)(z/res)+(*nBins)/2;
    if((bin<0)||(bin>=(*nBins))){
      fprintf(stderr,"Not enough bins\n");
      exit(1);
    }

    tot=0.0;
    z=0.0;
    do{
      y=(float)gaussian((double)z,(double)fSig,0.0);
      tot+=y*aRes;
      z+=aRes;
    }while(y>=TOL);
    tot*=2.0;
    wave[bin]+=tot*(float)gaussian((double)x,(double)fSig,0.0)*aRes;
    x+=aRes;
  }while(x<=footRad);

  /*smooth by pulse blurring*/
  smoothed=smooth(pSig,*nBins,wave,res);
  TIDY(wave);

  return(smoothed);
}/*makeWave*/


/*#############################*/

control *readCommands(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;

  /*allocate space*/
  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }

  /*switches*/
  dimage->writeWave=0;

  /*instrument parameters*/
  dimage->res=0.15;
  dimage->A=0.5;
  dimage->pSig=0.955414;
  dimage->fSig=5.5;

  /*slope range*/
  dimage->maxSlope=60.0;
  dimage->step=0.5;

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(dimage->outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-res",4)){
        checkArguments(1,i,argc,"-res");
        dimage->res=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-A",2)){
        checkArguments(1,i,argc,"-A");
        dimage->A=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-pSig",5)){
        checkArguments(1,i,argc,"-pSig");
        dimage->pSig=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-fSig",5)){
        checkArguments(1,i,argc,"-fSig");
        dimage->fSig=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxSlope",9)){
        checkArguments(1,i,argc,"-maxSlope");
        dimage->maxSlope=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-step",5)){
        checkArguments(1,i,argc,"-step");
        dimage->step=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-writeWave",10)){
        dimage->writeWave=1;
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n#####\nProgram to assess the impact of slope on lidar and test slope calculation\n#####\n\n-output name;        output filename\n-res res;            waveform resolution\n-A A;                ground ampltiude\n-pSig pSig;          pulse length (sigma)\n-fSig fSig;          footprint width (sigma)\n-maxSlope maxSlope;  maximum slope to go to\n-step step;          output slope resolution\n-writeWave;          write waveforms\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry gediRat -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }/*command parser*/



  /*convert to radians*/
  dimage->maxSlope*=M_PI/180.0;
  dimage->step*=M_PI/180.0;

  return(dimage);
}/*readCommands*/

/*the end*/
/*#############################*/

