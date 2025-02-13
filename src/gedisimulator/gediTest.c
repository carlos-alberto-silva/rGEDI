#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"



/*#########################*/
/*# Tests GEDI waveform   #*/ 
/*# generation   2015     #*/
/*# svenhancock@gmail.com #*/
/*#########################*/

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



/*####################################*/
/*pulse structure*/

typedef struct{
  int nBins;
  int centBin;  /*peak bin*/
  float *y;
  float *x;
}pulseStruct;

/*####################################*/
/*control structure*/

typedef struct{
  char inNamen[1000];
  char outNamen[1000];

  /*GEDI fotprint parameters*/
  double coord[3];
  float res;       /*range resolution*/
  float pFWHM;     /*pulse width in ns*/
  float pSigma;    /*pulse width in metres*/
  float fWidth;    /*footprint width*/
  float fSigma;    
  float pRes;
  pulseStruct *pulse;

  /*tolerances*/
  float iThresh;   /*intensity threshold*/
}control;


/*####################################*/
/*data structure*/

typedef struct{
  uint64_t nPoints;
  double *x;
  double *y;
  double *z;
  int *refl;
}dataStruct;


/*####################################*/
/*waveform structure*/

typedef struct{
  float **wave;  /*waveforms*/
  double minZ;   /*elevation bounds*/
  double maxZ;   /*elevation bounds*/
  int nBins;     /*number of wave bins*/
  int nWaves;    /*number of different ways*/
}waveStruct;


/*####################################*/
/*main*/

int main(int argc,char **argv)
{
  control *dimage=NULL;
  control *readCommands(int,char **);
  dataStruct *data=NULL;
  dataStruct *readAsciiData(char *);
  waveStruct *waves=NULL;
  waveStruct *makeGediWaves(control *,dataStruct *);
  void writeGEDIwave(control *,waveStruct *);
  void setGediPulse(control *);
 

  /*read command line*/
  dimage=readCommands(argc,argv);

  /*read data*/
  data=readAsciiData(dimage->inNamen);

  /*set up the pulse*/
  setGediPulse(dimage);

  /*make waveforms*/
  waves=makeGediWaves(dimage,data);

  /*output results*/
  writeGEDIwave(dimage,waves);

  /*tidy up*/
  if(data){
    TIDY(data->x);
    TIDY(data->y);
    TIDY(data->z);
    TIDY(data->refl);
    TIDY(data);
  }
  if(waves){
    TTIDY((void **)waves->wave,waves->nWaves);
    TIDY(waves);
  }
  if(dimage){
    if(dimage->pulse){
      TIDY(dimage->pulse->y);
      TIDY(dimage->pulse);
    }
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*####################################*/
/*write GEDI waveforms*/

void writeGEDIwave(control *dimage,waveStruct *waves)
{
  int i=0,j=0;
  float r=0;
  FILE *opoo=NULL;

 if((opoo=fopen(dimage->outNamen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",dimage->outNamen);
    exit(1);
  } 

  for(i=0;i<waves->nBins;i++){
    r=(float)waves->maxZ-(float)i*dimage->res;

    fprintf(opoo,"%f",r);
    for(j=0;j<waves->nWaves;j++)fprintf(opoo," %f",waves->wave[j][i]);
    fprintf(opoo,"\n");
  }

  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Written to %s\n",dimage->outNamen);
  return;
}/*writeGEDIwave*/


/*####################################*/
/*make GEDI waveforms*/

waveStruct *makeGediWaves(control *dimage,dataStruct *data)
{
  int bin=0,j=0,k=0;
  uint64_t i=0;
  double sep=0;
  double dX=0,dY=0;
  float refl=0,rScale=0;
  float tot=0;
  waveStruct *waves=NULL;
  waveStruct *allocateGEDIwaves(control *,dataStruct *);


  /*allocate*/
  waves=allocateGEDIwaves(dimage,data);


  /*make waves*/
  for(i=0;i<data->nPoints;i++){
    dX=data->x[i]-dimage->coord[0];
    dY=data->y[i]-dimage->coord[1];
    sep=sqrt(dX*dX+dY*dY);

    rScale=(float)gaussian(sep,(double)dimage->fSigma,0.0);

    if(rScale>dimage->iThresh){  /*if bright enough to matter*/
      refl=(float)data->refl[i]*rScale;
      for(j=0;j<dimage->pulse->nBins;j++){
        bin=(int)((waves->maxZ-data->z[i]+(double)dimage->pulse->x[j])/(double)dimage->res);
        if((bin>=0)&&(bin<waves->nBins)){
          waves->wave[0][bin]+=refl*dimage->pulse->y[j];
          waves->wave[1][bin]+=dimage->pulse->y[j];
        }
      }
    }
  }/*wave making*/

  /*normalise integral*/
  for(k=0;k<waves->nWaves;k++){
    tot=0.0;
    for(j=0;j<waves->nBins;j++)tot+=waves->wave[k][j];
    for(j=0;j<waves->nBins;j++)waves->wave[k][j]/=tot;
  }

  return(waves);
}/*makeGediWaves*/


/*####################################*/
/*allocate wave structure*/

waveStruct *allocateGEDIwaves(control *dimage,dataStruct *data)
{
  int j=0;
  uint64_t i=0;
  double maxZ=0,minZ=0;
  double buff=0;
  waveStruct *waves=NULL;

  if(!(waves=(waveStruct *)calloc(1,sizeof(waveStruct)))){
    fprintf(stderr,"error waveStruct allocation.\n");
    exit(1);
  }

  /*determine wave bounds*/
  buff=15.0;
  minZ=1000000000.0;
  maxZ=-1000000000.0;
  for(i=0;i<data->nPoints;i++){
    if(data->z[i]>maxZ)maxZ=data->z[i];
    if(data->z[i]<minZ)minZ=data->z[i];
  }/*bound finding*/

  waves->minZ=minZ-buff;
  waves->maxZ=maxZ+buff;

  waves->nBins=(int)((waves->maxZ-waves->minZ)/(double)dimage->res);
  waves->nWaves=2;
  waves->wave=fFalloc(waves->nWaves,"waveform",0);
  for(j=0;j<waves->nWaves;j++)waves->wave[j]=falloc((uint64_t)waves->nBins,"waveform",j+1);

  return(waves);
}/*allocateGEDIwaves*/


/*####################################*/
/*set GEDI pulse*/

void setGediPulse(control *dimage)
{
  int i=0;
  float fwhm=0;   /*FWHM in metres*/
  float x=0,y=0;
  float max=0,tot=0;

  if(!(dimage->pulse=(pulseStruct *)calloc(1,sizeof(pulseStruct)))){
    fprintf(stderr,"error pulseStruct allocation.\n");
    exit(1);
  }


  /*pulse length*/
  /*calculate sigma from FWHM*/
  fwhm=dimage->pFWHM*2.99792458/20.0;  /*time if for two way*/
  dimage->pSigma=(fwhm/2.0)/sqrt(log(2.0));

  /*determine number of bins*/
  dimage->pulse->nBins=0;
  x=0.0;
  do{
    y=(float)gaussian((double)x,(double)dimage->pSigma,0.0);
    x+=dimage->pRes;
    dimage->pulse->nBins+=2;  /*both sides of peak*/
  }while(y>=dimage->iThresh);

  dimage->pulse->x=falloc((uint64_t)dimage->pulse->nBins,"pulse x",0);
  dimage->pulse->y=falloc((uint64_t)dimage->pulse->nBins,"pulse y",0);
  dimage->pulse->centBin=(int)(dimage->pulse->nBins/2);

  max=-100.0;
  tot=0.0;
  x=-1.0*(float)dimage->pulse->centBin*dimage->pRes;
  for(i=0;i<dimage->pulse->nBins;i++){
    dimage->pulse->x[i]=x;
    dimage->pulse->y[i]=(float)gaussian((double)x,(float)dimage->pSigma,0.0);
    if(dimage->pulse->y[i]>max){
      max=dimage->pulse->y[i];
      dimage->pulse->centBin=i;
    }
    tot+=dimage->pulse->y[i];
    x+=dimage->pRes;
  }

  /*normalise to cope with rounding*/
  for(i=0;i<dimage->pulse->nBins;i++){
    dimage->pulse->y[i]/=tot;
  }


  /*footprint width*/
  dimage->fSigma=dimage->fWidth;  /*I think*/

  return;
}/*setGediPulse*/


/*####################################*/
/*read ASCII data*/

dataStruct *readAsciiData(char *inNamen)
{
  uint64_t i=0;
  dataStruct *data=NULL;
  char line[400],temp1[100],temp2[100];
  char temp3[100],temp4[100],temp5[100];
  char temp6[100],temp7[100],temp8[100];
  FILE *ipoo=NULL;


  if(!(data=(dataStruct *)calloc(1,sizeof(dataStruct)))){
    fprintf(stderr,"error dataStruct allocation.\n");
    exit(1);
  }


  if((ipoo=fopen(inNamen,"r"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",inNamen);
    exit(1);
  }

  /*count number of points*/
  data->nPoints=0;
  while(fgets(line,400,ipoo)!=NULL)if(strncasecmp(line,"#",1))data->nPoints++;

  /*allocate space*/
  data->x=dalloc(data->nPoints,"x",0);
  data->y=dalloc(data->nPoints,"y",0);
  data->z=dalloc(data->nPoints,"z",0);
  data->refl=ialloc(data->nPoints,"refl",0);


  /*rewind to start of file*/
  if(fseek(ipoo,(long)0,SEEK_SET)){ 
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  /*read data*/
  i=0;
  while(fgets(line,400,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      if(sscanf(line,"%s %s %s %s %s %s %s %s",temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8)==8){
        data->x[i]=atof(temp1);
        data->y[i]=atof(temp2);
        data->z[i]=atof(temp3);
        data->refl[i]=atoi(temp8);
        i++;
      }
    }
  }

  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return(data);
}


/*####################################*/
/*read command line*/

control *readCommands(int argc,char **argv)
{
  control *dimage=NULL;

  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }

  strcpy(dimage->inNamen,"/Users/dill/data/teast/maryland_play/teast.dat");
  strcpy(dimage->outNamen,"/Users/dill/data/teast/maryland_play/teast.wave");
  dimage->pFWHM=12.0;   /*12 ns FWHM*/
  dimage->fWidth=11.0;  /*22m diameter in e^-2*/
  dimage->res=0.15;
  dimage->pRes=0.01;
  dimage->coord[0]=624366.0;
  dimage->coord[1]=3.69810*pow(10.0,6.0);

  dimage->iThresh=0.00001;

  return(dimage);
}/*readCommands*/

/*the end*/
/*####################################*/

