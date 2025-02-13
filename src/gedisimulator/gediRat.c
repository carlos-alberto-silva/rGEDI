#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "hdf5.h"
#include "libLasRead.h"
#include "libDEMhandle.h"
#include "libLasProcess.h"
#include "libLidVoxel.h"
#include "libLidarHDF.h"
#include "libOctree.h"
#include "gediIO.h"


/*#########################*/
/*# Generated grids of    #*/
/*# waveforms    2015     #*/
/*# svenhancock@gmail.com #*/
/*#########################*/

/*#######################################*/
/*# Copyright 2015-2017, Steven Hancock #*/
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

#define TOL 0.000001   /*a tolerance*/


/*####################################*/
/*control structure*/

typedef struct{
  char **inList;
  char outNamen[1000];
  char waveNamen[400];

  /*IO structure*/
  gediIOstruct gediIO; /*generic IO options*/
  gediRatStruct gediRat; /*simulator options*/

  /*options*/
  char listFiles;      /*list waves only*/
  char overWrite;      /*overwrite old waveform switch*/
  uint64_t pBuffSize;  /*point buffer rading size in bytes*/
  char waveID[200];    /*wave ID if we are to use it*/
  char useID;          /*use wave ID*/
  char polyGr;         /*fit a polynomial to the ground*/
  char nnGr;           /*ground DEM from nearest neighbour*/

  /*HDF5 output*/
  char writeHDF;       /*write output as hdf5*/
  char writeL1B;       /*write L1B HDF5 output format*/
  int maxBins;         /*bins per wave for HDF5 output*/
  int hdfCount;        /*count used footprints*/
}control;


/*####################################*/
/*main*/

int main(int argc,char **argv)
{
  int i=0,j=0;
  control *dimage=NULL;
  control *readCommands(int,char **);
  lasFile *las=NULL;
  pCloudStruct **data=NULL;
  waveStruct *waves=NULL;
  gediHDF *hdfData=NULL;
  void writeGEDIwave(control *,waveStruct *,int);
  void tidySMoothPulse();
  void checkThisFile(lasFile *,control *,int);
  void groundFromDEM(pCloudStruct **,control *,waveStruct *);
  void checkWaveOverwrite(control *,int);

 
  /*read command line*/
  dimage=readCommands(argc,argv);

  /*set up the pulse*/
  setGediPulse(&dimage->gediIO,&dimage->gediRat);

  /*set up grid or batch if needed*/
  setGediGrid(&dimage->gediIO,&dimage->gediRat);

  /*loop over las files and read*/
  if(!(data=(pCloudStruct **)calloc(dimage->gediIO.nFiles,sizeof(pCloudStruct *)))){
    fprintf(stderr,"error waveStruct allocation.\n");
    exit(1);
  }
  for(i=0;i<dimage->gediIO.nFiles;i++){
    /*report progress if reading all data here*/
    if(dimage->gediRat.doGrid||dimage->gediRat.readALSonce)fprintf(stdout,"File %d of %d",i+1,dimage->gediIO.nFiles);
    /*read lasFile*/
    las=readLasHead(dimage->inList[i],dimage->pBuffSize);

    /*read data or write filename if needed*/
    if(dimage->listFiles==0)data[i]=readALSdata(las,&dimage->gediRat,i);
    else                    checkThisFile(las,dimage,i);
    if(dimage->gediRat.doGrid||dimage->gediRat.readALSonce)fprintf(stdout," nPoints %u\n",data[i]->nPoints);

    /*tidy lasFIle*/
    las=tidyLasFile(las);
  }/*file loop*/

  /*set up HDF5 if needed*/
  if(dimage->writeHDF)hdfData=setUpHDF(&dimage->gediIO,&dimage->gediRat,dimage->useID,dimage->waveID,&dimage->hdfCount,dimage->maxBins);

  /*make waveforms*/
  if(dimage->listFiles==0){
    /*loop over waveforms*/
    for(i=0;i<dimage->gediRat.gNx;i++){
      for(j=0;j<dimage->gediRat.gNy;j++){
        /*update centre coord*/
        updateGediCoord(&dimage->gediRat,i,j);

        /*see if that file already exists*/
        checkWaveOverwrite(dimage,i);

        /*if it is not to be overwritten*/
        if(dimage->gediRat.useFootprint){
          /*set up footprint*/
          setGediFootprint(&dimage->gediRat,&dimage->gediIO);

          /*make waveforms*/
          waves=makeGediWaves(&dimage->gediRat,&dimage->gediIO,data);
        }

        /*if it is usable*/
        if(dimage->gediRat.useFootprint){
          /*progress report*/
          if(dimage->writeHDF){
            if((dimage->hdfCount%dimage->gediIO.nMessages)==0){
              fprintf(stdout,"Wave %d of %d\n",i*dimage->gediRat.gNy+j,dimage->gediRat.gNx*dimage->gediRat.gNy);
            }
          }
          /*find the ground if needed*/
          if(dimage->gediIO.ground&&(dimage->polyGr||dimage->nnGr))groundFromDEM(data,dimage,waves);
  
          /*output results*/
          if(dimage->writeHDF)packGEDIhdf(waves,hdfData,i+j*dimage->gediRat.gNx,&dimage->gediIO,&dimage->gediRat,&dimage->hdfCount,dimage->useID,dimage->waveID);
          else                writeGEDIwave(dimage,waves,i+j*dimage->gediRat.gNx);
        }

        /*tidy up*/
        TIDY(dimage->gediRat.nGrid);
        TIDY(dimage->gediRat.lobe);
        waves=tidyWaveStruct(waves);
      }/*grid y loop*/
    }/*grid x loop*/

    /*write HDF if needed and not blank*/
    if(dimage->writeHDF&&(dimage->hdfCount>0)){
      hdfData->nWaves=dimage->hdfCount;  /*account for unusable footprints*/
      if(dimage->writeL1B)writeGEDIl1b(hdfData,dimage->outNamen,&(dimage->gediIO));
      else                writeGEDIhdf(hdfData,dimage->outNamen,&(dimage->gediIO));
    }
  }/*make and write a waveform if needed*/


  /*tidy up*/
  if(data){
    for(i=0;i<dimage->gediIO.nFiles;i++)data[i]=tidyPointCloud(data[i]);
    TIDY(data);
  }
  hdfData=tidyGediHDF(hdfData);
  if(dimage){
    if(dimage->gediIO.pulse){
      TIDY(dimage->gediIO.pulse->y);
      TIDY(dimage->gediIO.pulse->x);
      TIDY(dimage->gediIO.pulse);
    }
    TTIDY((void **)dimage->gediRat.coords,dimage->gediRat.gNx);
    TTIDY((void **)dimage->gediRat.geoCoords,dimage->gediRat.gNx);
    TTIDY((void **)dimage->gediRat.waveIDlist,dimage->gediRat.gNx);
    TTIDY((void **)dimage->inList,dimage->gediIO.nFiles);
    if(dimage->gediRat.wavefront){
      TTIDY((void **)dimage->gediRat.wavefront->front,dimage->gediRat.wavefront->nX);
      TIDY(dimage->gediRat.wavefront);
    }
    dimage->gediRat.octree=tidyOctree(dimage->gediRat.octree);
    TIDY(dimage);
  }
  tidySMoothPulse();
  return(0);
}/*main*/


/*##############################################*/
/*check whether we are to overwrite or not*/

void checkWaveOverwrite(control *dimage,int numb)
{
  FILE *opoo=NULL;

  /*set output filename*/
  if(dimage->gediRat.doGrid){
    sprintf(dimage->waveNamen,"%s.%d.%d.wave",dimage->outNamen,(int)dimage->gediRat.coord[0],(int)dimage->gediRat.coord[1]);
  }else if(dimage->gediRat.readALSonce){
    sprintf(dimage->waveNamen,"%s.%s.wave",dimage->outNamen,dimage->gediRat.waveIDlist[numb]);
  }else{
    strcpy(dimage->waveNamen,dimage->outNamen);
  }

  /*see if the file exists if we are checking for overwriting*/
  if(dimage->overWrite)dimage->gediRat.useFootprint=1;
  else{
    if((opoo=fopen(dimage->waveNamen,"r"))==NULL){
      dimage->gediRat.useFootprint=1;
    }else{
      dimage->gediRat.useFootprint=0;
    }
    if(opoo){
      fclose(opoo);
      opoo=NULL;
    }
  }

  return;
}/*checkWaveOverwrite*/


/*##############################################*/
/*determine ground properties with e DEM*/

void groundFromDEM(pCloudStruct **data,control *dimage,waveStruct *waves)
{
  int nX=0,nY=0,nBins=0;
  float res=0,rRes=0;
  float *gWave=NULL;
  float *waveFromDEM(double *,int,int,float,double,double,double,double,float,float,double *,int *);
  double minX=0,minY=0,maxX=0,maxY=0;
  double *gDEM=NULL,minZ=0;
  double *findGroundPoly(pCloudStruct **,int,double *,double *,double *,double *,float,int *,int *,double);
  double *findGroundNN(pCloudStruct **,int,double *,double *,float,int *,int *,double);
  void groundProperties(float *,int,double,float,waveStruct *,float);

  res=0.1;
  rRes=0.15;

  /*make DEM*/
  if(dimage->gediRat.doGrid){
    minX=dimage->gediRat.minX;
    maxX=dimage->gediRat.maxX;
    minY=dimage->gediRat.minY;
    maxY=dimage->gediRat.maxY;
  }else{
    minX=maxX=100000000000.0;
    maxX=maxY=-100000000000.0;
  }

  if(dimage->polyGr)   gDEM=findGroundPoly(data,dimage->gediIO.nFiles,&minX,&minY,&maxX,&maxY,res,&nX,&nY,waves->groundBreakElev);
  else if(dimage->nnGr)gDEM=findGroundNN(data,dimage->gediIO.nFiles,&minX,&minY,res,&nX,&nY,waves->groundBreakElev);

  if(gDEM){
    /*make gap filled ground waveform*/
    gWave=waveFromDEM(gDEM,nX,nY,res,minX,minY,dimage->gediRat.coord[0],dimage->gediRat.coord[1],dimage->gediIO.fSigma,rRes,&minZ,&nBins);
    TIDY(gDEM);

    /*ground properties*/
    groundProperties(gWave,nBins,minZ,rRes,waves,dimage->gediIO.fSigma);
  }else{
    waves->gElev=waves->gSlope=-10000.0;
  }

  TIDY(gWave);
  return;
}/*groundFromDEM*/


/*####################################*/
/*ground properties*/

void groundProperties(float *gWave,int nBins,double minZ,float rRes,waveStruct *waves,float fSigma)
{
  int i=0;
  float total=0;
  double z=0,gStdev=0;
  char hasGround=0;

  /*CofG*/
  hasGround=1;
  waves->gElev=0.0;
  total=0.0;
  for(i=0;i<nBins;i++){
    z=(double)i*(double)rRes+minZ;
    waves->gElev+=z*(double)gWave[i];
    total+=gWave[i];
  }
  if(total>0.0){
    waves->gElev/=total;
    for(i=0;i<nBins;i++)gWave[i]/=total;
  }else{
    hasGround=0;
  }

  if(hasGround){
    /*stdev*/
    gStdev=total=0.0;
    for(i=0;i<nBins;i++){
      z=(double)i*(double)rRes+minZ;
      gStdev+=pow(z-waves->gElev,2.0)*(double)gWave[i];
      total+=gWave[i];
    }
    gStdev=sqrt(gStdev/total);
    waves->gSlope=atan2(sqrt(gStdev*gStdev),fSigma)*180.0/M_PI;
  }else waves->gSlope=waves->gElev=-10000.0;

  return;
}/*groundProperties*/


/*####################################*/
/*make waveform from DEM*/

float *waveFromDEM(double *gDEM,int nX,int nY,float res,double minX,double minY,double x0,double y0,float fSigma,float rRes,double *minZ,int *nBins)
{
  int i=0,j=0;
  int bin=0,place=0;
  float *gWave=NULL;
  float weight=0;
  double maxZ=0,x=0,y=0;
  double sep=0;

  /*mind min and max*/
  maxZ=-10000000.0;
  *minZ=1000000000.0;


  for(i=nX*nY-1;i>=0;i--){
    if(gDEM[i]<*minZ)*minZ=gDEM[i];
    if(gDEM[i]>maxZ)maxZ=gDEM[i];
  }

  *nBins=(int)((maxZ-*minZ)/rRes+0.5);
  gWave=falloc((uint64_t)(*nBins),"ground wave",3);
  for(i=*nBins-1;i>=0;i--)gWave[i]=0.0;

  /*add up DEM*/
  for(i=0;i<nX;i++){
    x=((double)i+0.5)*(double)res+minX;
    for(j=0;j<nY;j++){
      y=((double)j+0.5)*(double)res+minY;
      place=j*nX+i;

      sep=sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));
      bin=(int)((gDEM[place]-(*minZ))/rRes);
      weight=(float)gaussian(sep,(double)fSigma,0.0);
      if((bin>=0)&&(bin<(*nBins)))gWave[bin]+=weight;
    }
  }

  return(gWave);
}/*waveFromDEM*/


/*####################################*/
/*write name if overlap*/

void checkThisFile(lasFile *las,control *dimage,int i)
{
  if(checkFileBounds(las,dimage->gediRat.minX,dimage->gediRat.maxX,dimage->gediRat.minY,dimage->gediRat.maxY)){
    fprintf(stdout,"Need %s\n",dimage->inList[i]);
  }
  return;
}/*checkThisFile*/


/*####################################*/
/*write GEDI waveforms*/

void writeGEDIwave(control *dimage,waveStruct *waves,int numb)
{
  int i=0,j=0;
  float r=0;
  char waveID[300];
  FILE *opoo=NULL;


  /*make waveID if needed*/
  if(dimage->gediRat.doGrid){
    if(dimage->useID){
      if(dimage->gediRat.geoCoords){
        sprintf(waveID,"%s.%d.%d",dimage->waveID,(int)dimage->gediRat.geoCoords[numb][0],(int)dimage->gediRat.geoCoords[numb][1]);
      }else{
        sprintf(waveID,"%s.%d.%d",dimage->waveID,(int)dimage->gediRat.coord[0],(int)dimage->gediRat.coord[1]);
      }
    }
  }else if(dimage->gediRat.readALSonce){
    strcpy(waveID,dimage->gediRat.waveIDlist[numb]);
  }else{
    if(dimage->useID)strcpy(waveID,dimage->waveID);
  }


  if((opoo=fopen(dimage->waveNamen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",dimage->waveNamen);
    exit(1);
  }

  /*write header*/
  if(dimage->gediIO.ground==0)fprintf(opoo,"# 1 elevation, 2 discrete intensity, 3 discrete count, 4 discrete fraction, 5 ALS pulse, 6 ALS and GEDI pulse, 7 ind decon, 8 ind decon GEDI, 9 decon GEDI, 10 ind decon\n");
  else                 fprintf(opoo,"# 1 elevation, 2 discrete intensity, 3 int canopy, 4 int ground, 5 discrete count, 6 count canopy, 7 count ground, 8 discrete fraction, 9 fraction canopy, 10 fraction ground, 11 ALS pulse, 12 ALS and GEDI pulse, 13 ind decon, 14 ind decon GEDI, 15 decon GEDI, 16 ind decon\n");
  fprintf(opoo,"# fSigma %f pSigma %f res %f sideLobes %d\n",dimage->gediIO.fSigma,dimage->gediIO.pSigma,dimage->gediIO.res,dimage->gediRat.sideLobe);
  if(dimage->gediRat.geoCoords)fprintf(opoo,"# coord %.2f %.2f\n",dimage->gediRat.geoCoords[numb][0],dimage->gediRat.geoCoords[numb][1]);
  else                         fprintf(opoo,"# coord %.2f %.2f\n",dimage->gediRat.coord[0],dimage->gediRat.coord[1]);
  fprintf(opoo,"# density point %f beam %f\n",dimage->gediRat.pointDense,dimage->gediRat.beamDense);
  fprintf(opoo,"# meanScanAng %f\n",waves->meanScanAng);
  if(dimage->useID)fprintf(opoo,"# waveID %s\n",waveID);
  if(dimage->gediIO.ground&&(dimage->polyGr||dimage->nnGr))fprintf(opoo,"# ground %f %f\n",waves->gElev,waves->gSlope);
  if(dimage->gediIO.ground)fprintf(opoo,"# simpleGround %f\n",waves->gElevSimp);

  /*write data*/
  for(i=0;i<waves->nBins;i++){
    r=(float)waves->maxZ-(float)i*dimage->gediIO.res;

    fprintf(opoo,"%f",r);
    for(j=0;j<dimage->gediIO.nTypeWaves;j++){
      fprintf(opoo," %.10f",waves->wave[j][i]);
       if(dimage->gediIO.ground&&(j<3))fprintf(opoo," %f %f",waves->canopy[j][i],waves->ground[j][i]);
    }
    fprintf(opoo,"\n");
  }

  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Written to %s\n",dimage->waveNamen);
  return;
}/*writeGEDIwave*/


/*####################################*/
/*read ASCII data*/

pCloudStruct *readAsciiData(char *inNamen)
{
  uint64_t i=0;
  pCloudStruct *data=NULL;
  char line[400],temp1[100],temp2[100];
  char temp3[100],temp4[100],temp5[100];
  char temp6[100],temp7[100],temp8[100];
  FILE *ipoo=NULL;


  if(!(data=(pCloudStruct *)calloc(1,sizeof(pCloudStruct)))){
    fprintf(stderr,"error pCloudStruct allocation.\n");
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
}/*readAsciiData*/


/*##############################################*/
/*read command line*/

control *readCommands(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;
  void writeGediRatHelpMessage();

  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }

  dimage->gediIO.nFiles=1;
  dimage->inList=chChalloc(dimage->gediIO.nFiles,"inList",0);
  dimage->inList[0]=challoc(200,"inList",0);
  strcpy(&(dimage->inList[0][0]),"/Users/dill/data/teast/maryland_play/sc_79_112_1.las");
  strcpy(dimage->outNamen,"teast.wave");
  dimage->gediIO.pFWHM=15.0;   /*12 ns FWHM*/
  dimage->gediIO.fSigma=5.5;  /*86% of energy within a diameter of 20-25m*/
  dimage->gediIO.res=0.15;
  dimage->gediIO.pRes=dimage->gediIO.res/4.0;
  dimage->gediRat.coord[0]=624366.0;
  dimage->gediRat.coord[1]=3.69810*pow(10.0,6.0);
  dimage->gediRat.decon=NULL;
  dimage->gediIO.aEPSG=4326;      /*default is not to reproject*/

  /*switches*/
  dimage->gediRat.readWave=0;
  dimage->gediIO.ground=0;
  dimage->gediRat.sideLobe=0;   /*no side lobes*/
  dimage->gediRat.lobeAng=0.0;
  dimage->gediRat.pulseAfter=1;  /*smooth after for speed*/
  dimage->listFiles=0;
  dimage->pBuffSize=(uint64_t)200000000;
  dimage->gediRat.checkCover=0;
  dimage->gediRat.normCover=1;
  dimage->gediRat.cleanOut=0;
  dimage->gediRat.topHat=0;
  dimage->useID=0;
  dimage->gediIO.readPulse=0;
  dimage->gediRat.useShadow=0;
  dimage->gediRat.vRes[0]=dimage->gediRat.vRes[1]=dimage->gediRat.vRes[2]=1.0;
  dimage->gediRat.beamRad=0.165;    /*33 cm*/
  dimage->gediRat.maxScanAng=1000000.0;   /*maximum scan angle*/
  dimage->polyGr=0;             /*don't fit a polynomial through the ground*/
  dimage->nnGr=0;               /*don't make a DEM from nearest neighbour*/
  dimage->overWrite=1;          /*over write any files with the same name if they exist*/
  dimage->gediRat.readALSonce=0;/*read each footprint separately*/
  dimage->writeHDF=0;           /*write output as ascii*/
  dimage->writeL1B=0;           /*write output as ascii*/
  dimage->gediRat.defWfront=0;  /*Gaussian footprint*/
  dimage->gediRat.wavefront=NULL;
  dimage->gediIO.pcl=0;         /*do not use PCL*/
  dimage->gediIO.pclPhoton=0;         /*do not use PCL*/

  /*beams*/
  dimage->gediIO.useCount=dimage->gediIO.useFrac=dimage->gediIO.useInt=1;

  /*octree*/
  dimage->gediRat.useOctree=1;
  dimage->gediRat.octLevels=0;
  dimage->gediRat.nOctTop=40;   
  dimage->gediRat.octree=NULL;

  /*gridding options*/
  dimage->gediRat.doGrid=0;           /*gridded switch*/
  dimage->gediRat.gRes=30.0;          /*grid resolution*/
  dimage->gediRat.gNx=dimage->gediRat.gNy=0;
  /*batch*/
  dimage->gediRat.coords=NULL;       /*list of coordinates*/
  dimage->gediRat.waveIDlist=NULL;   /*list of waveform IDs*/
  dimage->maxBins=1024;
  dimage->gediIO.nMessages=200;
  strcpy(dimage->waveID,"gediWave");

  dimage->gediRat.iThresh=0.0006;
  dimage->gediRat.meanN=12.0;
  dimage->gediRat.doDecon=0;
  dimage->gediRat.indDecon=0;
  dimage->gediRat.decimate=1.0;  /*accept all ALS beams*/
  dimage->gediIO.pSigma=-1.0;  /*leave blannk for default GEDI*/

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        TTIDY((void **)dimage->inList,dimage->gediIO.nFiles);
        dimage->gediIO.nFiles=1;
        dimage->inList=chChalloc(dimage->gediIO.nFiles,"input name list",0);
        dimage->inList[0]=challoc((uint64_t)strlen(argv[++i])+1,"input name list",0);
        strcpy(dimage->inList[0],argv[i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(dimage->outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-inList",7)){
        checkArguments(1,i,argc,"-inList");
        TTIDY((void **)dimage->inList,dimage->gediIO.nFiles);
        dimage->inList=readInList(&dimage->gediIO.nFiles,argv[++i]);
      }else if(!strncasecmp(argv[i],"-coord",6)){
        checkArguments(2,i,argc,"-coord");
        dimage->gediRat.coord[0]=atof(argv[++i]);
        dimage->gediRat.coord[1]=atof(argv[++i]);
        dimage->gediRat.useOctree=0;    /*no point using octree for single*/
      }else if(!strncasecmp(argv[i],"-decon",6)){
        dimage->gediRat.doDecon=1;
      }else if(!strncasecmp(argv[i],"-indDecon",9)){
        dimage->gediRat.indDecon=1;
        dimage->gediRat.doDecon=1;
      }else if(!strncasecmp(argv[i],"-LVIS",5)){
        dimage->gediIO.pSigma=0.6893;  /*two way trip*/
        dimage->gediIO.fSigma=6.25;
      }else if(!strncasecmp(argv[i],"-pSigma",7)){
        checkArguments(1,i,argc,"-pSigma");
        dimage->gediIO.pSigma=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-pFWHM",6)){
        checkArguments(1,i,argc,"-pFWHM");
        dimage->gediIO.pFWHM=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-fSigma",7)){
        checkArguments(1,i,argc,"-fSigma");
        dimage->gediIO.fSigma=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-readWave",9)){
        dimage->gediRat.readWave=1;
      }else if(!strncasecmp(argv[i],"-ground",8)){
        dimage->gediIO.ground=1;
        dimage->gediRat.cleanOut=1;
      }else if(!strncasecmp(argv[i],"-sideLobe",9)){
        dimage->gediRat.sideLobe=1;
      }else if(!strncasecmp(argv[i],"-lobeAng",8)){
        checkArguments(1,i,argc,"-lobeAng");
        dimage->gediRat.lobeAng=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-listFiles",10)){
        dimage->listFiles=1;
      }else if(!strncasecmp(argv[i],"-pBuff",6)){
        checkArguments(1,i,argc,"-pBuff");
        dimage->pBuffSize=(uint64_t)(atof(argv[++i])*1000000000.0);
      }else if(!strncasecmp(argv[i],"-noNorm",7)){
        dimage->gediRat.normCover=0;
      }else if(!strncasecmp(argv[i],"-checkCover",11)){
        dimage->gediRat.checkCover=1;
      }else if(!strncasecmp(argv[i],"-topHat",7)){
        dimage->gediRat.topHat=1;
      }else if(!strncasecmp(argv[i],"-waveID",7)){
        checkArguments(1,i,argc,"-waveID");
        dimage->useID=1;
        strcpy(dimage->waveID,argv[++i]);
      }else if(!strncasecmp(argv[i],"-readPulse",10)){
        checkArguments(1,i,argc,"-readPulse");
        dimage->gediIO.readPulse=1;
        strcpy(dimage->gediIO.pulseFile,argv[++i]);
      }else if(!strncasecmp(argv[i],"-pulseAfter",11)){
        dimage->gediRat.pulseAfter=1;
      }else if(!strncasecmp(argv[i],"-pulseBefore",12)){
        dimage->gediRat.pulseAfter=0;
      }else if(!strncasecmp(argv[i],"-maxScanAng",11)){
        checkArguments(1,i,argc,"-maxScanAng");
        dimage->gediRat.maxScanAng=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-useShadow",10)){
        dimage->gediRat.useShadow=1;
      }else if(!strncasecmp(argv[i],"-polyGround",11)){
        dimage->polyGr=1;
      }else if(!strncasecmp(argv[i],"-nnGround",99)){
        dimage->nnGr=1;
      }else if(!strncasecmp(argv[i],"-res",4)){
        checkArguments(1,i,argc,"-res");
        dimage->gediIO.res=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-gridBound",10)){
        checkArguments(4,i,argc,"-gridBound");
        dimage->gediRat.doGrid=1;
        dimage->useID=1;
        dimage->gediRat.gMinX=atof(argv[++i]);
        dimage->gediRat.gMaxX=atof(argv[++i]);
        dimage->gediRat.gMinY=atof(argv[++i]);
        dimage->gediRat.gMaxY=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-gridStep",9)){
        checkArguments(1,i,argc,"-gridStep");
        dimage->gediRat.gRes=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-keepOld",8)){
        dimage->overWrite=0;
      }else if(!strncasecmp(argv[i],"-listCoord",10)){
        checkArguments(1,i,argc,"-listCoord");
        dimage->gediRat.readALSonce=1;
        dimage->useID=1;
        strcpy(dimage->gediRat.coordList,argv[++i]);
      }else if(!strncasecmp(argv[i],"-hdf",4)){
        dimage->writeHDF=1;
      }else if(!strncasecmp(argv[i],"-l1b",4)){
        dimage->writeHDF=1;
        dimage->writeL1B=1;
      }else if(!strncasecmp(argv[i],"-ascii",6)){
        dimage->writeHDF=0;
      }else if(!strncasecmp(argv[i],"-maxBins",8)){
        checkArguments(1,i,argc,"-maxBins");
        dimage->maxBins=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-wavefront",10)){
        checkArguments(1,i,argc,"-wavefront");
        dimage->gediRat.defWfront=1;
        dimage->gediRat.wavefront=copyFrontFilename(argv[++i]);
      }else if(!strncasecmp(argv[i],"-noOctree",9)){
        dimage->gediRat.useOctree=0;
      }else if(!strncasecmp(argv[i],"-octLevels",9)){
        checkArguments(1,i,argc,"-octLevels");
        dimage->gediRat.octLevels=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-nOctPix",8)){
        checkArguments(1,i,argc,"-nOctPix");
        dimage->gediRat.nOctTop=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-countOnly",10)){
        dimage->gediIO.useCount=1;
        dimage->gediIO.useInt=0;
        dimage->gediIO.useFrac=0;
      }else if(!strncasecmp(argv[i],"-seed",5)){
        checkArguments(1,i,argc,"-seed");
        srand(atoi(argv[++i]));
      }else if(!strncasecmp(argv[i],"-decimate",9)){
        checkArguments(1,i,argc,"-decimate");
        dimage->gediRat.decimate=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-pcl",4)){
        dimage->gediIO.pcl=1;
      }else if(!strncasecmp(argv[i],"-aEPSG",6)){
        checkArguments(1,i,argc,"-aEPSG");
        dimage->gediIO.aEPSG=atoi(argv[++i]);;
      }else if(!strncasecmp(argv[i],"-help",5)){
        writeGediRatHelpMessage();
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry gediRat -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }/*command parser*/


  /*total number of beams*/
  dimage->gediIO.nTypeWaves=dimage->gediIO.useCount+dimage->gediIO.useFrac+dimage->gediIO.useInt;

  /*ensure pulse is sampled at same rate as waveform*/
  if(dimage->gediIO.readPulse==0)dimage->gediIO.pRes=dimage->gediIO.res;

  return(dimage);
}/*readCommands*/


/*####################################*/
/*Help message*/

void writeGediRatHelpMessage()
{
        fprintf(stdout,"\n#####\nProgram to create GEDI waveforms from ALS las or pts files. laz not yet supported\n#####\n\n\
\n# Input output filenames and format\n\
-input name;     lasfile input filename\n\
-inList list;    input file list (ASCII file) for multiple files\n\
-output name;    output filename\n\
-ground;         record separate ground and canopy waveforms\n\
-hdf;            write output as HDF5. Best with gridded or list of coords\n\
-l1b;            write output in the GEDI L1B HDF5 format. Best with gridded or list of coords\n\
-aEPSG epsg;     input EPSG code if the L1B output needs to be in degrees\n\
-ascii;          write output as ASCII (default). Good for quick tests\n\
-waveID id;      supply a waveID to pass to the output (only for single footprints)\n\
\n# Single footprint, list of footprints, or grid of footprints\n\
-coord lon lat;  footprint coordinate in same system as lasfile\n\
-listCoord name; list of coordinates\n\
-gridBound minX maxX minY maxY;     make a grid of waveforms in this box\n\
-gridStep res;   grid step size\n\
\n# Lidar characteristics. Defaults are expected GEDI values.\n\
-pSigma sig;     set Gaussian pulse width as 1 sigma\n\
-pFWHM fhwm;     set Gaussian pulse width as FWHM in ns. This is the outgoing laser pulse, which will be twice the detected pulse width\n\
-readPulse file; read pulse shape and width from a file insteda of making Gaussian\n\
-fSigma sig;     set footprint width\n\
-wavefront file; read wavefront shape from file instead of setting Gaussian. Note that footprint width is still set by fSigma\n\
-res res;        range resolution of waveform digitisation to output, in units of ALS data\n\
-LVIS;           use LVIS pulse length, sigma=6.25m\n\
-topHat;         use a top hat wavefront\n\
-sideLobe;       use side lobes\n\
-lobeAng ang;    lobe axis azimuth\n\
-pcl;            pulse comression lidar. Do not pad waveform\n\
\n# Input data quality filters\n\
-checkCover;     check that the footprint is covered by ALS data. Do not output if not\n\
-maxScanAng ang; maximum scan angle, degrees\n\
-decimate x;     probability of accepting an ALS beam\n\
\n# Computational speed options\n\
-pBuff s;        point reading buffer size in Gbytes\n\
-maxBins;        Optional: for HDF5, limit number of bins to save trimming.\n\
-countOnly;      only use count method\n\
-pulseAfter;     apply the pulse smoothing after binning for computational speed, at the risk of aliasing (default)\n\
-pulseBefore;    apply the pulse smoothing before binning to avoid the risk of aliasing, at the expense of computational speed\n\
-noNorm;         don't normalise for ALS density\n\
\n# Octree\n\
-noOctree;       do not use an octree\n\
-octLevels n;    number of octree levels to use\n\
-nOctPix n;      number of octree pixels along a side for the top level\n\
\n# Using full-waveform input data (not tested)\n\
-decon;          deconvolve\n\
-indDecon;       deconvolve individual beams\n\
-readWave;       read full-waveform where available\n\
\n# Miscellaneous\n\
-listFiles;      list files. Do not read them\n\
-keepOld;        do not overwrite old files, if they exist\n\
-useShadow;      account for shadowing in discrete return data through voxelisation\n\
-polyGround;     find mean ground elevation and slope through fitting a polynomial\n\
-nnGround;       find mean ground elevation and slope through nearest neighbour\n\
-seed n;         random number seed\n\n\nQuestions to svenhancock@gmail.com\n\n");
  return;
}/*writeGediRatHelpMessage*/


/*the end*/
/*####################################*/

