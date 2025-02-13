#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "libLasRead.h"
#include "libLasProcess.h"


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
/*lobe structure*/

typedef struct{
  double coord[2];  /*central coordinate*/
  float E;          /*fraction of energy here*/
  float fSigma;     /*footprint sigma*/
  double maxSepSq;  /*maximum distance from footprint needed*/
}lobeStruct;

/*####################################*/
/*control structure*/

typedef struct{
  char **inList;
  int nFiles;
  char outNamen[1000];

  double footRes;     /*footprint grid resolution*/

  /*options*/
  char ground;         /*only use ground points*/
  char readWave;       /*read waveform switch*/
  char checkCover;     /*check that the whole footprit is covered by data*/
  char cleanOut;       /*clean subterranean outliers*/
  uint64_t pBuffSize;  /*point buffer rading size in bytes*/

  /*GEDI fotprint parameters*/
  char sideLobe;     /*side lobe switch*/
  float lobeAng;     /*lobe major axis, degrees*/
  int nLobes;        /*number of side lobes*/
  lobeStruct *lobe;  /*lobe structure*/
  double minX;       /*minimum latitude of interest*/
  double maxX;       /*maximum latitude of interest*/
  double minY;       /*minimum longitude of interest*/
  double maxY;       /*maximum longitude of interest*/
  
  double coord[2];
  float res;       /*range resolution*/
  float pFWHM;     /*pulse width in ns*/
  float pSigma;    /*pulse width in metres*/
  float fWidth;    /*footprint width*/
  float fSigma;    
  float pRes;
  pulseStruct *pulse;

  /*tolerances*/
  float iThresh;   /*intensity threshold*/
  float meanN;
  char doDecon;    /*deconolution switch*/
  char indDecon;   /*deconolve individual ALS waves*/
  denPar *decon;   /*denoising parameters*/

  /*grid to normalise sampling density*/
  int *nGrid;
  int gX;
  int gY;
  float gridRes;
  double g0[2];  /*grid origin*/
}control;


/*####################################*/
/*grid structure*/

typedef struct{
  /*global grid*/
  uint64_t nX;
  uint64_t nY;
  uint64_t nDone;
  char *done;
  double minX;
  double minY;
  double res;   /*grid resolution*/
  float rRes;   /*range resolution*/

  /*file map*/
  int nWithin;   /*number of files intersecting*/
  int *map;      /*file map*/

  /*within file*/
  int nBins;
  double minZ;
  double maxZ;
  int wX;   /*number of waves in x*/
  int wY;   /*number of waves in y*/
  float **canopy;   /*canopy waveform*/
  float **ground;   /*ground waveform*/
}gridStruct;


/*####################################*/
/*data structure*/

typedef struct{
  uint32_t nPoints;
  double *x;
  double *y;
  double *z;
  unsigned char *class;     /*point classification*/
  int *refl;
  char *retNumb;            /*discrete return number*/
  unsigned char *packetDes; /*waveform or not*/
  float **grad;             /*Poynting vector*/
  float *time;              /*time in picoseconds of this wave*/
  uint64_t *waveMap;        /*pointer to waveform in file*/
  uint32_t *waveLen;        /*length of waveform in bins*/
  uint64_t waveStart;       /*offset to waveform data*/
  FILE *ipoo;               /*file pointer*/
  char hasWave;             /*waveform included*/
}dataStruct;


/*####################################*/
/*waveform structure*/

typedef struct{
  float **wave;  /*waveforms*/
  float **canopy;/*canopy waveform*/
  float **ground;/*ground waveform*/
  double minZ;   /*elevation bounds*/
  double maxZ;   /*elevation bounds*/
  int nBins;     /*number of wave bins*/
  int nWaves;    /*number of different ways*/
}waveStruct;


/*####################################*/
/*main*/

int main(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;
  control *readCommands(int,char **);
  lasFile **las=NULL;
  dataStruct **data=NULL;
  dataStruct *readALSdata(lasFile *,control *);
  waveStruct *waves=NULL;
  waveStruct *makeGediWaves(control *,dataStruct **);
  gridStruct *grid=NULL;
  gridStruct *setUpGrid(lasFile **,control *);
  void writeGEDIwave(control *,waveStruct *);
  void setGediFootprint(control *);
  void setGediPulse(control *);
  void tidySMoothPulse();
  void checkThisFile(lasFile *,control *,int);
  void setGediGrid(control *,lasFile **,int,gridStruct *);
 
  /*read command line*/
  dimage=readCommands(argc,argv);


  if(!(las=(lasFile **)calloc(dimage->nFiles,sizeof(lasFile *)))){
    fprintf(stderr,"error waveStruct allocation.\n");
    exit(1);
  }
  /*read headers and bounds*/
  for(i=0;i<dimage->nFiles;i++){
    las[i]=readLasHead(dimage->inList[i],dimage->pBuffSize);
    if(las[i]->ipoo){
      fclose(las[i]->ipoo);
      las[i]->ipoo=NULL;
    }
  }


  /*set up grid parameters*/
  grid=setUpGrid(las,dimage);

  /*set up the pulse*/
  setGediPulse(dimage);

  /*allocate space for footprint*/
  setGediFootprint(dimage);


  /*loop over las files*/
  for(i=0;i<dimage->nFiles;i++){
    /*set up footprint*/
    setGediGrid(dimage,las,i,grid);

    /*read data or write filename if needed*/
    //if(dimage->listFiles==0)data[i]=readALSdata(las,dimage);
    //else                    checkThisFile(las,dimage,i);

    /*tidy up*/
    TTIDY((void **)grid->canopy,grid->wX*grid->wY);
    TTIDY((void **)grid->ground,grid->wX*grid->wY);
  }/*file loop*/

  /*tidy up*/
  if(las){
    for(i=0;i<dimage->nFiles;i++)las[i]=tidyLasFile(las[i]);
    TIDY(las);
  }
  if(grid){
    TIDY(grid->done);
    TIDY(grid);
  }
  if(data){
    for(i=0;i<dimage->nFiles;i++){
      if(data[i]){
        TIDY(data[i]->x);
        TIDY(data[i]->y);
        TIDY(data[i]->z);
        TIDY(data[i]->refl);
        TIDY(data[i]->class);
        TIDY(data[i]->retNumb);
        TIDY(data[i]->packetDes);
        TTIDY((void **)data[i]->grad,3);
        data[i]->grad=NULL;
        TIDY(data[i]->time);
        TIDY(data[i]->waveMap);
        TIDY(data[i]->waveLen);
        if(data[i]->ipoo){
          fclose(data[i]->ipoo);
          data[i]->ipoo=NULL;
        }
      }
      TIDY(data[i]);
    }
    TIDY(data);
  }
  if(waves){
    TTIDY((void **)waves->wave,waves->nWaves);
    TTIDY((void **)waves->canopy,2);
    TTIDY((void **)waves->ground,2);
    TIDY(waves);
  }
  if(dimage){
    if(dimage->pulse){
      TIDY(dimage->pulse->y);
      TIDY(dimage->pulse->x);
      TIDY(dimage->pulse);
    }
    TIDY(dimage->lobe);
    TTIDY((void **)dimage->inList,dimage->nFiles);
    TIDY(dimage->nGrid);
    TIDY(dimage);
  }
  tidySMoothPulse();
  return(0);
}/*main*/


/*####################################*/
/*set footprint grid for this file*/

void setGediGrid(control *dimage,lasFile **las,int i,gridStruct *grid)
{
  int j=0;
  uint32_t k=0;
  double x=0,y=0,z=0;
  double dx=0,dy=0;


  /*determine how many files overlap. Need to include a buffer*/
/*include a buffer from the footprints*/
  grid->map=markInt(grid->nWithin,grid->map,i);
  grid->nWithin=1;
  for(j=0;j<dimage->nFiles;j++){
    if(i==j)continue;
    if((las[j]->minB[0]<=las[i]->maxB[0])&&(las[j]->minB[1]<=las[i]->maxB[1])&&(las[j]->maxB[0]>=las[i]->minB[0])&&(las[j]->maxB[1]>=las[i]->minB[1])){
      grid->map=markInt(grid->nWithin,grid->map,i);
      grid->nWithin++;
    }
  }/*determine how many files overlap*/

  /*determine bounds*/
  grid->maxZ=-1000000000.0;
  grid->minZ=1000000000.0;
  for(j=0;j<grid->nWithin;j++){
    if(las[j]->minB[2]<grid->minZ)grid->minZ=las[j]->minB[2];
    if(las[j]->maxB[2]<grid->maxZ)grid->maxZ=las[j]->maxB[2];
  }/*determine bounds*/
  

  /*allocate waveforms*/
  grid->nBins=(int)((grid->maxZ-grid->minZ)/(double)grid->rRes);
  grid->wX=(int)((las[i]->maxB[0]-las[i]->minB[0])/grid->res);
  grid->wY=(int)((las[i]->maxB[1]-las[i]->minB[1])/grid->res);
  grid->canopy=fFalloc(grid->wX*grid->wY,"canopy",0);
  grid->ground=fFalloc(grid->wX*grid->wY,"ground",0);
  for(j=grid->wX*grid->wY-1;j>=0;j--){
    grid->canopy[i]=falloc((uint64_t)grid->nBins,"canopy",j+1);
    grid->ground[i]=falloc((uint64_t)grid->nBins,"ground",j+1);
  }


/*make a grid of coordinates here*/



  /*loop over points and assign to waveforms*/
  for(j=0;j<grid->nWithin;j++){
    for(k=0;k<las[j]->nPoints;k++){
      readLasPoint(las[k],k);
      setCoords(&x,&y,&z,las[j]);


      /*see which footprints it is within*/

      //dX=dimage->coord[0]-x;
      //dY=dimage->coord[1]-y;
      //sepSq=dX*dX+dY*dY;

/*mark as done if it is done*/


    }/*point loop*/
  }/*file loop*/


  return;
}/*setGediGrid*/


/*####################################*/
/*set up grid*/

gridStruct *setUpGrid(lasFile **las,control *dimage)
{
  int i=0;
  double maxX=0,maxY=0;
  gridStruct *grid=NULL;

  if(!(grid=(gridStruct *)calloc(dimage->nFiles,sizeof(gridStruct)))){
    fprintf(stderr,"error waveStruct allocation.\n");
    exit(1);
  }

  grid->minX=grid->minY=1000000000000.0;
  maxX=maxY=-10000000000.0;
  for(i=0;i<dimage->nFiles;i++){
    if(las[i]->minB[0]<grid->minX)grid->minX=las[i]->minB[0];
    if(las[i]->minB[1]<grid->minY)grid->minY=las[i]->minB[1];
    if(las[i]->maxB[0]>maxX)maxX=las[i]->maxB[0];
    if(las[i]->maxB[1]>maxY)maxY=las[i]->maxB[1];
  }
  grid->done=NULL;
  grid->nDone=0;
  grid->res=dimage->footRes;
  grid->rRes=0.15;

  grid->nX=(uint64_t)((maxX-grid->minX)/grid->res);
  grid->nY=(uint64_t)((maxX-grid->minX)/grid->res);

  return(grid);
}/*setUpGrid*/


/*####################################*/
/*write name if overlap*/

void checkThisFile(lasFile *las,control *dimage,int i)
{
  if(checkFileBounds(las,dimage->minX,dimage->maxX,dimage->minY,dimage->maxY)){
    fprintf(stdout,"Need %s\n",dimage->inList[i]);
  }

  return;
}/*checkThisFile*/


/*####################################*/
/*read lasFile and save relevant data*/

dataStruct *readALSdata(lasFile *las,control *dimage)
{
  int j=0;
  int gX=0,gY=0;
  uint32_t i=0;
  uint32_t pUsed=0;    /*number of points used*/
  double x=0,y=0,z=0;
  double dX=0,dY=0,sepSq=0;
  dataStruct *data=NULL;
  char hasWave=0;   /*has waveform data, to save RAM*/

  /*allocate maximum number of points*/
  if(!(data=(dataStruct *)calloc(1,sizeof(dataStruct)))){
    fprintf(stderr,"error dataStruct allocation.\n");
    exit(1);
  }

  if(checkFileBounds(las,dimage->minX,dimage->maxX,dimage->minY,dimage->maxY)){
    data->x=dalloc(las->nPoints,"x",0);
    data->y=dalloc(las->nPoints,"y",0);
    data->z=dalloc(las->nPoints,"z",0);
    data->refl=ialloc(las->nPoints,"refl",0);
    data->class=uchalloc((uint64_t)las->nPoints,"class",0);
    data->retNumb=challoc((uint64_t)las->nPoints,"retNumb",0);
    data->packetDes=uchalloc((uint64_t)las->nPoints,"packetDes",0);
    data->grad=fFalloc(3,"grad",0);
    for(i=0;i<3;i++)data->grad[i]=falloc((uint64_t)las->nPoints,"grad",i+1);
    data->time=falloc((uint64_t)las->nPoints,"time",0);              /*time in picoseconds of this wave*/
    if(!(data->waveMap=(uint64_t *)calloc(las->nPoints,sizeof(uint64_t)))){
      fprintf(stderr,"error in input filename structure.\n");
      exit(1);
    }
    if(!(data->waveLen=(uint32_t *)calloc(las->nPoints,sizeof(uint32_t)))){
      fprintf(stderr,"error in input filename structure.\n");
      exit(1);
    }

    /*loop over points*/
    pUsed=0;
    hasWave=0;
    for(i=0;i<las->nPoints;i++){
      /*read one point*/
      readLasPoint(las,i);
      setCoords(&x,&y,&z,las);
      dX=dimage->coord[0]-x;
      dY=dimage->coord[1]-y;
      sepSq=dX*dX+dY*dY;

      if((x>=dimage->minX)&&(x<=dimage->maxX)&&(y>=dimage->minY)&&(y<=dimage->maxY)&&(z>-10000.0)&&(z<10000.0)&&((sepSq<=dimage->lobe[0].maxSepSq)||(dimage->nLobes>1))){  /*the point is of use*/
        data->x[pUsed]=x;
        data->y[pUsed]=y;
        data->z[pUsed]=z;
        if(las->refl>0)data->refl[pUsed]=(int)las->refl;
        else           data->refl[pUsed]=1;
        data->class[pUsed]=las->classif;

        if(checkOneWave(las)){
          hasWave=1;
          data->retNumb[pUsed]=(char)las->retNumb;
          data->packetDes[pUsed]=las->packetDes;
          for(j=0;j<3;j++)data->grad[j][pUsed]=las->grad[j];
          data->time[pUsed]=las->time;
          data->waveMap[pUsed]=las->waveMap;
          data->waveLen[pUsed]=las->waveLen;
        }else data->packetDes[pUsed]=0;
        pUsed++;
      }

      if(las->retNumb==las->nRet){  /*only once per beam*/
        /*mark sampling desnity for normalisation*/
        gX=(int)((x-dimage->g0[0])/(double)dimage->gridRes);
        gY=(int)((y-dimage->g0[1])/(double)dimage->gridRes);
        if((gX>=0)&&(gX<dimage->gX)&&(gY>=0)&&(gY<dimage->gY)){
          dimage->nGrid[gY*dimage->gX+gX]++;
        }
      }
    }/*point loop*/

    /*trim data arrays*/
    data->nPoints=pUsed;
    if(pUsed>0){
      if(!(data->x=(double *)realloc(data->x,data->nPoints*sizeof(double)))){
        fprintf(stderr,"Balls\n");
        exit(1);
      }
      if(!(data->y=(double *)realloc(data->y,data->nPoints*sizeof(double)))){
        fprintf(stderr,"Balls\n");
        exit(1);
      }
      if(!(data->z=(double *)realloc(data->z,data->nPoints*sizeof(double)))){
        fprintf(stderr,"Balls\n");
        exit(1);
      }
      if(!(data->refl=(int *)realloc(data->refl,data->nPoints*sizeof(int)))){
        fprintf(stderr,"Balls\n");
        exit(1);
      }
      if(!(data->class=(unsigned char *)realloc(data->class,data->nPoints*sizeof(unsigned char)))){
        fprintf(stderr,"Balls\n");
        exit(1);
      }
    }else{
      TIDY(data->x);
      TIDY(data->y);
      TIDY(data->z);
      TIDY(data->refl);
      TIDY(data->class);
    }
    if(hasWave==1){
      data->waveStart=las->waveStart;
      if(!(data->retNumb=(char *)realloc(data->retNumb,data->nPoints*sizeof(char)))){
        fprintf(stderr,"Balls\n");
        exit(1);
      }
      if(!(data->packetDes=(unsigned char *)realloc(data->packetDes,data->nPoints*sizeof(unsigned char)))){
        fprintf(stderr,"Balls\n");
        exit(1);
      }
      for(i=0;i<3;i++){
        if(!(data->grad[i]=(float *)realloc(data->grad[i],data->nPoints*sizeof(float)))){
          fprintf(stderr,"Balls\n");
          exit(1);
        }
      }
      if(!(data->time=(float *)realloc(data->time,data->nPoints*sizeof(float)))){
        fprintf(stderr,"Balls\n");
        exit(1);
      }
      if(!(data->waveMap=(uint64_t *)realloc(data->waveMap,data->nPoints*sizeof(uint64_t)))){
        fprintf(stderr,"Balls\n");
        exit(1);
      }
      if(!(data->waveLen=(uint32_t *)realloc(data->waveLen,data->nPoints*sizeof(uint32_t)))){
        fprintf(stderr,"Balls\n");
        exit(1);
      }
    }else{  /*clear out all the waveform bits*/
      TIDY(data->retNumb);
      TIDY(data->packetDes);
      TTIDY((void **)data->grad,3);
      data->grad=NULL;
      TIDY(data->time);
      TIDY(data->waveMap);
      TIDY(data->waveLen);
    }
  }else{/*file bounds check*/
    data->nPoints=0;
    data->retNumb=NULL;
    data->packetDes=NULL;
    data->grad=NULL;
    data->time=NULL;
    data->waveMap=NULL;
    data->waveLen=NULL;
    data->x=NULL;
    data->y=NULL;
    data->z=NULL;
    data->refl=NULL;
  }

  data->hasWave=hasWave;
  if(dimage->readWave&&hasWave){  /*only leave files open for waveform*/
    data->ipoo=las->ipoo;
    las->ipoo=NULL;
  }else{
    data->ipoo=NULL;
  }

  return(data);
}/*readALSdata*/


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

  if(dimage->ground==0)fprintf(opoo,"# 1 elevation, 2 discrete intensity, 3 discrete count, 4 ALS pulse, 5 ALS and GEDI pulse, 6 ind decon, 7 ind decon GEDI, 8 decon GEDI, 9 ind decon\n");
  else                 fprintf(opoo,"# 1 elevation, 2 discrete intensity, 3 int canopy, 4 int ground, 5 discrete count, 6 count canopy, 7 count ground, 8 ALS pulse, 9 ALS and GEDI pulse, 10 ind decon, 11 ind decon GEDI, 12 decon GEDI, 13 ind decon\n");


  for(i=0;i<waves->nBins;i++){
    r=(float)waves->maxZ-(float)i*dimage->res;

    fprintf(opoo,"%f",r);
    for(j=0;j<waves->nWaves;j++){
      fprintf(opoo," %f",waves->wave[j][i]);
       if(dimage->ground&&(j<2))fprintf(opoo," %f %f",waves->canopy[j][i],waves->ground[j][i]);
    }
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

waveStruct *makeGediWaves(control *dimage,dataStruct **data)
{
  int numb=0,bin=0,j=0,k=0;
  int gX=0,gY=0,n=0;
  uint32_t i=0;
  double sep=0;
  double dX=0,dY=0;
  float refl=0,rScale=0;
  float tot=0;
  waveStruct *waves=NULL;
  waveStruct *allocateGEDIwaves(control *,dataStruct **);
  void gediFromWaveform(dataStruct *,uint32_t,float,waveStruct *,control *);
  void processAggragate(control *,waveStruct *);
  void checkFootCovered(control *);
  void cleanOutliers(waveStruct *,control *);
  denPar *setDeconForGEDI(control *);


  /*check that whole footprint is covered*/
  if(dimage->checkCover)checkFootCovered(dimage);

  /*allocate*/
  waves=allocateGEDIwaves(dimage,data);

  /*set up denoising if using*/
  if(dimage->doDecon)dimage->decon=setDeconForGEDI(dimage);

  /*make waves*/
  for(n=0;n<dimage->nLobes;n++){
    for(numb=0;numb<dimage->nFiles;numb++){
      for(i=0;i<data[numb]->nPoints;i++){
        dX=data[numb]->x[i]-dimage->lobe[n].coord[0];
        dY=data[numb]->y[i]-dimage->lobe[n].coord[1];
        sep=sqrt(dX*dX+dY*dY);

        rScale=(float)gaussian(sep,(double)dimage->lobe[n].fSigma,0.0);

        if(rScale>dimage->iThresh){  /*if bright enough to matter*/
          /*scale by sampling density*/
          gX=(int)((data[numb]->x[i]-dimage->g0[0])/(double)dimage->gridRes);
          gY=(int)((data[numb]->y[i]-dimage->g0[1])/(double)dimage->gridRes);
          if((gX>=0)&&(gX<dimage->gX)&&(gY>=0)&&(gY<dimage->gY)){
            if(dimage->nGrid[gY*dimage->gX+gX]>0)rScale/=(float)dimage->nGrid[gY*dimage->gX+gX];
          }

          /*discrete return*/
          refl=(float)data[numb]->refl[i]*rScale;
          for(j=0;j<dimage->pulse->nBins;j++){
            bin=(int)((waves->maxZ-data[numb]->z[i]+(double)dimage->pulse->x[j])/(double)dimage->res);
            if((bin>=0)&&(bin<waves->nBins)){
              waves->wave[0][bin]+=refl*dimage->pulse->y[j];
                waves->wave[1][bin]+=rScale*dimage->pulse->y[j];
                if(dimage->ground){
                if(data[numb]->class[i]==2){
                  waves->ground[0][bin]+=refl*dimage->pulse->y[j];
                  waves->ground[1][bin]+=rScale*dimage->pulse->y[j];
                }else{
                  waves->canopy[0][bin]+=refl*dimage->pulse->y[j];
                  waves->canopy[1][bin]+=rScale*dimage->pulse->y[j];     
                }
              }/*ground recording if needed*/
            }/*bin bound check*/
          }/*pulse bin loop*/

          /*full-waveform*/
          if(dimage->readWave&&data[numb]->hasWave){
            if(data[numb]->packetDes[i]){  /*test for waveform*/
              gediFromWaveform(data[numb],i,rScale,waves,dimage);
            }
          }/*waveform test*/
        }
      }/*point loop*/
    }/*file loop*/
  }/*lobe loop*/

  /*clean outliers if needed*/
  if(dimage->cleanOut)cleanOutliers(waves,dimage);

  /*deconvolve aggragated waveform*/
  if(dimage->doDecon)processAggragate(dimage,waves);

  /*normalise integral*/
  for(k=0;k<waves->nWaves;k++){
    tot=0.0;
    for(j=0;j<waves->nBins;j++)tot+=waves->wave[k][j]*dimage->res;
    if(tot>0.0){
      for(j=0;j<waves->nBins;j++)waves->wave[k][j]/=tot;
      if(dimage->ground&&(k<2)){
        for(j=0;j<waves->nBins;j++){
          waves->canopy[k][j]/=tot;
          waves->ground[k][j]/=tot;
        }
      }
    }
  }

  if(dimage->decon){
    TTIDY((void **)dimage->decon->pulse,2);
    TIDY(dimage->decon);
  }
  return(waves);
}/*makeGediWaves*/


/*####################################*/
/*clean outlier points from waveform*/

void cleanOutliers(waveStruct *waves,control *dimage)
{
  int i=0,j=0,gStart=0;
  char pastGround=0;
  float gGap=0;  /*gap in ground return*/
  float maxGap=0;
  float maxGround=0,gThresh=0;
  float max=0,thresh=0;

  if(!dimage->ground){
    fprintf(stderr,"No need to clean without ground\n");
    exit(1);
  }

  maxGap=10.0;  /*maximum permittable gap in the ground return*/
  gGap=0.0;
  pastGround=0;

  /*determine max ground and max return*/
  maxGround=max=0.0;
  for(i=0;i<waves->nBins;i++){
    if(waves->ground[1][i]>maxGround)maxGround=waves->ground[1][i];
    if(waves->wave[1][i]>max)max=waves->wave[1][i];
  }
  gThresh=maxGround*0.01;
  thresh=max*0.001;

  for(i=0;i<waves->nBins;i++){
    if(waves->ground[1][i]>=gThresh){
      if(pastGround==0)gStart=i;
      pastGround=1;
      gGap=0.0;
    }else{
      if(pastGround)gGap+=dimage->res;
    }
    if(gGap>maxGap){  /*too big a break, delete*/
      for(;i<waves->nBins;i++){
        waves->canopy[0][i]=waves->canopy[1][i]=waves->ground[0][i]=waves->ground[1][i]=0.0;
        for(j=0;j<waves->nWaves;j++)waves->wave[j][i]=0.0;
      }
    }/*too big a break, delete*/
  }

  /*look for above canopy outliers*/
  gGap=0.0;
  maxGap=50.0;
  for(i=gStart;i>=0;i--){
    if(waves->wave[1][i]>=thresh)gGap=0.0;
    gGap+=dimage->res;

    if(gGap>=maxGap){
      for(;i>=0;i--){
        waves->canopy[0][i]=waves->canopy[1][i]=waves->ground[0][i]=waves->ground[1][i]=0.0;
        for(j=0;j<waves->nWaves;j++)waves->wave[j][i]=0.0;
      }
    }
  }

  return;
}/*cleanOutliers*/


/*####################################*/
/*check footprint is covered by ALS*/

void checkFootCovered(control *dimage)
{
  int i=0,j=0,nWithin=0;
  int thresh=0,nMissed=0;
  double dX=0,dY=0,sepSq=0;
  double useRad=0,radSq=0;

  useRad=10.0;
  if(useRad>(double)dimage->fSigma)useRad=(double)dimage->fSigma;
  radSq=useRad*useRad;

  for(i=0;i<dimage->gX;i++){
    dX=(double)i*(double)dimage->gridRes-dimage->coord[0];
    for(j=0;j<dimage->gY;j++){
      dY=(double)j*(double)dimage->gridRes-dimage->coord[1];
      sepSq=dX*dX+dY*dY;

      if(sepSq<radSq){
        if(dimage->nGrid[j*dimage->gX+i]==0)nMissed++;
        nWithin++;
      }
    }/*y loop*/
  }/*x loop*/

  thresh=(int)((float)nWithin*2.0/3.0);
  if(nMissed>thresh){
    fprintf(stderr,"Too many missed %d of %d\n",nMissed,nWithin);
    exit(1);
  }

  return;
}/*checkFootCovered*/


/*####################################*/
/*deconvolve aggragated wave*/

void processAggragate(control *dimage,waveStruct *waves)
{
  int i=0;
  float *processed=NULL,*smooPro=NULL;
  float *smooth(float,int,float *,float);

  /*Add background noise back*/
  for(i=0;i<waves->nBins;i++)waves->wave[6][i]+=dimage->meanN;

  /*deconvolve and reconvolve*/
  processed=processFloWave(&(waves->wave[6][0]),(int)waves->nBins,dimage->decon,1.0);
  for(i=0;i<waves->nBins;i++)waves->wave[7][i]=processed[i];
  smooPro=smooth(dimage->pSigma,(int)waves->nBins,processed,dimage->res);
  TIDY(processed);

  TIDY(waves->wave[6]);
  waves->wave[6]=smooPro;
  smooPro=NULL;

  return;
}/*processAggragate*/


/*####################################*/
/*GEDI wave from ALS waveforms*/

void gediFromWaveform(dataStruct *data,uint32_t i,float rScale,waveStruct *waves,control *dimage)
{
  int j=0,bin=0;
  int buffBins=0;
  uint32_t waveLen=0;
  float grad[3],*smoothed=NULL,*floWave=NULL;
  float *processed=NULL,*smooPro=NULL;
  float *smooth(float,int,float *,float);
  double x=0,y=0,z=0;
  unsigned char *wave=NULL,*temp=NULL;

  for(j=0;j<3;j++)grad[j]=data->grad[j][i];
  wave=readLasWave(data->waveMap[i],data->waveLen[i],data->ipoo,data->waveStart);

  /*buffer to give space for smoothing*/
  buffBins=80;
  waveLen=data->waveLen[i]+(uint32_t)(2*buffBins);
  temp=uchalloc((uint64_t)waveLen,"temp waveform",0);
  for(j=0;j<buffBins;j++)temp[j]=(unsigned char)dimage->meanN;
  for(j=0;j<(int)data->waveLen[i];j++)temp[j+buffBins]=wave[j];
  for(j=(int)data->waveLen[i]+buffBins;j<(int)waveLen;j++)temp[j]=(unsigned char)dimage->meanN;
  TIDY(wave);
  wave=temp;
  temp=NULL;

  /*deconvolve and reconvolve*/
  if(dimage->indDecon){
    processed=processWave(wave,(int)waveLen,dimage->decon,1.0);
    smooPro=smooth(dimage->pSigma,(int)waveLen,processed,dimage->res);
  }

  /*convolve with GEDI pulse*/
  floWave=falloc((uint64_t)waveLen,"",0);
  for(j=0;j<(int)waveLen;j++)floWave[j]=(float)wave[j]-dimage->meanN;
  smoothed=smooth(dimage->pSigma,(int)waveLen,floWave,dimage->res);
  TIDY(floWave);

  /*add up*/
  for(j=0;j<(int)waveLen;j++){
    binPosition(&x,&y,&z,j-buffBins,data->x[i],data->y[i],data->z[i],data->time[i],&(grad[0]));
    bin=(int)((waves->maxZ-z)/(double)dimage->res);
    if((bin>=0)&&(bin<waves->nBins)){
      waves->wave[2][bin]+=((float)wave[j]-dimage->meanN)*rScale;  /*with ALS pulse*/
      waves->wave[3][bin]+=smoothed[j]*rScale;
      if(dimage->doDecon){
        if(dimage->indDecon){
          waves->wave[4][bin]+=smooPro[j]*rScale;
          waves->wave[5][bin]+=processed[j]*rScale;
        }
        waves->wave[6][bin]+=((float)wave[j]-dimage->meanN)*rScale;
      }
    }
  }/*wave bin loop*/

  TIDY(wave);
  TIDY(smooPro);
  TIDY(smoothed);
  TIDY(processed);
  return;
}/*gediFromWaveform*/


/*####################################*/
/*set deconvolution parameters for GEDI*/

denPar *setDeconForGEDI(control *dimage)
{
  void setDenoiseDefault(denPar *);
  void readPulse(denPar *);
  denPar *decon=NULL;

  /*set defaults*/
  if(!(decon=(denPar *)calloc(1,sizeof(denPar)))){
    fprintf(stderr,"error decon structure allocation.\n");
    exit(1);
  }
  setDenoiseDefault(decon);

  /*particular values for here*/
  decon->deChang=pow(10.0,-8.0);  /*change between decon iterations to stop*/
  decon->thresh=17.0;
  decon->meanN=13.0;
  strcpy(decon->pNamen,"/Users/stevenhancock/data/bess/leica_shape/leicaPulse.dat");  /*pulse filename*/
  decon->deconMeth=0;     /*Gold's method*/
  decon->pScale=1.2;
  decon->noiseTrack=0;

  /*read system pulse*/
  readPulse(decon);

  return(decon);
}/*setDeconForGEDI*/


/*####################################*/
/*allocate wave structure*/

waveStruct *allocateGEDIwaves(control *dimage,dataStruct **data)
{
  int j=0,numb=0,k=0;
  uint64_t i=0;
  double maxZ=0,minZ=0;
  double buff=0;
  waveStruct *waves=NULL;
  char hasPoints=0;

  if(!(waves=(waveStruct *)calloc(1,sizeof(waveStruct)))){
    fprintf(stderr,"error waveStruct allocation.\n");
    exit(1);
  }

  /*determine wave bounds*/
  buff=35.0;
  minZ=100000000000.0;
  maxZ=-100000000000.0;
  hasPoints=0;
  for(numb=0;numb<dimage->nFiles;numb++){
    if(data[numb]->nPoints>0)hasPoints=1;
    for(i=0;i<data[numb]->nPoints;i++){
      if(data[numb]->z[i]>maxZ)maxZ=data[numb]->z[i];
      if(data[numb]->z[i]<minZ)minZ=data[numb]->z[i];
    }
  }/*bound finding*/

  if(hasPoints==0){
    fprintf(stderr,"No points included\n");
    exit(1);
  }

  waves->minZ=minZ-buff;
  waves->maxZ=maxZ+buff;
  /*fprintf(stdout,"Bounds are %f %f\n",waves->minZ,waves->maxZ);*/

  waves->nBins=(int)((waves->maxZ-waves->minZ)/(double)dimage->res);
  waves->nWaves=8;
  waves->wave=fFalloc(waves->nWaves,"result waveform",0);
  for(j=0;j<waves->nWaves;j++){
    waves->wave[j]=falloc((uint64_t)waves->nBins,"result waveform",j+1);
    for(k=0;k<waves->nBins;k++)waves->wave[j][k]=0.0;
  }
  if(dimage->ground){
    waves->canopy=fFalloc(2,"canopy",0);
    waves->ground=fFalloc(2,"ground",0);
    for(j=0;j<2;j++){
      waves->canopy[j]=falloc((uint64_t)waves->nBins,"canopy waveform",j+1);
      waves->ground[j]=falloc((uint64_t)waves->nBins,"ground waveform",j+1);
      for(k=0;k<waves->nBins;k++)waves->canopy[j][k]=waves->ground[j][k]=0.0;
    }
  }

  return(waves);
}/*allocateGEDIwaves*/


/*####################################*/
/*set GEDI footprint*/

void setGediFootprint(control *dimage)
{
  int i=0;
  float maxSep=0;
  float determineGaussSep(float,float);
  float totE=0;
  float az=0;
  double tX=0,tY=0;


  /*footprint width*/
  if(dimage->fSigma<0.0)dimage->fSigma=dimage->fWidth/2.0;
  az=dimage->lobeAng*M_PI/180.0;  /*convert anlge to radians*/

  /*number of lobes and allocate*/
  if(dimage->sideLobe==0)dimage->nLobes=1;
  else                   dimage->nLobes=7;
  if(!(dimage->lobe=(lobeStruct *)calloc(dimage->nLobes,sizeof(lobeStruct)))){
    fprintf(stderr,"error lobeStruct allocation.\n");
    exit(1);
  }

  /*central footprint*/
  i=0;
  dimage->lobe[i].coord[0]=dimage->coord[0];
  dimage->lobe[i].coord[1]=dimage->coord[1];
  dimage->lobe[i].fSigma=dimage->fSigma;
  maxSep=determineGaussSep(dimage->lobe[i].fSigma,dimage->iThresh*1.0);
  dimage->lobe[i].maxSepSq=(double)(maxSep*maxSep);
  if(dimage->sideLobe==0)dimage->lobe[i].E=1.0;
  else{  /*include side lobes*/
    totE=1.0+0.0599+0.0731+0.0317+0.0319+0.0167+0.0163;
    i=0;
    dimage->lobe[i].E=1.0/totE;

    /*first southern lobe*/
    i=1;
    dimage->lobe[i].E=0.0731/totE;
    dimage->lobe[i].coord[0]=dimage->coord[0]-20.0*sin(az);
    dimage->lobe[i].coord[1]=dimage->coord[1]-20.0*cos(az);
    dimage->lobe[i].fSigma=dimage->fSigma;
    maxSep=determineGaussSep(dimage->lobe[i].fSigma,dimage->iThresh*dimage->lobe[i].E);
    dimage->lobe[i].maxSepSq=(double)(maxSep*maxSep);

    /*first nothern lobe*/
    i=2;
    dimage->lobe[i].E=0.0599/totE;
    dimage->lobe[i].coord[0]=dimage->coord[0]+20.0*sin(az);
    dimage->lobe[i].coord[1]=dimage->coord[1]+20.0*cos(az);
    dimage->lobe[i].fSigma=dimage->fSigma;
    maxSep=determineGaussSep(dimage->lobe[i].fSigma,dimage->iThresh*dimage->lobe[i].E);
    dimage->lobe[i].maxSepSq=(double)(maxSep*maxSep);

    /*western lobe*/
    i=3;
    dimage->lobe[i].E=0.0319/totE;
    dimage->lobe[i].coord[0]=dimage->coord[0]-20.0*cos(az);
    dimage->lobe[i].coord[1]=dimage->coord[1]-20.0*sin(az);
    dimage->lobe[i].fSigma=dimage->fSigma;
    maxSep=determineGaussSep(dimage->lobe[i].fSigma,dimage->iThresh*dimage->lobe[i].E);
    dimage->lobe[i].maxSepSq=(double)(maxSep*maxSep);

    /*eastern lobe*/
    i=4;
    dimage->lobe[i].E=0.0317/totE;
    dimage->lobe[i].coord[0]=dimage->coord[0]+20.0*cos(az);
    dimage->lobe[i].coord[1]=dimage->coord[1]+20.0*sin(az);
    dimage->lobe[i].fSigma=dimage->fSigma;
    maxSep=determineGaussSep(dimage->lobe[i].fSigma,dimage->iThresh*dimage->lobe[i].E);
    dimage->lobe[i].maxSepSq=(double)(maxSep*maxSep);

    /*second southern lobe*/
    i=5;
    dimage->lobe[i].E=0.0167/totE;
    dimage->lobe[i].coord[0]=dimage->coord[0]-30.0*sin(az);
    dimage->lobe[i].coord[1]=dimage->coord[1]-30.0*cos(az);
    dimage->lobe[i].fSigma=dimage->fSigma;
    maxSep=determineGaussSep(dimage->lobe[i].fSigma,dimage->iThresh*dimage->lobe[i].E);
    dimage->lobe[i].maxSepSq=(double)(maxSep*maxSep);

    /*second northern lobe*/
    i=6;
    dimage->lobe[i].E=0.0163/totE;
    dimage->lobe[i].coord[0]=dimage->coord[0]+30.0*cos(az);
    dimage->lobe[i].coord[1]=dimage->coord[1]+30.0*sin(az);
    dimage->lobe[i].fSigma=dimage->fSigma;
    maxSep=determineGaussSep(dimage->lobe[i].fSigma,dimage->iThresh*dimage->lobe[i].E);
    dimage->lobe[i].maxSepSq=(double)(maxSep*maxSep);
  }/*side lobe test*/


  /*determine min and max bounds*/
  dimage->minX=dimage->minY=100000000000.0;
  dimage->maxX=dimage->maxY=-1000000000000.0;
  for(i=0;i<dimage->nLobes;i++){
    tX=dimage->lobe[i].coord[0]-dimage->lobe[i].maxSepSq;
    if(tX<dimage->minX)dimage->minX=tX;
    tX=dimage->lobe[i].coord[0]+dimage->lobe[i].maxSepSq;
    if(tX>dimage->maxX)dimage->maxX=tX;
    tY=dimage->lobe[i].coord[1]-dimage->lobe[i].maxSepSq;
    if(tY<dimage->minY)dimage->minY=tY;
    tY=dimage->lobe[i].coord[1]+dimage->lobe[i].maxSepSq;
    if(tY>dimage->maxY)dimage->maxY=tY;
  }

  /*grid for normalising sampling*/
  dimage->gridRes=1.5;
  dimage->gX=(int)((float)(dimage->maxX-dimage->minX)/dimage->gridRes)+2;
  dimage->gY=(int)((float)(dimage->maxY-dimage->minY)/dimage->gridRes)+2;
  dimage->g0[0]=dimage->minX+(double)dimage->gridRes;
  dimage->g0[1]=dimage->minY+(double)dimage->gridRes;
  dimage->nGrid=ialloc(dimage->gX*dimage->gY,"nGrid",0);
  for(i=dimage->gX*dimage->gY-1;i>=0;i--)dimage->nGrid[i]=0;

  return;
}/*setGediFootprint*/


/*####################################*/
/*maximum significant distance*/

float determineGaussSep(float fSigma,float thresh)
{
  float x=0,y=0;

  x=0.0;
  do{
    y=(float)gaussian((double)x,(double)fSigma,0.0);
    x+=0.2;
  }while(y>=thresh);

  return(x);
}/*determineGaussSep*/


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
  if(dimage->pSigma<0.0){  /*GEDI unless specificed*/
    fwhm=dimage->pFWHM*0.15;  /*time if for two way*/
    dimage->pSigma=fwhm/2.355;
  }

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
}/*readAsciiData*/


/*####################################*/
/*read command line*/

control *readCommands(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;

  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }

  dimage->nFiles=1;
  dimage->inList=chChalloc(dimage->nFiles,"inList",0);
  dimage->inList[0]=challoc(200,"inList",0);
  strcpy(&(dimage->inList[0][0]),"/Users/dill/data/teast/maryland_play/sc_79_112_1.las");
  strcpy(dimage->outNamen,"teast.wave");
  dimage->pFWHM=12.0;   /*12 ns FWHM*/
  dimage->fWidth=11.0;  /*22m diameter in e^-2*/
  dimage->res=0.15;
  dimage->pRes=0.01;
  dimage->coord[0]=624366.0;
  dimage->coord[1]=3.69810*pow(10.0,6.0);
  dimage->footRes=20.0;

  /*switches*/
  dimage->readWave=0;
  dimage->ground=0;
  dimage->sideLobe=0;   /*no side lobes*/
  dimage->lobeAng=0.0;
  dimage->pBuffSize=(uint64_t)200000000;
  dimage->checkCover=0;
  dimage->cleanOut=0;

  dimage->iThresh=0.0006;
  dimage->meanN=12.0;
  dimage->doDecon=0;
  dimage->indDecon=0;
  dimage->pSigma=-1.0;  /*leave blannk for default GEDI*/
  dimage->fSigma=-1.0;

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        TTIDY((void **)dimage->inList,dimage->nFiles);
        dimage->nFiles=1;
        dimage->inList=chChalloc(dimage->nFiles,"input name list",0);
        dimage->inList[0]=challoc((uint64_t)strlen(argv[++i])+1,"input name list",0);
        strcpy(dimage->inList[0],argv[i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(dimage->outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-inList",7)){
        checkArguments(1,i,argc,"-inList");
        TTIDY((void **)dimage->inList,dimage->nFiles);
        dimage->inList=readInList(&dimage->nFiles,argv[++i]);
      }else if(!strncasecmp(argv[i],"-coord",6)){
        checkArguments(2,i,argc,"-coord");
        dimage->coord[0]=atof(argv[++i]);
        dimage->coord[1]=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-decon",6)){
        dimage->doDecon=1;
      }else if(!strncasecmp(argv[i],"-indDecon",9)){
        dimage->indDecon=1;
        dimage->doDecon=1;
      }else if(!strncasecmp(argv[i],"-LVIS",5)){
        dimage->pSigma=0.6893/2.0;  /*two way trip*/
        dimage->fSigma=6.25;
      }else if(!strncasecmp(argv[i],"-pSigma",7)){
        checkArguments(1,i,argc,"-pSigma");
        dimage->pSigma=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-fSigma",7)){
        checkArguments(1,i,argc,"-fSigma");
        dimage->fSigma=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-readWave",9)){
        dimage->readWave=1;
      }else if(!strncasecmp(argv[i],"-ground",8)){
        dimage->ground=1;
        dimage->cleanOut=1;
      }else if(!strncasecmp(argv[i],"-sideLobe",9)){
        dimage->sideLobe=1;
      }else if(!strncasecmp(argv[i],"-lobeAng",8)){
        checkArguments(1,i,argc,"-lobeAng");
        dimage->lobeAng=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-pBuff",6)){
        checkArguments(1,i,argc,"-pBuff");
        dimage->pBuffSize=(uint64_t)(atof(argv[++i])*1000000000.0);
      }else if(!strncasecmp(argv[i],"-checkCover",11)){
        dimage->checkCover=1;
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n#####\nProgram to create GEDI waveforms from ALS las files\n#####\n\n-input name;     lasfile input filename\n-output name;    output filename\n-inList list;    input file list for multiple files\n-coord lon lat;  footprint coordinate in same system as lasfile\n-decon;          deconvolve\n-indDecon;       deconvolve individual beams\n-LVIS;           use LVIS pulse length, sigma=6.25m\n-pSigma sig;     set pulse width\n-fSigma sig;     set footprint width\n-readWave;       read full-waveform where available\n-ground;         split ground and canopy  points\n-sideLobe;       use side lobes\n-lobeAng ang;    lobe axis azimuth\n-pBuff s;        point reading buffer size in Gbytes\n-checkCover;     check that the footprint is covered by ALS data. Exit if not\n\nQuestions to svenhancock@gmail.com\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry gediRat -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }/*command parser*/


  return(dimage);
}/*readCommands*/

/*the end*/
/*####################################*/

