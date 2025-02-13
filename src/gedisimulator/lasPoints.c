#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "libLasRead.h"


/*########################*/
/*# Outputs point clouds #*/
/*# from las files       #*/
/*########################*/

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
/*control structure*/

typedef struct{
  char **inList;
  int nFiles;
  char outRoot[200];
  char canNamen[1000];  /*canopy filename*/
  char grNamen[1000];   /*ground filename*/
  float thresh;        /*maximum intensity of area of interest*/
  uint64_t pBuffSize;  /*point buffer rading size in bytes*/
  char writeAll;       /*write all points*/

  /*output files*/
  FILE *opoo;  /*non-ground*/
  FILE *gPoo;  /*ground points*/

  /*options*/
  char ground;   /*separate ground*/
  uint32_t decimate;  /*decimate the output cloud*/

  /*footprint parameters*/
  float fSigma;      /*footprint sigma*/
  double coord[2];   /*footprint centre*/
  double maxSepSq;
  double maxSep;
  char useBound;     /*output a rectangle*/
  double bound[4];   /*minX maxX minY maxY*/
}control;


/*####################################*/
/*main*/
  
int main(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;
  control *readCommands(int,char **);
  lasFile *las=NULL;
  void openOutput(control *);
  void readWritePoints(control *,lasFile *);
  void setArea(control *);


  /*read command line*/
  dimage=readCommands(argc,argv);

  /*determine maximum area of interest*/
  setArea(dimage);

  /*open output*/
  openOutput(dimage);

  /*loop over files*/
  for(i=0;i<dimage->nFiles;i++){
    /*read lasFile*/
    las=readLasHead(dimage->inList[i],dimage->pBuffSize);

    /*read and output data*/
    readWritePoints(dimage,las);

    /*tidy lasFIle*/
    las=tidyLasFile(las);
  }/*file loop*/


  fprintf(stdout,"Points to %s\n",dimage->canNamen);
  if(dimage->ground)fprintf(stdout,"Ground to %s\n",dimage->grNamen);


  /*tidy up*/
  if(dimage){
    if(dimage->opoo){
      fclose(dimage->opoo);
      dimage->opoo=NULL;
    }
    if(dimage->gPoo){
      fclose(dimage->gPoo);
      dimage->gPoo=NULL;
    }
    TTIDY((void **)dimage->inList,dimage->nFiles);
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*####################################*/
/*read and write points to file*/

void readWritePoints(control *dimage,lasFile *las)
{
  uint32_t i=0;
  double x=0,y=0,z=0;
  double dX=0,dY=0,sepSq=0;
  FILE *opoo=NULL;

  /*set bounds if not already done*/
  if(!dimage->useBound){
    dimage->bound[0]=dimage->coord[0]-dimage->maxSep;
    dimage->bound[1]=dimage->coord[0]+dimage->maxSep;
    dimage->bound[2]=dimage->coord[1]-dimage->maxSep;
    dimage->bound[3]=dimage->coord[1]+dimage->maxSep;
  }

  /*check file bounds*/
  if(checkFileBounds(las,dimage->bound[0],dimage->bound[1],dimage->bound[2],dimage->bound[3])||dimage->writeAll){
    for(i=0;i<las->nPoints;i+=dimage->decimate){
      /*read one point*/
      readLasPoint(las,i);
      setCoords(&x,&y,&z,las);

      dX=x-dimage->coord[0];
      dY=y-dimage->coord[1];
      sepSq=dX*dX+dY*dY;
      if((sepSq<=dimage->maxSepSq)||dimage->writeAll||(dimage->useBound&&((x>=dimage->bound[0])&&(x<=dimage->bound[1])&&(y>=dimage->bound[2])&&(y<=dimage->bound[3])))){
        if(dimage->ground&&(las->classif==2))opoo=dimage->gPoo;
        else                                 opoo=dimage->opoo;
        fprintf(opoo,"%f %f %f %d %d",x,y,z,(int)(las->refl),(int)las->scanAng);
        if((las->pointFormat==3)||(las->pointFormat==10)||(las->pointFormat==8)||(las->pointFormat==7)||(las->pointFormat==5)){
          fprintf(opoo," %.4f %.4f %.4f",(float)las->RGB[0]/pow(2.0,16),(float)las->RGB[1]/pow(2.0,16),(float)las->RGB[2]/pow(2.0,16));
        }   /*RGB scaled to be between 0 and 1*/
        fprintf(opoo,"\n");
        opoo=NULL;
      }/*separation check*/
    }/*point loop*/
  }/*file bounds check*/

  return;
}/*readWritePoints*/


/*####################################*/
/*set area of interest*/

void setArea(control *dimage)
{
  float x=0,y=0,res=0;

  /*only if we are not using a rectangular bound*/
  if(!dimage->useBound){
    if(dimage->fSigma>=0.0){  /*unless the radius has already been declared*/
      res=0.05;
      x=0.0;
      do{
        y=(float)gaussian((double)x,(double)dimage->fSigma,0.0);
        x+=res;
      }while(y>=dimage->thresh);
      dimage->maxSep=(double)x;
    }
    dimage->maxSepSq=dimage->maxSep*dimage->maxSep;
  }

  return;
}/*setArea*/


/*####################################*/
/*open output files*/

void openOutput(control *dimage)
{
  sprintf(dimage->canNamen,"%s.can.pts",dimage->outRoot);
  if((dimage->opoo=fopen(dimage->canNamen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",dimage->canNamen);
    exit(1);
  }
  fprintf(dimage->opoo,"# 1 x, 2 y, 3 z, 4 intensity, 5 scanAng, 6+ RGB if available\n");

  if(dimage->ground){
    sprintf(dimage->grNamen,"%s.ground.pts",dimage->outRoot);
    if((dimage->gPoo=fopen(dimage->grNamen,"w"))==NULL){
      fprintf(stderr,"Error opening output file %s\n",dimage->grNamen);
      exit(1);
    }
    fprintf(dimage->gPoo,"# 1 x, 2 y, 3 z, 4 intensity, 5 scanAng\n");
  }

  return;
}/*openOutput*/


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

  /*defaults*/
  dimage->nFiles=1;
  dimage->inList=chChalloc(dimage->nFiles,"inList",0);
  dimage->inList[0]=challoc(200,"inList",0);
  strcpy(&(dimage->inList[0][0]),"/Users/stevenhancock/data/gedi/laselva/ALS/tile197.las");
  strcpy(dimage->outRoot,"gTest");
  dimage->opoo=NULL;
  dimage->gPoo=NULL;
  dimage->writeAll=0;
  dimage->useBound=0;
  dimage->decimate=1;

  dimage->fSigma=5.5;  /*GEDI*/
  dimage->coord[0]=826367.0;
  dimage->coord[1]=1152271.0;
  dimage->ground=0;
  dimage->thresh=0.01;  /*1%*/

  dimage->pBuffSize=(uint64_t)200000000;


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
      }else if(!strncasecmp(argv[i],"-outRoot",8)){
        checkArguments(1,i,argc,"-outRoot");
        strcpy(dimage->outRoot,argv[++i]);
      }else if(!strncasecmp(argv[i],"-inList",7)){
        checkArguments(1,i,argc,"-inList");
        TTIDY((void **)dimage->inList,dimage->nFiles);
        dimage->inList=readInList(&dimage->nFiles,argv[++i]);
      }else if(!strncasecmp(argv[i],"-coord",6)){
        checkArguments(2,i,argc,"-coord");
        dimage->coord[0]=atof(argv[++i]);
        dimage->coord[1]=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-fSigma",7)){
        checkArguments(1,i,argc,"-fSigma");
        dimage->fSigma=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-rad",4)){
        checkArguments(1,i,argc,"-rad");
        dimage->fSigma=-1.0;
        dimage->maxSep=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-allPoints",10)){
        dimage->writeAll=1;
      }else if(!strncasecmp(argv[i],"-decimate",9)){
        checkArguments(1,i,argc,"-decimate");
        dimage->decimate=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-LVIS",5)){
        dimage->fSigma=6.25;  /*LVIS*/
      }else if(!strncasecmp(argv[i],"-ground",7)){
        dimage->ground=1;
      }else if(!strncasecmp(argv[i],"-pBuff",6)){
        checkArguments(1,i,argc,"-pBuff");
        dimage->pBuffSize=(uint64_t)(atof(argv[++i])*1000000000.0);
      }else if(!strncasecmp(argv[i],"-thresh",6)){
        checkArguments(2,i,argc,"-thresh");
        dimage->thresh=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-bounds",7)){
        checkArguments(4,i,argc,"-bounds");
        dimage->useBound=1;
        dimage->bound[0]=atof(argv[++i]);
        dimage->bound[1]=atof(argv[++i]);
        dimage->bound[2]=atof(argv[++i]);
        dimage->bound[3]=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n#####\nProgram to output ALS points within GEDI footprints\n#####\n\n-input name;     lasfile input filename\n-outRoot name;   output filename\n-inList list;    input file list for multiple files\n-coord lon lat;  footprint coordinate in same system as lasfile\n-fSigma sigma;   footprint width\n-rad rad;        radius to output\n-bounds minX maxX minY maxY;  output points within a rectangle\n-decimate n;     decimate output point cloud by a factor\n-allPoints;      write all points\n-LVIS;           use LVIS pulse length, sigma=6.25m\n-ground;         output canopy and ground separately\n-thresh t;       energy threshold to accept points\n-pBuff s;        point reading buffer size in Gbytes\n\nQuestions to svenhancock@gmail.com\n\n");
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

