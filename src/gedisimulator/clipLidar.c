#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "hdf5.h"
#include "libLasRead.h"
#include "libLidVoxel.h"
#include "libLasProcess.h"
#include "libLidarHDF.h"
#include "libOctree.h"
#include "gediIO.h"


/*#########################*/
/*# Clips out areas of    #*/
/*# lidar daya    2018    #*/
/*# svenhancock@gmail.com #*/
/*#########################*/

/*#######################################*/
/*# Copyright 2015-2018, Steven Hancock #*/
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
  char **inList;      /*list of input filenames*/
  char outNamen[1000]; /*output filename*/
  int nFiles;         /*number of input files*/
  double bounds[4];   /*bounds to clip*/
  uint64_t pBuffSize;  /*point buffer rading size in bytes*/
  uint64_t oBuffSize;  /*output buffer size in bytes*/
}control;


/*####################################*/
/*main*/

int main(int argc,char **argv)
{
  control *dimage=NULL;
  control *readCommands(int,char **);
  void clipLidar(control *);

  /*read command line*/
  dimage=readCommands(argc,argv);

  /*read and write files*/
  clipLidar(dimage);

  /*tidy up*/
  if(dimage){
    TTIDY((void **)dimage->inList,dimage->nFiles);
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*####################################*/
/*clip out lidar data*/

void clipLidar(control *dimage)
{
  int i=0;
  uint32_t j=0;
  uint64_t memDone=0;  /*keep track of how much memory has been done*/
  double x=0,y=0,z=0;
  lasFile *newLas=NULL,*las=NULL;
  FILE *opoo=NULL;


  /*open output*/
  if((opoo=fopen(dimage->outNamen,"wb"))==NULL){
    fprintf(stderr,"Error opening output file \"%s\"\n",dimage->outNamen);
    exit(1);
  }

  /*allocate output space*/
  newLas=readLasHead(dimage->inList[0],dimage->oBuffSize);
  /*set nonesense bounds*/
  for(i=0;i<3;i++){
    newLas->minB[i]=100000000.0;
    newLas->maxB[i]=-100000000.0;
  }


  /*loop over files*/
  for(i=0;i<dimage->nFiles;i++){
    las=readLasHead(dimage->inList[i],dimage->pBuffSize);

    /*is this file needed?*/
    if(checkFileBounds(las,dimage->bounds[0],dimage->bounds[2],dimage->bounds[1],dimage->bounds[3])){
      /*loop over data and copy new data*/
      for(j=0;j<las->nPoints;j++){
        /*read one point*/
        readLasPoint(las,j);
        setCoords(&x,&y,&z,las);
        /*is it within bounds?*/
        if((x>=dimage->bounds[0])&&(x<=dimage->bounds[2])&&(y>=dimage->bounds[1])&&(y<=dimage->bounds[3])){

          /*copy to write array*/

          /*write out array and reset once we have it full*/

          /*keep track of new global bounds*/
          if(x<newLas->minB[0])newLas->minB[0]=x;
          if(x>newLas->maxB[0])newLas->maxB[0]=x;
          if(y<newLas->minB[1])newLas->minB[1]=y;
          if(y>newLas->maxB[1])newLas->maxB[1]=y;
          if(z<newLas->minB[2])newLas->minB[2]=z;
          if(z>newLas->maxB[2])newLas->maxB[2]=z;
        }/*point bounds check*/
      }/*point loop*/
    }/*file bounds check*/

    /* at last file, copy all variable header records*/
    if(i==(dimage->nFiles-1)){

    }
    /*tidy up*/

  }/*file loop*/


  /*tidy up output*/
  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Written to %s\n",dimage->outNamen);
  return;
}/*clipLidar*/


/*##############################################*/
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
  dimage->inList=NULL;
  dimage->nFiles=0;
  dimage->pBuffSize=(uint64_t)200000000;
  dimage->oBuffSize=(uint64_t)200000000;
  strcpy(dimage->outNamen,"teast.las");

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
      }else if(!strncasecmp(argv[i],"-inList",7)){
        checkArguments(1,i,argc,"-inList");
        TTIDY((void **)dimage->inList,dimage->nFiles);
        dimage->inList=readInList(&dimage->nFiles,argv[++i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(dimage->outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-pBuff",6)){
        checkArguments(1,i,argc,"-pBuff");
        dimage->pBuffSize=(uint64_t)(atof(argv[++i])*1000000000.0)/2;
        dimage->oBuffSize=(uint64_t)(atof(argv[++i])*1000000000.0)/2;
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n#####\nProgram to clip out section of lidar data\n#####\n\n-input name;     lasfile input filename\n-output name;    output filename\n-inList list;    input file list for multiple files\n-bounds minX minY maxX maxX;  bounds of area to clip out\n-pBuff s;        point reading buffer size in Gbytes\n\n");
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

