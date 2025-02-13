#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tiffRead.h"
#include "libDEMhandle.h"


/*#########################*/
/*# Functions to handle  #*/
/*# DEMs to process lidar #*/
/*#########################*/


/*#######################################*/
/*# Copyright 2006-2016, Steven Hancock #*/
/*# The program is distributed under    #*/
/*# the terms of the GNU General Public #*/
/*# License.    svenhancock@gmail.com   #*/
/*#######################################*/


/*########################################################################*/
/*# This file is part of libCLidar.                                      #*/
/*#                                                                      #*/
/*# libCLidar is free software: you can redistribute it and/or modify    #*/
/*# it under the terms of the GNU General Public License as published by #*/
/*# the Free Software Foundation, either version 3 of the License, or    #*/
/*#  (at your option) any later version.                                 #*/
/*#                                                                      #*/
/*# libCLidar is distributed in the hope that it will be useful,         #*/
/*# but WITHOUT ANY WARRANTY; without even the implied warranty of       #*/
/*#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #*/
/*#   GNU General Public License for more details.                       #*/
/*#                                                                      #*/
/*#    You should have received a copy of the GNU General Public License #*/
/*#    along with libClidar.  If not, see <http://www.gnu.org/licenses/>.#*/
/*########################################################################*/



/*########################################################################################*/
/*read an ASCII DEM*/

demStruct *readAscDEM(char *namen,double minX,double minY,double maxX,double maxY)
{
  int i=0,j=0,place=0;
  int maxLen=0,lineN=0;
  demStruct *dem=NULL;
  char *line=NULL,temp1[100],temp2[100];
  char *token=NULL;
  FILE *ipoo=NULL;

  maxLen=100000;
  line=challoc((uint64_t)maxLen,"line",0);

  if(!(dem=(demStruct *)calloc(1,sizeof(demStruct)))){
    fprintf(stderr,"error demStruct allocation.\n");
    exit(1);
  }
  if((ipoo=fopen(namen,"r"))==NULL){
    fprintf(stderr,"Error opening dem file %s\n",namen);
    exit(1);
  }

  lineN=0;
  while(fgets(line,maxLen,ipoo)!=NULL){
    if(lineN<6){ /*read the header*/
      if(sscanf(line,"%s %s",temp1,temp2)==2){
        if(!strncasecmp(line,"ncols",5))dem->nX=atoi(temp2);
        else if(!strncasecmp(line,"nrows",5))dem->nY=atoi(temp2);
        else if(!strncasecmp(line,"xllcorner",9))dem->minX=atof(temp2);
        else if(!strncasecmp(line,"yllcorner",9))dem->minY=atof(temp2);
        else if(!strncasecmp(line,"cellsize",8))dem->res=atof(temp2);
        else if(!strncasecmp(line,"NODATA_value",12)){
          dem->noData=atof(temp2);
          j=dem->nY-1;
          dem->z=dalloc(dem->nX*dem->nY,"dem",0);
        }
      }
    }else{  /*read data*/
      token=strtok(line," ");
      i=0;
      while(token!=NULL) {
        place=j*dem->nX+i;
        dem->z[place]=atof(token);
        i++;
      }
      j--;
    }
    lineN++;
  }

  /*allocate space*/


  TIDY(line);
  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return(dem);
}/*readAscDEM*/


/*########################################################################################*/
/*read a geotiff DEM*/

demStruct *readTifDEM(char *namen,double minX,double minY,double maxX,double maxY)
{
  int i=0,j=0,i0=0,j0=0,i1=0,j1=0;
  int demI=0,demJ=0;
  uint64_t demPlace=0,tifPlace=0;
  demStruct *dem=NULL;
  geot *geotiff=NULL;

  /*allocate space*/
  if(!(dem=(demStruct *)calloc(1,sizeof(demStruct)))){
    fprintf(stderr,"error demStruct allocation.\n");
    exit(1);
  }
  if(!(geotiff=(geot *)calloc(1,sizeof(geot)))){
    fprintf(stderr,"error geotiff allocation.\n");
    exit(1);
  }

  /*read the geotiff*/
  readGeotiff(geotiff,namen,1);

  /*copy header*/
  dem->minX=minX;
  dem->minY=minY;
  dem->maxX=maxX;
  dem->maxY=maxY;
  dem->res=geotiff->scale[0];
  dem->noData=-9999.0;
  dem->nX=(int)((dem->maxX-dem->minX)/dem->res+1);
  dem->nY=(int)((dem->maxY-dem->minY)/dem->res+1);
  dem->z=dalloc(dem->nX*dem->nY,"z array",0);

  /*extract area of interest*/
  i0=(int)((minX-geotiff->tiepoints[3])/dem->res);
  i1=(int)((maxX-geotiff->tiepoints[3])/dem->res);
  j0=(int)((geotiff->tiepoints[4]-minY)/dem->res);
  j1=(int)((geotiff->tiepoints[4]-maxY)/dem->res);
  if(i0<0)i0=0;
  if(i1>=geotiff->nX)i1=geotiff->nX-1;
  if(j0<0)j0=0;
  if(j1<0)j1=0;
  if(j0>=geotiff->nY)j0=geotiff->nY-1;
  if(j1>=geotiff->nY)j1=geotiff->nY-1;

  dem->maxZ=-1000000.0;
  dem->minZ=1000000.0;
  for(i=i0;i<=i1;i++){
    for(j=j0;j>=j1;j--){
      /*check we are inside the dem*/
      demI=i-i0;
      demJ=dem->nY-(j-j1);  /*we are flipping the DEM over*/
      if((demI<0)||(demI>=dem->nX)||(demJ<0)||(demJ>=dem->nY))continue;

      demPlace=(uint64_t)demI+(uint64_t)demJ*(uint64_t)dem->nX;
      tifPlace=(uint64_t)i+(uint64_t)j*(uint64_t)geotiff->nX;

      if(geotiff->fImage){
        if(geotiff->fImage[tifPlace]>0.0)dem->z[demPlace]=(double)geotiff->fImage[tifPlace];
        else                             dem->z[demPlace]=(double)dem->noData;
      }else if(geotiff->dImage){
        if(geotiff->dImage[tifPlace]>0.0)dem->z[demPlace]=geotiff->dImage[tifPlace];
        else                             dem->z[demPlace]=dem->noData;
      }else{
        if(geotiff->image[tifPlace]<255)dem->z[demPlace]=(double)geotiff->image[tifPlace]*geotiff->scale[2];
        else                            dem->z[demPlace]=dem->noData;
      }

      /*update bounds*/
      if(dem->z[demPlace]<dem->minZ)dem->minZ=dem->z[demPlace];
      if(dem->z[demPlace]>dem->maxZ)dem->maxZ=dem->z[demPlace];
    }/*y loop*/
  }/*x loop*/

  /*tidy up*/
  geotiff=tidyTiff(geotiff);

  return(dem);
}/*readTifDEM*/


/*##################################################*/
/*find an elevationnon the DEM*/

double findDEMelev(double x,double y,demStruct *dem)
{
  int dI=0,dJ=0,dPlace=0;

  dI=(int)((x-dem->minX)/(double)dem->res);
  dJ=(int)((y-dem->minY)/(double)dem->res);

  /*make sure we don't go outside the DEM*/
  if(dI<0)dI=0;
  if(dJ<0)dJ=0;
  if(dI>=dem->nX)dI=dem->nX-1;
  if(dJ>=dem->nY)dJ=dem->nY-1;
  dPlace=dI+(dem->nY-(dJ+1))*dem->nX;

  return(dem->z[dPlace]);
}/*findDEMelev*/


/*#################################################*/
/*tidy up a DEM structuire*/

demStruct *tidyDEMstruct(demStruct *dem)
{
  if(dem){
    TIDY(dem->z);
    TIDY(dem);
  }

  return(dem);
}/*tidyDEMstruct*/

/*the end*/
/*#######################################*/

