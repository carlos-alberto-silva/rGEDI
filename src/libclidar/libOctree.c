#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "libOctree.h"

#if defined(_WIN32) && defined(_DLL)
#    ifndef DLL_EXPORT
#        define DLL_EXPORT __declspec(dllexport)
#    endif
// Windows or Linux static library, or Linux so
#else
#    ifndef DLL_EXPORT
#        define DLL_EXPORT
#    endif
#endif

/*########################*/
/*# Functions for octrees #*/
/*# in lidar programs     #*/
/*#########################*/

/*#######################################*/
/*# Copyright 2006-2017, Steven Hancock #*/
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


/*#######################################*/
/*allocate space for top level*/

octreeStruct *allocateOctree(int nLevels,int topN,double minX,double maxX,double minY,double maxY)
{
  int i=0;
  octreeStruct *octree=NULL;
  float dX=0,dY=0;

  /*allocate space*/
  if(!(octree=(octreeStruct *)calloc(1,sizeof(octreeStruct)))){
    fprintf(stderr,"error octree allocation.\n");
    exit(1);
  }

  /*set number of levels*/
  octree->nLevel=nLevels;

  /*copy bounds*/
  octree->minX=minX;
  octree->maxX=maxX;
  octree->minY=minY;
  octree->maxY=maxY;

  /*determine resolution*/
  dX=(float)(maxX-minX);
  dY=(float)(maxY-minY);
  octree->res=(dX>dY)?(float)(dX/(double)topN):(float)(dY/(double)topN);
  octree->nX=(int)(dX/octree->res)+1;
  octree->nY=(int)(dY/octree->res)+1;

  /*set pointers blank*/
  if(!(octree->tree=(treeStruct **)calloc(octree->nX*octree->nY,sizeof(treeStruct *)))){
    fprintf(stderr,"error octree allocation.\n");
    exit(1);
  }
  for(i=octree->nX*octree->nY-1;i>=0;i--)octree->tree[i]=NULL;
  octree->mapFile=NULL;
  octree->mapPoint=NULL;
  octree->nIn=NULL;
  octree->nMaps=0;

  return(octree);
}/*allocateOctree*/


/*#######################################*/
/*fill in octree*/

void fillOctree(double x,double y,double z,int nFile,uint32_t nPoint,octreeStruct *octree)
{
  int xBin=0,yBin=0;
  int place=0;
  double x0=0,y0=0;
  void mapOctree(int,octreeStruct *,treeStruct **,double,double,double,float,double,double,int,uint32_t);

  /*determine top level*/
  xBin=(int)((x-octree->minX)/(double)octree->res);
  yBin=(int)((y-octree->minY)/(double)octree->res);
  place=yBin*octree->nX+xBin;

  /*bounds check*/
  if((xBin>=0)&&(xBin<octree->nX)&&(yBin>=0)&&(yBin<octree->nY)){
    x0=(double)xBin*(double)octree->res+octree->minX;
    y0=(double)yBin*(double)octree->res+octree->minY;
    mapOctree(0,octree,&octree->tree[place],x,y,z,octree->res,x0,y0,nFile,nPoint);
  }/*bounds check*/

  return;
}/*fillOctree*/


/*#######################################*/
/*build octree map*/

void mapOctree(int level,octreeStruct *octree,treeStruct **tree,double x,double y,double z,float res,double x0,double y0,int nFile,uint32_t nPoint)
{
  int xBin=0,yBin=0,place=0,i=0,mapInd=0;
  void mapOctree(int,octreeStruct *,treeStruct **,double,double,double,float,double,double,int,uint32_t);
  uint32_t *markUint32(int,uint32_t *,uint32_t);
  int *markInt(int,int *,int);

  /*child coordinates*/
  xBin=(int)((x-x0)/(double)res);
  yBin=(int)((y-y0)/(double)res);
  place=xBin+yBin*2;

  /*if not allocated, allocate*/
  if(*tree==NULL){
    if(!(*tree=(treeStruct *)calloc(1,sizeof(treeStruct)))){
      fprintf(stderr,"error octree allocation.\n");
      exit(1);
    }
    if(level<(octree->nLevel-1)){  /*only allocate children for higher levels*/
      if(!(tree[0]->tree=(void **)calloc(4,sizeof(treeStruct *)))){
        fprintf(stderr,"error octree allocation.\n");
        exit(1);
      }
      for(i=0;i<4;i++)tree[0]->tree[i]=NULL;
    }else{  /*this is the lowest level, allocate map space*/
      if(octree->mapFile){
        if(!(octree->mapFile=(int **)realloc(octree->mapFile,(octree->nMaps+1)*sizeof(int *)))){
          fprintf(stderr,"Error in octree reallocation, %lu\n",(octree->nMaps+1)*sizeof(int *));
          exit(1);
        }
      }else octree->mapFile=iIalloc(octree->nMaps+1,"mapFile",0);
      if(octree->mapPoint){
        if(!(octree->mapPoint=(uint32_t **)realloc(octree->mapPoint,(octree->nMaps+1)*sizeof(uint32_t *)))){
          fprintf(stderr,"octree map allocation error\n");
          exit(1);
        }
      }else{
        if(!(octree->mapPoint=(uint32_t **)calloc(octree->nMaps+1,sizeof(uint32_t *)))){
          fprintf(stderr,"error octree map allocation.\n");
          exit(1);
        }
      }
      octree->mapPoint[octree->nMaps]=NULL;
      octree->mapFile[octree->nMaps]=NULL;
      octree->nIn=markUint32(octree->nMaps,octree->nIn,0);
      tree[0]->mapInd=octree->nMaps;
      octree->nMaps++;
    }
  }/*allocation*/

  /*which level are we in*/
  if(level<(octree->nLevel-1)){  /*keep recurssing*/
    res/=2.0;
    mapOctree(level+1,octree,(treeStruct **)(&tree[0]->tree[place]),x,y,z,res,x0+(double)xBin*(double)res,y0+(double)yBin*(double)res,nFile,nPoint);
  }else{  /*mark the points*/
    mapInd=tree[0]->mapInd;
    octree->mapPoint[mapInd]=markUint32(octree->nIn[mapInd],octree->mapPoint[mapInd],nPoint);
    octree->mapFile[mapInd]=markInt(octree->nIn[mapInd],octree->mapFile[mapInd],nFile);
    octree->nIn[mapInd]++;
  }

  return;
}/*mapOctree*/


/*####################################################*/
/*return list of files and point intersecting*/

pointMapStruct *mapFromOctree(int *octList,int nOct,octreeStruct *octree,double minX,double maxX,double minY,double maxY)
{
  int i=0,ind=0,xBin=0,yBin=0;
  double x0=0,y0=0;
  pointMapStruct *pointmap=NULL;
  void readOctree(treeStruct *,pointMapStruct *,int,octreeStruct *,double,double,float,double,double,double,double);

  /*allocate space*/
  if(!(pointmap=(pointMapStruct *)calloc(1,sizeof(pointMapStruct)))){
    fprintf(stderr,"error pointMapStruct allocation.\n");
    exit(1);
  }
  pointmap->nPoints=0;
  pointmap->fList=NULL;
  pointmap->pList=NULL;

  /*loop over top level*/
  for(i=0;i<nOct;i++){
    ind=octList[i];
    if(octree->tree[ind]){
      xBin=ind%octree->nX;
      yBin=ind/octree->nX;
      x0=(double)xBin*(double)octree->res+octree->minX;
      y0=(double)yBin*(double)octree->res+octree->minY;
      readOctree(octree->tree[ind],pointmap,0,octree,x0,y0,octree->res/2.0,minX,maxX,minY,maxY);
    }
  }/*top level loop*/

  return(pointmap);
}/*intersectOctree*/


/*#######################################*/
/*read the octree*/

void readOctree(treeStruct *tree,pointMapStruct *pointmap,int level,octreeStruct *octree,double x0,double y0,float res,double minX,double maxX,double minY,double maxY)
{
  int i=0,k=0,place=0;
  uint32_t j=0,mapInd=0,ind=0;
  double newX0=0,newY0=0,newX1=0,newY1=0;
  void readOctree(treeStruct *,pointMapStruct *,int,octreeStruct *,double,double,float,double,double,double,double);


  /*check if bottom level or not*/
  if(level<(octree->nLevel-1)){  /*recursively loop into octree*/
    /*loop over octree children*/
    for(i=0;i<2;i++){
      newX0=x0+(double)i*(double)res;
      for(k=0;k<2;k++){
        place=k*2+i;
        if(tree->tree[place]){
          newY0=y0+(double)k*(double)res;
          newX1=newX0+(double)res;
          newY1=newY0+(double)res;

          /*check bounds*/
          if((newX1<minX)||(newX0>maxX)||(newY0>maxY)||(newY1<minY))continue;
          readOctree((treeStruct *)tree->tree[place],pointmap,level+1,octree,newX0,newY0,res/2.0,minX,maxX,minY,maxY);
        }
      }
    }
  }else{  /*final level, count points*/
    mapInd=tree->mapInd;
    /*adjust array sizes*/
    if(pointmap->fList!=NULL){
      if(!(pointmap->fList=(int *)realloc(pointmap->fList,(pointmap->nPoints+octree->nIn[mapInd])*sizeof(int)))){
        fprintf(stderr,"Error allocating memory\n");
        exit(1);
      }
    }else pointmap->fList=ialloc(octree->nIn[mapInd],"fList",0);
    if(pointmap->pList!=NULL){
      if(!(pointmap->pList=(uint32_t *)realloc(pointmap->pList,(pointmap->nPoints+octree->nIn[mapInd])*sizeof(uint32_t)))){
        fprintf(stderr,"Error allocating memory\n");
        exit(1);
      }
    }else{
      if(!(pointmap->pList=(uint32_t *)calloc(octree->nIn[mapInd],sizeof(uint32_t)))){
        fprintf(stderr,"error tls allocation.\n");
        exit(1);
      }
    }
    /*copy data*/
    for(j=0;j<octree->nIn[mapInd];j++){
      ind=j+pointmap->nPoints;
      pointmap->fList[ind]=octree->mapFile[mapInd][j];
      pointmap->pList[ind]=octree->mapPoint[mapInd][j];
    }
    pointmap->nPoints+=octree->nIn[mapInd];
  }/*bottom level check*/

  return;
}/*readOctree*/


/*#######################################*/
/*deallocate the octree*/

DLL_EXPORT octreeStruct *tidyOctree(octreeStruct *octree)
{
  int i=0;

  if(octree){
    /*clean the map*/
    TTIDY((void **)octree->mapFile,octree->nMaps);
    TTIDY((void **)octree->mapPoint,octree->nMaps);       
    TIDY(octree->nIn);

    /*clean the octree*/
    for(i=octree->nX*octree->nY-1;i>=0;i--){
      TIDY(octree->tree[i]);
    }
    TIDY(octree);
  }

  return(octree);
}/*tidyOctree*/

/*the end*/
/*#######################################*/

