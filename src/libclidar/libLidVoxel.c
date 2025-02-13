#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "hdf5.h"
#include "tools.h"
#include "libLasRead.h"
#include "libDEMhandle.h"
#include "libLasProcess.h"
#include "libLidVoxel.h"
#include "gsl/gsl_sort.h"


/*#########################*/
/*# Functions for dealing #*/
/*# with voxels           #*/
/*# S Hancock             #*/
/*# 6th November 2014     #*/
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



#define TOLERANCE 0.000001   /*tolerance for intersection tests*/
#define VTOL 0.01            /*tolerance for voxel finding*/

/*global to pass between functions*/
double tanZen=0,cosZen=0,sinZen=0;  /*to save calculations*/
double tanAz=0,cosAz=0,sinAz=0;


/*#############################################*/
/*make silhouette image from point cloud*/

void silhouetteImage(int nFiles,pCloudStruct **alsData,tlsScan *tlsData,rImageStruct *rImage,lidVoxPar *lidPar,int *voxList,int nIn,tlsVoxMap *map)
{
  int i=0,k=0,bin=0;
  int vInd=0,pInd=0,fInd=0;
  uint32_t j=0;
  float zen=0,az=0;
  double vect[3];
  void markPointSilhouette(double *,rImageStruct *,int,lidVoxPar *,float,uint16_t,double);

  /*angles for rotation*/
  zen=(float)atan2(sqrt((double)rImage->grad[0]*(double)rImage->grad[0]+(double)rImage->grad[1]*(double)rImage->grad[1]),(double)rImage->grad[2]);
  az=(float)atan2((double)rImage->grad[0],(double)rImage->grad[1]);


  if(alsData){   /*use ALS data*/
    for(i=0;i<nFiles;i++){
      for(j=0;j<alsData[i]->nPoints;j++){
        vect[0]=alsData[i]->x[j]-rImage->x0;
        vect[1]=alsData[i]->y[j]-rImage->y0;
        vect[2]=alsData[i]->z[j]-rImage->z0;
        /*rotate to x-y plane*/
        rotateZ(vect,(double)(-1.0*az));
        rotateX(vect,(double)(-1.0*zen));
        bin=(int)(vect[2]/rImage->rRes+0.5);

        if((bin>=0)&&(bin<rImage->nBins)){
          /*black out the points*/
          markPointSilhouette(&(vect[0]),rImage,bin,lidPar,alsData[i]->gap[j],alsData[i]->refl[j],0.0);
        }
      }/*point loop*/
    }/*file loop*/
  }else if(tlsData){   /*use TLS data*/
    for(k=0;k<nIn;k++){
      vInd=voxList[k];
      for(i=0;i<map->nIn[vInd];i++){
        fInd=map->mapFile[vInd][i];
        pInd=map->mapPoint[vInd][i];

        vect[0]=(double)tlsData[fInd].point[pInd].x+tlsData[fInd].xOff-rImage->x0;
        vect[1]=(double)tlsData[fInd].point[pInd].y+tlsData[fInd].yOff-rImage->y0;
        vect[2]=(double)tlsData[fInd].point[pInd].z+tlsData[fInd].zOff-rImage->z0;

        /*rotate to x-y plane*/
        rotateZ(vect,(double)(-1.0*az));
        rotateX(vect,(double)(-1.0*zen));
        bin=(int)(vect[2]/rImage->rRes+0.5);

        if((bin>=0)&&(bin<rImage->nBins)){
          /*black out the points*/
          markPointSilhouette(&(vect[0]),rImage,bin,lidPar,tlsData[fInd].point[pInd].gap,tlsData[fInd].point[pInd].refl,tlsData[fInd].point[pInd].r);
        } 
      }/*point in voxel loop*/
    }/*voxel loop*/
  }else{  /*no data. Something is wrong*/
    fprintf(stderr,"No data provided\n");
    exit(1);
  }

  return;
}/*silhouetteImage*/


/*############################################*/
/*mark lidar point in range image*/

void markPointSilhouette(double *coord,rImageStruct *rImage,int bin,lidVoxPar *lidPar,float gap,uint16_t refl,double r)
{
  int xInd=0,yInd=0,rPlace=0;
  int xStart=0,xEnd=0;
  int yStart=0,yEnd=0;
  int xIcent=0,yIcent=0;
  float rad=0;
  float maxRsepSq=0,rSepSq=0;

  if(gap<lidPar->minGap)gap=lidPar->minGap;
  rad=tlsPointSize(r,refl,lidPar->beamTanDiv,lidPar->beamRad,lidPar->minRefl,lidPar->maxRefl,lidPar->appRefl,gap); //)*lidPar->appRefl/gap;

  /*range image*/
  xIcent=(int)((coord[0]/(double)rImage->iRes)+0.5*(double)rImage->nX);
  yIcent=(int)((coord[1]/(double)rImage->iRes)+0.5*(double)rImage->nY);
  xStart=xIcent-(int)(rad/rImage->iRes+0.5);
  xEnd=xIcent+(int)(rad/rImage->iRes+0.5);
  yStart=yIcent-(int)(rad/rImage->iRes+0.5);
  yEnd=yIcent+(int)(rad/rImage->iRes+0.5);
  maxRsepSq=rad*rad;

  if(xStart<0)xStart=-1;      /*enforce bounds*/
  else if(xStart>=rImage->nX)xStart=rImage->nX;
  if(xEnd<0)xEnd=-1;      /*enforce bounds*/
  else if(xEnd>=rImage->nX)xEnd=rImage->nX;
  if(yStart<0)yStart=-1;
  else if(yStart>=rImage->nY)yStart=rImage->nY;
  if(yEnd<0)yEnd=-1;
  else if(yEnd>=rImage->nY)yEnd=rImage->nY;   /*enforce bounds*/

  /*mark centre point*/
  if((xIcent>=0)&&(xIcent<rImage->nX)&&(yIcent>=0)&&(yIcent<rImage->nY)){
    rPlace=yIcent*rImage->nX+xIcent;
    rImage->image[bin][rPlace]=1;
  }

  /*mark other points*/
  for(xInd=xStart;xInd<=xEnd;xInd++){
    if((xInd<0)||(xInd>=rImage->nX))continue;
    for(yInd=yStart;yInd<=yEnd;yInd++){
      if((yInd<0)||(yInd>=rImage->nY))continue;
      rSepSq=(float)((xInd-xIcent)*(xInd-xIcent)+(yInd-yIcent)*(yInd-yIcent))*rImage->iRes*rImage->iRes;
      if(rSepSq<=maxRsepSq){
        rPlace=yInd*rImage->nX+xInd;
        rImage->image[bin][rPlace]=1;
      }/*check within point*/
    }/*loop around point*/
  }/*loop around point*/

  return;
}/*markPointSilhouette*/


/*############################################*/
/*determine hit size*/

float tlsPointSize(double range,uint16_t refl,float tanDiv,float beamRad,float min,float max,float rhoApp,float gap)
{
  float d=0;
  float appRefl=0;
  float reflScale=0;

  appRefl=(float)max-(float)min;   /*scale from DN to size, takes phase func and albedo into account*/

  d=range*tanDiv+beamRad;                /*beam radius*/
  reflScale=((float)refl-(float)min)/appRefl;
  if(reflScale<0.0)     reflScale=0.0;   /*keep to bounds*/
  else if(reflScale>1.0)reflScale=1.0;   /*keep to bounds*/
  if(rhoApp<0.0)rhoApp=0.0;
  if(gap>TOLERANCE)d*=sqrt(reflScale*rhoApp/gap);   /*take optics into account*/
  else             d*=sqrt(reflScale*rhoApp/TOLERANCE);
  return(d);
}/*tlsPointSize*/


/*#############################################*/
/*allocate structure for range image*/

rImageStruct *allocateRangeImage(float beamRad,float rRes,float iRes,float *grad,double *origin,double *bounds)
{
  int i=0,k=0;
  float zen=0;
  rImageStruct *rImage=NULL;

  if(!(rImage=(rImageStruct *)calloc(1,sizeof(rImageStruct)))){
    fprintf(stderr,"error range image structure allocation.\n");
    exit(1);
  }

  rImage->x0=origin[0];
  rImage->y0=origin[1];
  rImage->z0=origin[2];
  if(grad){
    if(fabs(grad[0]+grad[1]+grad[2])>TOLERANCE){
      for(i=0;i<3;i++)rImage->grad[i]=grad[i];
    }else{
      rImage->grad[0]=rImage->grad[1]=0.0;
      rImage->grad[2]=-1.0;
    }
  }else{
    rImage->grad[0]=rImage->grad[1]=0.0;
    rImage->grad[2]=-1.0;
  }

  /*angles for rotation*/
  zen=(float)atan2(sqrt((double)rImage->grad[0]*(double)rImage->grad[0]+\
       (double)rImage->grad[1]*(double)rImage->grad[1]),(double)rImage->grad[2]);

  rImage->bounds[0]=-1.0*(double)beamRad;
  rImage->bounds[1]=-1.0*(double)beamRad;
  rImage->bounds[2]=0.0;
  rImage->bounds[3]=(double)beamRad;
  rImage->bounds[4]=(double)beamRad;
  rImage->bounds[5]=fabs(bounds[5]-bounds[2])*-1.0*(double)cos(zen);


  rImage->rRes=rRes;
  rImage->iRes=iRes;
  rImage->nBins=(int)((rImage->bounds[5]-rImage->bounds[2])/rImage->rRes+0.5);
  rImage->nX=(int)((rImage->bounds[3]-rImage->bounds[0])/rImage->iRes+0.5);
  rImage->nY=(int)((rImage->bounds[4]-rImage->bounds[1])/rImage->iRes+0.5);

  /*allocate image and set blank*/
  rImage->image=chChalloc(rImage->nBins,"range image",0);
  for(i=0;i<rImage->nBins;i++){
    rImage->image[i]=challoc((uint64_t)rImage->nX*(uint64_t)rImage->nY,"range image",i+1);
    for(k=rImage->nX*rImage->nY-1;k>=0;k--)rImage->image[i][k]=0;
  }
  return(rImage);
}/*allocateRangeImage*/


/*#############################################*/
/*add up hits and misses for a single beam*/

void countVoxGap(double x,double y,double z,float *grad,voxStruct *vox,int retNumb,int nRet,float beamRad,int numb)
{
  int i=0;
  int *voxList=NULL,nTot=0;
  double *rangeList=NULL;

  /*only do this for last returns per beam*/
  if(retNumb<nRet)return;

  /*check that a vector is there*/
  if(grad){
    if(fabs(grad[0]+grad[1]+grad[2])>TOLERANCE){
      /*determine which voxels are intersected*/
      voxList=beamVoxels(&(grad[0]),x,y,z,&(vox->bounds[0]),&(vox->res[0]),vox->nX,vox->nY,vox->nZ,&nTot,beamRad,&rangeList,-1.0);

      /*loop along intersected voxels*/
      for(i=0;i<nTot;i++){
        if(rangeList[i]<=0.0)vox->hits[numb][voxList[i]]+=1.0;
        else                 vox->miss[numb][voxList[i]]+=1.0;
      }/*intersecting voxel loop*/
      TIDY(rangeList);
      TIDY(voxList);
    }
  }

  return;
}/*countVoxGap*/


/*#######################################*/
/*voxels intersecting beam with width*/

int *beamVoxels(float *gradIn,double x0,double y0,double z0,double *bounds,double *res,int nX,int nY,int nZ,int *nPix,double beamRad,double **rangeList,float vRes)
{
  int i=0,j=0,k=0;
  int tempPix=0;
  int *pixList=NULL;
  int *tempList=NULL;
  int *findVoxels(double *,double,double,double,double *,double *,int *,int,int,int,double **);
  int *markInt(int,int *,int);
  double *markDo(int,double *,double);
  double *tempRange=NULL;
  double grad[3];
  float ang=0,angStep=0;  /*angular steps around edge of beam*/
  float rad=0,radRes=0;   /*radius to step along radial lines*/
  double x=0,y=0,z=0;
  char foundNew=0;

  /*determine angular resolution*/
  if(vRes>0.0)angStep=atan2(vRes/3.0,beamRad);
  else        angStep=2.0*M_PI/90.0;

  /*determine radial resolution*/
  if(vRes>=beamRad)radRes=beamRad;
  else             radRes=vRes/2.0;

  /*central beam*/
  for(i=0;i<3;i++)grad[i]=(double)gradIn[i];
  pixList=findVoxels(&(grad[0]),x0,y0,z0,bounds,res,nPix,nX,nY,nZ,rangeList);

  /*loop around rim of the beam*/
  ang=0.0;
  while(ang<2.0*M_PI){  /*angular loop*/
    rad=0.0;
    while(rad<=(float)beamRad){   /*radial loop*/
      x=rad*sin(ang)+x0;  /*new start along edge of beam*/
      y=rad*cos(ang)+y0;  /*new start along edge of beam*/
      z=z0;     /*this should take into account the zenith angle of the beam*/

      /*find voxels intersected by the beam along that edge*/
      for(j=0;j<3;j++)grad[j]=(double)gradIn[j];
      tempList=findVoxels(&(grad[0]),x,y,z,bounds,res,&tempPix,nX,nY,nZ,&tempRange);
      /*now sort through*/
      for(j=0;j<tempPix;j++){
        foundNew=1;
        for(k=0;k<(*nPix);k++){
          if(pixList[k]==tempList[j]){
            foundNew=0;
            break;
          }
        }/*final list loop*/
        if(foundNew==1){  /*if new, mark it*/
          pixList=markInt(*nPix,pixList,tempList[j]);
          rangeList[0]=markDo(*nPix,rangeList[0],tempRange[j]);
          (*nPix)++;
        } /*if new, mark it*/
      }/*temporary list loop*/
      TIDY(tempList);
      TIDY(tempRange);

      rad+=radRes;
    }/*radial loop*/
    ang+=angStep;
  }/*sub step loop*/

  return(pixList);
}/*beamVoxels*/


/*#######################################*/
/*if outside bounds find closest facet*/

void findClosestFacet(double *coords,double *bounds,double *vect,int xDir,int yDir,int zDir)
{
  int i=0;
  char found=0;
  double xR=0,yR=0,zR=0;
  double minR=0;
  double transCoord[3];
  void checkOuterFacet(double,double *,double *,double *,double *,double *,char *);

  if(xDir<0)xR=(bounds[3]-coords[0])/vect[0];  /*check right side*/
  else      xR=(bounds[0]-coords[0])/vect[0];  /*check left side*/
  if(yDir<0)yR=(bounds[4]-coords[1])/vect[1];  /*check front side*/
  else      yR=(bounds[1]-coords[1])/vect[1];  /*check back side*/
  if(zDir<0)zR=(bounds[5]-coords[2])/vect[2];  /*check top side*/
  else      zR=(bounds[2]-coords[2])/vect[2];  /*check bottom side*/
  minR=1000000000000.0;

  /*check all the possible sides*/
  checkOuterFacet(zR,vect,coords,bounds,&minR,transCoord,&found);   /*top/bottom*/
  checkOuterFacet(xR,vect,coords,bounds,&minR,transCoord,&found);   /*left/right*/
  checkOuterFacet(yR,vect,coords,bounds,&minR,transCoord,&found);   /*front/back*/

  /*copy if found*/
  if(found){
    for(i=0;i<3;i++)coords[i]=transCoord[i];
  }

  return;
}/*findClosestFacet*/


/*#######################################*/
/*check outer facet*/

void checkOuterFacet(double range,double *vect,double *coords,double *bounds,double *minR,double *transCoord,char *found)
{
  double x=0,y=0,z=0;

  x=range*vect[0]+coords[0];
  y=range*vect[1]+coords[1];
  z=range*vect[2]+coords[2];
  if((x<=bounds[3]+TOLERANCE)&&(x>=bounds[0]-TOLERANCE)&&(y<=bounds[4]+TOLERANCE)&&\
     (y>=bounds[1]+TOLERANCE)&&(z<=bounds[5]+TOLERANCE)&&(z>=bounds[2]+TOLERANCE)){
    if(range<*minR){
      *minR=range;
      transCoord[0]=x;
      transCoord[1]=y;
      transCoord[2]=z;
      *found=1;
    }
  }
  return;
}/*checkOuterFacet*/


/*###############################################*/
/*make ground return slice solid*/

void fillInRimageGround(rImageStruct *rImage)
{
  int i=0,j=0,bin=0;
  char brEak=0;

  /*find lowest bin*/
  brEak=0;
  for(i=rImage->nBins-1;i>=0;i--){
    for(j=rImage->nX*rImage->nY-1;j>=0;j--){
      if(rImage->image[i][j]>0){
        brEak=1;
        bin=i;
        break;
      }
    }
    if(brEak)break;
  }
  if(brEak){
    for(j=rImage->nX*rImage->nY-1;j>=0;j--)rImage->image[bin][j]=1;
  }

  return;
}/*fillInRimageGround*/


/*###############################################*/
/*make a waveform from a point cloud image*/

void waveFromImage(rImageStruct *rImage,float **wave,char gaussFoot,float fSigma)
{
  int i=0,j=0,k=0;
  int place=0;
  float dx=0,dy=0;
  float weight=0,sep=0;
  float total=0,totWeight=0;
  char doneIt=0;

  for(i=0;i<rImage->nBins;i++){
    wave[0][i]=wave[1][i]=0.0;
  }
  if(gaussFoot==1){
    for(i=0;i<rImage->nX;i++){
      dx=(float)(i-rImage->nX/2)*rImage->iRes;
      for(j=0;j<rImage->nY;j++){
        dy=(float)(j-rImage->nY/2)*rImage->iRes;
        sep=sqrt(dx*dx+dy*dy);
        totWeight+=gaussian((double)sep,(double)fSigma,0.0);
      }
    }
  }else totWeight=(float)(rImage->nX*rImage->nY);

  /*turn range images into a waveforms*/
  for(i=0;i<rImage->nX;i++){
    dx=(float)(i-rImage->nX/2)*rImage->iRes;
    for(j=0;j<rImage->nY;j++){
      place=j*rImage->nX+i;
      dy=(float)(j-rImage->nY/2)*rImage->iRes;
      sep=sqrt(dx*dx+dy*dy);
      doneIt=0;

      for(k=0;k<rImage->nBins;k++){ /*image bin loop*/
        if(rImage->image[k][place]>0){
          if(gaussFoot==0){
            if(sep<=fSigma)weight=1.0;
            else           weight=0.0;
          }else if(gaussFoot==-1){
            weight=1.0;
          }else if(gaussFoot==1){
            weight=gaussian((double)sep,(double)fSigma,0.0);
          }

          if(doneIt==0){
            wave[0][k]+=weight; /*appRefl*(rimRes*rimRes)/(M_PI*beamRad*beamRad);*/
            total+=weight;
            doneIt=1;
          }/*first hit only*/
          wave[1][k]+=weight; /*appRefl*(rimRes*rimRes)/(M_PI*beamRad*beamRad);*/  /*hits per bin*/
        }
      }/*y loop*/
    }/*x loop*/
  }/*range image bin loop*/

  /*normalise waveforms*/
  total=0.0;
  for(j=0;j<rImage->nBins;j++)total+=wave[0][j];
  if(total>0.0){
    for(j=0;j<rImage->nBins;j++){
      wave[0][j]/=total;
      wave[1][j]/=totWeight;
    }
  }

  return;
}/*waveFromImage*/


/*############################################*/
/*allocate voxel structure*/

voxStruct *voxAllocate(int nFiles,float *vRes,double *bounds,char useRMSE)
{
  int i=0,j=0;
  voxStruct *vox=NULL;

  /*allocate sctructure*/
  if(!(vox=(voxStruct *)calloc(1,sizeof(voxStruct)))){
    fprintf(stderr,"error voxel structure allocation.\n");
    exit(1);
  }

  /*leave derived values as NULL, as not always needed*/
  vox->gap=NULL;
  vox->gapTo=NULL;
  vox->PAIb=NULL;
  vox->PAIrad=NULL;

  /*note that findVoxels() needs minX maxX etc, different to dimage's minX minY etc*/
  for(i=0;i<3;i++)vox->res[i]=(double)vRes[i];
  for(i=0;i<6;i++)vox->bounds[i]=bounds[i];
  vox->nX=(int)((vox->bounds[3]-vox->bounds[0])/vox->res[0]+0.99);  /*add 0.99 to avoid rounding*/
  vox->nY=(int)((vox->bounds[4]-vox->bounds[1])/vox->res[1]+0.99);  /*add 0.99 to avoid rounding*/
  vox->nZ=(int)((vox->bounds[5]-vox->bounds[2])/vox->res[2]+0.99);  /*add 0.99 to avoid rounding*/
  vox->volume=(float)(vox->res[0]*vox->res[1]*vox->res[2]);

  vox->savePts=1;   /*defaults*/
  vox->maxZen=1000000.0;  /*use all points*/

  /*check for memory wrapping*/
  if(((uint64_t)vox->nX*(uint64_t)vox->nY*(uint64_t)vox->nZ)>=2147483647){
    fprintf(stderr,"Voxel bounds are too big to handle. Reduce %d %d %d\n",vox->nX,vox->nY,vox->nZ);
    exit(1);
  }

  vox->nVox=vox->nX*vox->nY*vox->nZ;
  vox->nScans=nFiles;

  vox->hits=fFalloc(vox->nScans,"voxel beam hits",0);
  vox->miss=fFalloc(vox->nScans,"voxel beam miss",0);
  vox->inHit=fFalloc(vox->nScans,"voxel point hits",0);
  vox->inMiss=fFalloc(vox->nScans,"voxel point miss",0);
  vox->sampVol=fFalloc(vox->nScans,"voxel volume sampled",0);
  vox->totVol=fFalloc(vox->nScans,"voxel volume total",0);
  vox->sumRsq=fFalloc(vox->nScans,"sum of radius of TLS points, squared",0);
  vox->meanRefl=fFalloc(vox->nScans,"mean reflectance of intersecting beams",0);
  vox->meanZen=fFalloc(vox->nScans,"mean zenith angle of intersecting beams",0);
  vox->meanRange=fFalloc(vox->nScans,"mean range of intersecting beams",0);
  vox->contN=ialloc(vox->nVox,"voxel contribution",0);
  for(j=0;j<vox->nVox;j++)vox->contN[j]=0;

  vox->useRMSE=useRMSE;
  if(useRMSE){
    vox->rmse=falloc((uint64_t)vox->nVox,"voxel error",0);
    for(j=0;j<vox->nVox;j++)vox->rmse[j]=0.0;
  }

  for(i=0;i<vox->nScans;i++){
    vox->hits[i]=falloc((uint64_t)vox->nVox,"voxel hits",i+1);
    vox->miss[i]=falloc((uint64_t)vox->nVox,"voxel miss",i+1);
    vox->inHit[i]=falloc((uint64_t)vox->nVox,"voxel hits",i+1);
    vox->inMiss[i]=falloc((uint64_t)vox->nVox,"voxel miss",i+1);
    vox->sampVol[i]=falloc((uint64_t)vox->nVox,"voxel volume sampled",i+1);
    vox->totVol[i]=falloc((uint64_t)vox->nVox,"voxel volume total",i+1);
    vox->meanRefl[i]=falloc((uint64_t)vox->nVox,"mean reflectance of intersecting beams",i+1);
    vox->meanZen[i]=falloc((uint64_t)vox->nVox,"mean zenith angle of intersecting beams",i+1);
    vox->meanRange[i]=falloc((uint64_t)vox->nVox,"mean range of intersecting beams",i+1);
    vox->sumRsq[i]=falloc((uint64_t)vox->nVox,"sum of radius of TLS points, squared",i+1);
    for(j=0;j<vox->nVox;j++){
      vox->hits[i][j]=vox->miss[i][j]=vox->inHit[i][j]=vox->inMiss[i][j]=0.0;
      vox->sampVol[i][j]=vox->totVol[i][j]=vox->sumRsq[i][j]=0.0;
      vox->meanRefl[i][j]=vox->meanZen[i][j]=vox->meanRange[i][j]=0.0;
    }
  }/*file loop*/

  /*leave the DTM out for now*/
  vox->dem=NULL;

  return(vox);
}/*voxAllocate*/


/*#######################################*/
/*tidy voxel structure*/

voxStruct *tidyVox(voxStruct *vox)
{

  if(vox){
    TTIDY((void **)vox->hits,vox->nScans);
    TTIDY((void **)vox->miss,vox->nScans);
    TTIDY((void **)vox->inHit,vox->nScans);
    TTIDY((void **)vox->inMiss,vox->nScans);
    TTIDY((void **)vox->sampVol,vox->nScans);
    TTIDY((void **)vox->totVol,vox->nScans);
    TTIDY((void **)vox->meanRefl,vox->nScans);
    TTIDY((void **)vox->meanZen,vox->nScans);
    TTIDY((void **)vox->meanRange,vox->nScans);
    TTIDY((void **)vox->sumRsq,vox->nScans);
    TTIDY((void **)vox->gap,vox->nScans);
    TTIDY((void **)vox->gapTo,vox->nScans);
    TTIDY((void **)vox->PAIb,vox->nScans);
    TTIDY((void **)vox->PAIrad,vox->nScans);
    TIDY(vox->rmse);
    TIDY(vox->contN);
    TIDY(vox);
  }

  return(NULL);
}/*tidyVox*/


/*###########################################################################*/
/*tidy up voxel map*/

void tidyVoxelMap(tlsVoxMap *map,int nVox)
{
  TTIDY((void **)map->mapFile,nVox);
  TTIDY((void **)map->mapPoint,nVox);
  TIDY(map->nIn);

  return;
}/*tidyVoxelMap*/


/*###########################################################################*/
/*Find all voxels in a list of voxels and sort*/

int *findAllVoxels(double *vect,double xCent,double yCent,double zCent,voxStruct **vox,int nFiles,int **fileList,double **rangeList,int *nIn)
{
  int i=0,j=0;
  int *voxList=NULL;
  int *tempList=NULL;
  int tempInt=0;
  int *tempVL=NULL,*tempFL=NULL;  /*temp voxel and file indices*/
  double *tempRanges=NULL;
  double *tempR=NULL;

  /*reset counter*/
  *nIn=0;
  TIDY(*fileList);
  TIDY(*rangeList);

  /*loop over files*/
  for(i=0;i<nFiles;i++){
    /*find voxels for this file*/
    tempList=findVoxels(vect,xCent,yCent,zCent,&vox[i]->bounds[0],&vox[i]->res[0],&tempInt,vox[i]->nX,vox[i]->nY,vox[i]->nZ,&tempRanges);

    /*have we found any?*/
    if(tempInt==0){
      TIDY(tempList);
      TIDY(tempRanges);
      continue;
    }

    /*add space to main array*/
    if(*nIn>0){
      if(tempRanges[0]<*rangeList[0]){ /*insert before*/

        tempVL=ialloc(*nIn+tempInt,"",0);
        for(j=0;j<tempInt;j++)tempVL[j]=tempList[j];
        for(j=0;j<*nIn;j++)tempVL[j+tempInt]=voxList[j];
        //memcpy(&(tempVL[0]),tempList,tempInt*sizeof(int));
        //memcpy(&(tempVL[tempInt]),voxList,*nIn*sizeof(int));
        TIDY(voxList);
        voxList=tempVL;
        tempVL=NULL;

        tempFL=ialloc(*nIn+tempInt,"",0);
        for(j=0;j<tempInt;j++)tempFL[j]=i;
        for(j=0;j<*nIn;j++)tempFL[j+tempInt]=fileList[0][j];
        //memcpy(&(tempFL[tempInt]),*fileList,*nIn*sizeof(int));
        TIDY(*fileList);
        *fileList=tempFL;
        tempFL=NULL;

        tempR=dalloc(*nIn+tempInt,"",0);
        for(j=0;j<tempInt;j++)tempR[j]=tempRanges[j];
        for(j=0;j<*nIn;j++)tempR[j+tempInt]=rangeList[0][j];
        //memcpy(&(tempR[0]),tempRanges,tempInt*sizeof(double));
        //memcpy(&(tempR[tempInt]),*rangeList,*nIn*sizeof(double));
        TIDY(*rangeList);
        *rangeList=tempR;
        tempR=NULL;

      }else{   /*insert at end*/
        if(!(voxList=(int *)realloc(voxList,(*nIn+tempInt)*sizeof(int)))){
          fprintf(stderr,"Error in voxel index reallocation within multi-voxel tracing, allocating %lu\n",(*nIn+tempInt)*sizeof(int));
          exit(1);
        }
        if(!(*fileList=(int *)realloc(*fileList,(*nIn+tempInt)*sizeof(int)))){
          fprintf(stderr,"Error in file index reallocation within multi-voxel tracing, allocating %lu\n",(*nIn+tempInt)*sizeof(int));
          exit(1);
        }
        if(!(*rangeList=(double *)realloc(*rangeList,(*nIn+tempInt)*sizeof(double)))){
          fprintf(stderr,"Error in range reallocation within multi-voxel tracing, allocating %lu\n",(*nIn+tempInt)*sizeof(double));
          exit(1);
        }
        for(j=0;j<tempInt;j++)voxList[j+*nIn]=tempList[j];
        for(j=0;j<tempInt;j++)rangeList[0][j+*nIn]=tempRanges[j];
        //memcpy(&(voxList[*nIn]),tempList,tempInt*sizeof(int));
        //memcpy(&(*rangeList[*nIn]),tempRanges,tempInt*sizeof(double));
        for(j=0;j<tempInt;j++)fileList[0][*nIn+j]=i;
      }
      TIDY(tempList);
      TIDY(tempRanges);
    }else{  /*first allocation*/
      voxList=tempList;
      *rangeList=tempRanges;
      *fileList=ialloc(tempInt,"voxel file index",0);
      for(j=0;j<tempInt;j++)fileList[0][j]=i;
      tempList=NULL;
      tempRanges=NULL;
    }

    *nIn+=tempInt;
  }/*file loop*/

  return(voxList);
}/*findAllVoxels*/


/*###########################################################################*/
/*find intersecting voxels*/

int *findVoxels(double *grad,double xCent,double yCent,double zCent,double *bounds,double *vRes,int *nPix,int vX,int vY,int vZ,double **rangeList)
{
  int xBin=0,yBin=0,zBin=0;
  int xDir=0,yDir=0,zDir=0;
  int nextXbin=0,nextYbin=0,nextZbin=0;
  int *pixList=NULL;
  int *markInt(int,int *,int);
  float vectX=0,vectY=0,vectZ=0;
  double *markDo(int,double *,double);
  double nextX=0,nextY=0,nextZ=0;
  double rX=0,rY=0,rZ=0;
  double zen=0,az=0;
  double coord[3],vect[3];
  void findClosestFacet(double *,double *,double *,int,int,int);

  /*set vector*/
  if(grad[2]<-10.0){   /*grad is a polar vector*/
    zen=grad[0];
    az=grad[1];
    vectX=(float)(sin(zen)*sin(az));
    vectY=(float)(sin(zen)*cos(az));
    vectZ=(float)cos(zen);
  }else{  /*grad is a cartesian vector*/
    vectX=(float)grad[0];
    vectY=(float)grad[1];
    vectZ=(float)grad[2];
    zen=atan2(sqrt(grad[0]*grad[0]+grad[1]*grad[1]),grad[2]);
    az=atan2(grad[0],grad[1]);
  }

  /*set origin*/
  coord[0]=xCent;
  coord[1]=yCent;
  coord[2]=zCent;

  /*determine vector directions*/
  if(vectX!=0.0)xDir=(vectX>0.0)?1:-1;
  else          xDir=0;
  if(vectY!=0.0)yDir=(vectY>0.0)?1:-1;
  else          yDir=0;
  if(vectZ!=0.0)zDir=(vectZ>0.0)?1:-1;
  else          zDir=0;

  /*if we are outside test for intersection and move point to start of voxels*/
  if((coord[0]>bounds[3])||(coord[1]>bounds[4])||(coord[0]<bounds[0])||\
     (coord[1]<bounds[1])||(coord[2]>bounds[5])||(coord[2]<bounds[2])){
    /*for each voxel bound facet determine the range to. Reset coord as bound intersection*/
    vect[0]=vectX;vect[1]=vectY;vect[2]=vectZ;
    findClosestFacet(&(coord[0]),&(bounds[0]),&(vect[0]),xDir,yDir,zDir);
  }/*outside but heading towards voxel space check*/

  /*mark first voxel*/
  xBin=(int)((coord[0]-bounds[0])/vRes[0]);
  yBin=(int)((coord[1]-bounds[1])/vRes[1]);
  zBin=(int)((coord[2]-bounds[2])/vRes[2]);
  *nPix=0;
  if((xBin>=0)&&(xBin<vX)&&(yBin>=0)&&(yBin<vY)&&(zBin>=0)&&(zBin<vZ)){
    pixList=markInt(*nPix,pixList,xBin+vX*yBin+vX*vY*zBin);
    rangeList[0]=markDo(*nPix,rangeList[0],sqrt((coord[0]-xCent)*\
        (coord[0]-xCent)+(coord[1]-yCent)*(coord[1]-yCent)+(coord[2]-zCent)*(coord[2]-zCent)));
    (*nPix)++;
  }

  /*loop along intersections*/
  while((coord[0]>=bounds[0])&&(coord[0]<=bounds[3])&&(coord[1]>=bounds[1])&&(coord[1]<=bounds[4])&&(coord[2]>=bounds[2])&&(coord[2]<=bounds[5])){
    /*where would it go next*/
    nextXbin=xBin+xDir;
    nextYbin=yBin+yDir;
    nextZbin=zBin+zDir;
    nextX=(double)nextXbin*vRes[0]+bounds[0];
    nextY=(double)nextYbin*vRes[1]+bounds[1];
    nextZ=(double)nextZbin*vRes[2]+bounds[2];

    /*range to next hit*/
    if(vectX!=0.0)rX=(nextX-coord[0])/vectX;
    else          rX=0.0;
    if(vectY!=0.0)rY=(nextY-coord[1])/vectY;
    else          rY=0.0;
    if(vectZ!=0.0)rZ=(nextZ-coord[2])/vectZ;
    else          rZ=0.0;

    /*if nothing has changed*/
    if(rX==0.0)rX=(bounds[3]-bounds[0])*1000.0;
    if(rY==0.0)rY=(bounds[4]-bounds[1])*1000.0;
    if(rZ==0.0)rZ=(bounds[5]-bounds[2])*1000.0;

    if((rX<rY)&&(rX<rZ)){       /*x change*/
      coord[0]=nextX;
      coord[1]+=rX*vectY;
      coord[2]+=rX*vectZ;
      xBin+=xDir;
    }else if((rY<rZ)&&(rY<rX)){  /*y change*/
      coord[0]+=rY*vectX;
      coord[1]=nextY;
      coord[2]+=rY*vectZ;
      yBin+=yDir;
    }else if((rZ<rY)&&(rZ<rX)){  /*z change*/
      coord[0]+=rZ*vectX;
      coord[1]+=rZ*vectY;
      coord[2]=nextZ;
      zBin+=zDir;
    }else if((rX<rZ)&&(rX==rY)){ /*x and y change*/
      coord[0]=nextX;
      coord[1]+=rX*vectY;
      coord[2]+=rX*vectZ;
      xBin+=xDir;
      yBin+=yDir;
    }else if((rX<rY)&&(rX==rZ)){ /*x and z change*/
      coord[0]=nextX;
      coord[1]+=rX*vectY;
      coord[2]+=rX*vectZ;
      xBin+=xDir;
      zBin+=zDir;
    }else if((rY<rX)&&(rY==rZ)){ /*y and z change*/
      coord[0]+=rY*vectX;
      coord[1]=nextY;
      coord[2]+=rY*vectZ;
      yBin+=yDir;
      zBin+=zDir;
    }else if((rX==rY)&&(rX==rZ)){ /*x and y and z change*/
      coord[0]=nextX;
      coord[1]+=rX*vectY;
      coord[2]+=rX*vectZ;
      xBin+=xDir;
      yBin+=yDir;
      zBin+=zDir;
    }else{
      fprintf(stderr,"Intersection test issue: %f %f %f\n",rX,rY,rZ);
      exit(1);
    }

    /*Is it within bounds? If so record the hit*/
    if((xBin>=0)&&(xBin<vX)&&(yBin>=0)&&(yBin<vY)&&(zBin>=0)&&(zBin<vZ)){
      /*record*/
      pixList=markInt(*nPix,pixList,xBin+vX*yBin+vX*vY*zBin);
      rangeList[0]=markDo(*nPix,rangeList[0],sqrt((coord[0]-xCent)*(coord[0]-xCent)+\
                 (coord[1]-yCent)*(coord[1]-yCent)+(coord[2]-zCent)*(coord[2]-zCent)));
      (*nPix)++;
    }
  }/*intersection loop*/

  /*record the range to the very end*/
  rangeList[0]=markDo(*nPix,rangeList[0],sqrt((coord[0]-xCent)*\
      (coord[0]-xCent)+(coord[1]-yCent)*(coord[1]-yCent)+(coord[2]-zCent)*(coord[2]-zCent)));

  return(pixList);
}/*findVoxels*/


/*###########################################################################*/
/*find bounds of filled voxels*/

double *findVoxelBounds(int *voxList,int nIn,voxStruct *vox,tlsVoxMap *map,float *grad,double x0,double y0,double z0)
{
  int i=0,k=0,vInd=0;
  int ii=0,jj=0,kk=0;
  int xBin=0,yBin=0,zBin=0;
  float zen=0,az=0;
  double *bounds=NULL,vect[3];

  bounds=dalloc(6,"bounds",0);
  bounds[0]=bounds[1]=bounds[2]=10000000000.0;
  bounds[3]=bounds[4]=bounds[5]=-10000000000.0;

  zen=(float)atan2(sqrt((double)grad[0]*(double)grad[0]+(double)grad[1]*(double)grad[1]),(double)grad[2]);
  az=(float)atan2((double)grad[0],(double)grad[1]);

  /*loop over intersected voxels*/
  for(i=0;i<nIn;i++){
    vInd=voxList[i];
    if(map->nIn[vInd]>0){
      xBin=vInd%(vox->nX*vox->nY);
      yBin=(vInd-xBin)%vox->nX;
      zBin=((vInd-xBin)-yBin*vox->nX)/(vox->nX*vox->nY);

      /*each corner in turn*/
      for(ii=0;ii<2;ii++){
        for(jj=0;jj<2;jj++){
          for(kk=0;kk<2;kk++){
            /*rotate to beam vector*/
            vect[0]=(double)(xBin+ii)*vox->res[0]+vox->bounds[0]-x0;
            vect[1]=(double)(yBin+jj)*vox->res[1]+vox->bounds[1]-y0;
            vect[2]=(double)(zBin+kk)*vox->res[2]+vox->bounds[2]-z0;
            rotateZ(vect,(double)(-1.0*az));
            rotateX(vect,(double)(-1.0*zen));
            for(k=0;k<3;k++){   /*bound check*/
              if(vect[k]<bounds[k])bounds[k]=vect[k];
              if(vect[k]>bounds[k+3])bounds[k+3]=vect[k];
            }/*bound check*/
          }
        }
      }/*eight corner loop*/
    }/*filled voxel check*/
  }/*intersected voxel loop*/



  return(bounds);
}/*findVoxelBounds*/


/*###########################################################################*/
/*set waveform elevation along vector*/

void setWaveformRange(float *range,double z0,float *grad,int nBins,float res)
{
  int i=0;
  float r=0;
  double zen=0;

  zen=(float)atan2(sqrt((double)grad[0]*(double)grad[0]+(double)grad[1]*(double)grad[1]),(double)grad[2]);

  for(i=0;i<nBins;i++){
    r=(float)i*res;
    range[i]=z0+r*(float)cos(zen);
  }

  return;
}/*setWaveformRange*/


/*###########################################################################*/
/*read bounds for voxels from TLS*/

void readBoundsFromTLS(double *bounds,char **inList,int nScans)
{
  int i=0,k=0;
  uint32_t j=0,tInd=0;
  double x=0,y=0,z=0;
  double xCent=0,yCent=0,zCent=0;
  tlsScan *tempTLS=NULL;

  bounds[0]=bounds[1]=bounds[2]=10000000000.0;
  bounds[3]=bounds[4]=bounds[5]=-10000000000.0;

  for(i=0;i<nScans;i++){  /*file loop*/
    readTLSpolarBinary(inList[i],0,&tempTLS);
    for(j=0;j<tempTLS->nBeams;j++){/*point loop*/
      /*update TLS beams if needed*/
      readTLSpolarBinary(inList[i],j,&tempTLS);
      tInd=j-tempTLS->pOffset;   /*update index to account for buffered memory*/
      /*beam origin*/
      xCent=(double)tempTLS->beam[tInd].x+tempTLS->xOff;
      yCent=(double)tempTLS->beam[tInd].y+tempTLS->yOff;
      zCent=(double)tempTLS->beam[tInd].z+tempTLS->zOff;

      for(k=0;k<tempTLS->beam[tInd].nHits;k++){  /*hit loop*/
        /*point coordinate*/
        x=xCent+tempTLS->beam[tInd].r[k]*sin(tempTLS->beam[tInd].az)*sin(tempTLS->beam[tInd].zen);
        y=yCent+tempTLS->beam[tInd].r[k]*cos(tempTLS->beam[tInd].az)*sin(tempTLS->beam[tInd].zen);
        z=zCent+tempTLS->beam[tInd].r[k]*cos(tempTLS->beam[tInd].zen);

        /*determine bounds*/
        if(x<bounds[0])bounds[0]=x;
        if(y<bounds[1])bounds[1]=y;
        if(z<bounds[2])bounds[2]=z;
        if(x>bounds[3])bounds[3]=x;
        if(y>bounds[4])bounds[4]=y;
        if(z>bounds[5])bounds[5]=z;
      }/*hit loop*/
    }/*point loop*/
    tempTLS=tidyTLScan(tempTLS);
  }/*file loop*/

  return;
}/*readBoundsFromTLS*/


/*###########################################################################*/
/*clip x and y bounds to a beam*/

void beamVoxelBounds(double *origin,float *grad,float fSigma,char gaussFoot,double *bounds)
{
  int i=0;
  float rad=0;
  double x=0,y=0;

  if(gaussFoot)rad=determineGaussSep(fSigma,0.001);
  else         rad=fSigma;

  /*put beam start at top of voxel space*/
  origin[2]=bounds[5];

  bounds[0]=bounds[1]=100000000000.0;
  bounds[3]=bounds[4]=-100000000000.0;

  for(i=-1;i<=1;i+=2){  /*loop over edges*/
    /*top*/
    x=origin[0]+(float)i*rad;
    y=origin[1]+(float)i*rad;
    if(x<bounds[0])bounds[0]=x;
    if(x>bounds[3])bounds[3]=x;
    if(y<bounds[1])bounds[1]=y;
    if(y>bounds[4])bounds[4]=y;

    /*bottom*/
    x=origin[0]+(float)i*rad+((float)origin[2]-(float)bounds[2])*grad[0];
    y=origin[1]+(float)i*rad+((float)origin[2]-(float)bounds[2])*grad[1];
    if(x<bounds[0])bounds[0]=x;
    if(x>bounds[3])bounds[3]=x;
    if(y<bounds[1])bounds[1]=y;
    if(y>bounds[4])bounds[4]=y;
  }/*edges loop*/

  return;
}/*beamVoxelBounds*/


/*############################################*/
/*the below are TLS specific things*/

/*###############################################*/
/*mark an image with a point and add up area*/

void makeBinImage(double *coord,float *area,float *count,char *rImage,int numb,float gap,double r,uint16_t refl,float beamRad,int rNx,int rNy,float rimRes,lidVoxPar *tlsPar)
{
  int xInd=0,yInd=0,rPlace=0;
  int xStart=0,xEnd=0;
  int yStart=0,yEnd=0;
  double xIcent=0,yIcent=0;
  float rad=0;
  float maxRsepSq=0,rSepSq=0;

  if(gap<tlsPar->minGap)gap=tlsPar->minGap;
  rad=tlsPointSize(r,refl,tlsPar->beamTanDiv,tlsPar->beamRad,tlsPar->minRefl,tlsPar->maxRefl,tlsPar->appRefl,gap);

  (*area)+=rad*rad*tlsPar->appRefl/(beamRad*beamRad);
  (*count)+=1.0;

  /*range image*/
  xIcent=(int)((coord[0]+beamRad)/rimRes);
  yIcent=(int)((coord[1]+beamRad)/rimRes);
  xStart=xIcent-(int)(rad/rimRes);
  xEnd=xIcent+(int)(rad/rimRes);
  yStart=yIcent-(int)(rad/rimRes);
  yEnd=yIcent+(int)(rad/rimRes);
  maxRsepSq=rad*rad;

  if(xStart<0)xStart=0;      /*enforce bounds*/
  if(xEnd>=rNx)xEnd=rNx-1;
  if(yStart<0)yStart=0;
  if(yEnd>=rNy)yEnd=rNy-1;   /*enforce bounds*/

  for(xInd=xStart;xInd<=xEnd;xInd++){
    for(yInd=yStart;yInd<=yEnd;yInd++){
      rSepSq=(float)((xInd-xIcent)*(xInd-xIcent)+(yInd-yIcent)*(yInd-yIcent))*rimRes*rimRes;
      if(rSepSq<=maxRsepSq){
        rPlace=yInd*rNx+xInd;
        if(rImage[rPlace]==0)rImage[rPlace]=1;
      }/*check within point*/
    }/*loop around point*/
  }/*loop around point*/

  return;
}/*makeBinImage*/


/*############################################*/
/*read cnaopy bounds*/

void readCanBounds(canBstruct *canB,char *canNamen,double *bounds)
{
  int xBin=0,yBin=0;
  int totN=0,place=0;
  double x=0,y=0;
  double sX=0,eX=0;     /*first and second points to determine res*/
  double buff=0;
  char line[400],temp1[100],temp2[100];
  char temp3[100],temp4[100];
  FILE *ipoo=NULL;

  if((ipoo=fopen(canNamen,"r"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",canNamen);
    exit(1);
  }

  buff=20.0;

  canB->cUbound[0]=canB->cUbound[1]=100000000.0;
  canB->cUbound[2]=canB->cUbound[3]=-100000.0;

  sX=eX=-1.0;

  /*read number of pixels*/
  xBin=0;
  while(fgets(line,400,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      if(sscanf(line,"%s %s",temp1,temp2)==2){
        x=atof(temp1);
        y=atof(temp2);

        /*first and second point for res*/
        if(sX<0.0)sX=x;
        else if(eX<0.0){
          if(sX!=x)eX=x;
        }
        /*if(sY<0.0)sY=y;
        else if(eY<0.0){
          if(sY!=y)eY=y;
        }*/

        if((x>=(bounds[0]-buff))&&(x<=(bounds[3]+buff))&&(y>=(bounds[1]-buff))&&(y<=(bounds[4]+buff))){
          if(x<canB->cUbound[0])canB->cUbound[0]=x;
          if(y<canB->cUbound[1])canB->cUbound[1]=y;
          if(x>canB->cUbound[2])canB->cUbound[2]=x;
          if(y>canB->cUbound[3])canB->cUbound[3]=y;
        }

        xBin++;
      }
    }
  }

  if(canB->cUbound[2]<0.0){  /*then there is no data*/
    fprintf(stderr,"Canopy bound issue %f in %s from %d xBound %.2f %.2f yBound %.2f %.2f\n",canB->cUbound[2],canNamen,xBin,bounds[0],bounds[3],bounds[1],bounds[4]);
    exit(1);
  }
  canB->cUbound[0]-=buff;
  canB->cUbound[1]-=buff;
  canB->cUbound[2]+=buff;
  canB->cUbound[3]+=buff;

  canB->cRes=eX-sX;
  canB->cNx=(int)((canB->cUbound[2]-canB->cUbound[0])/canB->cRes)+1;
  canB->cNy=(int)((canB->cUbound[3]-canB->cUbound[1])/canB->cRes)+1;

  if((canB->cNx<=0)||(canB->cNy<=0)){
    fprintf(stderr,"canopy pixel error\n");
    fprintf(stderr,"x %f %f y %f %f res %f\n",canB->cUbound[2],canB->cUbound[0],canB->cUbound[3],canB->cUbound[1],canB->cRes);
    exit(1);
  }

  totN=canB->cNx*canB->cNy;
  canB->canMax=falloc((uint64_t)totN,"canMax",0);

  if(fseek(ipoo,(long)0,SEEK_SET)){ /*rewind to start of file*/
    fprintf(stderr,"fseek error\n");
    exit(1);
  }
  while(fgets(line,400,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      if(sscanf(line,"%s %s %s %s",temp1,temp2,temp3,temp4)==4){  /*read top and bottom*/
        x=atof(temp1);
        y=atof(temp2);

        xBin=(int)((x-canB->cUbound[0])/canB->cRes+0.5);
        yBin=(int)((y-canB->cUbound[1])/canB->cRes+0.5);

        /*within bounds check*/
        if((xBin>=0)&&(xBin<canB->cNx)&&(yBin>=0)&&(yBin<canB->cNy)){
          place=yBin*canB->cNx+xBin;
          canB->canMin[place]=atof(temp3);
          canB->canMax[place]=atof(temp4);
        }
      }else if(sscanf(line,"%s %s %s",temp1,temp2,temp3)==3){  /*only read bottom*/
        x=atof(temp1);
        y=atof(temp2);

        xBin=(int)((x-canB->cUbound[0])/canB->cRes+0.5);
        yBin=(int)((y-canB->cUbound[1])/canB->cRes+0.5);

        if((xBin>=0)&&(xBin<canB->cNx)&&(yBin>=0)&&(yBin<canB->cNy)){
          place=yBin*canB->cNx+xBin;
          canB->canMin[place]=atof(temp3);
          canB->canMax[place]=10000.0;   /*a large number*/
        }
      }
    }
  }

  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }

  return;
}/*readCanBounds*/


/*############################################*/
/*determine indices of bounds*/

void setCanInd(int *wStart,int *wEnd,double xCent,double yCent,double zCent,lasFile *lasIn,canBstruct *canB)
{
  int i=0;
  int xBin=0,yBin=0,place=0;
  double x=0,y=0,z=0;
  double sSep=0;

  (*wStart)=-1;
  (*wEnd)=lasIn->waveLen+2;

  for(i=0;i<lasIn->waveLen;i++){
    binPosition(&x,&y,&z,i,xCent,yCent,zCent,lasIn->time,lasIn->grad);

    xBin=(int)((x-canB->cUbound[0])/(double)canB->cRes+0.5);
    yBin=(int)((y-canB->cUbound[1])/(double)canB->cRes+0.5);

    if((xBin<0)||(xBin>=canB->cNx)||(yBin<0)||(yBin>=canB->cNy)){
      //fprintf(stderr,"Canopy bounds not wide enough for ind %d %d %d %d point %f %f can bound %f %f\n",xBin,yBin,canB->cNx,canB->cNy,x,y,canB->cUbound[0],canB->cUbound[1]);
      (*wStart)=0;
      (*wEnd)=lasIn->waveLen;
      break;
    }

    place=yBin*canB->cNx+xBin;

    if(((*wStart)<0)&&(z>canB->canMax[place])){
      sSep=sqrt((z-canB->canMax[place])*(z-canB->canMax[place]));
      if(sSep<0.15)(*wStart)=i;
    }

    if(z<canB->canMin[place]){
      (*wEnd)=i+1;  /*add one as we loop up to but not over here*/
      break;
    }
  }
  if((*wStart)<0)(*wStart)=0;
  if((*wEnd)>lasIn->waveLen)(*wEnd)=lasIn->waveLen;

  return;
}/*setCanInd*/


/*#######################################*/
/*impose canopy bounds*/

void imposeCanBound(float *denoised,int wStart,int wEnd,int numb)
{
  int i=0;

  for(i=wStart;i>=0;i--)denoised[i]=0.0;
  for(i=wEnd;i<numb;i++)denoised[i]=0.0;

  return;
}/*imposeCanBound*/


/*############################################*/
/*set min and max for this pixel*/

double setCanGround(double x,double y,canBstruct *canB)
{
  int xBin=0,yBin=0,place=0;
  double z=0;

  xBin=(int)((x-canB->cUbound[0])/(double)canB->cRes+0.5);
  yBin=(int)((y-canB->cUbound[1])/(double)canB->cRes+0.5);

  if((xBin>=0)&&(xBin<canB->cNx)&&(yBin>=0)&&(yBin<canB->cNy)){
    place=yBin*canB->cNx+xBin;
    z=canB->canMin[place];
  }else{
    //fprintf(stderr,"Canopy bounds not wide enough for z %d %d %d %d\n",xBin,yBin,canB->cNx,canB->cNy);
    z=-10000.0;
  }

  return(z);
}/*setCanBounds*/



/*#######################################*/
/*set higher voxels bank*/

void setTopVoxBlank(voxStruct *vox)
{
  int i=0,j=0,k=0;
  int numb=0,place=0;
  float total=0;

  for(i=0;i<vox->nX;i++){
    for(j=0;j<vox->nY;j++){
      /*loop down to find first filled voxel*/
      total=0.0;
      for(k=vox->nZ-1;k>=0;k--){
        place=k*vox->nX*vox->nY+j*vox->nX+i;
        for(numb=0;numb<vox->nScans;numb++)total+=vox->hits[numb][place]+vox->miss[numb][place];
        if(total>0.0){
          break;
        }
      }
      /*set all above here to blank*/
      for(k++;k<vox->nZ;k++){
        place=k*vox->nX*vox->nY+j*vox->nX+i;
        total=0.0;
        for(numb=0;numb<vox->nScans;numb++)total+=vox->hits[numb][place]+vox->miss[numb][place];
        if(total>0.0){
          fprintf(stderr,"Balls\n");
          exit(1);
        }
        for(numb=0;numb<vox->nScans;numb++)vox->miss[numb][place]=1.0;
      }
    }/*y loop*/
  }/*x loop*/

  return;
}/*setTopVoxBlank*/


/*#######################################*/
/*voxelise from waveform lidar*/

void voxelate(voxStruct *vox,float *wave,lasFile *lasIn,double xCent,double yCent,double zCent,float beamRad)
{
  int i=0,j=0,nTot=0;
  int lastBin=0,zInd=0;
  int *voxList=NULL;
  int *nHitCont=NULL;
  double x0=0,y0=0,z0=0;
  double zTop=0,zBot=0;
  float *meanHit=NULL;
  char *passList=NULL;

  /*determine when the beam is blocked*/
  lastBin=lasIn->waveLen+1;
  for(i=0;i<lasIn->waveLen;i++){
    if(wave[i]>0.0)lastBin=i;
  }

  /*map for all beams*/
  binPosition(&x0,&y0,&z0,0,xCent,yCent,zCent,lasIn->time,lasIn->grad);
  voxList=beamVoxels(lasIn->grad,x0,y0,z0,&(vox->bounds[0]),&(vox->res[0]),vox->nX,vox->nY,vox->nZ,&nTot,beamRad,NULL,-1.0);
  passList=challoc((uint64_t)nTot,"pass list",0);
  for(j=0;j<nTot;j++)passList[j]=0;

  /*gap fraction for each voxel for this beam*/
  meanHit=falloc((uint64_t)nTot,"meanHit",0);
  nHitCont=ialloc(nTot,"nHitCont",0);
  for(i=0;i<nTot;i++){
    meanHit[i]=0.0;
    nHitCont[i]=0;
  }

  /*loop along until last beam*/
  for(i=0;i<lastBin;i++){
    binPosition(&x0,&y0,&z0,i,xCent,yCent,zCent,lasIn->time,lasIn->grad);
    for(j=0;j<nTot;j++){
      zInd=(int)((float)voxList[j]/(float)(vox->nX*vox->nY));
      zBot=(double)(zInd)*(double)vox->res[2]+vox->bounds[2];
      zTop=(double)(zInd+1)*(double)vox->res[2]+vox->bounds[2];
      if((z0<zTop)&&(z0>zBot)){  /*check that this voxel is within bin*/
        vox->hits[0][voxList[j]]+=wave[i]/0.15;
        vox->miss[0][voxList[j]]+=(1.0-wave[i])/0.15;
        vox->contN[voxList[j]]++;
      }
    }/*voxel loop*/
  }/*bin loop*/

  TIDY(meanHit);
  TIDY(nHitCont);
  TIDY(passList);
  TIDY(voxList);
  return;
}/*voxelate*/


/*############################################*/
/*write out ASCII voxels*/

void writeAsciiVox(voxStruct *vox,char *outRoot)
{
  int i=0,j=0,k=0;
  int place=0;
  char namen[200];
  FILE *opoo=NULL;

  sprintf(namen,"%s.vox",outRoot);
  if((opoo=fopen(namen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",namen);
    exit(1);
  }
  fprintf(opoo,"# 1 x, 2 y, 3 z, 4 cover, 5 nALS\n");
  /*write out the voxel map*/
  for(i=0;i<vox->nX;i++){
    for(j=0;j<vox->nY;j++){
      for(k=0;k<vox->nZ;k++){
        place=k*vox->nX*vox->nY+j*vox->nX+i;
        if(vox->contN[place]>0){
          fprintf(opoo,"%f %f %f %f %d\n",((float)i+0.5)*vox->res[0]+vox->bounds[0],\
            ((float)j+0.5)*vox->res[1]+vox->bounds[1],((float)k+0.5)*vox->res[2]+\
            vox->bounds[2],vox->hits[0][place]/(vox->hits[0][place]+vox->miss[0][place]),vox->contN[place]);
        }else{
          fprintf(opoo,"%f %f %f %f %d\n",((float)i+0.5)*vox->res[0]+vox->bounds[0],\
            ((float)j+0.5)*vox->res[1]+vox->bounds[1],((float)k+0.5)*vox->res[2]+\
            vox->bounds[2],-1.0,vox->contN[place]);
        }
      }
    }
  }
  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Written to %s\n",namen);

  return;
}/*writeAsciiVox*/

/*the end*/
/*#######################################*/

