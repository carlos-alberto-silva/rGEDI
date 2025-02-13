#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "libLasRead.h"
#include "libLasProcess.h"
#include "tiffWrite.h"


/*#########################*/
/*# Makes geotiffs from   #*/
/*# las files   2015      #*/
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


/*##################################################*/
/*control structure*/

typedef struct{
  /*input output*/
  char **inList;
  int nFiles;
  char outNamen[1000];
  char bNamen[1000];
  FILE *bFile;        /*bounds file output*/

  /*switches*/
  char drawInt;       /*intensity image switch*/
  char drawHeight;    /*height.elevation image switch*/
  char drawDens;      /*draw ensity images*/
  char onlyGround;    /*use only ground points switch*/
  char gapFill;       /*fill gap switch*/
  char findDens;      /*find point and footprint density*/
  char drawVegVol;    /*draw volume of vegetation, can set height thesholds*/
  char drawCov;       /*draw canopy cover switch*/
  char writeBounds;   /*write out file bounds*/
  uint64_t pBuffSize; /*point buffer rading size in bytes*/
  char printNpoint;   /*write number of points to the screen*/
  char charImage;     /*char or float image*/

  /*input filters*/
  int16_t maxScanAng; /*maximum scan angle*/

  /*enforced bounds*/
  char findBounds;    /*find bounds from data switch*/
  double bounds[4];   /*minX, minY, maxX, maxY*/

  /*geotiff parts*/
  float res;
  float maxDN;
  uint16_t epsg;

  /*CHM settings*/
  float hThresh;  /*threshold to avoid noise on CHM*/
  float hRes;     /*resolution to make psuedo-waveform for CHM*/
  float hRange;   /*maximum expected extent of points*/

  /*height bounds for volume*/
  float minVolH;      /*minimum height to do volume over*/
  float maxVolH;      /*maximum height to do volume over*/
  uint16_t maxVint; /*maximum vegetatiob intensity*/
}control;


/*##################################################*/
/*image structure*/

typedef struct{
  int nX;
  int nY;
  uint64_t *nIn;  /*number of points in pixel*/
  uint64_t *nCan; /*number of canopy returns*/
  int *nFoot;     /*footprint density*/
  float *jimlad;  /*FLOATING POINT IMAGE*/
  float min;      /*min intensity*/
  float max;      /*max intensity*/
  unsigned char *image;  /*IMAGE TO WRITE*/
  double minX;
  double minY;
  double maxX;
  double maxY;
  double geoL[2];
  int geoI[2];
  uint16_t epsg;   /*EPSG code*/
  int maxPoint;    /*max number of points per pixels*/
  int maxFoot;     /*max number of footprints per pixels*/
  float **heightStack;  /*if making a height image, a stack of data*/
  float *minH;          /*bottom of height array*/
  float hRes;           /*height resolution to use*/
  float hThresh;        /*threshold to avoid noise in height*/
  int nHeight;
  float hRange;
}imageStruct;


/*##################################################*/
/*main*/

int main(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;
  control *readCommands(int,char **);
  lasFile *las=NULL;
  imageStruct *image=NULL;
  imageStruct *allocateImage(control *);
  imageStruct *tidyImage(imageStruct *);
  void collateImage(control *,lasFile *,imageStruct *);
  void finishImage(control *,imageStruct *);
  void writeImage(control *,imageStruct *);
  void fillGaps(control *,imageStruct *);
  void writeFileBounds(lasFile *,char *,control *);
  void updateBounds(double *,lasFile *);

  /*read command line*/
  dimage=readCommands(argc,argv);

  /*read file bounds if needed*/
  if(dimage->writeBounds||dimage->findBounds||dimage->printNpoint){
    for(i=0;i<dimage->nFiles;i++){
      las=readLasHead(dimage->inList[i],dimage->pBuffSize);
      if(dimage->epsg==0)readLasGeo(las);
      else               las->epsg=dimage->epsg;
      if(dimage->writeBounds)writeFileBounds(las,dimage->inList[i],dimage);
      if(dimage->printNpoint)fprintf(stdout,"nPoints %s %u\n",dimage->inList[i],las->nPoints);
      if(dimage->findBounds)updateBounds(dimage->bounds,las);
      las=tidyLasFile(las);
    }
  }

  /*create imge if needed*/
  if(dimage->drawInt||dimage->drawHeight||dimage->findDens||dimage->drawDens||dimage->drawCov||dimage->drawVegVol){
    /*allocate image array*/
    image=allocateImage(dimage);

    /*loop over las files*/
    for(i=0;i<dimage->nFiles;i++){
      /*re-read data*/
      las=readLasHead(dimage->inList[i],dimage->pBuffSize);

      /*is this file needed?*/
      if(checkFileBounds(las,dimage->bounds[0],dimage->bounds[2],dimage->bounds[1],dimage->bounds[3]))collateImage(dimage,las,image);
      las=tidyLasFile(las);
    }
    /*finish off image*/
    finishImage(dimage,image);

    /*fill gaps if needed*/
    if(dimage->gapFill)fillGaps(dimage,image);
  }/*image drawing check*/

  /*write image*/
  if(dimage->drawInt||dimage->drawHeight||dimage->drawCov||dimage->drawVegVol)writeImage(dimage,image);
  if(dimage->writeBounds)fprintf(stdout,"Written to %s\n",dimage->bNamen);


  /*tidy up arrays*/
  image=tidyImage(image);
  if(dimage){
    TTIDY((void **)dimage->inList,dimage->nFiles);
    if(dimage->bFile){
      fclose(dimage->bFile);
      dimage->bFile=NULL;
    }
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*##################################################*/
/*tidy image structure*/

imageStruct *tidyImage(imageStruct *image)
{
  if(image){
    TIDY(image->jimlad);
    TIDY(image->nIn);
    TIDY(image->image);
    TIDY(image->nFoot);
    TIDY(image->nCan);
    TIDY(image->minH);
    TTIDY((void **)image->heightStack,image->nHeight);
    TIDY(image);
  }
  return(image);
}/*tidyImage*/


/*##################################################*/

void writeFileBounds(lasFile *las,char *namen,control *dimage)
{
  fprintf(dimage->bFile,"%s %.2f %.2f %.2f %.2f %.2f %.2f\n",namen,las->minB[0],las->minB[1],las->minB[2],las->maxB[0],las->maxB[1],las->maxB[2]);
  
  return;
}/*writeFileBounds*/


/*##################################################*/
/*fill gaps in an image*/

void fillGaps(control *dimage,imageStruct *image)
{
  int i=0,j=0;
  uint64_t place=0;
  float *newJimlad=NULL,newFloat=0;
  unsigned char *newImage=NULL,newChar=0;
  void fillDTMhole(int,int,uint64_t,imageStruct *,control *,char,float *,unsigned char *);

  if(dimage->charImage)newImage=uchalloc(image->nX*image->nY,"newImage",0);
  else                 newJimlad=falloc(image->nX*image->nY,"newImage",0);

  /*loop over image*/
  for(j=0;j<image->nY;j++){
    for(i=0;i<image->nX;i++){
      place=(uint64_t)j*(uint64_t)image->nX+(uint64_t)i;
      /*print progress*/
      if(((place*100)%((uint64_t)image->nX*(uint64_t)image->nY))<100){
        fprintf(stdout,"Gap filling %d%%\n",(int)((place*100)/((uint64_t)image->nX*(uint64_t)image->nY)));
      }

      /*is the data missing?*/
      if(image->nIn[place]==0){
        if(dimage->charImage){
          fillDTMhole(i,j,place,image,dimage,1,&newFloat,&newChar);
          newImage[place]=newChar;
        }else{
          fillDTMhole(i,j,place,image,dimage,2,&newFloat,&newChar);
          newJimlad[place]=newFloat;
        }
      }else{
        if(dimage->charImage)newImage[place]=image->image[place];
        else                 newJimlad[place]=image->jimlad[place];
       }
    }/*x loop*/
  }/*y loop*/

  if(dimage->charImage)image->image=newImage;
  else                 image->jimlad=newJimlad;

  newImage=NULL;
  newJimlad=NULL;

  return;
}/*fillGaps*/


/*##################################################*/
/*fill a char image*/

void fillDTMhole(int i,int j,uint64_t place0,imageStruct *image,control *dimage,char mode,float *newFloat,unsigned char *newChar)
{
  int ii=0,ti=0,tj=0;
  int nIn=0,window=0;
  int maxWindow=0;
  int minIn=0;
  int nTest=0;
  uint64_t *indList=NULL,place=0;
  uint64_t *setFillList(int,int,int,int *,imageStruct *);
  float dx=0,dy=0;
  float fill=0,weight=0;
  float totWeight=0,dist=0;

  /*set window size*/
  window=1;
  maxWindow=(int)(20.0/dimage->res);
  minIn=8;
  totWeight=0.0;

  /*loop over window sizes until we find some*/
  nIn=0;
  fill=0.0;
  do{
    /*choose pixels to test*/
    indList=setFillList(i,j,window,&nTest,image);

    for(ii=0;ii<nTest;ii++){
      place=indList[ii];
      if(image->nIn[place]>0){
        /*find distance*/
        ti=place%image->nX;
        tj=place/image->nX;
        dx=(float)(ti-i);
        dy=(float)(tj-j);
        dist=sqrt(dx*dx+dy*dy);

        /*fill gap with distance and number weighted average*/
        weight=(float)image->nIn[place]/(dist*dist);
        if(mode==1)fill+=(float)image->image[place]*weight;
        else if(mode==2)fill+=image->jimlad[place]*weight;
        totWeight+=weight;
        nIn+=image->nIn[place];
      }
    }
    TIDY(indList);
    window++;
  }while((nIn<minIn)&&(window<maxWindow));

  if(mode==1){
    if(nIn>=minIn)*newChar=(unsigned char)(fill/totWeight);
    else          *newChar=255;
  }else if(mode==2){
    if(nIn>=minIn)*newFloat=fill/totWeight;
    else          *newFloat=0.0;
  }
}/*fillDTMhole*/


/*##################################################*/
/*choose indexes for gap filling*/

uint64_t *setFillList(int i,int j,int w,int *nTest,imageStruct *image)
{
  int ii=0,jj=0,tempN=0;
  uint64_t *indList=NULL;

  /*allocate maximum possible space*/
  tempN=(2*w+1)*(2*w+1)-((2*w-1)*(2*w-1));
  if(!(indList=(uint64_t *)calloc(tempN,sizeof(uint64_t)))){
    fprintf(stderr,"error fill index list allocation.\n");
    exit(1);
  }
  (*nTest)=0;

  /*x edges*/
  for(ii=i-w;ii<=(i+w);ii+=w*2){
    if((ii<0)||(ii>=image->nX))continue;
    for(jj=j-w;jj<=j+w;jj++){
      if((jj<0)||(jj>=image->nY))continue;
      indList[*nTest]=(uint64_t)jj*(uint64_t)image->nX+(uint64_t)ii;
      (*nTest)++;
    }
  }

  /*y edges*/
  for(jj=j-w;jj<=(j+w);jj+=w*2){
    if((jj<0)||(jj>=image->nY))continue;
    for(ii=(i-w)+1;ii<i+w;ii++){
      if((ii<0)||(ii>=image->nX))continue;
      indList[*nTest]=(uint64_t)jj*(uint64_t)image->nX+(uint64_t)ii;
      (*nTest)++;
    }
  }

  return(indList);
}/*setFillList*/


/*##################################################*/
/*write image to geotiff*/

void writeImage(control *dimage,imageStruct *image)
{
  if(dimage->charImage){
    drawTiff(dimage->outNamen,&(image->geoL[0]),&(image->geoI[0]),(double)dimage->res,image->image,image->nX,image->nY,255.0/(image->max-image->min),image->epsg);
  }else{
    drawTiffFlo(dimage->outNamen,&(image->geoL[0]),&(image->geoI[0]),(double)dimage->res,image->jimlad,image->nX,image->nY,1.0,image->epsg);
  }

  return;
}/*writeImage*/


/*##################################################*/
/*collate image*/

void collateImage(control *dimage,lasFile *las,imageStruct *image)
{
  int place=0,hBin=0;
  int xBin=0,yBin=0;
  uint32_t j=0;
  double x=0,y=0,z=0;
  void testVegVol(float *,float,uint64_t,control *,uint64_t *,unsigned char);

  if(las->epsg==0)las->epsg=dimage->epsg;

  /*check EPSG*/
  if(las->epsg!=image->epsg){
    fprintf(stderr,"EPSG mismatch %d %d\n",(int)image->epsg,las->epsg);
    exit(1);
  }

  /*loop over points*/
  for(j=0;j<las->nPoints;j++){
    readLasPoint(las,j);
    setCoords(&x,&y,&z,las);

    /*check against filters*/
    if(abs(las->scanAng)<=dimage->maxScanAng){
      xBin=(int)((x-image->minX)/(double)dimage->res+0.5);
      yBin=(int)((image->maxY-y)/(double)dimage->res+0.5);

      if((xBin>=0)&&(xBin<image->nX)&&(yBin>=0)&&(yBin<image->nY)){
        place=yBin*image->nX+xBin;
        if(dimage->drawInt)image->jimlad[place]+=(float)las->refl;
        else if(dimage->drawHeight){
          if((dimage->onlyGround)&&(las->classif==2))image->jimlad[place]+=(float)z;
          else if(dimage->onlyGround==0){
            if(image->minH[place]<-999.0)image->minH[place]=(float)z-image->hRange/2.0;
            hBin=(int)(((float)z-image->minH[place])/image->hRes+0.5);
            if((hBin>=0)&&(hBin<image->nHeight))image->heightStack[place][hBin]+=1.0;
            else fprintf(stderr,"Height bounds not quite wide enough. Point %f bounds %f %f\n",z,image->minH[place],image->minH[place]+image->hRange);
          }
        }else if(dimage->drawVegVol)testVegVol(&image->jimlad[place],(float)z,las->refl,dimage,&image->nIn[place],las->classif);
        if(dimage->findDens&&(las->retNumb==las->nRet))image->nFoot[place]++;
        if(dimage->drawCov&&(las->classif!=2))image->nCan[place]++;

        if(dimage->drawVegVol==0){
          if((dimage->onlyGround==0)||(las->classif==2))image->nIn[place]++;
        }
      }
    }/*filter check*/
  }/*point loop*/

  return;
}/*collateImage*/


/*##################################################*/
/*build up vegetation volume*/

void testVegVol(float *value,float z,uint64_t refl,control *dimage,uint64_t *nIn,unsigned char classif)
{

  if((refl<dimage->maxVint)&&(classif!=2)){
    if(z>(*value))(*value)=z;
    (*nIn)=1;
  }

  return;
}/*testVegVol*/


/*##################################################*/
/*finish off image*/

void finishImage(control *dimage,imageStruct *image)
{
  int i=0;
  int nContP=0,nContF=0;
  float meanPoint=0,meanFoot=0;
  float findTop(imageStruct *,int);
  void processHedges(imageStruct *,control *);


  /*normalise and find bounds*/
  meanPoint=meanFoot=0.0;
  nContP=nContF=0;
  if(!dimage->drawVegVol){
    for(i=image->nX*image->nY-1;i>=0;i--){
      if(image->nIn[i]>0){
        /*normalise elevationsor find top  if needed*/
        if(dimage->drawInt||dimage->drawHeight){
          /*normalise or find top*/
          if(dimage->drawInt||dimage->onlyGround)             image->jimlad[i]/=(float)image->nIn[i];
          else if(dimage->drawHeight&&(dimage->onlyGround==0))image->jimlad[i]=findTop(image,i);
          /*bounds for scaling geotiff*/
          if(image->jimlad[i]<image->min)image->min=image->jimlad[i];
          if(image->jimlad[i]>image->max)image->max=image->jimlad[i];
        }

        /*scale to canopy cover percent if needed*/
        if(dimage->drawCov)image->jimlad[i]=((float)image->nCan[i]/(float)image->nIn[i])*100.0;

        /*find density of needed*/
        if(dimage->findDens){
          if(image->nIn[i]>(uint64_t)image->maxPoint)image->maxPoint=image->nIn[i];
          if(image->nFoot[i]>(uint64_t)image->maxFoot)image->maxFoot=image->nFoot[i];
          if(image->nIn[i]>0){
            meanPoint+=(float)image->nIn[i];
            nContP++;
          }
          if(image->nFoot[i]>0){
            meanFoot+=(float)image->nFoot[i];
            nContF++;
          }
        }
      }else{  /*mark as missing data*/
        if(dimage->drawCov)image->jimlad[i]=255.0;
      }
    }
  }

  /*process for hedges only*/
  if(dimage->drawVegVol)processHedges(image,dimage);

  if((!dimage->drawDens)&&(dimage->gapFill==0))TIDY(image->nIn);
  if(dimage->findDens){
    if(nContF>0)meanFoot/=(float)nContF;
    if(nContP>0)meanPoint/=(float)nContP;
    fprintf(stdout,"Mean point density %f per m2\n",meanPoint/(dimage->res*dimage->res));
    fprintf(stdout,"Mean footprint density %f per m2\n",meanFoot/(dimage->res*dimage->res));
  }


  if(dimage->maxDN>0.0){
    if(image->max>dimage->maxDN)dimage->maxDN=image->max;
  }

  /*copy to uchar array*/
  if(dimage->charImage){
    if(dimage->drawInt||dimage->drawHeight){
      image->image=uchalloc((uint64_t)image->nX*(uint64_t)image->nY,"image",0);
      for(i=image->nX*image->nY-1;i>=0;i--){
        if(image->jimlad[i]<image->max)image->image[i]=(unsigned char)((image->jimlad[i]-image->min)*255.0/(image->max-image->min));
        else                           image->image[i]=255;
      }
    }
    if(dimage->drawCov){
      image->image=uchalloc((uint64_t)image->nX*(uint64_t)image->nY,"image",0);
      for(i=image->nX*image->nY-1;i>=0;i--)image->image[i]=(unsigned char)image->jimlad[i];
      image->max=255.0;
      image->min=0.0;
    }
    TIDY(image->jimlad);
  }
  return;
}/*finishImage*/


/*##################################################*/
/*find the tree top*/

float findTop(imageStruct *image,int i)
{
  int j=0;
  float top=0,tot=0;
  float *cumul=NULL;
  float x1=0,x2=0;
  float y1=0,y2=0;
  float m=0,c=0,thresh=0;

  /*find total and cumulative*/
  tot=0.0;
  cumul=falloc(image->nHeight,"cumulative wave",0);
  for(j=0;j<image->nHeight;j++){
    tot+=image->heightStack[i][j];
    cumul[j]=tot;
  }

  /*find crossing point and interpolate*/
  thresh=tot*image->hThresh;
  for(j=image->nHeight-2;j>=0;j--){
    if(cumul[j]<=thresh){
      /*points either side of crossing*/
      x1=(float)(j+1)*image->hRes+image->minH[i];
      x2=(float)j*image->hRes+image->minH[i];
      y1=cumul[j+1];
      y2=cumul[j];

      /*interpolate with a line*/
      m=(y2-y1)/(x2-x1);
      c=y1-m*x1;
      top=(thresh-c)/m;

      break;
    }
  }/*height bin loop*/

  TIDY(cumul);
  return(top);
}/*findTop*/


/*##################################################*/
/*processing to try and extract hedges*/

void processHedges(imageStruct *image,control *dimage)
{
  int i=0,j=0,place=0;
  int ii=0,jj=0,w=0;
  int nAbove=0,nBelow=0,nRight=0;
  float *newJim=NULL;

  newJim=falloc(image->nX*image->nY,"new jim",0);

  w=2;

  /*loop over image*/
  for(i=0;i<image->nX;i++){
    for(j=0;j<image->nY;j++){
      place=i+j*image->nX;
      nAbove=nBelow=nRight=0;
      if((image->jimlad[place]>=dimage->minVolH)&&(image->jimlad[place]<=dimage->maxVolH)&&(image->nIn[place]>0)){
        /*loop over focal area*/
        for(ii=i-w;ii<=i+w;ii++){
          if((ii<0)||(ii>=image->nX))continue;
          for(jj=j-w;jj<=j+w;jj++){
            if((jj<0)||(jj>=image->nY))continue;
            place=ii+jj*image->nX;
            if(image->nIn[place]>0){
              if(image->jimlad[place]>dimage->maxVolH)nAbove++;
              else if(image->jimlad[place]<dimage->minVolH)nBelow++;
              else nRight++;
            }
          }
        }

        place=i+j*image->nX;
        if((nRight>0)&&(nAbove<18)&&(nBelow>=0))newJim[place]=image->jimlad[place]*dimage->res*dimage->res;
        else                                   newJim[place]=0.0;
//newJim[place]=image->jimlad[place]*dimage->res*dimage->res;
      }else newJim[place]=0.0;
    }
  }


  for(i=image->nX*image->nY-1;i>=0;i--)image->jimlad[i]=newJim[i]*1.78;
  TIDY(newJim);

  return;
}/*processHedges*/


/*##################################################*/
/*update image bounds*/

void updateBounds(double *bounds,lasFile *las)
{
  if(las->minB[0]<bounds[0])bounds[0]=las->minB[0];
  if(las->minB[1]<bounds[1])bounds[1]=las->minB[1];
  if(las->maxB[0]>bounds[2])bounds[2]=las->maxB[0];
  if(las->maxB[1]>bounds[3])bounds[3]=las->maxB[1];
  return;
}/*updateBounds*/

/*##################################################*/
/*allocate image struture*/

imageStruct *allocateImage(control *dimage)
{
  int i=0,j=0;
  imageStruct *image=NULL;

  if(!(image=(imageStruct *)calloc(1,sizeof(imageStruct)))){
    fprintf(stderr,"error imageStruct allocation.\n");
    exit(1);
  }

  /*set image bounds*/
  image->minX=dimage->bounds[0];
  image->minY=dimage->bounds[1];
  image->maxX=dimage->bounds[2];
  image->maxY=dimage->bounds[3];

  /*size of image*/
  image->nX=(int)((image->maxX-image->minX)/(double)dimage->res+1);
  image->nY=(int)((image->maxY-image->minY)/(double)dimage->res+1);
  fprintf(stdout,"Image will be %d by %d\n",image->nX,image->nY);

  /*allocate data arrays*/
  if(dimage->drawInt||dimage->drawHeight||dimage->drawCov||dimage->drawVegVol)image->jimlad=falloc((uint64_t)image->nX*(uint64_t)image->nY,"jimlad",0);
  else                                                                        image->jimlad=NULL;
  if(dimage->drawCov){
    if(!(image->nCan=(uint64_t *)calloc(image->nX*image->nY,sizeof(uint64_t)))){
      fprintf(stderr,"error in canopy allocation\n");
      exit(1);
    }
  }
  if(!(image->nIn=(uint64_t *)calloc(image->nX*image->nY,sizeof(uint64_t)))){
    fprintf(stderr,"error in canopy allocation\n");
    exit(1);
  }
  if(dimage->findDens)image->nFoot=ialloc(image->nX*image->nY,"nFoot",0);
  else                image->nFoot=NULL;
  /*height image if needed*/
  if(dimage->drawHeight&&(dimage->onlyGround==0)){
    image->heightStack=fFalloc(image->nX*image->nY,"height stack",0);
    image->minH=falloc(image->nX*image->nY,"height array start",0);
    image->hRes=dimage->hRes;
    image->hThresh=dimage->hThresh;
    image->hRange=dimage->hRange;
    image->nHeight=(int)(image->hRange/image->hRes+1.0);
  }else image->heightStack=NULL;

  /*set all arrays to zero*/
  for(i=image->nX*image->nY-1;i>=0;i--){
    if(dimage->drawInt||dimage->drawHeight)image->jimlad[i]=0.0;
    if(dimage->findDens)image->nFoot[i]=0;
    if(dimage->drawCov)image->nCan[i]=0;
    if(dimage->drawHeight&&(dimage->onlyGround==0)){
      image->minH[i]=-9999.0;
      image->heightStack[i]=falloc(image->nHeight,"height stack",i+1);
      for(j=0;j<image->nHeight;j++)image->heightStack[i][j]=0.0;
    }
    image->nIn[i]=0;
  }
  image->min=1000000.0;
  image->max=-1000000.0;
  image->maxFoot=image->maxPoint=0;


  /*geolocation*/
  image->geoI[0]=image->geoI[1]=0;
  image->geoL[0]=image->minX+0.5*(double)dimage->res;
  image->geoL[1]=image->maxY-0.5*(double)dimage->res;
  image->epsg=dimage->epsg;

  return(image);
}/*allocateImage*/


/*##################################################*/
/*read command line*/

control *readCommands(int argc,char **argv)
{
  int i=0,j=0;
  control *dimage=NULL;
  char **readInList(int *,char *);

  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error contN allocation.\n");
    exit(1);
  }

  dimage->nFiles=1;
  dimage->inList=chChalloc(dimage->nFiles,"inList",0);
  dimage->inList[0]=challoc(200,"inList",0);
  strcpy(&(dimage->inList[0][0]),"/Users/stevenhancock/data/gedi/USDA_pilot/waveform/Lift02/WF_V13_-Riegl680i-HRZ-140701_130919_1-originalpoints.las");
  strcpy(dimage->outNamen,"teast.tif");
  dimage->drawInt=1;
  dimage->drawHeight=0;
  dimage->drawDens=0;
  dimage->findDens=0;
  dimage->drawCov=0;
  dimage->writeBounds=0;
  dimage->printNpoint=0;
  dimage->bFile=NULL;
  dimage->epsg=0;    /*leave bank*/
  dimage->pBuffSize=(uint64_t)200000000;
  dimage->charImage=1;
  dimage->onlyGround=0;
  dimage->gapFill=0;
  /*chm settings*/
  dimage->hThresh=0.99;  /*threshold to avoid noise on CHM*/
  dimage->hRes=0.25;     /*resolution to make psuedo-waveform for CHM*/
  dimage->hRange=100.0;  /*maximum expected extent of points*/
  /*bounds*/
  dimage->findBounds=1;
  dimage->bounds[0]=dimage->bounds[1]=1000000000.0;
  dimage->bounds[2]=dimage->bounds[3]=-1000000000.0;

  /*filters*/
  dimage->maxScanAng=120;

  dimage->res=100.0;
  dimage->maxDN=-1.0;


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
      }else if(!strncasecmp(argv[i],"-res",4)){
        checkArguments(1,i,argc,"-res");
        dimage->res=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxDN",6)){
        checkArguments(1,i,argc,"-maxDN");
        dimage->maxDN=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-epsg",5)){
        checkArguments(1,i,argc,"-epsg");
        dimage->epsg=(uint16_t)atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-noInt",6)){
        dimage->drawInt=dimage->drawCov=0;
      }else if(!strncasecmp(argv[i],"-height",7)){
        dimage->drawHeight=1;
        dimage->drawInt=dimage->drawCov=0;
      }else if(!strncasecmp(argv[i],"-DTM",4)){
        dimage->drawHeight=1;
        dimage->drawInt=dimage->drawCov=0;
        dimage->onlyGround=1;
        dimage->gapFill=1;
      }else if(!strncasecmp(argv[i],"-findDens",9)){
        dimage->findDens=1;
      }else if(!strncasecmp(argv[i],"-cover",6)){
        dimage->drawCov=1;
        dimage->drawInt=dimage->drawHeight=0;
      }else if(!strncasecmp(argv[i],"-hThresh",8)){
        checkArguments(1,i,argc,"-hThres");
        dimage->hThresh=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-hRes",5)){
        checkArguments(1,i,argc,"-hRes");
        dimage->hRes=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-hRange",7)){
        checkArguments(1,i,argc,"-hRange");
        dimage->hRange=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-writeBound",11)){
        checkArguments(1,i,argc,"-writeBound");
        dimage->writeBounds=1;
        strcpy(dimage->bNamen,argv[++i]);
        if((dimage->bFile=fopen(dimage->bNamen,"w"))==NULL){
          fprintf(stderr,"Error opening output file %s\n",dimage->bNamen);
          exit(1);
        }
      }else if(!strncasecmp(argv[i],"-pBuff",6)){
        checkArguments(1,i,argc,"-pBuff");
        dimage->pBuffSize=(uint64_t)(atof(argv[++i])*1000000000.0);
      }else if(!strncasecmp(argv[i],"-printNpoint",12)){
        dimage->printNpoint=1;
      }else if(!strncasecmp(argv[i],"-float",6)){
        dimage->charImage=0;
      }else if(!strncasecmp(argv[i],"-bounds",7)){
        checkArguments(4,i,argc,"-bounds");
        dimage->findBounds=0;
        for(j=0;j<4;j++)dimage->bounds[j]=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-vegVol",7)){
        dimage->drawVegVol=1;
        dimage->charImage=0;
        dimage->drawInt=dimage->drawCov=dimage->drawHeight=0;
      }else if(!strncasecmp(argv[i],"-maxVh",6)){
        checkArguments(1,i,argc,"-maxVh");
        dimage->maxVolH=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-minVh",6)){
        checkArguments(1,i,argc,"-minVh");
        dimage->minVolH=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxVint",8)){
        checkArguments(1,i,argc,"-maxVint");
        dimage->maxVint=(uint16_t)atoi(argv[++i]);   
      }else if(!strncasecmp(argv[i],"-maxScanAng",11)){
        checkArguments(1,i,argc,"-maxScanAng");
        dimage->maxScanAng=(int16_t)atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n#####\nProgram to create GEDI waveforms from ALS las files\n#####\n\n-input name;     lasfile input filename\n-output name;    output filename\n-inList list;    input file list for multiple files\n-res res;        image resolution, in metres\n-bounds minX minY maxX maxY;     user defined image bounds\n-float;          output as float\n-height;         draw height image\n-DTM;            make a bare Earth DEM\n-cover;          draw canopy cover map\n-noInt;          no image\n-findDens;       find point and footprint density\n-epsg n;         geolocation code if not read from file\n-hRange x;       range to expect points over for CHM\n-hThresh x;      percentile threshold to use for CHM in presence of noise\n-writeBound n;   write file bounds to a file\n-pBuff s;        point reading buffer size in Gbytes\n-printNpoint;    print number of points in each file\n\n-vegVol;     draw hedge volume\n-minVh h;\n-maxVh h;\n-maxVint dn;\n-maxScanAng ang;    maximum ALS scan angle to use, in degrees\nQuestions to svenhancock@gmail.com\n\n");
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
/*##################################################*/

