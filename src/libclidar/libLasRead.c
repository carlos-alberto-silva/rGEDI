#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "libLasRead.h" 


/*######################*/
/*# A library for      #*/
/*# handling las files #*/
/*# S Hancock, 2015    #*/
/*######################*/


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


int nFilesOpen;    /*record the number of files open*/
#define MAX_OPEN 250


/*##################################################################*/
/*read the jheader*/

lasFile *readLasHead(char *namen,uint64_t pBuffSize)
{
  int i=0;
  int offset=0;         /*to step around header byte arrays*/
  int tempLen=0;
  uint64_t temp64=0;
  char *pubHead=NULL;   /*public header*/
  lasFile *las=NULL;

  if(!(las=(lasFile *)calloc(1,sizeof(lasFile)))){
    fprintf(stderr,"error lasFile allocation.\n");
    exit(1);
  }

  /*open file*/
  if((las->ipoo=fopen(namen,"rb"))==NULL){
    fprintf(stderr,"Error opening input file \"%s\"\n",namen);
    exit(1);
  }
  nFilesOpen++;
  strcpy(las->namen,namen);


  /*determine type*/
  tempLen=96;
  pubHead=challoc((uint64_t)tempLen,"pubHead",0);
  if(fread(&(pubHead[0]),sizeof(char),tempLen,las->ipoo)!=tempLen){
    fprintf(stderr,"error reading data from %s\n",namen);
    exit(1);
  }
  if(strncasecmp(pubHead,"LASF",4)){  /*check is a las file*/
    fprintf(stderr,"Incorrect filetype for %s\n",namen);
    exit(1);
  }/*check it is a LASF file*/

  offset=24;  /*version major*/
  memcpy(&las->vMajor,&pubHead[offset],1);
  offset=25;  /*version minor*/
  memcpy(&las->vMinor,&pubHead[offset],1);
  /*fprintf(stdout,"Version %d.%d\n",las->vMajor,las->vMinor);*/
  offset=94;  /*header size*/
  memcpy(&las->headSize,&pubHead[offset],2);
  TIDY(pubHead);

  if((las->vMajor!=1)||(las->vMinor>4)){
    fprintf(stderr,"Version too new for this program\n");
    fprintf(stderr,"Version %d.%d\n",las->vMajor,las->vMinor);
    exit(1);
  }

  /*rewind*/
  if(fseek(las->ipoo,(long)0,SEEK_SET)){ /*rewind to start of file*/
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  /*set header length depending on version*/
  pubHead=challoc((uint64_t)las->headSize,"pubHead",0);

  if(fread(&(pubHead[0]),sizeof(char),las->headSize,las->ipoo)!=las->headSize){
    fprintf(stderr,"error reading data from %s\n",namen);
    exit(1);
  }

  /*copy memory bits over*/
  offset=90;   /*day of year*/
  memcpy(&las->doy,&pubHead[offset],2);
  offset=92;  /*year*/
  memcpy(&las->year,&pubHead[offset],2);
  offset=94;  /*header size*/
  memcpy(&las->headSize,&pubHead[offset],2);
  offset=96;  /*offset to point data*/
  memcpy(&las->offsetToP,&pubHead[offset],4);
  offset=100;  /*number of variable records*/
  memcpy(&las->nVarRec,&pubHead[offset],4);
  offset=104;  /*point data format ID*/
  memcpy(&las->pointFormat,&pubHead[offset],1);
  offset=105;  /*point data format ID*/
  memcpy(&las->pRecLen,&pubHead[offset],1);
  offset=107;  /*number of point records*/
  memcpy(&las->nPoints,&pubHead[offset],4);
  offset=111;  /*number of points by return*/
  for(i=0;i<7;i++)memcpy(&las->nPbyRet[i],&pubHead[offset+4*i],4);

  /*we have lost 8 bytes somewhere here, because the above is 5 byte, fool*/
  offset=139-8;
  for(i=0;i<3;i++){  /*point scaling factors*/
    memcpy(&las->posScale[i],&pubHead[offset],8);
    offset+=8;
  }/*read scaling factors*/
  offset=163-8;
  for(i=0;i<3;i++){  /*point offsets*/
    memcpy(&las->posOffset[i],&pubHead[offset+i*8],8);
  }/*read offsets*/

  offset=187-8;
  for(i=0;i<3;i++){  /*point bounds*/
    memcpy(&las->maxB[i],&pubHead[offset+i*2*8],8);
    memcpy(&las->minB[i],&pubHead[offset+i*2*8+8],8);
    //fprintf(stdout,"Bounds %d %f %f\n",i,las->minB[i],las->maxB[i]);
  }/*read bounds*/

  if((las->vMajor==1)&&(las->vMinor>2)){
    offset=235-8;   /*Waveform packet start*/
    memcpy(&(las->waveStart),&(pubHead[offset]),sizeof(uint64_t));
  }

  if((las->vMajor==1)&&(las->vMinor==4)){  /*8 byte point record*/
    offset=255-8;
    memcpy(&temp64,&(pubHead[offset]),sizeof(uint64_t));
    if((las->nPoints==0)&&(temp64>0)){
      if(temp64>=4294967296){
        fprintf(stderr,"Currently the library cannot handle more than 4294967296 points per file\n");
        exit(1);
      }
      las->nPoints=(uint32_t)temp64;
    }
  }

  TIDY(pubHead);

  /*make sure we don't open too many files*/
  if(nFilesOpen>=MAX_OPEN){
    if(las->ipoo){
      fclose(las->ipoo);
      las->ipoo=NULL;
    }
    nFilesOpen--;
  }

  /*set up point reading buffer size*/
  if(pBuffSize>0)las->maxBuffLen=(uint32_t)(pBuffSize/(uint64_t)las->pRecLen);
  else           las->maxBuffLen=1;
  las->maxBuffSize=(uint64_t)las->maxBuffLen*(uint64_t)las->pRecLen;
  las->buffStart=0;    /*buffer start in number of points*/
  las->buffLen=0;      /*buffer length in number of points*/
  las->buffByte=0;     /*buffer length in number of bytes*/

  return(las);
}/*readLasHead*/


/*##############################################*/
/*read geolocation*/

void readLasGeo(lasFile *las)
{
  int i=0,offset=0;
  int headLen=0,varLen=0;
  uint16_t recID=0,reserve=0;
  char **varHead=NULL,*varData=NULL;
  char namen[16];

  /*open file*/
  if((las->ipoo=fopen(las->namen,"rb"))==NULL){
    fprintf(stderr,"Error opening input file \"%s\"\n",las->namen);
    exit(1);
  }
  nFilesOpen++;

  if(fseeko(las->ipoo,(long)las->headSize,SEEK_SET)){
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  /*read the variable header - needed for waveform information*/
  headLen=54;
  las->epsg=0;
  varHead=chChalloc(las->nVarRec,"variable headers",0);
  for(i=0;i<las->nVarRec;i++){
    varHead[i]=challoc((uint64_t)headLen,"variable headers",i+1);
    if(fread(&(varHead[i][0]),sizeof(char),headLen,las->ipoo)!=headLen){
      fprintf(stderr,"Error reading variable header\n");
      exit(1);
    }
    offset=0;
    memcpy(&reserve,&varHead[i][offset],2);
    offset=2;
    memcpy(&(namen[0]),&varHead[i][offset],16);
    offset=18;
    memcpy(&recID,&varHead[i][offset],2);
    offset=20;
    memcpy(&varLen,&varHead[i][offset],2);

    /*read variable part*/
    varData=challoc((uint64_t)varLen,"variable data",0);
    if(fread(&(varData[0]),sizeof(char),varLen,las->ipoo)!=varLen){
      fprintf(stderr,"Error reading variable header\n");
      exit(1);
    }
    if(!strncasecmp(namen,"LASF_Projection",16)){  /*geo projection*/
      if(recID==34735){
fprintf(stdout,"vrLen %d\n",(int)varLen);
        offset=102;
        memcpy(&las->epsg,&varData[offset],2);
        fprintf(stdout,"EPSG %d\n",las->epsg);
      }
    }else{/*geo projection*/
      fprintf(stdout,"%s\n",namen);
    }
    TIDY(varData);
  }/*variable header loop*/

  if(las->epsg==0){
    fprintf(stderr,"No geolocation information. Setting default\n");
    las->epsg=32619;
  }
  TTIDY((void **)varHead,las->nVarRec);

  /*close of too many files*/
  if(nFilesOpen>=MAX_OPEN){
    if(las->ipoo){
      fclose(las->ipoo);
      las->ipoo=NULL;
    }
    nFilesOpen--;
  }
  return;
}/*readLasGeo*/


/*##############################################*/
/*read a single point*/

void readLasPoint(lasFile *las,uint32_t j)
{
  uint64_t offset=0;
  uint64_t offTo=0;
  char tempByte=0;

  /*if not already, open file*/
  if(las->ipoo==NULL){
    if((las->ipoo=fopen(las->namen,"rb"))==NULL){
      fprintf(stderr,"Error opening input file \"%s\"\n",las->namen);
      exit(1);
    }
    nFilesOpen++;
  }


  /*do we need to read new buffer*/
  if((las->pointBuff==NULL)||(j>=(las->buffStart+las->buffLen))||(j<las->buffStart)){
    TIDY(las->pointBuff);
    las->buffStart=j;
    if((j+(uint32_t)las->maxBuffLen)<las->nPoints)las->buffByte=las->maxBuffSize;
    else                                          las->buffByte=(uint64_t)(las->nPoints%las->maxBuffLen)*(uint64_t)las->pRecLen;
    las->buffLen=(uint32_t)(las->buffByte/(uint64_t)las->pRecLen);

    if(!(las->pointBuff=(char *)calloc(las->buffByte,sizeof(char)))){
      fprintf(stderr,"error point buffer allocation.\n");
      exit(1);
    } 

    if(fseeko(las->ipoo,(long)las->offsetToP+(long)((uint64_t)las->buffStart*(uint64_t)las->pRecLen),SEEK_SET)){
      fprintf(stderr,"fseek error\n");
      exit(1);
    }
    if(fread(&(las->pointBuff[0]),sizeof(char),las->buffByte,las->ipoo)!=las->buffByte){
      fprintf(stderr,"Error reading point data, size %d\n",(int)las->buffByte);
      exit(1);
    }
  }/*read buffer*/


  if(nFilesOpen>=MAX_OPEN){
    if(las->ipoo){
      fclose(las->ipoo);
      las->ipoo=NULL;
    }
    nFilesOpen--;
  }

  /*offset to start of this point record*/
  offTo=(uint64_t)(j-las->buffStart)*(uint64_t)las->pRecLen;

  /*point format 3 and 4*/
  offset=offTo;
  memcpy(&las->x,&las->pointBuff[offset],4);
  offset+=4;
  memcpy(&las->y,&las->pointBuff[offset],4);
  offset+=4;
  memcpy(&las->z,&las->pointBuff[offset],4);
  offset+=4;
  memcpy(&las->refl,&las->pointBuff[offset],2);
  offset+=2;

  if(las->pointFormat<6){
    memcpy(&las->field,&las->pointBuff[offset],1);
    offset+=1;
    las->retNumb=las->field.retNumb;
    las->nRet=las->field.nRet;
    las->sDir=las->field.sDir;
    las->edge=las->field.edge;
  }else{
    memcpy(&las->newField,&las->pointBuff[offset],2);
    offset+=2;
    las->retNumb=las->newField.retNumb;
    las->nRet=las->newField.nRet;
    las->classF=las->newField.classF;
    las->sChan=las->newField.sChan;
    las->sDir=las->newField.sDir;
    las->edge=las->newField.edge;
  }

  memcpy(&las->classif,&las->pointBuff[offset],1);
  offset+=1;

  if(las->pointFormat<6){
    memcpy(&tempByte,&las->pointBuff[offset],1);
    las->scanAng=(int16_t)tempByte;
    offset+=1;
    memcpy(&las->userData,&las->pointBuff[offset],1);
    offset+=1;
  }else{
    memcpy(&las->userData,&las->pointBuff[offset],1);
    offset+=1;
    memcpy(&las->scanAng,&las->pointBuff[offset],2);
    offset+=2;
  }

  memcpy(&las->psID,&las->pointBuff[offset],2);
  offset+=2;
  /*GPS time*/
  memcpy(&las->gpsTime,&las->pointBuff[offset],8);
  offset+=8;

  if((las->pointFormat==3)||(las->pointFormat==10)||(las->pointFormat==8)||(las->pointFormat==7)||(las->pointFormat==5)){  /*there is RGB*/
    for(j=0;j<3;j++){
      memcpy(&(las->RGB[j]),&las->pointBuff[offset],2);
      offset+=2;
    }
  }/*there is RGB*/

  if((las->pointFormat==4)||(las->pointFormat==5)||(las->pointFormat==9)||(las->pointFormat==10)){   /*full waveform data*/
    memcpy(&las->packetDes,&las->pointBuff[offset],1);
    offset+=1;
    memcpy(&las->waveMap,&las->pointBuff[offset],8);
    offset+=8;
    memcpy(&las->waveLen,&las->pointBuff[offset],4);
    offset+=4;
    memcpy(&las->time,&las->pointBuff[offset],4);
    offset+=4;
    memcpy(&las->grad[0],&las->pointBuff[offset],4);
    offset+=4;
    memcpy(&las->grad[1],&las->pointBuff[offset],4);
    offset+=4;
    memcpy(&las->grad[2],&las->pointBuff[offset],4);
    offset+=4;
  }else{
    las->packetDes=0;
    las->grad[0]=las->grad[1]=las->grad[2]=0.0;  /*leave the grad bits blank*/
  }

  return;
}/*readLasPoint*/


/*#########################################################################*/
/*read a waveform*/

unsigned char *readLasWave(uint64_t waveMap,int32_t waveLen,FILE *ipoo,uint64_t waveStart)
{
  unsigned char *wave=NULL;

  wave=uchalloc((uint64_t)waveLen,"waveform",(int)waveMap);
  if(fseeko(ipoo,(off_t)((uint64_t)waveMap+(uint64_t)waveStart),SEEK_SET)){
    printf("Error seeking through las file\n");
    exit(1);
  }
  if(fread(&(wave[0]),sizeof(unsigned char),waveLen,ipoo)!=waveLen){
    fprintf(stderr,"Error reading waveform at %ld for %ld\n",(long int)((uint64_t)waveMap+(uint64_t)waveStart),(long int)waveLen);
    exit(1);
  }
  return(wave);
}/*readLasWaveform*/


/*##############################################*/
/*read input file list*/

char **readInList(int *nFiles,char *inList)
{
  int i=0;
  char line[400];
  char **namen=NULL;
  FILE *ipoo=NULL;

  if((ipoo=fopen(inList,"r"))==NULL){
    fprintf(stderr,"Error opening input file list \"%s\"\n",inList);
    exit(1);
  }
  i=0;   /*count up files*/
  while(fgets(line,399,ipoo)!=NULL){
    if(strncasecmp(line,"#",1))i++;
  }/*count up files*/

  if(fseek(ipoo,(long)0,SEEK_SET)){ /*rewind to start of file*/
    fprintf(stderr,"fseek error\n");
    exit(1);
  }
  *nFiles=i;
  namen=chChalloc(*nFiles,"file names",0);
  i=0;
  while(fgets(line,399,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      namen[i]=challoc((uint64_t)strlen(line)+1,"file names",i+1);
      sscanf(line,"%s",namen[i]);
      i++;
    }
  }

  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return(namen);
}/*readInList*/



/*######################################*/
/*tidy las file*/

lasFile *tidyLasFile(lasFile *las)
{
  if(las){
    if(las->ipoo){
      fclose(las->ipoo);
      las->ipoo=NULL;
    }
    TIDY(las->wave);
    TIDY(las->pointBuff);
    TIDY(las);
  }

  return(las);
}/*tidyLasFile*/


/*#########################################################################*/
/*set coordinates*/

void setCoords(double *lon,double *lat,double *height,lasFile *lasIn)
{
  *lon=(double)lasIn->x*lasIn->posScale[0]+lasIn->posOffset[0];
  *lat=(double)lasIn->y*lasIn->posScale[1]+lasIn->posOffset[1];
  *height=(double)lasIn->z*lasIn->posScale[2]+lasIn->posOffset[2];

  return;
}/*setCoords*/


/*#########################################################################*/
/*determine coordinate of a bin*/

void binPosition(double *x,double *y,double *z,int bin,double xCent,double yCent,double zCent,float time,float *grad)
{
  double r=0;     /*distance along beam froma anchor point*/

  r=(double)bin*1000.0-(double)time;

  *x=xCent+(double)grad[0]*r;
  *y=yCent+(double)grad[1]*r;
  *z=zCent+(double)grad[2]*r;

  return;
}/*binPosition*/


/*############################################*/
/*check there is a waveform and just one per beam*/

char checkOneWave(lasFile *lasIn)
{
  if((lasIn->packetDes)&&(lasIn->nRet==lasIn->retNumb)&&(lasIn->waveLen>0))return(1);
  else                                                                                 return(0);
}/*checkOneWave*/


/*############################################*/
/*check file bounds*/

char checkFileBounds(lasFile *lasIn,double minX,double maxX,double minY,double maxY)
{
  if((lasIn->minB[0]<=maxX)&&(lasIn->minB[1]<=maxY)&&\
     (lasIn->maxB[0]>=minX)&&(lasIn->maxB[1]>=minY))return(1);
  else                                              return(0);
}/*checkFileBounds*/


/*############################################*/
/*read GBIC table for Leica instruments*/

void readGBIC(char appGBIC,char balFlights,lasFile **lasIn,listLas *lasList)
{
  int i=0,j=0;
  char line[200];
  char temp1[100];
  char temp2[100];
  FILE *ipoo=NULL;

  for(i=0;i<lasList->nFiles;i++){
    lasIn[i]->gbLen=257;
    lasIn[i]->gbic=falloc((uint64_t)lasIn[i]->gbLen,"GBIC",i+1);
    if(balFlights==0)lasIn[i]->flightBal=1.0;
    else             lasIn[i]->flightBal=lasList->flightBal[i];
    if(appGBIC){  /*open file if needed*/
      if((ipoo=fopen(lasList->gbicNamen[i],"r"))==NULL){
        fprintf(stderr,"Error opening GBIC file \"%s\"\n",lasList->gbicNamen[i]);
        exit(1);
      }
      while(fgets(line,200,ipoo)!=NULL){
        if(strncasecmp(line,"#",1)){
          sscanf(line,"%s %s",temp1,temp2);
          j=atoi(temp1);
          if((j>=0)&&(j<lasIn[i]->gbLen))lasIn[i]->gbic[j]=atof(temp2)/lasIn[i]->flightBal;
        }
      }/*read names*/
      if(ipoo){
        fclose(ipoo);
        ipoo=NULL;
      }
    }else{  /*fill with ones*/
      for(j=0;j<lasIn[i]->gbLen;j++)lasIn[i]->gbic[j]=1.0;
    }
  }/*file loop*/
  return;
}/*readGBIC*/


/*#########################################################################*/
/*allocate an array of lasfiles*/

lasFile **lfalloc(int numb)
{
  lasFile **lasIn=NULL;

  if(!(lasIn=(lasFile **)calloc(numb,sizeof(lasFile *)))){
    fprintf(stderr,"error in las file structure\n");
    fprintf(stderr,"Allocating %d\n",numb);
    exit(1);
  }

  return(lasIn);
}/*lfalloc*/


/*#########################################################################*/
/*read list of input files names*/

listLas *readLasList(char *namen)
{
  int i=0;
  listLas *lasList=NULL;
  char line[400],temp3[200];
  char temp1[200],temp2[200];
  FILE *ipoo=NULL;

  if(!(lasList=(listLas *)calloc(1,sizeof(listLas)))){
    fprintf(stderr,"error in input filename structure.\n");
    exit(1);
  }

  if((ipoo=fopen(namen,"r"))==NULL){
    fprintf(stderr,"Error opening list file \"%s\"\n",namen);
    exit(1);
  }


  /*count number of files*/
  lasList->nFiles=0;
  while(fgets(line,400,ipoo)!=NULL)if(strncasecmp(line,"#",1))lasList->nFiles++;
  if(lasList->nFiles==0){
    fprintf(stderr,"No files in %s\n",namen);
    exit(1);
  }
  lasList->nameList=chChalloc(lasList->nFiles,"name list",0);
  lasList->gbicNamen=chChalloc(lasList->nFiles,"GBIC file list",0);
  /*TLS only parameters, left blank for ALS*/
  lasList->scanCent=dDalloc(lasList->nFiles,"scan centre",0);
  lasList->align=dDalloc(lasList->nFiles,"scan alignment",0);
  lasList->flightBal=falloc((uint64_t)lasList->nFiles,"flight balance",0);
  for(i=0;i<lasList->nFiles;i++){
    lasList->scanCent[i]=dalloc(3,"scan centre",i+1);
    lasList->align[i]=dalloc(3,"alignment",i+1);
    lasList->flightBal[i]=1.0;
  }/*TLS only parameters*/


  if(fseek(ipoo,(long)0,SEEK_SET)){ /*rewind to start of file*/
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  /*read*/
  i=0;
  while(fgets(line,400,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      /*really I should read the number of spaces here*/

      if(sscanf(line,"%s %s %s",temp1,temp2,temp3)==3){ /*read las, GBIC file and flight intensity balance*/
        lasList->nameList[i]=challoc((uint64_t)strlen(temp1)+1,"name list",i+1);
        strcpy(&(lasList->nameList[i][0]),temp1);
        lasList->gbicNamen[i]=challoc((uint64_t)strlen(temp2)+1,"GBIC file list",i+1);
        strcpy(&(lasList->gbicNamen[i][0]),temp2);
        lasList->flightBal[i]=atof(temp3);
      }else if(sscanf(line,"%s %s",temp1,temp2)==2){ /*read las and GBIC file*/
        lasList->nameList[i]=challoc((uint64_t)strlen(temp1)+1,"name list",i+1);
        strcpy(&(lasList->nameList[i][0]),temp1);
        lasList->gbicNamen[i]=challoc((uint64_t)strlen(temp2)+1,"GBIC file list",i+1);
        strcpy(&(lasList->gbicNamen[i][0]),temp2);
      }else{
        lasList->nameList[i]=challoc((uint64_t)strlen(line)+1,"name list",i+1);
        sscanf(line,"%s",lasList->nameList[i]);
      }
      i++;
    }
  }/*read names*/

  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }

  return(lasList);
}/*readLasList*/


/*#########################################################################*/
/*tidy lasList*/

void tidyListLas(listLas *lasList)
{
  if(lasList){/*input list*/
    TTIDY((void **)lasList->nameList,lasList->nFiles);
    TTIDY((void **)lasList->scanCent,lasList->nFiles);
    TTIDY((void **)lasList->align,lasList->nFiles);
    TTIDY((void **)lasList->gbicNamen,lasList->nFiles);
    TIDY(lasList->flightBal);
    TIDY(lasList);
  }/*input list*/
  return;
}/*tidyListLas*/


/*#########################################################################*/
/*read a file of coordinates*/

double **readCoordList(int nFiles,char **nameList,char *coordNamen)
{
  int i=0;
  double **coords=NULL;
  char line[1000],temp1[300],temp2[100];
  char temp3[100],temp4[100];
  FILE *ipoo=NULL;

  if((ipoo=fopen(coordNamen,"r"))==NULL){
    fprintf(stderr,"Error opening coords input file \"%s\"\n",coordNamen);
    exit(1);
  }
  coords=dDalloc(nFiles,"coords",0);
  for(i=0;i<nFiles;i++)coords[i]=dalloc(3,"coords",i+1);

  while(fgets(line,200,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      if(sscanf(line,"%s %s %s %s",temp1,temp2,temp3,temp4)==4){
        for(i=0;i<nFiles;i++){
          if(!strncasecmp(temp1,nameList[i],strlen(nameList[i]))){
            coords[i][0]=atof(temp2);
            coords[i][1]=atof(temp3);
            coords[i][2]=atof(temp4);
            break;
          }/*if file names match*/
        }/*file list loop*/
      }/*read data*/
    }/*comment check*/
  }/*line loop*/

  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return(coords);
}/*readCentres*/


/*#########################################################################*/
/*tidy up point cloud structure*/

pCloudStruct *tidyPointCloud(pCloudStruct *data)
{
  if(data){
    TIDY(data->x);
    TIDY(data->y);
    TIDY(data->z);
    TIDY(data->refl);
    TIDY(data->class);
    TIDY(data->nRet);
    TIDY(data->retNumb);
    TIDY(data->packetDes);
    TTIDY((void **)data->grad,3);
    data->grad=NULL;
    TIDY(data->time);
    TIDY(data->waveMap);
    TIDY(data->waveLen);
    TIDY(data->scanAng);
    if(data->ipoo){
      fclose(data->ipoo);
      data->ipoo=NULL;
    }
  }
  TIDY(data->gap);
  TIDY(data->range);
  TIDY(data);
  return(data);
}/*tidyPointCloud*/


/*the end*/
/*######################################*/

