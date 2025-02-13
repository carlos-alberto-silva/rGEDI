#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "libLasRead.h"
#include "libLasProcess.h"


/*##########################*/
/*# Fits a single Gaussian #*/
/*# to a waveform  2017    #*/
/*# svenhancock@gmail.com  #*/
/*##########################*/


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



/*############################################*/
/*control structure*/

typedef struct{
  char inNamen[1000];  /*waveform filename*/
  int useCol;         /*column to read the waveform from*/
  denPar den;         /*denoising parameters*/
}control;

int main(int argc,char **argv)
{
  int nBins=0;
  control *dimage=NULL;
  control *readCommands(int,char **);
  float **data=NULL,sigma=0;
  float **readData(char *,int *,int);
  float findGround(float **,int,control *,float *);
  float ground=0;

  /*read command line*/
  dimage=readCommands(argc,argv);

  /*read data*/
  data=readData(dimage->inNamen,&nBins,dimage->useCol);

  /*find ground*/
  ground=findGround(data,nBins,dimage,&sigma);

  /*print results*/
  fprintf(stdout,"%.4f sigma %f\n",ground,sigma);

  /*tidy up*/
  TTIDY((void **)data,2);
  if(dimage){
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*############################################*/
/*fit single gaussian for ground*/

float findGround(float **data,int nBins,control *dimage,float *sigma)
{
  int nGauss=0;
  float ground=0;
  float *denoised=NULL,*fitWave=NULL,*gaussPar=NULL;
  float *fitSingleGauss(float *,float *,int,float,int *,float **);

  /*denoise*/
  denoised=processFloWave(data[1],nBins,&dimage->den,1.0);

  /*fit a Gaussian*/
  fitWave=fitSingleGauss(data[0],denoised,nBins,0.5,&nGauss,&gaussPar);
  TIDY(denoised);

  /*copy result and tidy up*/
  ground=gaussPar[0];
  *sigma=gaussPar[2];
  TIDY(gaussPar);
  TIDY(fitWave);
  return(ground);
}/*findGround*/


/*############################################*/
/*read data*/

float **readData(char *inNamen,int *nBins,int useCol)
{
  int i=0;
  float **data=NULL;
  char line[10000],temp1[50],temp2[50],temp3[50];
  char temp4[50],temp5[50];
  FILE *ipoo=NULL;

  /*open file*/
  if((ipoo=fopen(inNamen,"r"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",inNamen);
    exit(1);
  }

  /*first count the number of lines*/
  (*nBins)=0;
  while(fgets(line,10000,ipoo)!=NULL){
    if(strncasecmp(line,"#",1))(*nBins)++;
  }/*count number of wavefoms*/

  /*rewind to start of file*/
  if(fseek(ipoo,(long)0,SEEK_SET)){
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  /*allocate*/
  data=fFalloc(2,"data",0);
  for(i=0;i<2;i++)data[i]=falloc((uint64_t)(*nBins),"data",i+1);

  i=0;
  while(fgets(line,10000,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      if(useCol==5){
        if(sscanf(line,"%s %s %s %s %s",temp1,temp2,temp3,temp4,temp5)==5){
          data[0][i]=atof(temp1);
          data[1][i]=atof(temp5);
        }
      }else if(useCol==2){
        if(sscanf(line,"%s %s",temp1,temp2)==2){
          data[0][i]=atof(temp1);
          data[1][i]=atof(temp2);
        }
      }else{
        fprintf(stderr,"Are you sure you want to read column %d?\n",useCol);
        exit(1);
      }
      i++;
    }
  }/*read data*/

  /*close file*/
  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return(data);
}/*readData*/


/*############################################*/
/*read command line*/

control *readCommands(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;
  void setDenoiseDefault(denPar *);


  /*allocate structures*/
  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }

  dimage->useCol=2;

  /*denoising*/
  setDenoiseDefault(&dimage->den);
  dimage->den.statsLen=20.0;
  dimage->den.varNoise=1;
  dimage->den.threshScale=3.0;
  dimage->den.noiseTrack=1;
  dimage->den.sWidth=0.0;


  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        strcpy(dimage->inNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-col",4)){
        checkArguments(1,i,argc,"-col");
        dimage->useCol=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n#####\nProgram to fit single Gaussian to wave\n#####\n\n-input name;   input filaname\n-col n;      column to read waveform from\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry gediRat -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }
 
  return(dimage);
}/*readCommands*/

/*the end*/
/*############################################*/

