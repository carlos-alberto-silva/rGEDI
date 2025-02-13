#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "libLasProcess.h"


/*##############################*/
/*# Tests adding noise to GEDI #*/
/*# waveforms                  #*/
/*# 2016 svenhancock@gmail.com #*/
/*##############################*/

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


/*tolerances*/
#define TOL 0.00001
#define XRES 0.00005
#define YTOL 0.0000001
#define MINERR 0.0000001


/*element reflectance*/
float rhoG;
float rhoC;


/*##############################*/
/*control structure*/

typedef struct{
  /*input/output*/
  int nFiles;   /*number of waveforms*/
  char **inList;
  char listNamen[1000];
  char outRoot[200];

  /*options*/
  char ground;
  char useInt;
  char linkNoise;  /*relate noise to link margin*/
  char physNoise;  /*physical noise*/
  char writeWave;  /*write waveforms switch*/

  /*noise parameters*/
  char bitRate;   /*digitiser bit rate*/
  float rho;      /*surface reflectance*/
  float maxDN;    /*maximum DN we need to digitise*/
  float offset;   /*waveform DN offset*/
  float lPower;   /*aser power in Joules*/
  float tran;     /*atmospheric transmission*/
  float solarI;   /*solar intensity at 1064nm*/
  float linkM;    /*link margin*/
  float linkCov;  /*cover at which link margin is true*/
  float trueSig;  /*true noise sigma. 5 Night, 10 Day*/
  float linkSig;  /*width of link margin noise*/

  /*laser parmaeters*/
  float pSigma;   /*pulse width*/
  float fSigma;   /*footprint width*/
  float deSig;    /*detector response*/
  float res;      /*range resolution*/
}control;


/*###########################################################*/
/*data structure*/

typedef struct{
  int nBins;     /*number of wavefor bins*/
  float *wave;   /*original waveform*/
  float *ground; /*ground wave if we are to read it*/
  float *noised; /*noised waveform, if needed*/
  double *z;     /*wave bin elevation*/
  double gElev;  /*mean ground elevation*/
  double tElev;  /*top elevation*/
  float gStdev;  /*measure of ground width*/
  float cov;     /*ALS canopy cover*/
}dataStruct;


/*##############################*/
/*main*/

int main(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;
  control *readCommands(int,char **);
  dataStruct *data=NULL;
  dataStruct *readData(char *,control *);
  void addPhysNoise(dataStruct *,control *);
  void addLinkNoise(dataStruct *,control *);
  void writeWaves(dataStruct *,control *,int);
  float setNoiseSigma(float,float,float,float);
  void addNoise(dataStruct *,control *);


  /*read command Line*/
  dimage=readCommands(argc,argv);

  /*set up link noise if needed*/
  dimage->linkSig=setNoiseSigma(dimage->linkM,dimage->linkCov,dimage->pSigma,dimage->fSigma);

  /*loop over files*/
  for(i=0;i<dimage->nFiles;i++){
    //fprintf(stdout,"Wave %d of %d\n",i+1,dimage->nFiles);

    /*read waveform*/
    data=readData(dimage->inList[i],dimage);

    /*add noise*/
    addNoise(data,dimage);

    /*write results*/
    if(dimage->writeWave)writeWaves(data,dimage,i);


    /*tidy up*/
    if(data){
      TIDY(data->noised);
      TIDY(data->ground);
      TIDY(data->wave);
      TIDY(data->z);
      TIDY(data);
    }
  }/*file loop*/


  if(dimage){
    TIDY(dimage);
  }
  return(0);
}/*main*/



/*####################################################*/
/*add noise to waveform*/

void addNoise(dataStruct *data,control *dimage)
{
  int i=0;
  float tot=0.0;
  float *tempNoise=NULL;
  float GaussNoise();
  float *smooNoise=NULL;
  float *digitiseWave(float *,int,char,float,float);
  float reflScale=0;
  void scaleNoiseDN(float *,int,float,float,float);

  if(dimage->linkNoise){   /*link margin basED NOISE*/
    /*Gaussian noise*/
    tempNoise=falloc((uint64_t)data->nBins,"temp noised",0);
    tot=0.0;
    for(i=0;i<data->nBins;i++)tot+=data->wave[i]*dimage->res;
    reflScale=(data->cov*rhoC+(1.0-data->cov)*rhoG)*tot/(dimage->linkCov*rhoC+(1.0-dimage->linkCov)*rhoG);
    for(i=0;i<data->nBins;i++)tempNoise[i]=dimage->linkSig*GaussNoise()*reflScale;
    /*smooth noise by detector response*/
    smooNoise=smooth(dimage->deSig,data->nBins,tempNoise,dimage->res);
    for(i=0;i<data->nBins;i++)tempNoise[i]=data->wave[i]+smooNoise[i];
    TIDY(smooNoise);
    /*scale to match sigma*/
    scaleNoiseDN(tempNoise,data->nBins,dimage->linkSig*reflScale,dimage->trueSig,dimage->offset);
    /*digitise*/
    data->noised=digitiseWave(tempNoise,data->nBins,dimage->bitRate,dimage->maxDN,tot);
    TIDY(tempNoise);
  }else{
    fprintf(stderr,"Not ready for that yet\n");
  }
  return;
}/*addNoise*/


/*####################################################*/
/*scale noise to match Bryan's numbers*/

void scaleNoiseDN(float *noised,int nBins,float noiseSig,float trueSig,float offset)
{
  int i=0;
  float sigScale=0;

  sigScale=trueSig/noiseSig;

  for(i=0;i<nBins;i++)noised[i]=noised[i]*sigScale+offset;

  return;
}/*scaleNoiseDN*/


/*####################################################*/
/*digitse*/

float *digitiseWave(float *wave,int nBins,char bitRate,float maxDN,float tot)
{
  int i=0;
  int nDN=0;
  float *sampled=NULL;
  float resDN=0;

  sampled=falloc((uint64_t)nBins,"sampled wave",0);

  /*number of bins*/
  nDN=1;
  for(i=0;i<bitRate;i++)nDN*=2;

  resDN=maxDN/(float)nDN;
  for(i=0;i<nBins;i++)sampled[i]=floor(wave[i]/resDN)*resDN;

  return(sampled);
}/*digitiseWave*/


/*####################################################*/
/*noise from a Gaussian*/

float GaussNoise()
{
  float noise=0,max=0;
  float x1=0,x2=0,w=0;

  if(RAND_MAX>0)max=(float)RAND_MAX;
  else          max=-1.0*(float)RAND_MAX;

  /*Box approximation to Gaussian random number*/
  w=0.0;
  do{
    x1=2.0*(float)rand()/max-1.0;
    x2=2.0*(float)rand()/max-1.0;
    w=x1*x1+x2*x2;
  }while(w>=1.0);
  w=sqrt((-2.0*log(w))/w);

  noise=x1*w;

  return(noise);
}/*GaussNoise*/


/*####################################################*/
/*calculate sigma for link noise*/

float setNoiseSigma(float linkM,float cov,float pSigma,float fSigma)
{
  float sig=0;
  float groundAmp=0;
  float slope=0;
  float gRefl=0;
  float tanSlope=0;
  float sigEff=0;          /*effective ground return width*/
  float probNoise=0,probMiss=0;
  float findSigma(float,float,float,float);

  slope=2.0*M_PI/180.0;

  gRefl=(1.0-cov)*rhoG;

  tanSlope=sin(slope)/cos(slope);
  sigEff=sqrt(pSigma*pSigma+fSigma*fSigma*tanSlope*tanSlope);
  groundAmp=(gRefl/(gRefl+rhoC*cov))/(sigEff*sqrt(2.0*M_PI));  /*normalise by total waveform reflectance*/

  probNoise=0.05;
  probMiss=0.1;
  sig=findSigma(probNoise,probMiss,groundAmp,linkM);

  return(sig);
}/*setNoiseSigma*/


/*####################################################*/
/*find sigma for this combination, from Xioali's notes*/

float findSigma(float probNoise,float probMiss,float groundAmp,float linkM)
{
  float sig=0;
  float step=0;
  float err=0,minErr=0;
  float thisLink=0;
  float nNsig=0,nSsig=0;
  float threshS=0,threshN=0;
  void gaussThresholds(float,float,float,float,float *,float *);

  char direction=0;

  /*initial guess*/
  sig=groundAmp/2.0;
  step=groundAmp/10.0;

  /*determine threshold in terms of standard deviations*/
  gaussThresholds(1.0,XRES,probNoise,probMiss,&nNsig,&nSsig);

  direction=0;
  minErr=0.00015;
  do{
    /*scale thresholds by sigma*/
    threshN=nNsig*sig;
    threshS=nSsig*sig;

    thisLink=10.0*log((groundAmp-threshS)/threshN)/log(10.0);
    if(thisLink<linkM){
      if(direction==1)step/=2.0;
      sig-=step;
      direction=-1;
      if(sig<=0.0){
        step/=2.0;
        sig+=step;
        direction=0;
      }
    }else if(thisLink>linkM){
      if(direction==-1)step/=2.0;
      sig+=step;
      direction=1;
    }
    err=fabs(thisLink-linkM);
  }while(err>minErr);

  return(sig);
}/*findSigma*/

/*####################################################*/
/*calculate Gaussian hresholds*/

void gaussThresholds(float sig,float res,float probNoise,float probMiss,float *threshN,float *threshS)
{
  float x=0,y=0;
  float cumul=0;
  char foundS=0,foundN=0;
  float probN=0,probS=0;

  probN=1.0-probNoise/(30.0/0.15);
  probS=1.0-probMiss;

  /*determine start*/
  x=0.0;
  do{
    y=(float)gaussian((double)x,(double)sig,0.0);
    x-=res;
  }while(y>=YTOL);

  do{
    y=(float)gaussian((double)x,(double)sig,0.0);
    cumul+=y*res;

    if(foundS==0){
      if(cumul>=probS){
        foundS=1;
        *threshS=x;
      }
    }

    if(foundN==0){
      if(cumul>=probN){
        foundN=1;
        *threshN=x;
      }
    }

    x+=res;
  }while((foundS==0)||(foundN==0));

  return;
}/*gaussThresholds*/


/*####################################################*/
/*write waves*/

void writeWaves(dataStruct *data,control *dimage,int numb)
{
  int i=0;
  char namen[200];
  FILE *opoo=NULL;

  sprintf(namen,"%s.%d.wave",dimage->outRoot,numb);
  if((opoo=fopen(namen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",namen);
    exit(1);
  }


  fprintf(opoo,"# 1 range, 2 signal, 3 noised\n");
  for(i=0;i<data->nBins;i++)fprintf(opoo,"%f %f %f\n",data->z[i],data->wave[i],data->noised[i]);

  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Written to %s\n",namen);
  return;
}/*writeWaves*/


/*####################################################*/
/*add physically based noise*/

void addPhysNoise(dataStruct *data,control *dimage)
{
  float *sampled=NULL;
  float *digitiseWave(float *,int,char,float,float);
  float *scaleReflectance(float *,int,float,float,float);
  float *reflWave=NULL,*backWave=NULL;
  float *addBackground(float *,int,float,float);

  /*scale to reflectance and digitse*/
  reflWave=scaleReflectance(data->wave,data->nBins,dimage->rho,dimage->tran,dimage->lPower);

  /*add background light*/
  backWave=addBackground(reflWave,data->nBins,dimage->rho,dimage->solarI);
  TIDY(reflWave);

  /*scale to detected electrons*/


  /*add shot noise*/


  /*digitise*/
  sampled=digitiseWave(backWave,data->nBins,dimage->bitRate,dimage->maxDN,1.0);
  TIDY(backWave);



  /*add electronic noise*/


  return;
}/*addPhysNoise*/


/*####################################################*/
/*add background light*/

float *addBackground(float *reflWave,int nBins,float rho,float solarI)
{
  int i=0;
  float *backWave=NULL;

  backWave=falloc((uint64_t)nBins,"backWave",0);

  for(i=0;i<nBins;i++){
    backWave[i]=reflWave[i]+solarI;   /*I AM NOT CONVINCED BY THIS*/
  }

  return(backWave);
}/*addBackground*/


/*####################################################*/
/*scale to joules leaving the ground*/

float *scaleReflectance(float *wave,int nBins,float rho,float tran,float lPower)
{
  int i=0;
  float tot=0;
  float *refl=NULL;

  /*total energy*/
  tot=0.0;
  for(i=0;i<nBins;i++)tot+=wave[i];

  /*scale*/
  refl=falloc((uint64_t)nBins,"refl",0);
  for(i=0;i<nBins;i++)refl[i]=wave[i]*(rho/tot)*tran*lPower;


  return(refl);
}/*scaleReflectance*/


/*####################################################*/
/*read data*/

dataStruct *readData(char *namen,control *dimage)
{
  int i=0;
  dataStruct *data=NULL;
  char line[400],temp1[100],temp2[100];
  char temp3[100],temp4[100],temp5[100];
  char temp6[100],temp7[100];
  FILE *ipoo=NULL;

  /*open input*/
  if((ipoo=fopen(namen,"r"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",namen);
    exit(1);
  }

  if(!(data=(dataStruct *)calloc(1,sizeof(dataStruct)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }

  /*count number of wavebins*/
  data->nBins=0;
  while(fgets(line,400,ipoo)!=NULL)if(strncasecmp(line,"#",1))data->nBins++;

  data->wave=falloc((uint64_t)data->nBins,"waveform",0);
  data->z=dalloc(data->nBins,"z",0);
  if(dimage->ground)data->ground=falloc((uint64_t)data->nBins,"ground",0);

  /*rewind to start of file*/
  if(fseek(ipoo,(long)0,SEEK_SET)){
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  /*read data*/
  i=0;
  while(fgets(line,400,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      if(dimage->useInt){  /*read intensity*/
        if(dimage->ground==0){   /*don't read ground*/
          if(sscanf(line,"%s %s",temp1,temp2)==2){
            data->z[i]=atof(temp1);
            data->wave[i]=atof(temp2);
            i++;
          }
        }else{                   /*read ground*/
          if(sscanf(line,"%s %s %s %s",temp1,temp2,temp3,temp4)==4){
            data->z[i]=atof(temp1);
            data->wave[i]=atof(temp2);
            data->ground[i]=atof(temp4);
            i++;
          }
        }/*ground switch*/
      }else{ /*read count*/
        if(dimage->ground==0){   /*don't read ground*/
          if(sscanf(line,"%s %s %s",temp1,temp2,temp3)==3){
            data->z[i]=atof(temp1);
            data->wave[i]=atof(temp3);
            i++;
          }
        }else{                   /*read ground*/
          if(sscanf(line,"%s %s %s %s %s %s %s",temp1,temp2,temp3,temp4,temp5,temp6,temp7)==7){
            data->z[i]=atof(temp1);
            data->wave[i]=atof(temp5);
            data->ground[i]=atof(temp7);
            i++;
          }
        }/*ground switch*/
      }/*intensity or count switch*/
    }/*line breaking*/
  }/*line loop*/

  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return(data);
}/*readData*/


/*####################################################*/
/*read command line*/

control *readCommands(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;


  /*allocate structures*/
  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  strcpy(dimage->outRoot,"teast");
  dimage->writeWave=1;

  dimage->ground=1;
  dimage->useInt=0;
  dimage->pSigma=0.764331; /*pulse length*/
  dimage->fSigma=5.5;
  dimage->deSig=0.1; //4.0*0.15/2.355;
  dimage->res=0.15;

  dimage->lPower=9.0/1000.0;
  dimage->tran=0.8;
  dimage->rho=0.5;
  dimage->bitRate=12;
  dimage->maxDN=4096.0; //1.0/(dimage->pSigma*sqrt(2.0*M_PI));
  dimage->offset=94.0;
  dimage->solarI=668.0*0.7;
  dimage->linkM=3.0;
  dimage->trueSig=5.0;   /*night*/

  dimage->linkNoise=1;
  dimage->physNoise=0;

  /*others*/
  rhoG=0.4;
  rhoC=0.57;


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
      }else if(!strncasecmp(argv[i],"-outRoot",8)){
        checkArguments(1,i,argc,"-outRoot");
        strcpy(dimage->outRoot,argv[++i]);
      }else if(!strncasecmp(argv[i],"-linkM",6)){
        checkArguments(2,i,argc,"-linkM");
        dimage->linkM=atof(argv[++i]);
        dimage->linkCov=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-deSigma",8)){
        checkArguments(1,i,argc,"-deSigma");
        dimage->deSig=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-seed",5)){
        checkArguments(1,i,argc,"-seed");
        srand(atoi(argv[++i]));
      }else if(!strncasecmp(argv[i],"-bitRate",8)){
        checkArguments(1,i,argc,"-bitRate");
        dimage->bitRate=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-trueSig",8)){
        checkArguments(1,i,argc,"-trueSig");
        dimage->trueSig=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n#####\nProgram to add noise to GEDI waveforms\n#####\n\n-input name;     lasfile input filename\n-outRoot name;   output filename root\n-inList list;    input file list for multiple files\n-linkM m c;      link margin in db and cover in fraction\n-deSigma s;      detector response width\n-bitRate bits;   digitiser bit rate\n-seed n;         random number seed\n-trueSig sig;    true noise sigma for scaling\n\n");
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
/*##############################*/

