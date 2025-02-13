#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_fft_complex.h"
#include "mpfit.h"
#include "libLasRead.h"
#include "libLasProcess.h"


/*########################*/
/*# Functions to process #*/
/*# waveform lidar data  #*/
/*########################*/

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


#define TOL 0.00001      /*generic floating plint tolerance*/
smoothPulse smooPulse;   /*global structure to save reallocation*/


/*##################################################*/
/*process waveform to denoise and deconvolve*/

float *processWave(unsigned char *wave,int waveLen,denPar *decon,float gbic)
{
  int i=0;
  float *temp=NULL,*processed=NULL;
  float *processFloWave(float *,int,denPar *,float);

  /*convert to a float array for ease*/
  temp=falloc((uint64_t)waveLen,"presmoothed",0);
  for(i=0;i<waveLen;i++)temp[i]=(float)wave[i];
  decon->maxDN=255.0;
  decon->bitRate=8;
  processed=processFloWave(temp,waveLen,decon,gbic);
  TIDY(temp);

  return(processed);
}/*processWave*/


/*##################################################*/
/*process floating point waveform to denoise and deconvolve*/

float *processFloWave(float *wave,int waveLen,denPar *decon,float gbic)
{
  int i=0;
  float *temp=NULL;
  float *mediated=NULL;
  float *preSmoothed=NULL;
  float *mSmoothed=NULL;
  float *denoised=NULL;
  float *smoothed=NULL;
  float *processed=NULL;
  float *denoise(float,float,int,int,float *,float,char);
  float *medianFloat(float *,int,int);
  float *deconvolve(float *,int,float **,int,float,int,double,char);
  float thisTail=0;        /*tail threshold to use here*/
  float *gaussWave=NULL;
  float *fitGaussians(float *,int,denPar *);
  float *CofGhard(float *,uint32_t);
  float *hardHitWave(denPar *,int);
  float *sampled=NULL;
  float *digitise(float *,int,char,float);
  float *correctDrift(float *,int,int,denPar *);
  void medNoiseStats(float *,uint32_t,float *,float *,float *,float,float,char);
  void meanNoiseStats(float *,uint32_t,float *,float *,float *,float,float,int);
  char checkHardEnergy(int *,float *,float);
  char testHard(denPar *,float *,int,float);
  char hardTarget=0;

  /*smooth before denoising*/
  if(decon->psWidth>0.0){   /*Gaussian smoothing*/
    preSmoothed=smooth(decon->psWidth,waveLen,wave,decon->res);
  }else if(decon->preMatchF){  /*matched filter*/
    preSmoothed=matchedFilter(wave,waveLen,decon,decon->res);
  }else{
    preSmoothed=falloc((uint64_t)waveLen,"",0);
    for(i=0;i<waveLen;i++)preSmoothed[i]=wave[i];
  }

  /*determine noise statistics*/
  if(decon->varNoise){
    if(decon->medStats){
      /*convert to a set bit rate*/
      sampled=digitise(preSmoothed,waveLen,decon->bitRate,decon->maxDN);
      medNoiseStats(sampled,waveLen,&(decon->meanN),&(decon->thresh),&thisTail,decon->tailThresh,decon->threshScale,decon->bitRate);
      TIDY(sampled);
    }else meanNoiseStats(preSmoothed,waveLen,&(decon->meanN),&(decon->thresh),&thisTail,decon->tailThresh,decon->threshScale,(int)(decon->statsLen/decon->res));
  }else thisTail=decon->tailThresh;
  if(thisTail<0)thisTail=decon->thresh;

  /*smooth again if needed*/
  if(decon->msWidth>0.0){   /*Gaussian smoothing after noise stats*/
    mSmoothed=smooth(decon->msWidth,waveLen,preSmoothed,decon->res);
    TIDY(preSmoothed);
  }else{
    mSmoothed=preSmoothed;
    preSmoothed=NULL;
  }

  /*median filter if needed*/
  if(decon->medLen>0){
    mediated=medianFloat(mSmoothed,decon->medLen,waveLen);
    TIDY(mSmoothed);
  }else{
    mediated=mSmoothed;
    mSmoothed=NULL;
  }

  /*correct for detector drift if needed*/
  temp=correctDrift(mediated,waveLen,(int)(decon->statsLen/decon->res),decon);
  TIDY(mediated);

  /*remove background noise*/
  if((decon->meanN>0.0)||(decon->thresh>0.0)){
    denoised=denoise(decon->meanN,decon->thresh,decon->minWidth,waveLen,temp,thisTail,decon->noiseTrack);
    TIDY(temp);
  }else{
    denoised=temp;
    temp=NULL;
  }

  /*see if it a single return*/
  if(decon->matchHard)hardTarget=testHard(decon,denoised,waveLen,decon->res);
  else                hardTarget=0;
  if(hardTarget){  /*bunch up the energy*/
    processed=CofGhard(denoised,waveLen);
    TIDY(denoised);
    return(processed);
  }

  /*smooth if required. Note that pulse is smoothed in readPulse()*/
  if(decon->sWidth>0.0){
    smoothed=smooth(decon->sWidth,waveLen,denoised,decon->res);
    TIDY(denoised);
  }else if(decon->posMatchF){  /*matched filter*/
    smoothed=matchedFilter(denoised,waveLen,decon,decon->res);
    TIDY(denoised);
  }else{
    smoothed=denoised;
    denoised=NULL;
  }

  /*scale by GBIC*/
  if((gbic>0.0)&&(gbic!=1.0))for(i=0;i<waveLen;i++)smoothed[i]/=gbic;



  /*deconvolve if required*/
  if((decon->deconMeth>=0)&&(hardTarget==0)){
    processed=deconvolve(smoothed,waveLen,decon->pulse,decon->pBins,\
                  decon->res,decon->maxIter,decon->deChang,decon->deconMeth);
    TIDY(smoothed);
  }else{
    processed=smoothed;    /*otherwise just use the denoised array*/
    smoothed=NULL;
  }

  /*Gaussian fitting*/
  if(decon->fitGauss||decon->gaussFilt){
    gaussWave=fitGaussians(processed,waveLen,decon);
    /*test for hard target*/
    if(decon->gaussFilt){
      if((decon->nGauss==1)&&(decon->gPar[2]<=decon->hardWidth))hardTarget=1;
      else hardTarget=checkHardEnergy(&decon->nGauss,decon->gPar,decon->hardWidth);
    }else hardTarget=0;

    /*copy Gaussian if to be used*/
    if(hardTarget){ /*single hit*/
      TIDY(gaussWave);
      TIDY(processed);
      gaussWave=hardHitWave(decon,waveLen);
    }else if(decon->fitGauss){ /*pass on Gaussian waveform*/
      TIDY(processed);
    }else{                /*delete fitting and pass original*/
      TIDY(gaussWave);
      gaussWave=processed;
      processed=NULL;
    }
  }else{   /*don't fit Gaussians*/
    gaussWave=processed;
    processed=NULL;
    hardTarget=0;
  }

  return(gaussWave);
}/*processFloWave*/


/*################################################*/
/*if hard target, replace with CofG*/

float *CofGhard(float *denoised,uint32_t waveLen)
{
  uint32_t i=0,bin=0;
  float *wave=NULL;
  float CofG=0,contN=0;

  wave=falloc((uint64_t)waveLen,"CofG wave",0);
  for(i=0;i<waveLen;i++)wave[i]=0.0;

  CofG=contN=0.0;
  for(i=0;i<waveLen;i++){
    CofG+=denoised[i]*(float)i;
    contN+=denoised[i];
  }

  if(contN>0.0){
    CofG/=contN;
    bin=(int)CofG;
    wave[bin]=contN;
  }

  return(wave);
}/*CofGhard*/


/*################################################*/
/*test correlation*/

char testHard(denPar *denoise,float *floWave,int nBins,float res)
{
  int maxPulse=0,maxWave=0;
  int nFeat=0;
  int countSigFeats(int,float *,float,float);
  char isHard=0;
  char checkWidth();
  float *matchWave=NULL;
  float *smoothMatchedPulse(int,float *,float,int,float *,float);
  float pulseE=0,waveE=0;
  float RMSE=0,pRes=0;
  float matchRMSE(int,float *,float,int,float,int,float *,float,int,float);
  void waveStats(int *,float *,int,float *,float);


  pRes=denoise->pulse[0][1]-denoise->pulse[0][0];

  /*find peak and energy*/
  waveStats(&maxWave,&waveE,nBins,floWave,res);

  /*see if there is a single feature*/
  nFeat=countSigFeats(nBins,floWave,waveE,res);

  if(nFeat==1){
    /*smooth waveform by pulse*/
    //matchWave=smoothMatchedPulse(denoise->pBins,denoise->pulse[1],pRes,nBins,floWave,res);
    matchWave=smooth(0.3,nBins,floWave,0.15);

    /*align both by peak*/
    waveStats(&maxPulse,&pulseE,denoise->pBins,denoise->pulse[1],pRes);

    /*calculate RMSE*/
    RMSE=matchRMSE(denoise->pBins,denoise->hardPulse,pRes,maxPulse,pulseE,nBins,matchWave,res,maxWave,waveE);
    TIDY(matchWave);

    /*fprintf(stderr,"RMSE %f thresh %f\n",RMSE,denoise->hardThresh);*/
    if(RMSE<=denoise->hardThresh)isHard=1;
    else                         isHard=0;
  }else isHard=0;

  return(isHard);
}/*testHard*/


/*#############################*/
/*RMSE between wave and pulse*/

float matchRMSE(int pBins,float *pulse,float pRes,int maxPulse,float pulseE,int nBins,float *floWave,float res,int maxWave,float waveE)
{
  int i=0,j=0;
  int sBin=0,eBin=0;
  float RMSE=0,meanP=0;
  float eWithout=0,withoutThresh=0;

  withoutThresh=0.1;

  if(pRes>res){
    fprintf(stderr,"Pulse not sampled well enough\n");
    exit(1);
  }

  eWithout=0.0;
  for(i=0;i<nBins;i++){
    sBin=(int)(((i-maxWave)*res+maxPulse*pRes)/pRes);
    eBin=(int)(((i+1-maxWave)*res+maxPulse*pRes)/pRes);
    if(sBin<0)sBin=0;
    else if(sBin>pBins)sBin=pBins;
    if(eBin<0)eBin=0;
    else if(eBin>pBins)eBin=pBins;

    meanP=0.0;
    for(j=sBin;j<eBin;j++)meanP+=pulse[j];
    if((eBin-sBin)>0)meanP*=pRes/((float)(eBin-sBin)*pulseE*res);

    /*fprintf(stdout,"%d %f %f\n",i,meanP,floWave[i]*res/waveE);*/
    if((meanP>0.0)||(floWave[i]>0.0)){
      RMSE+=(meanP-floWave[i]*res/waveE)*(meanP-floWave[i]*res/waveE);
    }
    if((meanP<0.001)&&(floWave[i]>0.001))eWithout+=floWave[i]*res/waveE-meanP;
  }

  if(eWithout>withoutThresh)RMSE=1000.0;    /*significant energy outside*/
  else if(eWithout<0.01)    RMSE=0.0;       /*small hit*/
  else                      RMSE=sqrt(RMSE);/*use RMSE*/
  return(RMSE);
}/*matchRMSE*/


/*#############################*/
/*count significant features*/

int countSigFeats(int nBins,float *floWave,float waveE,float res)
{
  int i=0;
  int nFeat=0;
  float thresh=0;
  float cumul=0;
  char inFeat=0;

  thresh=waveE*0.002;

  cumul=0.0;
  nFeat=0;
  for(i=0;i<nBins;i++){
    if(floWave[i]>TOL){
      inFeat=1;
      cumul+=floWave[i]*res;
    }else if(floWave[i]<=TOL){
      if(inFeat&&(cumul>=thresh))nFeat++;
      cumul=0.0;
      inFeat=0;
    }
  }
  return(nFeat);
}/*countSigFeats*/


/*#############################*/
/*find energy and max location*/

void waveStats(int *maxBin,float *E,int nBins,float *wave,float res)
{
  int i=0;
  float max=0;

  (*E)=0.0;
  max=-1000.0;
  for(i=0;i<nBins;i++){
    if(wave[i]>max){
      max=wave[i];
      (*maxBin)=i;
    }
    (*E)+=wave[i]*res;
  }

  return;
}/*waveStats*/


/*####################################################*/
/*digitse*/

float *digitise(float *wave,int nBins,char bitRate,float maxDN)
{
  int i=0;
  int nDN=0;
  float *sampled=NULL;
  float resDN=0,max=0;
  float tot=0,newTot=0;

  sampled=falloc((uint64_t)nBins,"sampled wave",0);

  /*number of bins*/
  nDN=1;
  for(i=0;i<bitRate;i++)nDN*=2;

  /*find maximum if maxDN not defined*/
  if(maxDN<0.0){
    max=-10000.0;
    for(i=0;i<nBins;i++){
      if(wave[i]>max)max=wave[i];
    }
    if(max<=0.0)max=1.0;
  }else max=maxDN;

  /*find total*/
  tot=0.0;
  for(i=0;i<nBins;i++)tot+=wave[i];

  resDN=max/(float)nDN;
  newTot=0.0;
  for(i=0;i<nBins;i++){
    sampled[i]=floor(wave[i]/resDN+0.5)*resDN;
    newTot+=sampled[i];
  }

  /*rescale energy*/
  if(newTot>0.0){
    for(i=0;i<nBins;i++){
      sampled[i]*=tot/newTot;
    }
  }

  return(sampled);
}/*digitise*/


/*##########################################*/
/*see if one return dominates*/

char checkHardEnergy(int *nGauss,float *gPar,float hardWidth)
{
  int i=0,maxInd=0;
  float total=0,max=0;
  float energy=0,thresh=0;
  char hardHit=0;

  max=-100.0;
  total=0.0;
  for(i=0;i<(*nGauss);i++){
    energy=sqrt(2.0*M_PI)*gPar[3*i+1]*gPar[3*i+2];
    if(energy>max){
      max=energy;
      maxInd=i;
    }
    total+=energy;
  }

  thresh=0.97;   /*97% of energy in a single hit*/
  if((max>=(thresh*total))&&(gPar[3*maxInd+2]<=hardWidth)){
    hardHit=1;
    gPar[0]=gPar[3*maxInd];   /*overwrite Gaussian parameters with brightest*/
    gPar[1]=gPar[3*maxInd+1];
    gPar[2]=gPar[3*maxInd+2];
    (*nGauss)=1;
  }else hardHit=0;

  return(hardHit);
}/*checkHardEnergy*/


/*##########################################*/
/*single hard return*/

float *hardHitWave(denPar *decon,int numb)
{
  int i=0;
  float *hardWave=NULL;

  /*set to 0*/
  hardWave=falloc((uint64_t)numb,"hard waveform",0);
  for(i=0;i<numb;i++)hardWave[i]=0.0;

  /*the single return*/
  if(decon->gPar[0]>=0.0){
    i=(int)(decon->gPar[0]/decon->res);
    if((i>=0)&&(i<numb))hardWave[i]=1.0; /*sqrt(2.0*M_PI)*decon->gPar[1]*decon->gPar[2];*/
    else{
      fprintf(stderr,"Beyond bounds %f %f\n",decon->gPar[0],(float)i*decon->res);
    }
  }

  return(hardWave);
}/*hardHitWave*/


/*##########################################*/
/*correct for attenuation*/

float *correctAttenuation(float *denoised,int numb)
{
  int i=0;
  float tot=0,gap=0;
  float *trueArea=NULL;

  /*count up energy to nromalise visible to unity*/
  tot=0.0;
  for(i=0;i<numb;i++)tot+=denoised[i];

  /*count up gap and apply correction*/
  trueArea=falloc((uint64_t)numb,"",0);
  gap=1.0;
  for(i=0;i<numb;i++){
    if(gap>0.0)trueArea[i]=denoised[i]/gap;
    else       trueArea[i]=0.0;
    gap-=denoised[i]/tot;
  }

  for(i=0;i<numb;i++){
    trueArea[i]/=tot;
    if(trueArea[i]>1.0)trueArea[i]=1.0;
  }

  return(trueArea);
}/*correctAttenuation*/


/*##################################################*/
/*fit Gaussians*/

float *fitGaussians(float *wave,int waveLen,denPar *decon)
{
  int i=0;
  float *fitMultiGauss(float *,float *,int,float,int *,float);
  float *gaussWave=NULL;
  float *x=NULL;

  /*fit Gaussians. Parameters are packed at array end */
  x=falloc((uint64_t)waveLen,"x",0);
  for(i=0;i<waveLen;i++)x[i]=(float)i*decon->res;  /*Gaussian fitting is base 1*/
  decon->nGauss=0;
  gaussWave=fitMultiGauss(x,wave,waveLen,decon->gWidth,&(decon->nGauss),decon->minGsig);
  TIDY(x);

  /*transfer Gaussians parameters*/
  decon->gPar=falloc(3*(uint64_t)decon->nGauss,"Gaussian parameters",0);
  for(i=3*decon->nGauss-1;i>=0;i--)decon->gPar[i]=gaussWave[i+waveLen];

  /*trim the fitted wave*/
  if(!(gaussWave=(float *)realloc(gaussWave,waveLen*sizeof(float)))){
    fprintf(stderr,"Error reallocating %lu\n",waveLen*sizeof(float));
    exit(1);
  }

  return(gaussWave);
}/*fitGaussians*/


/*##################################################*/
/*waveform stats*/

void medNoiseStats(float *wave,uint32_t waveLen,float *meanN,float *thresh,float *tailThresh,float tailOff,float threshScale,char bitRate)
{
  uint32_t i=0;
  int *hist=NULL;
  int maxHist=0;
  int nDN=0,ind=0;
  float hRes=0;
  float min=0,max=0;
  float modQuart=0;
  void modalDeviation(float,float *,uint32_t,float *,int,float);

  nDN=1;
  for(i=0;i<bitRate;i++)nDN*=2;

  min=10000000.0;
  max=-10000000.0;
  for(i=0;i<waveLen;i++){
    if(wave[i]>max)max=wave[i];
    if(wave[i]<min)min=wave[i];
  }
  hRes=(max-min)/(float)(nDN-1);
  if(hRes<=0.0){
    min=0.0;
    max=1.0;
    hRes=1.0;   /*it is a blank wave*/
  }

  /*calculate a histogram and work out the modal noise value*/
  hist=ialloc(nDN+1,"noise histogram",0);
  maxHist=-1000;
  for(i=0;i<waveLen;i++){
    ind=(int)((wave[i]-min)/hRes);
    if((ind>=0)&&(ind<(nDN+1))){
      hist[ind]++;
      if(hist[ind]>maxHist){
        maxHist=hist[ind];
        (*meanN)=wave[i];  /*modal noise*/
      }
    }else{
      fprintf(stderr,"index issue %d %f\n",ind,wave[i]);
      exit(1);
    }
  }

  /*modal quartile*/
  modalDeviation(*meanN,wave,waveLen,&modQuart,nDN,hRes);
  if(modQuart<2)modQuart=2;

  (*thresh)=(*meanN)+threshScale*(float)modQuart;
  if((*thresh)<14.0)(*thresh)=14.0;
  (*tailThresh)=(*thresh)+tailOff;

  TIDY(hist);
  return;
}/*medNoiseStats*/


/*############################################*/
/*mean noise statistics*/

void meanNoiseStats(float *sampled,uint32_t waveLen,float *meanN,float *thresh,float *thisTail,float tailThresh,float threshScale,int statBins)
{
  int i=0,start=0;
  float stdev=0;

  start=2;
  if((uint32_t)statBins>waveLen){
    fprintf(stderr,"Not enough bins for this statistics length %d %d\n",statBins,(int)waveLen);
    exit(1);
  }else if((statBins-start)<=0){
    fprintf(stderr,"What are you doing? Start too soon %d %d\n",statBins,start);
    exit(1);
  }

  /*NOTE, subtract 2 from noise bins to avoid low smoothing at signal start*/
  (*meanN)=0.0;
  for(i=start;i<statBins;i++)(*meanN)+=sampled[i];
  (*meanN)/=(float)(statBins-start);
  stdev=0.0;
  for(i=start;i<statBins;i++)stdev+=((*meanN)-sampled[i])*((*meanN)-sampled[i]);
  stdev=sqrt(stdev/(float)(statBins-start));


  (*thresh)=(*meanN)+threshScale*stdev;
  if(tailThresh<=0.0)(*thisTail)=(*thresh);
  else               (*thisTail)=(*meanN)+tailThresh*stdev;

  return;
}/*meanNoiseStats*/



/*############################################*/
/*deviation from the mode histogram*/

void modalDeviation(float modal,float *wave,uint32_t waveLen,float *modQuart,int nDN,float hRes)
{
  int i=0;
  int *devHist=NULL;
  int bin=0,total=0;
  int quartThresh=0;
  int cumul=0,subTot=0;

  total=subTot=0;
  devHist=ialloc(nDN,"",0);
  for(i=0;i<nDN;i++)devHist[i]=0;
  for(i=0;i<waveLen;i++){
    bin=(int)((wave[i]-(int)modal)/hRes+0.5);
    if(bin>=0){
      devHist[bin]++;
      subTot++;
    }
    total++;
  }

  /*calculate mode and top quartile*/
  quartThresh=(int)((float)total*3.0/4.0);
  cumul=total-subTot;
  for(i=0;i<nDN;i++){
    cumul+=devHist[i];
    if(cumul>=quartThresh){
      (*modQuart)=(unsigned char)i;
      break;
    }
  }
  TIDY(devHist);
  return;
}/*modalDeviation*/


/*##################################################*/
/*matched filter*/

float *matchedFilter(float *wave,int nBins,denPar *denoise,float res)
{
  int i=0,j=0,bin=0;
  float *smoothed=NULL;
  float energy=0;

  /*allocate space*/
  smoothed=falloc((uint64_t)nBins,"matched filtered wave",0);

  /*smooth wave*/
  for(i=0;i<nBins;i++){
    smoothed[i]=0.0;
    energy=0.0;
    for(j=0;j<denoise->pBins;j++){
      bin=(int)ceil((denoise->pulse[0][j]-denoise->pulse[0][denoise->maxPbin])/res)+i;  /*round up to match standard*/
      if((bin>=0)&&(bin<nBins)){
        smoothed[i]+=denoise->matchPulse[j]*wave[bin];
        energy+=denoise->matchPulse[j];
      }
    }/*pulse loop*/
    if(energy>0.0)smoothed[i]/=energy;
  }/*wave loop*/

  return(smoothed);
}/*matchedFilter*/


/*##################################################*/
/*smooth waveform*/

float *smooth(float sWidth,int nBins,float *data,float res)
{
  int i=0,j=0;
  int tP=0;
  int p1=0,p2=0;         /*array places*/
  int nPulse=0;          /*number of pulse bins*/
  float *smoothed=NULL;
  float energy=0,newRes=0;
  float *setPulse(float,int *,float);
  float *markFloat(int,float *,float);
  int *markInt(int,int *,int);
  char newPulse=0;


  /*Nyquist sample*/
  newRes=res/3.0;

  /*smooth as required*/
  smoothed=falloc((uint64_t)nBins,"",0);
  if(sWidth<TOL)memcpy(smoothed,data,sizeof(float)*nBins);
  else{
    /*test to see if a new pulse array is needed*/
    newPulse=1;
    for(tP=0;tP<smooPulse.nPulses;tP++){
      if((fabs(sWidth-smooPulse.sWidth[tP])<TOL)&&((fabs((newRes)-smooPulse.res[tP])<TOL))){
        newPulse=0;
        break;
      }
    }/*new pulse test*/

    /*create new pulse if needed*/
    if(newPulse){
      tP=smooPulse.nPulses;
      if(!(smooPulse.pulse=(float **)realloc(smooPulse.pulse,(smooPulse.nPulses+1)*sizeof(float *)))){
        fprintf(stderr,"Error reallocating %lu\n",(smooPulse.nPulses+1)*sizeof(float *));
        exit(1);
      }
      smooPulse.res=markFloat(smooPulse.nPulses,smooPulse.res,newRes);
      smooPulse.pulse[tP]=setPulse(sWidth,&nPulse,smooPulse.res[smooPulse.nPulses]);
      smooPulse.sWidth=markFloat(smooPulse.nPulses,smooPulse.sWidth,sWidth);
      smooPulse.nBins=markInt(smooPulse.nPulses,smooPulse.nBins,nPulse);
      smooPulse.nPulses++;
    }

    for(i=0;i<nBins;i++){
      smoothed[i]=smooPulse.pulse[tP][0]*data[i];
      energy=smooPulse.pulse[tP][0];
      for(j=1;j<smooPulse.nBins[tP];j++){
        p1=i-(int)((float)j*smooPulse.res[tP]/res+0.5);
        p2=i+(int)((float)j*smooPulse.res[tP]/res+0.5);
        if((p1>=0)&&(p1<nBins)){
          smoothed[i]+=smooPulse.pulse[tP][j]*data[p1];
          energy+=smooPulse.pulse[tP][j];
        }
        if((p2>=0)&&(p2<nBins)){
          smoothed[i]+=smooPulse.pulse[tP][j]*data[p2];
          energy+=smooPulse.pulse[tP][j];
        }
      }
      smoothed[i]/=energy;
    }/*bin loop*/
  }/*smooth check*/

  return(smoothed);
}/*smooth*/


/*##################################################*/
/*tidy up smoothing pulse arrays*/

void tidySMoothPulse()
{
  TIDY(smooPulse.res);
  TIDY(smooPulse.nBins);
  TIDY(smooPulse.sWidth);
  TTIDY((void **)smooPulse.pulse,smooPulse.nPulses);

  return;
}/*tidySMoothPulse*/


/*##################################################*/
/*set smoothing pulse*/

float *setPulse(float sWidth,int *nPulse,float res)
{
  int i=0;
  float *pulse=NULL;
  float y=0,x=0;
  float minY=0;
  double gaussian(double,double,double);

  minY=0.0001;
  for(i=0;;i++){
    x=(float)i*res;
    y=(float)gaussian((double)x,(double)sWidth,0.0);
    if(y<=minY)break;
  }

  (*nPulse)=i;
  pulse=falloc((uint64_t)(*nPulse),"smoothing pulse",0);
  for(i=0;i<(*nPulse);i++){
    x=(float)i*res;
    pulse[i]=(float)gaussian((double)x,(double)sWidth,0.0);
  }

  return(pulse);
}/*setPulse*/


/*##################################################*/
/*deconvolve using Gold's method*/

float *deconvolve(float *data,int nBins,float **pulse,int pBins,float res,int maxIter,double minChange,char meth)
{
  int i=0,numb=0;
  float *decon=NULL;
  double *deconDo=NULL;
  double *pulseDo=NULL;
  double *dataDo=NULL;
  double *resamplePulse(int,float **,float,int);
  double *goldMeth(double *,double *,int,int,double);
  double *richLucy(double *,double *,int,int,double);
  double energy=0,newEn=0;   /*to balance energies*/

  /*arrays have to be of base 2 length*/
  numb=pow(2.0,(float)((int)(log((double)nBins)/log(2.0)+0.5)+1));
  dataDo=dalloc(numb,"dataDo",0);

  energy=0.0;
  for(i=0;i<nBins;i++){
    dataDo[i]=(double)data[i];
    energy+=dataDo[i];
  }

  /*pulse needs resampling to the correct resolution*/
  pulseDo=resamplePulse(numb,pulse,res,pBins);

  /*call deconvolution method*/
  if(meth==0)     deconDo=goldMeth(dataDo,pulseDo,numb,maxIter,minChange);
  else if(meth==1)deconDo=richLucy(dataDo,pulseDo,numb,maxIter,minChange);
  else{
    fprintf(stderr,"Deconvolution method not defined\n");
    exit(1);
  }
  TIDY(pulseDo);  /*tidy as we go along*/
  TIDY(dataDo);

  /*transfer data*/
  decon=falloc((uint64_t)nBins,"",0);
  newEn=0;
  for(i=0;i<nBins;i++)newEn+=deconDo[i];
  for(i=0;i<nBins;i++)decon[i]=(float)(deconDo[i]*energy/newEn);

  TIDY(deconDo);
  return(decon);
}/*deconvolve*/


/*##################################################*/
/*Richardson-Lucy deconvolution*/

double *richLucy(double *data,double *pulse,int numb,int maxIter,double minChange)
{
  int i=0,j=0;
  int real=0,imag=0;
  double *o=NULL,*work=NULL;
  double *smooth=NULL,*denom=NULL;
  double changeSq=0,minChangeSq=0; //tot=0
  double new=0;
  int gsl_fft_complex_radix2_forward(gsl_complex_packed_array,size_t,size_t);
  int gsl_fft_complex_radix2_backward(gsl_complex_packed_array, size_t,size_t);

  minChangeSq=minChange*minChange;

  /*perform FFT*/
  gsl_fft_complex_radix2_forward((gsl_complex_packed_array)pulse,1,numb);

  /*real arrays*/
  o=dalloc(numb,"o",0);
  denom=dalloc(numb,"denominator",0);
  /*complex arrays*/
  work=dalloc(2*numb,"workspace",0);
  smooth=dalloc(2*numb,"workspace",0);

  /*initial guess*/
  for(i=0;i<numb;i++){
    o[i]=data[i];
    //tot+=data[i];
  }

  do{  /*the iterative deconvolution*/
    /*convolve pulse and estimate*/
    for(i=0;i<numb;i++){
      work[2*i]=o[i];
      work[2*i+1]=0.0;
    }
    gsl_fft_complex_radix2_forward((gsl_complex_packed_array)work,1,numb);
    for(i=0;i<numb;i++){
      real=i*2;
      imag=i*2+1;
      smooth[real]=pulse[real]*work[real]-pulse[imag]*work[real];
      smooth[imag]=pulse[imag]*work[real]+pulse[real]*work[imag];
    }
    gsl_fft_complex_radix2_backward((gsl_complex_packed_array)smooth,1,numb);

    /*calculate new estimate*/
    for(i=0;i<numb;i++){
      real=i*2;
      imag=i*2+1;
      denom[i]=sqrt(smooth[real]*smooth[real]+smooth[imag]*smooth[imag]);
      if(denom[i]>0.0)work[real]=data[i]/denom[i];
      else            work[real]=0.0;
      work[imag]=0.0;
    }
    gsl_fft_complex_radix2_forward((gsl_complex_packed_array)work,1,numb);
    for(i=0;i<numb;i++){
      real=2*i;
      imag=2*i+1;
      smooth[real]=pulse[real]*work[real]-pulse[imag]*work[real];
      smooth[imag]=pulse[imag]*work[real]+pulse[real]*work[imag];
    }
    gsl_fft_complex_radix2_backward((gsl_complex_packed_array)smooth,1,numb);

    //tot=0.0;
    changeSq=0.0;
    for(i=0;i<numb;i++){
      real=2*i;
      imag=2*i+1;
      new=o[i]*sqrt(smooth[real]*smooth[real]+smooth[imag]*smooth[imag]);
      changeSq+=(o[i]-new)*(o[i]-new);
      o[i]=new;
      //tot+=o[i];
    }
    changeSq/=(double)numb;
    /*fprintf(stdout,"Iter %d change %.20f\n",j,changeSq);*/
    if((minChange>=0.0)&&(changeSq<=minChangeSq))break;
    j++;
  }while(j<maxIter);

  TIDY(work);
  TIDY(smooth);
  TIDY(denom);
  return(o);
}/*richLucy*/


/*##################################################*/
/*Gold's method*/

double *goldMeth(double *data,double *pulse,int numb,int maxIter,double minChange)
{
  int i=0,j=0;
  int real=0,imag=0;
  double *o=NULL,*work=NULL;
  double *smooth=NULL,*denom=NULL;
  double changeSq=0,minChangeSq=0;
  double new=0;
  int gsl_fft_complex_radix2_forward(gsl_complex_packed_array,size_t,size_t);
  int gsl_fft_complex_radix2_backward(gsl_complex_packed_array, size_t,size_t);

  minChangeSq=minChange*minChange;

  /*perform FFT*/
  gsl_fft_complex_radix2_forward((gsl_complex_packed_array)pulse,1,numb);

  /*real arrays*/
  o=dalloc(numb,"o",0);
  denom=dalloc(numb,"denominator",0);
  /*complex arrays*/
  work=dalloc(2*numb,"workspace",0);
  smooth=dalloc(2*numb,"workspace",0);

  /*initial guess*/
  for(i=0;i<numb;i++)o[i]=data[i];

  /*iterate over Gold's method*/
  do{
    /*fourier current estimate*/
    for(i=0;i<numb;i++){
      work[2*i]=o[i];
      work[2*i+1]=0.0;
    }
    gsl_fft_complex_radix2_forward((gsl_complex_packed_array)work,1,numb);

    /*blur with pulse*/
    for(i=0;i<numb;i++){
      real=i*2;
      imag=i*2+1;
      smooth[real]=pulse[real]*work[real]-pulse[imag]*work[real];
      smooth[imag]=pulse[imag]*work[real]+pulse[real]*work[imag];
    }
    gsl_fft_complex_radix2_backward((gsl_complex_packed_array)smooth,1,numb);

    /*reblur deoniminator*/
    changeSq=0.0;
    for(i=0;i<numb;i++){
      real=i*2;
      imag=i*2+1;
      denom[i]=sqrt(smooth[real]*smooth[real]+smooth[imag]*smooth[imag]);
      if(denom[i]>0.0){
        new=o[i]*data[i]/denom[i];
      }else new=0.0;
      changeSq+=(o[i]-new)*(o[i]-new);
      o[i]=new;
    }
    changeSq/=(double)numb;
    /*fprintf(stdout,"Iter %d change %.20f\n",j,changeSq);*/
    if((minChange>=0.0)&&(changeSq<=minChangeSq))break;
    j++;
  }while(j<maxIter);

  TIDY(work);
  TIDY(smooth);
  TIDY(denom);

  return(o);
}/*goldMeth*/


/*##################################################*/
/*denoise by thresholding*/

float *denoise(float meanN,float thresh,int minWidth,int nBins,float *data,float tailThresh,char noiseTrack)
{
  int i=0,j=0;
  int start=0;
  float *denoised=NULL;
  float thisMean=0;
  char waveStart=0;

  denoised=falloc((uint64_t)nBins,"denoised",0);
  for(i=0;i<nBins;i++)denoised[i]=0.0;
  waveStart=0;

  if(noiseTrack)thisMean=meanN;
  else          thisMean=thresh;

  /*threshold*/
  start=-1;
  for(i=0;i<nBins;i++){
    if((data[i]>=tailThresh)||((waveStart==0)&&(data[i]>=thresh))){        /*check whether we're within a feature*/
    //if((waveStart==0)&&(data[i]>=thresh)){        /*check whether we're within a feature*/
      if(start<0)start=i;       /*mark start*/
      waveStart=1;
    }else if(start>=0){         /*left a feature*/
      if((i-start)>=minWidth){  /*check feature width*/
        for(j=start;j>=0;j--){      /*noise tracking*/
          if(data[j]<=thisMean)break;
          denoised[j]=data[j]-meanN;
          if(denoised[j]<0.0)denoised[j]=0.0;  /*force positive to avoid nan in Gold*/
        }
        for(i=start+1;i<nBins;i++){
          if(data[i]<=thisMean)break;
          denoised[i]=data[i]-meanN;
          if(denoised[i]<0.0)denoised[i]=0.0;  /*force positive to avoid nan in Gold*/
        }
      }
      start=-1;  /*reset counter*/
    }/*within feature check*/
  }/*threshold*/

  return(denoised);
}/*denoise*/


/*##################################################*/
/*resample pulse to complex array at correct res*/

double *resamplePulse(int numb,float **pulse,float res,int pBins)
{
  int i=0,bin=0,step=0,maxBin=0;
  float max=0,maxRange=0;
  double *pulseDo=NULL,total=0;

  pulseDo=dalloc(2*numb,"pulseDo",0);
  for(i=2*numb-1;i>=0;i--)pulseDo[i]=0.0;

  /*find the peak to set at zero as pulse is aligned by CofG*/
  max=-1000.0;
  for(i=0;i<pBins;i++){
    if(pulse[1][i]>max){
      max=pulse[1][i];
      maxRange=pulse[0][i];
      maxBin=i;
    }
  }/*peak finding*/

  /*find nearest pulse point to bin, start from centre*/
  step=(int)(res/(pulse[0][1]-pulse[0][0]));
  if(step<1)step=1;

  total=0.0;
  for(i=maxBin;i<pBins;i+=step){
    bin=(int)((pulse[0][i]-maxRange)/res);
    if(bin<0)bin+=numb;
    pulseDo[2*bin]=pulse[1][i];
    total+=pulse[1][i];
  }/*from centre up*/
  for(i=maxBin-1;i>=0;i-=step){  /*then work from the centre backwards*/
    bin=(int)((pulse[0][i]-maxRange)/res);
    if(bin<0)bin+=numb;
    pulseDo[2*bin]=pulse[1][i];
    total+=pulse[1][i];
  }/*from centre back*/

  /*normalise pulse and prevent zeroes*/
  for(i=0;i<numb;i++){
    if(pulseDo[2*i]<0.0)pulseDo[2*i]=0.0;
    else                pulseDo[2*i]/=total;
  }/*normalisation*/

  return(pulseDo);
}/*reamplePulse*/


/*##################################################*/
/*read pulse shape*/

void readPulse(denPar *denoise)
{
  int i=0,nGauss=0;
  int nMax=0;
  char line[200];
  char temp[2][100];
  FILE *ipoo=NULL;
  /*smoothing*/
  float *smoothed=NULL;
  float *tempPulse=NULL;
  float *fitSingleGauss(float *,float *,int,float,int *,float **);
  float *gaussPar=NULL;
  float max=0,pRes=0,mu=0;
  float *setPulse(float,int *,float);

  /*is it an assymmetric or gaussian pulse?*/
  if(!denoise->deconGauss){
    if((ipoo=fopen(denoise->pNamen,"r"))==NULL){
      fprintf(stderr,"Error opening pulse file %s\n",denoise->pNamen);
      exit(1);
    }

    denoise->pBins=0;
    /*count number of bins*/
    while(fgets(line,200,ipoo)!=NULL){
      if(strncasecmp(line,"#",1))denoise->pBins++;
    }/*line counting*/

    /*rewind*/
    if(fseek(ipoo,(long)0,SEEK_SET)){ /*rewind to start of file*/
      fprintf(stderr,"fseek error\n");
      exit(1);
    }

    denoise->pulse=fFalloc(2,"",0);
    for(i=0;i<2;i++)denoise->pulse[i]=falloc((uint64_t)denoise->pBins,"",i+1);

    /*read data*/
    i=0;
    max=-1000.0;
    nMax=0;
    while(fgets(line,200,ipoo)!=NULL){
      if(strncasecmp(line,"#",1)){
        if(sscanf(line,"%s %s",temp[0],temp[1])==2){
          if(i>denoise->pBins){
            fprintf(stderr,"Error\n");
            exit(1);
          }
          denoise->pulse[0][i]=atof(&(temp[0][0]))*denoise->pScale;
          denoise->pulse[1][i]=atof(&(temp[1][0]));
          if(denoise->pulse[1][i]>=max){
            max=denoise->pulse[1][i];
            denoise->maxPbin=i;
            if(i>0){
              if(denoise->pulse[1][i-1]>=denoise->pulse[1][i])nMax++;
            }
          }
          i++;
        }
      }/*comment check*/
    }/*data reading*/

    /*if it is a pulse compressed chirp, middle is centre*/
    if(nMax>2)denoise->maxPbin=denoise->pBins/2;

    if(ipoo){
      fclose(ipoo);
      ipoo=NULL;
    }
  }else{  /*Gaussian pulse*/
    pRes=0.01;
    denoise->pulse=fFalloc(2,"",0);
    tempPulse=setPulse(denoise->pSigma*denoise->pScale,&denoise->pBins,pRes);
    for(i=0;i<2;i++)denoise->pulse[i]=falloc(2*(uint64_t)denoise->pBins,"",i+1);
    mu=(float)denoise->pBins*pRes;
    for(i=0;i<2*denoise->pBins;i++){
      denoise->pulse[0][i]=(float)i*pRes-mu;
      if(i>=denoise->pBins)denoise->pulse[1][i]=tempPulse[i-denoise->pBins];
      else if(i>0)         denoise->pulse[1][i]=tempPulse[denoise->pBins-i];
      else                 denoise->pulse[1][i]=0.0;
    }
    TIDY(tempPulse);
  }/*Gaussian or assymmetric pulse test*/

  /*if we want a pulse to detect hard targets*/
  if(denoise->matchHard){
    denoise->hardPulse=smooth(0.3,denoise->pBins,denoise->pulse[1],denoise->pulse[0][1]-denoise->pulse[0][0]);
  }

  /*matched filter if we need it*/
  if(denoise->preMatchF||denoise->posMatchF){
    denoise->matchPulse=falloc((uint64_t)denoise->pBins,"matched pulse",0);
    for(i=0;i<denoise->pBins;i++)denoise->matchPulse[i]=denoise->pulse[1][i];
  }

  /*smooth if required*/
  if(denoise->sWidth>0.0){
    smoothed=smooth(denoise->sWidth,denoise->pBins,&(denoise->pulse[1][0]),denoise->pulse[0][1]-denoise->pulse[0][0]);
    TIDY(denoise->pulse[1])
    denoise->pulse[1]=smoothed;
    smoothed=NULL;
  }/*smoothing*/
  if(denoise->psWidth>0.0){
    smoothed=smooth(denoise->psWidth,denoise->pBins,&(denoise->pulse[1][0]),denoise->pulse[0][1]-denoise->pulse[0][0]);
    TIDY(denoise->pulse[1])
    denoise->pulse[1]=smoothed;
    smoothed=NULL;
  }/*smoothing*/

  /*smooth deconvolution by matched flter if needed*/
  if((denoise->preMatchF||denoise->posMatchF)&&(denoise->deconMeth>=0)){
    smoothed=matchedFilter(&(denoise->pulse[1][0]),denoise->pBins,denoise,denoise->pulse[0][1]-denoise->pulse[0][0]);
    TIDY(denoise->pulse[1]);
    denoise->pulse[1]=smoothed;
    smoothed=NULL;
  }/*smoothing*/

  /*fit a Gaussian or determine width for hard limits if required*/
  if((denoise->gaussPulse)||(denoise->gaussFilt)){
    nGauss=0;
    tempPulse=fitSingleGauss(denoise->pulse[0],denoise->pulse[1],denoise->pBins,0.5,&(nGauss),&gaussPar);

    /*set width if needed*/
    if(denoise->gaussFilt){
      if(nGauss>1){
        fprintf(stderr,"Multiple Gaussians fitted to pulse\n");
        exit(1);
      }
      denoise->hardWidth=gaussPar[2]*denoise->hardTol;
    }
    TIDY(gaussPar);

    /*copy Gaussian over pulse if needed*/
    if(denoise->gaussPulse){
      TIDY(denoise->pulse[1]);
      denoise->pulse[1]=tempPulse;
      tempPulse=NULL;
    }else TIDY(tempPulse);
  }

  return;
}/*readPulse*/


/*########################################################################*/
/*set default denoising parameters*/

void setDenoiseDefault(denPar *denoise)
{
  /*denoising*/
  denoise->meanN=12.0;
  denoise->tailThresh=-1.0;
  denoise->thresh=15.0;
  denoise->minWidth=6;
  denoise->sWidth=0.0;
  denoise->msWidth=0.0;
  denoise->psWidth=0.0;
  denoise->medLen=0;
  denoise->varNoise=0;
  denoise->medStats=0;
  denoise->statsLen=30.0;
  denoise->noiseTrack=1;
  denoise->threshScale=1.0;
  denoise->bitRate=8;
  denoise->maxDN=-1.0;
  denoise->preMatchF=0;    /*no matched filter before denoising*/
  denoise->posMatchF=0;    /*no matched filter after denoising*/
  /*detector drift*/
  denoise->corrDrift=0;    /*do not correct for drift*/
  denoise->varDrift=1;     /*if we do correct, use a variable drift factor*/
  denoise->fixedDrift=0.0; /*if we are using a fixed drift factor, use this*/

  /*deconvolution*/
  denoise->deconMeth=-1;     /*do not deconvolve*/
  denoise->pScale=1.0;      /*scale pulse length by*/
  denoise->maxIter=2000;     /*maximum number of iterations*/
  denoise->deChang=pow(10,0-7.0);  /*change between decon iterations to stop*/
  denoise->deconGauss=1;     /*use Gaussian pulse by default*/
  denoise->pSigma=1.0;
  strcpy(denoise->pNamen,"/home/sh563/data/bess/leica_shape/leicaPulse.dat");  /*pulse filename*/
  denoise->pulse=NULL;       /*pulse to deconvolve by*/
  denoise->pBins=0;          /*number of pulse bins*/
  denoise->res=0.15;

  /*Gaussian fitting*/
  denoise->gWidth=1.5;
  denoise->fitGauss=0;       /*do not fit Gaussians*/
  denoise->gaussPulse=0;     /*do not turn pulse to Gaussian*/
  denoise->minGsig=0.00001;  /*minimum Gaussian fitting width*/

  /*Gaussian hard target identification*/
  denoise->gaussFilt=0;    /*switch*/
  denoise->hardWidth=0.0;  /*maxWidth of hard feature*/
  denoise->hardTol=1.0;    /*tolerance to scale width by*/

  /*correlation hard target finding*/
  denoise->matchHard=0;

  return;
}/*setDenoiseDefault*/


/*########################################################################*/
/*rotate about x axis*/

void rotateX(double *vect,double theta)
{
  int i=0;
  double temp[3];

  temp[0]=vect[0];
  temp[1]=vect[1]*cos(theta)+vect[2]*sin(theta);
  temp[2]=vect[2]*cos(theta)-vect[1]*sin(theta);

  for(i=0;i<3;i++)vect[i]=temp[i];
  return;
}/*rotateX*/


/*########################################################################*/
/*rotate about y axis*/

void rotateY(double *vect,double theta)
{
  int i=0;
  double temp[3];

  temp[0]=vect[0]*cos(theta)-vect[1]*sin(theta);
  temp[1]=vect[1];
  temp[2]=vect[0]*sin(theta)+vect[2]*cos(theta);

  for(i=0;i<3;i++)vect[i]=temp[i];
  return;
}/*rotateY*/


/*########################################################################*/
/*rotate about z axis*/

void rotateZ(double *vect,double theta)
{
  int i=0;
  double temp[3];

  temp[0]=vect[0]*cos(theta)+vect[1]*sin(theta);
  temp[1]=vect[1]*cos(theta)-vect[0]*sin(theta);
  temp[2]=vect[2];

  for(i=0;i<3;i++)vect[i]=temp[i];
  return;
}/*rotateZ*/


/*########################################################################*/
/*check point bounds*/

char boundsCheck(double x,double y,double z,double *bounds)
{
  if((x>=(bounds[0]))&&(y>=(bounds[1]))&&(z>=(bounds[2]))&&\
     (x<=(bounds[3]))&&(y<=(bounds[4]))&&(z<=(bounds[5])))return(1);
  else                                                    return(0);
}/*boundsCheck*/


/*########################################################################*/
/*DEM by nearest neighbour*/

double *findGroundNN(pCloudStruct **data,int nFiles,double *minX,double *minY,float res,int *nX,int *nY,double groundBreakElev)
{
  double *gDEM=NULL;

  fprintf(stderr,"The nearest neighbour DEM option is not ready yet\n");
  exit(1);

  return(gDEM);
}/*findGroundNN*/


/*########################################################################*/
/*fit a polynomial to the ground*/

double *findGroundPoly(pCloudStruct **data,int nFiles,double *minX,double *minY,double *maxX,double *maxY,float res,int *nX,int *nY,double groundBreakElev)
{
  int i=0;
  double *gDEM=NULL;
  double *polyDEM(groundDstruct *,double,double,float,int,int);
  groundDstruct *groundD=NULL;   /*ground data structure*/
  groundDstruct *arrangeGroundData(pCloudStruct **,int,double,double,double,double,double);
  void fitManyPlanes(groundDstruct *,int,int);


  /*allocate and load relevant data*/
  groundD=arrangeGroundData(data,nFiles,groundBreakElev,*minX,*minY,*maxX,*maxY);

  /*find bounds if needed*/
  if(*maxX<-1000.0){  /*if <1000, no vounds defined*/
    *minX=*minY=100000000000.0;
    *maxX=*maxY=-100000000000.0;
    for(i=0;i<nFiles;i++){
      if(data[i]->nPoints>0){
        if(data[i]->bounds[0]<*minX)*minX=data[i]->bounds[0];
        if(data[i]->bounds[1]<*minY)*minY=data[i]->bounds[1];
        if(data[i]->bounds[3]>*maxX)*maxX=data[i]->bounds[3];
        if(data[i]->bounds[4]>*maxY)*maxY=data[i]->bounds[4];
      }
    }
  }
  *nX=(int)((*maxX-*minX)/(double)res);
  *nY=(int)((*maxY-*minY)/(double)res);

  if(groundD->nPoints>0){  /*check that there are ground returns*/
    /*fit the ground*/
    fitManyPlanes(groundD,1,1);

    /*make a ground DEM*/
    gDEM=polyDEM(groundD,*minX,*minY,res,*nX,*nY);
  }else gDEM=NULL;

  /*tidy up*/
  if(groundD){
    TIDY(groundD->xUse);
    TIDY(groundD->yUse);
    TIDY(groundD->zUse);
    TIDY(groundD->par);
    TIDY(groundD);
  }
  return(gDEM);
}/*findGroundPoly*/


/*########################################*/
/*make a ground DEM*/

double *polyDEM(groundDstruct *groundD,double minX,double minY,float res,int nX,int nY)
{
  int i=0,j=0,place=0;
  double *gDEM=NULL;
  double x=0,y=0;
  double polyGround(double,double,groundDstruct *);

  if(nX*nY>0)gDEM=dalloc(nX*nY,"ground DEM",0);

  for(i=0;i<nX;i++){
    x=((double)i+0.5)*(double)res+minX;
    for(j=0;j<nY;j++){
      place=j*nX+i;
      y=((double)j+0.5)*(double)res+minY;
      gDEM[place]=polyGround(x,y,groundD);
    }
  }

  return(gDEM);
}/*polyDEM*/


/*########################################*/
/*calculate elevation from polynomial*/

double polyGround(double x,double y,groundDstruct *groundData)
{
  int k=0;
  double z=0;

  z=0.0;
  for(k=groundData->nPoly[0]-1;k>=0;k--)z+=groundData->par[k]*pow(x,(double)(groundData->nPoly[0]-k));
  for(k=groundData->nPoly[1]-1;k>=0;k--)z+=groundData->par[k+groundData->nPoly[0]]*pow(y,(double)(groundData->nPoly[1]-k-1));

  return(z);
}/*polyGround*/


/*########################################*/
/*fit many polynomial planes and average*/

void fitManyPlanes(groundDstruct *groundData,int cNx,int cNy)
{
  int i=0,j=0;
  void fitPolyPlane(groundDstruct *);

  for(i=cNx*cNy-1;i>=0;i--){
    groundData[i].nPoly[0]=6;
    groundData[i].nPoly[1]=7;
    if(groundData[i].nPoints>(groundData[i].nPoly[0]+groundData[i].nPoly[1])){
      fprintf(stdout,"Fitting %d of %d\n",i,cNx*cNy);
      fitPolyPlane(&(groundData[i]));
    }else{
      groundData[i].par=dalloc(groundData[i].nPoly[0]+groundData[i].nPoly[1],"pars",0);
      for(j=0;j<(groundData[i].nPoly[0]+groundData[i].nPoly[1]);j++)groundData[i].par[j]=0.0;
    }
  }

  return;
}/*fitManyPlanes*/


/*########################################*/
/*fit a polynomial plane to ground*/

void fitPolyPlane(groundDstruct *groundData)
{
  int i=0,nPar=0;
  int fitCheck=0;
  int mpfit(mp_func,int,int,double *, mp_par *,mp_config *,void *,mp_result *);
  int groundErr(int,int,double *,double *,double **,void *);
  double meanZ=0;
  mp_par *setGroundBounds(int *);
  mp_par *parStruct=NULL;
  mp_result *result=NULL;
  mp_config *config=NULL;

  /*ninth order for now. No real reason*/

  /*allocate and initial guess*/
  nPar=groundData->nPoly[0]+groundData->nPoly[1];
  groundData->par=dalloc(nPar,"parameters",0);
  for(i=0;i<nPar;i++)groundData->par[i]=0.0;
  meanZ=0.0;
  for(i=0;i<groundData->nPoints;i++)meanZ+=groundData->zUse[i];
  meanZ/=(double)groundData->nPoints;
  groundData->par[nPar-1]=meanZ;

  /*allocate arrays*/
  if(!(config=(mp_config *)calloc(1,sizeof(mp_config)))){
    fprintf(stderr,"error in control structure.\n");
    exit(1);
  }
  config->nofinitecheck=1;
  if(!(result=(mp_result *)calloc(1,sizeof(mp_result)))){
    fprintf(stderr,"error in mpfit structure.\n");
    exit(1);
  }
  result->resid=dalloc(groundData->nPoints,"",0);
  result->xerror=dalloc(nPar,"",0);
  result->covar=dalloc(nPar*nPar,"",0);

  /*set parameter bounds*/
  parStruct=setGroundBounds(groundData->nPoly);

  /*the fitting*/
  fitCheck=mpfit(groundErr,groundData->nPoints,nPar,groundData->par,parStruct,config,(void *)groundData,result);
  if(fitCheck<0){
    fprintf(stderr,"fitCheck %d\n",fitCheck);
    exit(1);
  }

  /*tidy up*/
  if(result){
    TIDY(result->resid);
    TIDY(result->xerror);
    TIDY(result->covar);
    TIDY(result);
  }
  TIDY(config);
  TIDY(parStruct);
  return;
}/*fitPolyPlane*/


/*################################################*/
/*allocate and load relevant data*/

groundDstruct *arrangeGroundData(pCloudStruct **data,int nFiles,double groundBreakElev,double minX,double minY,double maxX,double maxY)
{
  int numb=0;
  uint32_t i=0,j=0;
  groundDstruct *groundD=NULL;   /*ground data structure*/


  /*initialise*/
  if(!(groundD=(groundDstruct *)calloc(1,sizeof(groundDstruct)))){
    fprintf(stderr,"error in groundDstruct structure.\n");
    exit(1);
  }
  groundD->nPoints=0;
  groundD->xUse=NULL;
  groundD->yUse=NULL;
  groundD->zUse=NULL;
  groundD->par=NULL;

  /*count the number of points*/
  for(numb=0;numb<nFiles;numb++){
    for(i=0;i<data[numb]->nPoints;i++){
      /*are we wthin bounds, or have bounds not been defined yet? In the latter case, read everything*/
      if((maxX<-1000000.0)||((data[numb]->x[i]>=minX)&&(data[numb]->x[i]<=maxX)&&(data[numb]->y[i]>=minY)&&(data[numb]->y[i]<=maxY))){
        if(data[numb]->class[i]==2)groundD->nPoints++;
      }/*within bound check*/
    }/*point loop*/
  }/*file loop*/

  /*allocate space*/
  groundD->xUse=dalloc(groundD->nPoints,"ground x",0);
  groundD->yUse=dalloc(groundD->nPoints,"ground x",0);
  groundD->zUse=dalloc(groundD->nPoints,"ground x",0);

  /*copy over data*/
  j=0;
  for(numb=0;numb<nFiles;numb++){
    for(i=0;i<data[numb]->nPoints;i++){
      /*are we wthin bounds, or have bounds not been defined yet? In the latter case, read everything*/
      if((maxX<-1000000.0)||((data[numb]->x[i]>=minX)&&(data[numb]->x[i]<=maxX)&&(data[numb]->y[i]>=minY)&&(data[numb]->y[i]<=maxY))){
        if((data[numb]->class[i]==2)&&(data[numb]->z[i]>=groundBreakElev)){
          groundD->xUse[j]=data[numb]->x[i];
          groundD->yUse[j]=data[numb]->y[i];
          groundD->zUse[j]=data[numb]->z[i];
          j++;
        }/*is point ground*/
      }/*within bound check*/
    }/*point loop*/
  }/*file loop*/

  return(groundD);
}/*arrangeGroundData*/


/*###############################################*/
/*Ground error function*/

int groundErr(int numb, int npar, double *p, double *deviates,double **derivs, void *private)
{
  int i=0,j=0;
  double x=0,y=0,z=0;
  double xBit=0,yBit=0;
  groundDstruct *g=NULL;


  g=(groundDstruct *)private;

  /*point loop*/
  for(i=0;i<numb;i++){
    x=g->xUse[i];
    y=g->yUse[i];

    xBit=yBit=0.0;
    for(j=g->nPoly[0]-1;j>=0;j--)xBit+=p[j]*pow(x,(double)(g->nPoly[0]-j));
    for(j=g->nPoly[1]-1;j>=0;j--)yBit+=p[j+g->nPoly[0]]*pow(y,(double)(g->nPoly[1]-j-1));
    z=xBit+yBit;
    deviates[i]=z-g->zUse[i];

    /*derivatives*/
    /*if(derivs){
      for(j=g->nPoly[0]-1;j>=0;j--)derivs[j][i]=(double)(j+1)*pow(x,(double)(g->nPoly[0]-j));
      for(j=g->nPoly[0]-1;j>=0;j--)derivs[j+g->nPoly[0]][i]=(double)j*pow(y,(double)(g->nPoly[1]-j-1));
    }*/
  }/*point loop*/

  g=NULL;
  return(0);
}/*groundErr*/


/*########################################*/
/*set ground fit parameter controls*/

mp_par *setGroundBounds(int *nPar)
{
  int i=0,j=0,count=0;
  mp_par *parStruct=NULL;

  if(!(parStruct=(mp_par *)calloc(nPar[0]+nPar[1],sizeof(mp_par)))){
    fprintf(stderr,"error in bound structure.\n");
    exit(1);
  }

  count=0;
  for(j=0;j<2;j++){
    for(i=0;i<nPar[j];i++){
      parStruct[count].fixed=0;
      parStruct[count].side=2;
      parStruct[count].limited[0]=0;
      parStruct[count].limits[0]=-3000.0; //pow(-350.0,1.0/(float)(nPar[j]-i));
      parStruct[count].limited[1]=0;
      parStruct[count].limits[1]=3000.0; //pow(350.0,1.0/(float)(nPar[j]-i));
      parStruct[count].step=pow(20.0,1.0/(float)(nPar[j]-i));
      count++;
    }
  }

  return(parStruct);
}/*setGroundBounds*/


/*########################################*/
/*check waveform is usable*/

char checkWaveform(float *wave,uint32_t nBins)
{
  char usable=1;
  uint32_t i=0;

  for(i=0;i<nBins;i++){
    if(isnan(wave[i])){
      usable=0;
      break;
    }
  }
  return(usable);
}/*checkWaveform*/


/*####################################################*/
/*rh metrics*/

float *findRH(float *wave,double *z,int nBins,double gHeight,float rhRes,int *nRH)
{
  int i=0,j=0;
  float cumul=0;
  float totE=0,r=0;
  float *rh=NULL;
  float lastRet=0;
  char *done=NULL;
  char hasEnergy=0;

  /*total energy*/
  totE=0.0;
  for(i=0;i<nBins;i++)totE+=wave[i];

  /*allocate space*/
  *nRH=(int)(100.0/rhRes+1);
  rh=falloc((uint64_t)(*nRH),"rh metrics",0);
  done=challoc((uint64_t)(*nRH),"RH done flag",0);
  for(i=0;i<(*nRH);i++){
    done[i]=0;
    rh[i]=-9999.0;
  }

  /*which way is the waveform arranged?*/
  cumul=0.0;
  if(z[nBins-i]<z[0]){   /*wave is from from top to bottom*/
    for(i=nBins-1;i>=0;i--){
      cumul+=wave[i];
      if(wave[i]>0.0){
        for(j=0;j<(*nRH);j++){
          r=((float)j*(float)(rhRes/100.0))*totE;
          if(r>totE)r=totE;
          if((done[j]==0)&&(cumul>=r)){
            rh[j]=(float)(z[i]-gHeight);
            done[j]=hasEnergy=1;
          }
        }
        lastRet=(float)(z[i]-gHeight);
      }
    }
  }else{
    for(i=0;i<nBins;i++){  /*wave is from bottom to top*/
      cumul+=wave[i];
      if(wave[i]>0.0){
        for(j=0;j<(*nRH);j++){
          r=((float)j*(float)(rhRes/100.0))*totE;
          if(r>totE)r=totE;
          if((done[j]==0)&&(cumul>=r)){
            rh[j]=(float)(z[i]-gHeight);
            done[j]=hasEnergy=1;
          }
        }
        lastRet=(float)(z[i]-gHeight);
      }
    }
  }
  TIDY(done);

  /*in case there was a rounding error for RH100*/
  if(rh[(*nRH)-1]<-9000.0)rh[(*nRH)-1]=lastRet;

  /*was this waveform empty*/
  if(!hasEnergy){
    for(i=0;i<*nRH;i++)rh[i]=0.0;
  }

  return(rh);
}/*findRH*/


/*####################################################*/
/*foliage height diversity from waveform*/

float foliageHeightDiversity(float *wave,int nBins)
{
  int i=0;
  float FHD=0,thresh=0;
  float p=0,total=0;

  /*only if there is a valid wave*/
  if(wave==NULL)return(-100000.0);

  /*determine a threshold*/
  thresh=TOL;

  /*total*/
  total=0.0;
  for(i=0;i<nBins;i++)total+=wave[i];

  FHD=0.0;
  for(i=0;i<nBins;i++){
    if(wave[i]>thresh){
      p=wave[i]/total;
      FHD-=p*log(p);
    }
  }

  return(FHD);
}/*foliageHeightDiversity*/


/*####################################################*/
/*foliage height diversity from histogram*/

float foliageHeightDiversityHist(float *wave,int nBins,float res)
{
  int i=0,histBins=0;
  int bin=0,maxBins=200;
  float FHD=0,thresh=0;
  float max=0,min=0;
  float p=0,total=0;
  float *hist=NULL;

  /*exit if there is not a valid wave*/
  if(wave==NULL)return(0.0);

  /*determine a threshold*/
  thresh=TOL;

  /*find bounds*/
  max=-100000.0;
  min=10000000.0;
  total=0.0;
  for(i=0;i<nBins;i++){
    if(wave[i]>thresh){
      if(wave[i]>max)max=wave[i];
      if(wave[i]<min)min=wave[i];
      total+=1.0;
    }
  }
  histBins=(int)((max-min)/res+1.0);
  /*truncate to prevent daft numbers of computations*/
  if(histBins>maxBins){
    histBins=maxBins;
    res=(max-min)/(float)(histBins-1);
  }

  /*exit if the waveform is blank*/
  if((total<TOL)||(histBins<1))return(0.0);

  /*make histogram*/
  hist=falloc((uint64_t)histBins,"FHD histogram",0);
  for(i=0;i<nBins;i++){
    if(wave[i]>thresh){
      bin=(int)((wave[i]-min)/res);
      if(bin<0)bin=0;
      if(bin>=histBins)bin=histBins-1;
      hist[bin]+=1.0;
    }
  }

  FHD=0.0;
  for(i=0;i<histBins;i++){
    p=hist[i]/total;
    if(p>0.0)FHD-=p*log(p);
  }

  TIDY(hist);
  return(FHD);
}/*foliageHeightDiversityHist*/


/*####################################################*/
/*subtract lowest Gaussian from canopy and correct for attenuation*/

float *subtractGaussFromCan(float *wave,int nBins,float mu,float A,float sig,double *z)
{
  int i=0;
  float tot=0,gap=0;
  float *canWave=NULL;
  float delta=0;

  canWave=falloc((uint64_t)nBins,"canopy only wave",0);

  /*subtract ground*/
  tot=0.0;
  for(i=0;i<nBins;i++){
    if((float)z[i]>=mu){
      canWave[i]=wave[i]-A*(float)gaussian(z[i],(double)sig,(double)mu);
      if(canWave[i]<0.0)canWave[i]=0.0;
      tot+=canWave[i];
    }else  canWave[i]=0.0;
  }

  /*correct for attenuation*/
  gap=1.0;
  for(i=0;i<nBins;i++){
    delta=canWave[i]/tot;
    if(gap>0.0)canWave[i]/=gap;
    gap-=delta;
  }

  return(canWave);
}/*subtractGaussFromCan*/


/*####################################################*/
/*subtract ground from canopy and correct fors attenuation*/

float *canProfile(float *wave,float *ground,int nBins)
{
  int i=0;
  float tot=0,gap=0;
  float *canProf=NULL;

  /*allocate*/
  canProf=falloc((uint64_t)nBins,"canope wave",0);

  /*wave integram for attenuation correction*/
  tot=0.0;
  for(i=0;i<nBins;i++)tot+=wave[i];

  /*calculate*/
  gap=1.0;
  for(i=0;i<nBins;i++){
    if(gap>0.0)canProf[i]=(wave[i]-ground[i])/gap;
    if(canProf[i]<0.0)canProf[i]=0.0;
    gap-=wave[i]/tot;
  }

  return(canProf);
}/*canProfile*/


/*####################################################*/
/*make a canopy only waveform*/

float *subtractGroundFromCan(float *wave,float *ground,int nBins)
{
  int i=0;
  float *canWave=NULL;

  /*allocate*/
  canWave=falloc((uint64_t)nBins,"canope wave",0);

  ///*wave integra; for attenuation correction*/
  //tot=0.0;
  //for(i=0;i<nBins;i++)tot+=wave[i];

  /*calculate*/
  for(i=0;i<nBins;i++){
    canWave[i]=wave[i]-ground[i];
    if(canWave[i]<0.0)canWave[i]=0.0;
  }

  return(canWave);
}/*canProfile*/


/*####################################################*/
/*L moments*/

float *waveLmoments(float *rh,int nRH,float rhRes,int nLm)
{
  int i=0;
  float res=0,r=0;
  float *Lmoments=NULL;

  if(nLm==0)return(Lmoments);
  else if(nLm>4){
    fprintf(stderr,"Not set up to deal with more than 4 L-moments yet\n");
    exit(1);
  }

  res=rhRes/100.0;
  Lmoments=falloc((uint64_t)nLm,"L-moments",0);
  for(i=0;i<nLm;i++)Lmoments[i]=0.0;

  for(i=0;i<nRH;i++){
    r=(float)i*rhRes;
    Lmoments[0]+=rh[i]*res;
    if(nLm>1)Lmoments[1]+=rh[i]*(2.0*r-1)*res;
    if(nLm>2)Lmoments[2]+=rh[i]*(6.0*r*r-6.0*r+1.0)*res;
    if(nLm>3)Lmoments[3]+=rh[i]*(20.0*r*r*r-30.0*r*r+12.0*r-1.0)*res;
  }/*RH metric loop*/

  return(Lmoments);
}/*waveLmoments*/


/*####################################################*/
/*Wenge Ni-Meister's biomass metric, less the a factor*/

float niMetric(float *wave,double *z,int nBins,float res,double gElev,float c)
{
  int i=0;
  float niM=0,height=0;

  niM=0.0;
  for(i=0;i<nBins;i++){
    height=(float)(z[i]-gElev);
    if(height>0.0)niM+=wave[i]*pow(height,c)*res;
    else          niM+=wave[i]*(-1.0)*pow(fabs(height),c)*res;
  }

  return(niM);
}/*niMetric*/


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


/*####################################################*/
/*correct detector drift*/

float *correctDrift(float *wave,int nBins,int noiseBins,denPar *den)
{
  int i=0,sBin=0,eBin=0;
  int iters=0,maxIter=0;
  float xi=0,meanAft=0;
  float *drift=NULL;
  float step=0,diff=0,lastDiff=0;
  float initialDriftGuess(int,int,float *,float,float);
  void removeDrift(float *,float *,int,float,float,float);
  void endPointForDrift(int *,int *,float *,float *,int,int,float,int,float);
  char dir=0;


  /*allocate*/
  drift=falloc((uint64_t)nBins,"drift corrected wave",0);


  /*do we fix detector drift?*/
  if(den->corrDrift){
    if(!den->varDrift){  /*fixed drift factor*/
      xi=den->fixedDrift;
      removeDrift(drift,wave,nBins,den->meanN,den->res,xi);
    }else{               /*variable drift factor*/
      /*find end point and mean after noise level*/
      endPointForDrift(&sBin,&eBin,&meanAft,wave,nBins,noiseBins,den->threshScale,den->minWidth,den->meanN);

      /*initial guess*/
      xi=initialDriftGuess(sBin,eBin,wave,den->meanN,meanAft);
      if(den->meanN>meanAft)dir=1;
      else                  dir=-1;
      step=(xi>0.01)?xi:0.01;

      /*iterate over factors*/
      maxIter=1000;
      iters=0;
      diff=fabs(den->meanN-meanAft);
      lastDiff=1000.0;
      while((fabs(diff)>0.01)&&(fabs(diff-lastDiff)>0.000000001)&&(iters<maxIter)){
        /*try this value for xi*/
        removeDrift(drift,wave,nBins,den->meanN,den->res,xi);

        /*check for success*/
        meanAft=0.0;
        for(i=nBins-1;i>=eBin;i--)meanAft+=drift[i];
        meanAft/=(int)(nBins-eBin);

        /*adjust xi*/
        if(den->meanN>meanAft){
          if(dir<0)step/=2.0;
          xi+=step;
        }else{
          if(dir>0)step/=2.0;
          xi-=step;
        }
        lastDiff=diff;
        diff=fabs(den->meanN-meanAft);
        iters++;
      }
    }
  }else{   /*if not fixing, copy old wave*/
    for(i=0;i<nBins;i++)drift[i]=wave[i];
  }

  return(drift);
}/*correctDrift*/


/*####################################################*/
/*intial guess of the drift factor*/

float initialDriftGuess(int sBin,int eBin,float *wave,float meanN,float meanAft)
{
  int i=0;
  float tot=0;

  /*total energy above mean noise before*/
  tot=0.0;
  for(i=sBin;i<=eBin;i++)tot+=wave[i]-meanN;

  return((meanN-meanAft)/tot);
}/*initialDriftGuess*/


/*####################################################*/
/*remove drift with known factor*/

void removeDrift(float *drift,float *wave,int nBins,float meanN,float res,float xi)
{
  int i=0;
  float cumul=0.0;

  cumul=0.0;
  for(i=0;i<nBins;i++){
    drift[i]=wave[i]+cumul*xi;
    cumul+=(drift[i]-meanN)*res;
  }

  return;
}/*removeDrift*/


/*####################################################*/
/*find end point and after mean noise to use for detector drift*/

void endPointForDrift(int *sBin,int *eBin,float *meanAft,float *wave,int nBins,int noiseBins,float threshScale,int minWidth,float meanN)
{
  int i=0,nIn=0;
  float thresh=0,stdev=0;

  stdev=0.0;
  for(i=0;i<noiseBins;i++)stdev+=(wave[i]-meanN)*(wave[i]-meanN);
  stdev=sqrt(stdev/(float)noiseBins);

  /*find start*/
  thresh=meanN+stdev*threshScale;
  nIn=0;
  *sBin=-1;
  for(i=0;i<nBins;i++){
    if((wave[i]-meanN)>=thresh)nIn++;
    else                       nIn=0;
    if(nIn>=minWidth){
      for(;i>=0;i--){
        if(wave[i]<=meanN){
          *sBin=i;
          break;
        }
      }
      if(*sBin<0)*sBin=0;
      break;
    }
  }


  /*find stats at end*/
  *meanAft=stdev=0.0;
  for(i=nBins-noiseBins;i<nBins;i++)*meanAft+=wave[i];
  *meanAft/=(float)noiseBins;
  for(i=nBins-noiseBins;i<nBins;i++)stdev+=(wave[i]-*meanAft)*(wave[i]-*meanAft);
  stdev=sqrt(stdev/(float)noiseBins);

  /*find end*/
  thresh=*meanAft+stdev*threshScale;
  nIn=0;
  *eBin=-1;
  for(i=nBins-1;i>=0;i--){
    if((wave[i]-(*meanAft))>=thresh)nIn++;
    else                            nIn=0;
    if(nIn>=minWidth){
      for(;i<nBins;i++){
        if(wave[i]<=meanN){
          *eBin=i;
          break;
        }
      }
      if(*eBin<0)*eBin=nBins-1;
      break;
    }
  }

  return;
}/*startEndPoints*/

/*the end*/
/*################################################*/

