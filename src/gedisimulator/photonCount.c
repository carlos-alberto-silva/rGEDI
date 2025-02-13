#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "hdf5.h"
#include "libLasRead.h"
#include "libLasProcess.h"
#include "libLidarHDF.h"
#include "libOctree.h"
#include "tools.h"
#include "gediIO.h"
#include "gediNoise.h"
#include "photonCount.h"
#include "gsl/gsl_fft_complex.h"

//#define DEBUG


/*##############################*/
/*# Generates photon count from#*/
/*# simulated GEDI waveforms   #*/
/*# 2019 svenhancock@gmail.com #*/
/*##############################*/

/*#######################################*/
/*# Copyright 2015-2019, Steven Hancock #*/
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
#define TOL 0.0000001


/*####################################################*/
/*turn a waveform in to a photon-count pseudo-wave*/

float *uncompressPhotons(float *wave,dataStruct *data,photonStruct *photonCount,noisePar *noise,gediIOstruct *gediIO)
{
  float *photWave=NULL;
  float *corrWave=NULL,*tempSmoo=NULL;
  float *crossCorrelateWaves(float *,float,int,pulseStruct *,float);
  float *crossCorrelateTime(float *,float,int,pulseStruct *,float,float);
  float *applyHannFilter(float *,int,pulseStruct *);
  void writePCLwaves(dataStruct *,float *,float *,float *);
  pulseStruct *setHannFilter(float,float);


  #ifdef DEBUG
  int i=0;
  static int c=0;
  #endif

  /*do we have a usable pulse?*/
  if(gediIO->pulse==NULL){
    fprintf(stderr,"No pulse. Cannot use PCL\n");
    exit(1);
  }


  /*first perform photon counting, if needed*/
  if(gediIO->pclPhoton)photWave=countWaveform(wave,data,photonCount,gediIO->den,noise);
  else                 photWave=wave;

  /*perform cross-correlation*/
  corrWave=crossCorrelateTime(photWave,data->res,data->nBins,gediIO->pulse,gediIO->pRes,gediIO->pclSwidth);
  //corrWave=crossCorrelateWaves(photWave,data->res,data->nBins,gediIO->pulse,gediIO->pRes);

  /*Hann filter if needed*/
  if(gediIO->hannWidth>TOL){
    if(gediIO->hannFilt==NULL)gediIO->hannFilt=setHannFilter(gediIO->hannWidth,data->res);
    tempSmoo=applyHannFilter(corrWave,data->nBins,gediIO->hannFilt);
    TIDY(corrWave);
    corrWave=tempSmoo;
    tempSmoo=NULL;
  }

  #ifdef DEBUG
  for(i=0;i<data->nBins;i++)fprintf(stdout,"%d %f %f %f %f debug2\n",c,data->z[i],wave[i],photWave[i],corrWave[i]);
  c++;
  #endif
  if(gediIO->writePcl)writePCLwaves(data,wave,photWave,corrWave);


  /*tidy up*/
  if(photWave!=wave){
    TIDY(photWave);
  }else photWave=NULL;

  return(corrWave);
}/*uncompressPhotons*/


/*####################################################*/
/*apply a Hann filter*/

float *applyHannFilter(float *corrWave,int nBins,pulseStruct *hannFilt)
{
  int i=0,j=0,bin=0;
  float contN=0;
  float *smoo=NULL;

  /*allocate space*/
  smoo=falloc(nBins,"Hann filtered wave",0);

  /*loop over waveform bins*/
  for(i=0;i<nBins;i++){
    smoo[i]=contN=0.0;

    /*loop over pulse and smooth*/
    for(j=0;j<hannFilt->nBins;j++){
      bin=i+(j-hannFilt->centBin);   /*Hann filter resolution is matched to waveform*/

      if((bin>=0)&&(bin<nBins)){
        smoo[i]+=corrWave[bin]*hannFilt->y[j];
        contN+=hannFilt->y[j];
      }
    }
    if(contN>0.0)smoo[i]/=contN;
  }

  return(smoo);
}/*applyHannFilter*/


/*####################################################*/
/*set the Hann filter structure*/

pulseStruct *setHannFilter(float sWidth,float res)
{
  int i=0;
  pulseStruct *hannFilt=NULL;

  /*allocate space*/
  if(!(hannFilt=(pulseStruct *)calloc(1,sizeof(pulseStruct)))){
    fprintf(stderr,"error Hann filter allocation.\n");
    exit(1);
  }

  /*determine number of bins*/
  hannFilt->nBins=((int)(sWidth/res+1)/2)*2+1;  /*force to be odd*/
  hannFilt->centBin=hannFilt->nBins/2;

  /*allocate array space*/
  hannFilt->x=falloc(hannFilt->nBins,"Hann filter x",0);
  hannFilt->y=falloc(hannFilt->nBins,"Hann filter y",0);

  /*set values*/
  for(i=0;i<hannFilt->nBins;i++){
    hannFilt->x[i]=(float)(i-hannFilt->centBin)*res;
    hannFilt->y[i]=pow(cos(M_PI*(float)(i-hannFilt->centBin)/(float)hannFilt->nBins),2);
  }

  return(hannFilt);
}/*setHannFilter*/


/*####################################################*/
/*write PCL waveforms*/

void writePCLwaves(dataStruct *data,float *wave,float *photWave,float *corrWave)
{
  int i=0;
  static int c=0;
  char waveNamen[250];
  FILE *opoo=NULL;

  /*name and open file*/
  sprintf(waveNamen,"pclWave.%s.count.%d.txt",data->waveID,c);
  if((opoo=fopen(waveNamen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",waveNamen);
    exit(1);
  }

  /*write data*/
  fprintf(opoo,"# 1 z, 2 waveform, 3 photonWave, 4 correlated, 5 original\n");
  for(i=0;i<data->nBins;i++)fprintf(opoo,"%f %.10f %f %.10f %.10f\n",data->z[i],wave[i],photWave[i],corrWave[i],data->wave[0][i]);

  /*close up*/
  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"PCL waves written to %s\n",waveNamen);
  c++;

  return;
}/*writePCLwaves*/


/*####################################################*/
/*resample a pulse for PCL*/

void resamplePclPulse(pulseStruct *pulse,float res,float pRes)
{
  int i=0,*contN=NULL;
  int bin=0;

  /*allocate space and zero*/
  pulse->rBins=(int)((float)pulse->nBins*pRes/res);
  pulse->resamp=falloc(pulse->rBins,"",0);
  contN=ialloc(pulse->rBins,"",0);
  for(i=0;i<pulse->rBins;i++){
    pulse->resamp[i]=0.0;
    contN[i]=0;
  }


  /*resample pulse*/
  for(i=0;i<pulse->nBins;i++){
    bin=(int)floor((float)i*pRes/res);
    if((bin>=0)&&(bin<pulse->rBins)){
      pulse->resamp[bin]+=pulse->y[i];
      contN[bin]++;
    }
  }

  /*normalise resampled*/
  for(i=0;i<pulse->rBins;i++){
    if(contN[i]>0)pulse->resamp[i]/=(float)contN[i];
  }
  TIDY(contN);
  pulse->rCent=(int)floor((float)pulse->centBin*pRes/res);

  return;
}/*resamplePclPulse*/


/*####################################################*/
/*perform a cross-correlation in the time domain*/

float *crossCorrelateTime(float *photWave,float res,int nBins,pulseStruct *pulse,float pRes,float pclSwidth)
{
  int i=0,j=0,bin=0;
  int thisCont=0;
  float *compCorr=NULL;
  float *smooWave=NULL;
  float meanP=0,meanW=0;
  float stdevP=0,stdevW=0;
  void resamplePclPulse(pulseStruct *,float,float);

  /*allocate space*/
  compCorr=falloc(nBins,"compCorr",0);

  /*allocate resampled pulse if needed*/
  if(pulse->resamp==NULL)resamplePclPulse(pulse,res,pRes);

  /*if not already done, smooth the pulse*/
  if(pulse->pclSmoo==NULL){
    if(pclSwidth>TOL)pulse->pclSmoo=smooth(pclSwidth,pulse->nBins,pulse->resamp,res);
    else             pulse->pclSmoo=pulse->resamp;
  }

  /*smooth waveform if needed*/
  if(pclSwidth>TOL)smooWave=smooth(pclSwidth,nBins,photWave,res);
  else             smooWave=photWave;

  /*find the average of the pulse*/
  meanP=0.0;
  for(i=0;i<pulse->rBins;i++)meanP+=pulse->pclSmoo[i];
  meanP/=(float)pulse->rBins;
  //meanP=singleMedian(pulse->pclSmoo,pulse->rBins);

  /*find the stdev of the pulse*/
  stdevP=0.0;
  for(i=0;i<pulse->rBins;i++)stdevP+=(pulse->pclSmoo[i]-meanP)*(pulse->pclSmoo[i]-meanP);
  stdevP=sqrt(stdevP/(float)pulse->rBins);

  /*find the average of the wave*/
  meanW=0.0;
  for(i=0;i<nBins;i++)meanW+=smooWave[i];
  meanW/=(float)nBins;
  //meanW=singleMedian(smooWave,nBins);

  /*find the stdev of the wave*/
  stdevW=0.0;
  for(i=0;i<nBins;i++)stdevW+=(smooWave[i]-meanW)*(smooWave[i]-meanW);
  stdevW=sqrt(stdevW/(float)nBins);

  /*loop over bins in time domain*/
  //for(i=pulse->rCent;i<(nBins-pulse->rCent);i++){ /*step in by the cent bins of the pulse*/
  for(i=0;i<nBins;i++){
    compCorr[i]=0.0;
    thisCont=0;

    /*loop over pulse to convolve*/
    for(j=0;j<pulse->rBins;j++){
      bin=i+j-pulse->rCent;  /*bin on the pulse*/

      /*are we within the pulse array?*/
      if((bin>=0)&&(bin<nBins)){
        compCorr[i]+=(smooWave[bin]-meanW)*(pulse->pclSmoo[pulse->rBins-(j+1)]-meanP)/(stdevP*stdevW);
        //compCorr[i]+=smooWave[bin]*pulse->pclSmoo[pulse->rBins-(j+1)]/(stdevP*stdevW);
        thisCont++;
      }
    }/*pulse bin loop*/

    /*normalise*/
    //if(thisCont>0)compCorr[i]/=(float)thisCont;
    compCorr[i]/=(float)pulse->rBins;
  }/*wave bin loop*/

  if(smooWave!=photWave)TIDY(smooWave);

  return(compCorr);
}/*crossCorrelateTime*/


/*####################################################*/
/*perform a cross-correlation using Fourier*/

float *crossCorrelateWaves(float *photWave,float res,int nBins,pulseStruct *pulse,float pRes)
{
  int i=0,bin=0;
  int numb=0;
  int *contN=NULL;
  float *corrWave=NULL;
  double *compPulse=NULL;
  double *compWave=NULL;
  double *compCorr=NULL;
  int gsl_fft_complex_radix2_forward(gsl_complex_packed_array,size_t,size_t);
  int gsl_fft_complex_radix2_backward(gsl_complex_packed_array, size_t,size_t);

  /*FFT requires that array is a power of 2 long*/
  numb=pow(2.0,(float)(int)(log((double)nBins)/log(2.0)+1.0));

  /*allocate space for resampled complex pulse*/
  compPulse=dalloc(2*numb,"complex pulse",0);
  contN=ialloc(numb,"contribution counter",0);
  for(i=0;i<numb;i++){
    compPulse[2*i]=compPulse[2*i+1]=0.0;
    contN[i]=0;
  }

  /*split pulse over 0*/
  /*for(i=pulse->centBin;i<pulse->nBins;i++){
    bin=(int)((float)(i-pulse->centBin)*pRes/res+0.5);
    if((bin<0)||(bin>=numb)){
      fprintf(stderr,"Out of bounds when splitting pulse, %d of %d\n",bin,numb);
      continue;
    }
    compPulse[2*bin]+=(double)pulse->y[i];
    contN[bin]++;
  }
  for(i=0;i<pulse->centBin;i++){
    bin=numb-((int)((float)(pulse->centBin-i)*pRes/res+0.5)+1);
    if((bin<0)||(bin>=numb)){
      fprintf(stderr,"Out of bounds when splitting pulse, %d of %d\n",bin,numb);
      continue;
    }
    compPulse[2*bin]+=(double)pulse->y[i];
    contN[bin]++;
  }*/


  /*resample pulse to match waveform*/
  for(i=0;i<pulse->nBins;i++){
    bin=(int)((float)i*pRes/res+0.5);
    if((bin<0)||(bin>=numb))continue;
    compPulse[2*bin]+=(double)pulse->y[i];
    contN[bin]++;
  }

  /*normalise*/
  for(i=0;i<numb;i++){
    if(contN[i]>0)compPulse[2*i]/=(float)contN[i];
  }
  TIDY(contN);

  /*median of pulse*/
  //meanP=singleMedian(pulse->y,pulse->nBins);

  /*make waveform complex*/
  compWave=dalloc(2*numb,"complex wave",0);
  for(i=0;i<nBins;i++){
    compWave[2*i]=(double)photWave[i];
    compWave[2*i+1]=0.0;  /*imaginary part*/
  }
  //meanW=singleMedian(photWave,nBins);
  for(i=2*nBins;i<2*numb;i++)compWave[i]=0.0;

  /*subtract means*/
  /*for(i=numb-1;i>=0;i--){
    compWave[2*i]-=meanW;
    compPulse[2*i]-=meanP;
  }*/

  /*remove assymmetry of signal*/
  //removeAsymmetryPCL(compWave,numb);

  /*fourier transform both*/
  gsl_fft_complex_radix2_forward((gsl_complex_packed_array)compPulse,1,numb);
  gsl_fft_complex_radix2_forward((gsl_complex_packed_array)compWave,1,numb);

  /*correlate*/
  compCorr=dalloc(2*numb,"complex correlation",0);
  for(i=0;i<numb;i++){
    compCorr[2*i]=compPulse[2*i]*compWave[2*i]+compPulse[2*i+1]*compWave[2*i+1];
    compCorr[2*i+1]=compPulse[2*i+1]*compWave[2*i]-compPulse[2*i]*compWave[2*i+1];
  }

  /*inverse fourier*/
  gsl_fft_complex_radix2_backward((gsl_complex_packed_array)compCorr,1,numb);

  /*make real*/
  corrWave=falloc(nBins,"correlated wave",0);
  for(i=0;i<nBins;i++){
    corrWave[i]=(float)compCorr[2*i]; //sqrt(compCorr[2*i]*compCorr[2*i]+compCorr[2*i+1]*compCorr[2*i+1]);
  }

  /*tidy up*/
  TIDY(compPulse);
  TIDY(compWave);
  TIDY(compCorr);

  return(corrWave);
}/*crossCorrelateWaves*/


/*####################################################*/
/*truncate the assymmetry of a waveform for PCL*/

void removeAsymmetryPCL(float *wave,int numb)
{
  int i=0;
  float medianW=0;

  /*find the median*/
  medianW=singleMedian(wave,numb);

  /*find start point*/
  for(i=0;i<numb;i++){
    /*if less than mean, set to mean*/
    if(wave[i]>=medianW){
      for(;i>=0;i--)wave[i]=medianW;
      break;
    }
  }

  /*find end point*/
  for(i=numb-1;i>=0;i--){
    /*if less then mean, set to mean*/
    if(wave[i]>=medianW){
      for(;i<numb;i++)wave[i]=medianW;
      break;
    }
  }

  return;
}/*removeAsymmetryPCL*/


/*####################################################*/
/*produce a photon-counting pseudo-waveform*/

float *countWaveform(float *denoised,dataStruct *data,photonStruct *photonCount,denPar *den,noisePar *noise)
{
  int i=0,nPhot=0;
  int bin=0;
  float minI=0;
  float *temp=NULL;
  float **phots=NULL;
  void applyShotNoise(float *,int);

  #ifdef DEBUG
  static int count=0;
  fprintf(stdout,"Photon counting\n");
  fflush(stdout);
  #endif

  /*set minimum to zero*/
  minI=10000.0;
  for(i=0;i<data->nBins;i++){
    if(denoised[i]<minI)minI=denoised[i];
  }
  temp=falloc(data->nBins,"temp pcl waveform",0);
  for(i=0;i<data->nBins;i++)temp[i]=denoised[i]-minI;

  /*set window size for background noise photons*/
  photonCount->H=data->res*(float)data->nBins*2.0;  /*two way distance*/

  /*extract photon coords along with their flags*/
  phots=countPhotons(temp,data,photonCount,&nPhot,den,noise,1);

  /*reset temp array*/
  for(i=0;i<data->nBins;i++)temp[i]=0.0;

  /*bin up in to new wave*/
  for(i=0;i<nPhot;i++){
    bin=(int)(((float)data->z[0]-phots[0][i])/data->res);
    if((bin>=0)&&(bin<data->nBins))temp[bin]+=1.0;
  }
  TTIDY((void **)phots,3);

  /*apply shot noise if neeed*/
  if(noise->shotNoise)applyShotNoise(temp,data->nBins);

  #ifdef DEBUG
  for(i=0;i<data->nBins;i++)fprintf(stdout,"%d %f %f %f ta\n",count,data->z[i],temp[i],data->wave[0][i]);
  fflush(stdout);
  fprintf(stderr,"count %d\n",count);fflush(stderr);
  count++;
  #endif

  return(temp);
}/*countWaveform*/


/*####################################################*/
/*apply shot noise*/

void applyShotNoise(float *temp,int nBins)
{
  int i=0,bin=0;
  int lostPhots=0;
  int nCulled=0,totIn=0;
  float shotSig=0;
  float shotNoise=0;
  float photThresh=0;
  float *mask=NULL;

  /*loop over bins*/
  for(i=0;i<nBins;i++){
    /*only if there are photons*/
    if(temp[i]>TOL){
      /*set sigma*/
      shotSig=sqrt(temp[i]);

      /*draw Gaussian random number*/
      shotNoise=(float)round(GaussNoise(NULL)*shotSig);

      /*count truncated negative*/
      temp[i]+=shotNoise;
      if(temp[i]<0.0){
        lostPhots-=(int)(temp[i]+0.5);
        temp[i]=0.0;
      }

      totIn+=(int)(temp[i]+0.5);
    }/*return check*/
  }/*bin loop*/

  /*redeploy negative numbers*/
  if(lostPhots>0){
    if(lostPhots>=totIn){  /*have we lost so many that the signal is empty?*/
      for(i=0;i<nBins;i++)temp[i]=0.0;
    }else{                 /*otherwise redeploy negative numbers*/
      mask=falloc(nBins,"temporary mask",0);
      nCulled=0;

      while(nCulled<lostPhots){
        /*make the mask array*/
        for(i=0;i<nBins;i++){
          if(temp[i]>0.0)mask[i]=1.0;
          else           mask[i]=0.0;
        }

        photThresh=(float)rand()/(float)RAND_MAX;
        bin=(int)pickArrayElement(photThresh,mask,nBins,0);
        if(temp[bin]>0.0){
          temp[bin]-=1.0;
          nCulled++;
        }
      }
      TIDY(mask);
    }/*redploy lost photons*/
  }/*lost photon check*/

  return;
}/*applyShotNoise*/


/*####################################################*/
/*produce photon counting photons*/

float **countPhotons(float *denoised,dataStruct *data,photonStruct *photonCount,int *nPhot,denPar *den,noisePar *noise,char pcl)
{
  int i=0;
  int nPhotons=0,nNoise=0;
  int setNumberNoise(float,float,float);
  float **phots=NULL;  /*arrray with z, isSignal and isGround*/
  float photThresh=0,d=0,thisZ=0;
  float photonNoiseIntensity(float);
  float minZ=0,maxZ=0;
  float *thisGr=NULL;
  float *wave=NULL;
  float *adjustPhotonProb(float *,dataStruct *,denPar *,noisePar *,int,photonStruct *);
  void knockOffNegativeWaves(float *,dataStruct *);
  void adjustTotalPhotRate(photonStruct *,float);
  void setPhotonGround(float *,float *,float,double,float *,float *,double *,int);
  char testPhotonGround(dataStruct *,float);


  /*remove negatives if needed*/
  knockOffNegativeWaves(denoised,data);

  /*rescale waveform for reflectance*/
  wave=adjustPhotonProb(denoised,data,den,noise,data->useType,photonCount);

  /*do we need to set up the probability array?*/
  if(photonCount->prob==NULL){
    if(photonCount->reflDiff)adjustTotalPhotRate(photonCount,data->cov);
    setPhotonProb(photonCount);
  }

  /*choose a number of signal photons to use*/
  photThresh=(float)rand()/(float)RAND_MAX;
  nPhotons=(int)pickArrayElement(photThresh,photonCount->prob,photonCount->pBins,1);

  /*generate noise photons*/
  if(fabs(photonCount->noise_mult)>TOL)nNoise=setNumberNoise(data->cov,photonCount->noise_mult,photonCount->H);
  else                                 nNoise=0;
  *nPhot=nPhotons+nNoise;

  /*allocate space*/
  phots=fFalloc(3,"photon coordinates",0);
  for(i=0;i<3;i++)phots[i]=falloc(*nPhot,"photon coordinates",i+1);

  /*generate signal photons*/
  for(i=0;i<nPhotons;i++){
    /*pick a point along the waveform*/
    photThresh=(float)rand()/(float)RAND_MAX;
    d=pickArrayElement(photThresh,wave,data->nBins,1);

    phots[0][i]=(float)data->z[0]-d*data->res;   /*determine range*/
    phots[1][i]=1.0;                             /*is signal*/
    phots[2][i]=(float)testPhotonGround(data,d); /*is this ground or canopy?*/
  }/*signal photon loop*/

  /*Noise*/
  /*set bounds of search window*/
  if(!pcl){  /*for ICESat-2 mode, set bounds relative to ground*/
    if(data->ground)thisGr=data->ground[data->useType];
    else            thisGr=NULL;
    setPhotonGround(&minZ,&maxZ,photonCount->H,data->gElev,data->wave[data->useType],thisGr,data->z,data->nBins);
    thisGr=NULL;
  }else{     /*for PCL mode, set bounds from data*/
    if(data->z[0]>data->z[data->nBins-1]){
      maxZ=data->z[0];
      minZ=data->z[data->nBins-1];
    }else{
      minZ=data->z[0];
      maxZ=data->z[data->nBins-1];
    }
  }

  /*add noise photons*/
  #ifdef DEBUG
  fprintf(stdout,"Adding %d noise over %f signal %d\n",nNoise,photonCount->H,nPhotons);
  #endif

  for(i=0;i<nNoise;i++){
    thisZ=(maxZ-minZ)*((float)rand()/(float)RAND_MAX)+minZ;

    phots[0][i+nPhotons]=thisZ;   /*determine range*/
    phots[1][i+nPhotons]=0.0;     /*is signal*/
    phots[2][i+nPhotons]=0.0;     /*is this ground or canopy?*/
  }/*noise loop*/

  /*clear out prob if needed*/
  if(photonCount->reflDiff)TIDY(photonCount->prob);

  /*tidy up*/
  if(wave!=denoised){
    TIDY(wave);
  }else wave=NULL;

  return(phots);
}/*countPhotons*/


/*####################################################*/
/*select photons for photon counting*/

void photonCountCloud(float *denoised,dataStruct *data,photonStruct *photonCount,char *outRoot,int numb,denPar *den,noisePar *noise)
{
  int i=0,nRH=0,nPhot=0;
  float *rhReal=NULL,noiseInt=0;
  float photonNoiseIntensity(float);
  float **phots=NULL;

  /*open file if needed*/
  if(photonCount->opoo==NULL){
    sprintf(photonCount->outNamen,"%s.pts",outRoot);
    if((photonCount->opoo=fopen(photonCount->outNamen,"w"))==NULL){
      fprintf(stderr,"Error opening input file %s\n",photonCount->outNamen);
      exit(1);
    }
    fprintf(photonCount->opoo,"# 1 X, 2 Y, 3 Z, 4 minht, 5 WFGroundZ, 6 RH50, 7 RH60, 8 RH75, 9 RH90, 10 RH95, 11 CanopyZ, 12 canopycover, 13 shot#, 14 photon#, 15 iteration#, 16 refdem, 17 noiseInt, 18 signal, 19 ground\n");
  }

  /*generate photons*/
  phots=countPhotons(denoised,data,photonCount,&nPhot,den,noise,0);

  /*get true RH metrics*/
  rhReal=findRH(data->wave[data->useType],data->z,data->nBins,data->gElev,5.0,&nRH);

  /*determine reflectance for noise intensity*/
  noiseInt=photonNoiseIntensity(data->cov);

  /*write out photons*/
  for(i=0;i<nPhot;i++){
    fprintf(photonCount->opoo,"%.2f %.2f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %d %d 1 %.3f %.3f %d %d\n",data->lon,data->lat,phots[0][i],rhReal[0],data->gElev,rhReal[10],rhReal[12],rhReal[15],rhReal[18],rhReal[19],rhReal[nRH-1],data->cov,numb,i,data->gElev,noiseInt,(int)phots[1][i],(int)phots[2][i]);
  }/*mutiple photon loop*/

  TTIDY((void **)phots,3);
  TIDY(rhReal);
  return;
}/*photonCountCloud*/


/*########################################################*/
/*is this photon from the ground?*/

char testPhotonGround(dataStruct *data,float d)
{
  int bin=0;
  float gFrac=0;
  char isGround=0;

  if(data->ground==NULL)isGround=0;
  else{
    bin=(int)d;
    if(data->wave[data->useType][bin]>0.0){
      gFrac=data->ground[data->useType][bin]/data->wave[data->useType][bin];
      isGround=(gFrac>=0.5)?1:0;
    }else isGround=0;
  }

  return(isGround);
}/*testPhotonGround*/


/*########################################################*/
/*remove negative values for photon counting*/

void knockOffNegativeWaves(float *denoised,dataStruct *data)
{
  int i=0;
  float min=0;

  /*find the minimum*/
  min=100000.0;
  for(i=0;i<data->nBins;i++){
    if(denoised[i]<min)min=denoised[i];
  }


  /*translate waves if needed*/
  if(min<0.0){
    for(i=0;i<data->nBins;i++){
      denoised[i]-=min;
      if(data->ground)data->ground[data->useType][i]-=min;
    }
  }

  return;
}/*knockOffNegativeWaves*/


/*########################################################*/
/*adjust waveform to account for refl difference*/

float *adjustPhotonProb(float *denoised,dataStruct *data,denPar *den,noisePar *noise,int numb,photonStruct *phot)
{
  int i=0;
  float tot=0;
  float *wave=NULL,*canopy=NULL;
  float *smooGr=NULL,*smooCan=NULL;

  /*do we have a ground*/
  if(data->ground==NULL){  /*no ground*/
    wave=denoised;
  }else{   /*there is a ground*/
    /*is any adjustment needed*/
    if(fabs(1.0-phot->rhoVrhoG)<TOL)wave=denoised;
    else{  /*if it is needed, do it*/
      if(den->varNoise||noise->linkNoise){
        fprintf(stderr,"Not able to readjust denoised waveforms just yet\n");
        exit(1);
      }else{
        /*find canopy portion*/
        canopy=falloc((uint64_t)data->nBins,"canopy wave",0);
        for(i=0;i<data->nBins;i++)canopy[i]=data->wave[numb][i]-data->ground[numb][i];

        /*smooth ground if needed*/
        if((den->sWidth>0.001)||(den->psWidth>0.001)||(den->msWidth>0.001)){
          /*smooth if needed*/
          smooCan=processFloWave(canopy,data->nBins,den,1.0);
          smooGr=processFloWave(data->ground[numb],data->nBins,den,1.0);
        }else{
          smooCan=canopy;
          smooGr=data->ground[numb];
        }

        /*add up and normalise*/
        wave=falloc((uint64_t)data->nBins,"rescaled erflectance wave",0);
        tot=0.0;
        for(i=0;i<data->nBins;i++){
          wave[i]=smooCan[i]*phot->nPhotC+smooGr[i]*phot->nPhotG;
          tot+=wave[i];
        }
        if(fabs(1.0-tot)>TOL){
          for(i=0;i<data->nBins;i++)wave[i]/=tot;
        }
      }
    }

    /*tidy up*/
    if(smooCan!=canopy){
      TIDY(smooCan);
    }
    TIDY(canopy);
    if(smooGr!=data->ground[numb]){
    TIDY(smooGr);
    }
  }/*is there a ground check?*/

  return(wave);
}/*adjustPhotonProb*/


/*########################################################*/
/*determine ground for photon counting noise*/

void setPhotonGround(float *minZ,float *maxZ,float H,double gElev,float *wave,float *ground,double *z,int nBins)
{
  int i=0;
  float tot=0;
  float CofG=0;

  /*do we have any ground energy?*/
  tot=0.0;
  if(ground){
    for(i=0;i<nBins;i++)tot+=ground[i];
  }

  /*if no ground return, use CofG*/
  if(tot<=TOL){
    tot=0.0;
    for(i=0;i<nBins;i++){
      tot+=wave[i];
      CofG+=(float)z[i]*wave[i];
    }
    CofG/=tot;
    *maxZ=(float)gElev+H/4.0;  /*divided by 4 as H is the 2 way distance*/
    *minZ=(float)gElev-H/4.0;
  }else CofG=gElev;  /*otherwise use the ground elevation*/
  *maxZ=(float)CofG+H/4.0;
  *minZ=(float)CofG-H/4.0;

  return;
}/*setPhotonGround*/


/*########################################################*/
/*set number of noise photons for photon-counting*/

int setNumberNoise(float cov,float noise_mult,float H)
{
  int nNoise=0;
  /*float refl=0; OLD, for varying reflectance*/
  float photThresh=0;
  float noiseRate=0;
  float c=299792458.0;
  photonStruct tempPhot;

  /*is any noise being added?*/
  if(noise_mult>TOL){   /*noise defined as a rate*/
    /*surface reflectance*/
    /*if((cov<0.0)||(cov>1.0))cov=0.5;   OLD, for varying ground reflectance
    refl=cov*0.15+(1.0-cov)*0.22;*/  /*assuming ground and canopy reflectance values*/

    /*noise rate in photons per window*/
    noiseRate=noise_mult*1000000.0; //*refl     used to be *refl to account for xchanging surface reflectance
    /*nNoise=(int)(50.0*(H/c)*noiseRate+0.5);  This is to match Kaitlin's matlab code, but unsure where the 50 came from*/
    tempPhot.designval=(H/c)*noiseRate;

    /*pick from a Poisson*/
    setPhotonProb(&tempPhot);
    photThresh=(float)rand()/(float)RAND_MAX;
    nNoise=(int)pickArrayElement(photThresh,tempPhot.prob,tempPhot.pBins,1);
    TIDY(tempPhot.prob);
  }else if(noise_mult<(-1.0*TOL)){  /*noise defined as mean number of photons*/
    /*mean numb3er of photons*/
    tempPhot.designval=-1.0*noise_mult;   /*mean number of noise photons stored as a negative*/
    /*pick from a Poisson*/
    setPhotonProb(&tempPhot);
    photThresh=(float)rand()/(float)RAND_MAX;
    nNoise=(int)pickArrayElement(photThresh,tempPhot.prob,tempPhot.pBins,1);
    TIDY(tempPhot.prob);
  }else nNoise=0;

  return(nNoise);
}/*setNumberNoise*/


/*########################################################*/
/*determine solar background intensity for photon count*/

float photonNoiseIntensity(float cov)
{
  float noiseInt=0;

  /*it should be a backwards average of this, using the FOV (1D slice thereof)*/
  /*noiseInt=1.0+(1.0-cov)*1.5;*/

  /*this is what is in the ICEsat-2 code, just about*/
  if(cov<0.25)noiseInt=2.5;
  else if(cov>0.75)noiseInt=1.0;
  else noiseInt=1.5+(float)rand()/(float)RAND_MAX;

  return(noiseInt);
}/*photonNoiseIntensity*/


/*####################################################*/
/*adjust photon rate for changing reflectance*/

void adjustTotalPhotRate(photonStruct *photonCount,float cov)
{
  if(cov>=0.0)photonCount->designval=cov*photonCount->nPhotC+(1.0-cov)*photonCount->nPhotG;

  return;
}/*adjustTotalPhotRate*/


/*####################################################*/
/*set photon probability*/

void setPhotonProb(photonStruct *photonCount)
{
  int i=0;
  float y=0;
  float poissonPDF(float,float);

  /*determine number of steps*/
  photonCount->pBins=0;
  do{
    y=poissonPDF((float)photonCount->pBins,photonCount->designval);
    photonCount->pBins++;
  }while((y>0.00001)||((float)photonCount->pBins<photonCount->designval));  /*at least to mean and then to low prob*/

  /*allocate space*/
  photonCount->prob=falloc((uint64_t)photonCount->pBins,"photon prob",0);

  /*set probabilities*/
  for(i=0;i<photonCount->pBins;i++)photonCount->prob[i]=poissonPDF((float)i,photonCount->designval);

  return;
}/*setPhotonProb*/


/*####################################################*/
/*Amy's Poission function*/

float poissonPDF(float n,float lambda)
{
  float y=0;
  y=exp(-1.0*lambda+n*log(lambda)-lgamma(n+1.0));
  return(y);
}/*poissonPDF*/


/*####################################################*/
/*determine point at which array exceeds threshold*/

float pickArrayElement(float photThresh,float *jimlad,int nBins,char interpolate)
{
  int i=0;
  float x=0,y0=0;
  float tot=0,*cumul=NULL;

  /*determine total energy and adjust threshold*/
  tot=0.0;
  cumul=falloc((uint64_t)nBins,"cumul",0);
  for(i=0;i<nBins;i++){
    tot+=jimlad[i];
    if(i>0)cumul[i]=cumul[i-1]+jimlad[i];
    else   cumul[i]=jimlad[i];
  }
  photThresh*=tot;

  /*determine cross point with binary search*/
  for(i=0;i<nBins;i++)if(cumul[i]>=photThresh)break;
  //start=0;
  //end=nBins-1;
  //while((end-start)>1){
  //  i=(end+start)/2;
  //  if(cumul[i]>photThresh)end=i;
  //  else if(cumul[i]<photThresh)start=i;
  //  else break;
  //}
  //if(cumul[i]>photThresh)i--;  /*get to the right side of the divide*/
  //if(i<0)i=0;
  //else if(i>=nBins)i=nBins-1;

  /*extrapolate between two elements*/
  if(interpolate){
    if(fabs(cumul[i]-photThresh)<TOL)x=(float)i;
    else{
      if(i>0)y0=cumul[i-1];
      else   y0=0.0;
      if(fabs(cumul[i]-y0)<TOL)x=(float)i;
      else{
        if(i<(nBins-1))x=(cumul[i]-photThresh)/(cumul[i]-y0)+(float)i;
        else           x=(float)(nBins-1);
      }
    }
  }else x=(float)i;
  TIDY(cumul);

  return(x);
}/*pickArrayElement*/


/*########################################################*/
/*set photon rates*/

void setPhotonRates(photonStruct *photonCount)
{

  /*what mode has been used to define photon rates?*/
  if(fabs(photonCount->rhoVrhoG-1.0)>TOL){                   /*adjusted rhoV/rhoG*/
    photonCount->nPhotG=2.0*photonCount->designval/photonCount->rhoVrhoG-1.0;
    photonCount->nPhotC=photonCount->nPhotG*photonCount->rhoVrhoG;
    photonCount->reflDiff=1;
    /*prevent negative rates. only needed for black soil etc.*/
    if(photonCount->nPhotC<0.0)photonCount->nPhotC=0.0;
    if(photonCount->nPhotG<0.0)photonCount->nPhotG=0.0;

  }else if((photonCount->nPhotG+photonCount->nPhotC)>0.0){   /*separately defined photon rates*/
    photonCount->designval=(photonCount->nPhotC+photonCount->nPhotG)/2.0;
    photonCount->rhoVrhoG=photonCount->nPhotC/photonCount->nPhotG;
    photonCount->reflDiff=1;

  }else if(fabs(photonCount->designval)<TOL){                /*no reflectance correction. seperate rates*/
    photonCount->designval=(photonCount->nPhotC+photonCount->nPhotG)/2.0;
    photonCount->rhoVrhoG=1.0;
    photonCount->reflDiff=0;

  }else{                                                     /*no reflectance correction. total rate*/
    photonCount->nPhotG=photonCount->nPhotC=photonCount->designval;
    photonCount->rhoVrhoG=1.0;
    photonCount->reflDiff=0;

  }

  return;
}/*setPhotonRates*/


/*the end*/
/*########################################################*/

