#include <Rcpp.h>

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "hdf5.h"
extern "C" {
#include "tools/tools.h"
#include "libclidar/libLasRead.h"
#include "libclidar/libLasProcess.h"
#include "libclidar/libLidarHDF.h"
#include "libclidar/libOctree.h"
#include "gedisimulator/gediIO.h"
#include "gedisimulator/gediNoise.h"
#include "gedisimulator/gediMetric.h"

/*#define USEPHOTON*/

  #ifdef USEPHOTON
  #include "photonCount.h"
  #endif

}


#include <gediMetric.hpp>

  /*##############################*/
  /*# Generates metrics from     #*/
  /*# simulated GEDI waveforms   #*/
  /*# 2015 svenhancock@gmail.com #*/
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

using namespace Rcpp;

  /*element reflectance*/
float rhoG;
float rhoC;

control* metric_makeControl(
  const char*     input,
  bool            writeFit        = false,
  bool            writeGauss      = false,
  bool            readBinLVIS     = false,
  bool            readHDFlvis     = false,
  bool            readHDFgedi     = false,
  CharacterVector level2          = CharacterVector(0),
  NumericVector   bounds          = NumericVector(0),
  bool            ground          = false,
  bool            useInt          = false,
  bool            useFrac         = false,
  float           rhRes           = 5,
  float           laiRes          = 10.0,
  float           laiH            = 30.0,
  bool            noRHgauss       = false,
  float           gTol            = 0.0,
  float           fhdHistRes      = 0.001,
  bool            forcePsigma     = false,
  bool            bayesGround     = false,
  bool            dontTrustGround = false,
  bool            noRoundCoord    = false,
  NumericVector   dcBias          = NumericVector(0),
  float           nSig            = 0.0,
  int             seed            = 1,
  float           hNoise          = 0.0,
  NumericVector   linkNoise       = NumericVector(0),
  float           linkFsig        = 5.5,
  float           linkPsig        = 0.764331,
  float           trueSig         = 5.0,
  int             bitRate         = 12,
  float           maxDN           = 4096.0,
  bool            renoise         = false,
  float           newPsig         = 1.0,
  float           oldPsig         = 0.764331,
  float           addDrift        = 0.0,
  bool            missGround      = 0,
  NumericVector   minGap          = NumericVector(0),
  bool            photonCount     = false,
  float           nPhotons        = 2.1,
  float           photonWind      = 200,
  float           noiseMult       = 0.1,
  float           meanN           = 0.0,
  float           thresh          = 0.00000001,
  bool            varNoise        = false,
  float           varScale        = 1.5,
  float           statsLen        = 30.0,
  bool            noiseTrack      = false,
  float           sWidth          = 0.0,
  float           psWidth         = 0.0,
  float           msWidth         = 0.0,
  bool            preMatchF       = false,
  bool            postMatchF      = false,
  CharacterVector pFile           = CharacterVector(0),
  float           gWidth          = 1.2,
  float           minGsig         = 0.764331,
  float           minWidth        = 0,
  bool            medNoise        = false,
  bool            varDrift        = false,
  NumericVector   driftFac        = NumericVector(0),
  float           varRhoG         = 0.4,
  float           varRhoC         = 0.57,
  float           pSigma          = 1.0,
  bool            gold            = false,
  float           deconTol        = 0.0000001
){
  control *dimage=NULL;

  /*allocate structures*/
  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  if(!(dimage->gediIO.den=(denPar *)calloc(1,sizeof(denPar)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  if(!(dimage->gediIO.gFit=(denPar *)calloc(1,sizeof(denPar)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }


  /*defaults*/
  /*input/output*/
  dimage->gediIO.nFiles=1;
  dimage->gediIO.inList = chChalloc(dimage->gediIO.nFiles,"inList",0);
  dimage->gediIO.inList[0] = challoc(200,"inList",0);
  strcpy(&(dimage->gediIO.inList[0][0]),"/Users/stevenhancock/data/gedi/analysis/side_lobe/laselva/contLaselva.0.12.wave");
  strcpy(dimage->outRoot,"teastMetric");
  dimage->maxGauss=20;
  dimage->opooMet=NULL;
  dimage->opooGauss=NULL;
  dimage->hdfGedi=NULL;
  /*scan settings*/
  dimage->gediIO.pSigma=0.764331; /*pulse length*/
  dimage->gediIO.fSigma=5.5;      /*footprint width*/
  dimage->gediIO.res=0.15;

  /*switches*/
  dimage->writeFit=0;
  dimage->gediIO.ground=0;
  dimage->gediIO.useInt=0;
  dimage->gediIO.useCount=1;
  dimage->gediIO.useFrac=0;
  dimage->rhRes=5.0;
  dimage->laiRes=10.0;
  dimage->maxLAIh=30.0;
  dimage->bayesGround=0;
  dimage->noise.missGround=0;
  dimage->noise.linkNoise=0;
  dimage->noise.driftFact=0.0;
  dimage->gediIO.linkPsig=0.764331; /*pulse length*/
  dimage->gediIO.linkFsig=5.5;      /*footprint width*/
  dimage->noise.trueSig=5.0;
  dimage->noise.deSig=0.0; //0.1; //4.0*0.15/2.355;
  dimage->noise.bitRate=12;
  dimage->noise.maxDN=4096.0; //1.0/(dimage->pSigma*sqrt(2.0*M_PI));
  dimage->noise.minGap=0.0;
  dimage->noRHgauss=0;    /*do find RH metrics by Gaussian fitting*/
  dimage->renoiseWave=0;  /*do not denoise "truth"*/
  dimage->noise.newPsig=-1.0;   /*leave blank*/
  dimage->gediIO.dontTrustGround=0;  /*do trust ground in waveforms, if there*/
  dimage->readBinLVIS=0;      /*read ASCII rather than binary LVIS*/
  dimage->readHDFlvis=0;      /*read ASCII rather than HDF5 LVIS*/
  dimage->readHDFgedi=0;      /*read ASCII rather than HDF5 GEDI*/
  dimage->gediIO.readPsigma=1;       /*read pSigma from file*/
  dimage->coord2dp=1;         /*round up coords in output*/
  dimage->useBounds=0;        /*process all data provided*/
  dimage->writeGauss=0;       /*do not write Gaussian parameters*/

  /*set default denoising parameters*/
  setDenoiseDefault(dimage->gediIO.den);
  dimage->gediIO.den->meanN=0.0;  /*we haven't added noise yet*/
  dimage->gediIO.den->thresh=0.00000001;  /*tiny number as no noise yet*/
  dimage->gediIO.den->noiseTrack=0;
  dimage->gediIO.den->minWidth=0;
  dimage->gediIO.den->varNoise=0;
  dimage->gediIO.den->threshScale=1.5;
  dimage->gediIO.den->fitGauss=0;
  dimage->gediIO.den->psWidth=0.0;

  /*set default Gaussian fitting  parameters*/
  setDenoiseDefault(dimage->gediIO.gFit);
  dimage->gediIO.gFit->meanN=0.0;    /*no denoising here*/
  dimage->gediIO.gFit->thresh=0.000000005;    /*no denoising here*/
  dimage->gediIO.gFit->noiseTrack=0;    /*no denoising here*/
  dimage->gediIO.gFit->minWidth=0;    /*no denoising here*/
  dimage->gediIO.gFit->varNoise=0;    /*no denoising here*/
  dimage->gediIO.gFit->gWidth=1.2;
  dimage->gediIO.gFit->sWidth=0.0;
  dimage->gediIO.gFit->fitGauss=1;
  dimage->gediIO.gFit->minGsig=0.764331;
  /*noise parameters*/
  dimage->noise.meanN=0.0;
  dimage->noise.nSig=nSig;
  dimage->bThresh=0.001;
  dimage->noise.hNoise=hNoise;
  dimage->noise.offset=94.0;
  /*projection, not yet used*/
  dimage->gediIO.wEPSG=4326;  /*waveforms*/
  dimage->gediIO.bEPSG=4326;  /*bounds*/
  /*LVIS data*/
  dimage->lvis.data=NULL;
  dimage->hdfLvis=NULL;
  /*LVIS level2 data*/
  dimage->readL2=0;   /*do not read L2*/
  /*photon counting*/
  dimage->ice2=0;             /*GEDI mode, rather than ICESat-2*/


  // Change parameters based on user input



  // TODO if readHDFlvis
  dimage->readHDFlvis=1;
  dimage->gediIO.readPsigma=0;

  // TODO if varNoise
  dimage->gediIO.den->varNoise=1;

  // TODO if varScale
  dimage->gediIO.den->varNoise=1;
  dimage->gediIO.den->threshScale=3;

  // TODO minWidth
  dimage->gediIO.den->minWidth=2;

  // TODO sWidth
  dimage->gediIO.den->sWidth=8;

  // TODO statsLen
  dimage->gediIO.den->statsLen=10;




#ifdef USEPHOTON
  dimage->photonCount.designval=2.1;
  dimage->photonCount.prob=NULL;
  dimage->photonCount.pBins=0;
  dimage->photonCount.H=200.0;
  dimage->photonCount.noise_mult=0.1;
  dimage->photonCount.rhoVrhoG=1.0;
  dimage->photonCount.writeHDF=0;  /*write ASCII by default*/
  dimage->photonCount.hdf=NULL;
#endif
  /*others*/
  rhoG=0.4;
  rhoC=0.57;
  dimage->rhoRatio=rhoC/rhoG;
  dimage->gTol=0.0;
  dimage->gediIO.nMessages=200;
  dimage->fhdHistRes=0.001;

  /*read the command line*/
  TTIDY((void **)dimage->gediIO.inList,dimage->gediIO.nFiles);
  dimage->gediIO.inList=NULL;
  dimage->gediIO.nFiles=1;
  dimage->gediIO.inList=chChalloc(dimage->gediIO.nFiles,"input name list",0);
  dimage->gediIO.inList[0]=challoc((uint64_t)strlen(input)+1,"input name list",0);
  strcpy(dimage->gediIO.inList[0],input);
  if (level2.length() == 1) {
    dimage->readL2=1;
    strcpy(dimage->l2namen,level2[0]);
  }else if(writeFit){
    dimage->writeFit=1;
  }else if(writeGauss){
    dimage->writeGauss=1;
  }else if(dcBias.length() == 1){
    dimage->noise.meanN=dimage->noise.offset=dcBias[0];
  }else if(ground){
    dimage->gediIO.ground=1;
  }else if(varNoise){
    dimage->gediIO.den->varNoise=1;
  }else if(medNoise){
    dimage->gediIO.den->medStats=1;
  }else if(noiseTrack){
    dimage->gediIO.den->noiseTrack=1;
  }else if(pFile.length() == 1){
    strcpy(dimage->gediIO.den->pNamen,pFile[0]);
    dimage->gediIO.den->deconGauss=0;
  }else if(gold){
    dimage->gediIO.den->deconMeth=0;
  }else if(preMatchF){
    dimage->gediIO.den->preMatchF=1;
  }else if(postMatchF){
    dimage->gediIO.den->posMatchF=1;
  }else if(useInt){
    dimage->gediIO.useInt=1;
    dimage->gediIO.useCount=0;
    dimage->gediIO.useFrac=0;
  }else if(useFrac){
    dimage->gediIO.useInt=0;
    dimage->gediIO.useCount=0;
    dimage->gediIO.useFrac=1;
  }else if(linkNoise.length() == 2){
    dimage->noise.linkNoise=1;
    dimage->noise.linkM=linkNoise[0];
    dimage->noise.linkCov=linkNoise[1];
  }else if(missGround){
    dimage->noise.missGround=1;
  }else if(minGap.length() == 1){
    dimage->noise.missGround=1;
    dimage->noise.minGap=minGap[0];
  }else if(bayesGround){
    dimage->bayesGround=1;
  }else if(bitRate != 12){
    dimage->noise.bitRate=bitRate;
    dimage->noise.maxDN=pow(2.0,(float)dimage->noise.bitRate);
  }else if(noRHgauss){
    dimage->noRHgauss=1;
  }else if(renoise){
    dimage->renoiseWave=1;
  }else if(dontTrustGround){
    dimage->gediIO.dontTrustGround=1;
  }else if(readBinLVIS){
    dimage->readBinLVIS=1;
  }else if(readHDFlvis){
    dimage->readHDFlvis=1;
    dimage->gediIO.readPsigma=0;
  }else if(readHDFgedi){
    dimage->readHDFgedi=1;
  }else if(forcePsigma){
    dimage->gediIO.readPsigma=0;
  }else if(noRoundCoord){
    dimage->coord2dp=0;
  }else if(bounds.length() == 4){
    dimage->useBounds=1;
    dimage->minX=dimage->gediIO.bounds[0]=bounds[0];
    dimage->minY=dimage->gediIO.bounds[1]=bounds[1];
    dimage->maxX=dimage->gediIO.bounds[2]=bounds[2];
    dimage->maxY=dimage->gediIO.bounds[3]=bounds[3];
  }else if(varDrift){
    dimage->gediIO.den->corrDrift=1;
    dimage->gediIO.den->varDrift=1;
  }else if(driftFac.length() == 1) {
    dimage->gediIO.den->corrDrift=1;
    dimage->gediIO.den->varDrift=0;
    dimage->gediIO.den->fixedDrift=driftFac[0];
  #ifdef USEPHOTON
  }else if(!strncasecmp(argv[i],"-photonCount",12)){
    dimage->ice2=1;
  }else if(!strncasecmp(argv[i],"-nPhotons",9)){
    checkArguments(1,i,argc,"-nPhotons");
    dimage->photonCount.designval=atof(argv[++i]);
  }else if(!strncasecmp(argv[i],"-photonWind",11)){
    checkArguments(1,i,argc,"-photonWind");
    dimage->photonCount.H=atof(argv[++i]);
  }else if(!strncasecmp(argv[i],"-noiseMult",10)){
    checkArguments(1,i,argc,"-noiseMult");
    dimage->photonCount.noise_mult=atof(argv[++i]);
  }else if(!strncasecmp(argv[i],"-rhoVrhoG",9)){
    checkArguments(1,i,argc,"-rhoVrhoG");
    dimage->photonCount.rhoVrhoG=atof(argv[++i]);
  }else if(!strncasecmp(argv[i],"-photHDF",8)){
    dimage->photonCount.writeHDF=1;
  #endif
  }


  // Parameters default
  srand(seed);
  dimage->gediIO.den->meanN=meanN;
  dimage->gediIO.gFit->gWidth=gWidth;
  dimage->gediIO.den->thresh=thresh;
  dimage->gediIO.den->sWidth=sWidth;
  dimage->gediIO.den->psWidth=psWidth;
  dimage->gediIO.den->msWidth=msWidth;
  dimage->gediIO.gFit->minGsig=minGsig;
  dimage->gediIO.den->minWidth=minWidth;
  dimage->gediIO.den->statsLen=statsLen;
  dimage->gediIO.den->threshScale=varScale;
  dimage->gediIO.den->pSigma=pSigma;
  dimage->gediIO.den->deChang=deconTol;
  dimage->rhRes=rhRes;
  dimage->maxLAIh=laiH;
  dimage->noise.trueSig=trueSig;
  rhoG=varRhoG;
  rhoC=varRhoC;
  dimage->noise.maxDN=maxDN;
  dimage->gTol=gTol;
  dimage->noise.newPsig=newPsig;
  dimage->gediIO.pSigma=oldPsig;
  dimage->gediIO.linkPsig=linkPsig;
  dimage->gediIO.linkFsig=linkFsig;
  dimage->fhdHistRes=fhdHistRes;
  dimage->noise.driftFact=addDrift;

  /*read deconvolution pulse if needed*/
  if(dimage->gediIO.den->preMatchF||dimage->gediIO.den->preMatchF||dimage->gediIO.den->deconMeth>=0)readPulse(dimage->gediIO.den);
  if((!dimage->gediIO.ground)&&(dimage->noise.missGround)){
    fprintf(stderr,"Noise option conflict. Cannot use missGround without ground\n");
    Rcpp::stop("1");
  }

  return(dimage);
}/*readCommands*/


// [[Rcpp::export]]
Rcpp::DataFrame processFloWave2(
  CharacterVector input,
  bool            writeFit        = false,
  bool            writeGauss      = false,
  bool            readBinLVIS     = false,
  bool            readHDFlvis     = false,
  bool            readHDFgedi     = false,
  CharacterVector level2          = CharacterVector(0),
  NumericVector   bounds          = NumericVector(0),
  bool            ground          = false,
  bool            useInt          = false,
  bool            useFrac         = false,
  float           rhRes           = 5,
  float           laiRes          = 10.0,
  float           laiH            = 30.0,
  bool            noRHgauss       = false,
  float           gTol            = 0.0,
  float           fhdHistRes      = 0.001,
  bool            forcePsigma     = false,
  bool            bayesGround     = false,
  bool            dontTrustGround = false,
  bool            noRoundCoord    = false,
  NumericVector   dcBias          = NumericVector(0),
  float           nSig            = 0.0,
  int             seed            = 1,
  float           hNoise          = 0.0,
  NumericVector   linkNoise       = NumericVector(0),
  float           linkFsig        = 5.5,
  float           linkPsig        = 0.764331,
  float           trueSig         = 5.0,
  int             bitRate         = 12,
  float           maxDN           = 4096.0,
  bool            renoise         = false,
  float           newPsig         = 1.0,
  float           oldPsig         = 0.764331,
  float           addDrift        = 0.0,
  bool            missGround      = 0,
  NumericVector   minGap          = NumericVector(0),
  bool            photonCount     = false,
  float           nPhotons        = 2.1,
  float           photonWind      = 200,
  float           noiseMult       = 0.1,
  float           meanN           = 0.0,
  float           thresh          = 0.00000001,
  bool            varNoise        = false,
  float           varScale        = 1.5,
  float           statsLen        = 30.0,
  bool            noiseTrack      = false,
  float           sWidth          = 0.0,
  float           psWidth         = 0.0,
  float           msWidth         = 0.0,
  bool            preMatchF       = false,
  bool            postMatchF      = false,
  CharacterVector pFile           = CharacterVector(0),
  float           gWidth          = 1.2,
  float           minGsig         = 0.764331,
  float           minWidth        = 0,
  bool            medNoise        = false,
  bool            varDrift        = false,
  NumericVector   driftFac        = NumericVector(0),
  float           varRhoG         = 0.4,
  float           varRhoC         = 0.57,
  float           pSigma          = 1.0,
  bool            gold            = false,
  float           deconTol        = 0.0000001
)
{
  int i=0;
  control *dimage=NULL;
  dataStruct *data=NULL;
  metStruct *metric=NULL;
  Rcpp::DataFrame resultDF=DataFrame::create();
  std::string input_str = Rcpp::as<std::string>(input[0]);
  float *processed=NULL,*denoised=NULL;;

  /*read command Line*/
    dimage=metric_makeControl(
      input_str.c_str(),
      writeFit,
      writeGauss,
      readBinLVIS,
      readHDFlvis,
      readHDFgedi,
      level2,
      bounds,
      ground,
      useInt,
      useFrac,
      rhRes,
      laiRes,
      laiH,
      noRHgauss,
      gTol,
      fhdHistRes,
      forcePsigma,
      bayesGround,
      dontTrustGround,
      noRoundCoord,
      dcBias,
      nSig,
      seed,
      hNoise,
      linkNoise,
      linkFsig,
      linkPsig,
      trueSig,
      bitRate,
      maxDN,
      renoise,
      newPsig,
      oldPsig,
      addDrift,
      missGround,
      minGap,
      photonCount,
      nPhotons,
      photonWind,
      noiseMult,
      meanN,
      thresh,
      varNoise,
      varScale,
      statsLen,
      noiseTrack,
      sWidth,
      psWidth,
      msWidth,
      preMatchF,
      postMatchF,
      pFile,
      gWidth,
      minGsig,
      minWidth,
      medNoise,
      varDrift,
      driftFac,
      varRhoG,
      varRhoC,
      pSigma,
      gold,
      deconTol
    );


  /*set link noise if needed*/
    dimage->noise.linkSig=setNoiseSigma(dimage->noise.linkM,dimage->noise.linkCov,dimage->gediIO.linkPsig,dimage->gediIO.linkFsig,rhoC,rhoG);

  /*allocate metric array*/
    if(!(metric=(metStruct *)calloc(1,sizeof(metStruct)))){
      fprintf(stderr,"error metric structure allocation.\n");
      Rcpp::stop("1");
    }

  /*loop over files*/
    for(i=0;i<dimage->gediIO.nFiles;i++){
      Rcout << "Wave " <<  i+1 << " of " << dimage->gediIO.nFiles << "          \r";

    /*read waveform*/
    if(dimage->readBinLVIS)     data=readBinaryLVIS(dimage->gediIO.inList[0],&dimage->lvis,i,&dimage->gediIO);
    else if(dimage->readHDFlvis)data=unpackHDFlvis(dimage->gediIO.inList[0],&dimage->hdfLvis,&dimage->gediIO,i);
    else if(dimage->readHDFgedi)data=unpackHDFgedi(dimage->gediIO.inList[0],&dimage->gediIO,&dimage->hdfGedi,i);
    else                        data=readASCIIdata(dimage->gediIO.inList[i],&(dimage->gediIO));
    if(dimage->readL2)setL2ground(data,i,dimage);

    /*check bounds if needed*/
    if(dimage->useBounds)checkWaveformBounds(data,dimage);

    if (resultDF.ncol() == 0) {
      resultDF = createMetricsDataFrame(dimage);
    }
    /*is the data usable*/
    if(data->usable){
      /*denoise and change pulse if needed*/
      if(dimage->renoiseWave)modifyTruth(data,&dimage->noise);

      /*determine truths before noising*/
      determineTruth(data,dimage);

/*set link noise if needed*/

    if (  dimage->noise.driftFact ||
          dimage->noise.missGround ||
          dimage->noise.linkNoise ||
          dimage->noise.nSig>0.0 ||
          dimage->noise.meanN>0.0 ||
          dimage->noise.hNoise>0.0 )
    {
      addNoise(data,&dimage->noise,dimage->gediIO.fSigma,dimage->gediIO.pSigma,dimage->gediIO.res,rhoC,rhoG);
      denoised=processFloWave(data->noised,data->nBins,dimage->gediIO.den,1.0);
    } else {
      /*process waveform*/
      /*denoise*/
      denoised=processFloWave(data->wave[data->useType],data->nBins,dimage->gediIO.den,1.0);
    }

      /*are we in GEDI mode?*/
      if(!dimage->ice2){
        /*Gaussian fit*/
        if(dimage->noRHgauss==0)processed=processFloWave(denoised,data->nBins,dimage->gediIO.gFit,1.0);

        /*shift Gaussian centres to align to absolute elevation*/
        alignElevation(data->z[0],data->z[data->nBins-1],dimage->gediIO.gFit->gPar,dimage->gediIO.gFit->nGauss);

        /*determine metrics*/
        findMetrics(metric,dimage->gediIO.gFit->gPar,dimage->gediIO.gFit->nGauss,denoised,data->noised,data->nBins,data->z,dimage,data);

        /*write results*/
        if(dimage->readBinLVIS||dimage->readHDFlvis||dimage->readHDFgedi)writeMetricsDataFrame(data,dimage,metric,i,denoised,processed,resultDF);
        else                                                             writeMetricsDataFrame(data,dimage,metric,i,denoised,processed,resultDF);
      }else{  /*ICESat-2 mode*/
        photonCountCloud(denoised,data,&dimage->photonCount,dimage->outRoot,i,dimage->gediIO.den,&dimage->noise);
      }/*operation mode switch*/
    }/*is the data usable*/


    /*tidy as we go along*/
    TIDY(processed);
    TIDY(denoised);
    if(data){
      TIDY(data->noised);
      if(dimage->readHDFgedi){  /*pointer to array. do not free*/
        data->wave[0]=NULL;
        if(data->ground)data->ground[0]=NULL;
      }
      TTIDY((void **)data->ground,data->nWaveTypes);
      TTIDY((void **)data->wave,data->nWaveTypes);
      TIDY(data->totE);
      TIDY(data->z);
      TIDY(data);
    }
    TIDY(dimage->gediIO.gFit->gPar);
    TIDY(dimage->gediIO.den->gPar);
    dimage->gediIO.den->nGauss=0;
    dimage->gediIO.gFit->nGauss=0;
    TIDY(metric->rhMax);
    TIDY(metric->rhInfl);
    TIDY(metric->rhReal);
    TIDY(metric->rh);
    TIDY(metric->bGr);
    TIDY(metric->tLAI);
    TIDY(metric->gLAI);
    TIDY(metric->hgLAI);
    TIDY(metric->hiLAI);
    TIDY(metric->hmLAI);
    //TIDY(metric->LmomGau);
    //TIDY(metric->LmomRea);
    //TIDY(metric->LmomInf);
    //TIDY(metric->LmomMax);
  }/*file loop*/

  /*TIDY LVIS data if it was read*/
  if(dimage->readBinLVIS)TIDY(dimage->lvis.data);
  if(dimage->readHDFgedi)dimage->hdfGedi=tidyGediHDF(dimage->hdfGedi);


  #ifdef USEPHOTON
  else             Rprintf("Written to %s\n",dimage->photonCount.outNamen);
  #endif


  /*tidy up arrays*/
  tidySMoothPulse();
  TIDY(metric);
  if(dimage){
    if(dimage->lvisL2){
      TIDY(dimage->lvisL2->lfid);
      TIDY(dimage->lvisL2->shotN);
      TIDY(dimage->lvisL2->zG);
      TIDY(dimage->lvisL2);
    }
    if(dimage->readBinLVIS||dimage->readHDFlvis||dimage->readHDFgedi)TTIDY((void **)dimage->gediIO.inList,1);
    else                                        TTIDY((void **)dimage->gediIO.inList,dimage->gediIO.nFiles);
    dimage->gediIO.inList=NULL;
    if(dimage->opooMet){
      fclose(dimage->opooMet);
      dimage->opooMet=NULL;
    }
    if(dimage->opooGauss){
      fclose(dimage->opooGauss);
      dimage->opooGauss=NULL;
    }
    #ifdef USEPHOTON
    if(dimage->photonCount.opoo){
      fclose(dimage->photonCount.opoo);
      dimage->photonCount.opoo=NULL;
    }
    TIDY(dimage->photonCount.prob);
    #endif
    if(dimage->gediIO.den){
      TTIDY((void **)dimage->gediIO.den->pulse,2);
      TIDY(dimage->gediIO.den->matchPulse);
      TIDY(dimage->gediIO.den->hardPulse);
      TIDY(dimage->gediIO.den);
    }
    if(dimage->gediIO.gFit){
      TTIDY((void **)dimage->gediIO.gFit->pulse,2);
      TIDY(dimage->gediIO.gFit);
    }
    dimage->hdfLvis=tidyLVISstruct(dimage->hdfLvis);
    TIDY(dimage);
  }
  return(resultDF);
}/*main*/

DataFrame createMetricsDataFrame(control* dimage) {
  char name[22];
  DataFrame df = DataFrame::create();


  df.push_back(CharacterVector(dimage->gediIO.nFiles), "waveID");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"true ground");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"true top");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"ground slope");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"ALS cover");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"gHeight");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"maxGround");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"inflGround");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"signal top");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"signal bottom");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"cover");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"leading edge ext");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"trailing edge extent");

  // RHs
  int nRH = (int)(100.0/dimage->rhRes+1);
  for (int i=0;i<nRH;i++) {
    sprintf(name, "rhGauss %g", (float)i*dimage->rhRes);
    df.push_back(NumericVector(dimage->gediIO.nFiles), name);
  }
  for (int i=0;i<nRH;i++) {
    sprintf(name, "rhMax %g", (float)i*dimage->rhRes);
    df.push_back(NumericVector(dimage->gediIO.nFiles), name);
  }
  for (int i=0;i<nRH;i++) {
    sprintf(name, "rhInfl %g", (float)i*dimage->rhRes);
    df.push_back(NumericVector(dimage->gediIO.nFiles), name);
  }
  for (int i=0;i<nRH;i++) {
    sprintf(name, "rhReal %g", (float)i*dimage->rhRes);
    df.push_back(NumericVector(dimage->gediIO.nFiles), name);
  }

  df.push_back(NumericVector(dimage->gediIO.nFiles),"gaussHalfCov");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"maxHalfCov");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"infHalfCov");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"bayHalfCov");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"pSigma");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"fSigma");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"linkM");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"linkCov");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"lon");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"lat");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"groundOverlap");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"groundMin");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"groundInfl");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"waveEnergy");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"blairSense");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"pointDense");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"beamDense");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"zenith");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"FHD");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"niM2");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"niM2.1");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"meanNoise");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"noiseStdev");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"noiseThresh");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"FHDhist");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"FHDcan");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"FHDcanHist");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"FHDcanGauss");
  df.push_back(NumericVector(dimage->gediIO.nFiles),"FHDcanGhist");

  // LAI
  int laiBins = (int)(dimage->maxLAIh/dimage->laiRes+0.5)+1;
  for (int i=0;i<laiBins;i++) {
    sprintf(name, "tLAI%dt%d", i*(int)dimage->laiRes,(i+1)*(int)dimage->laiRes);
    df.push_back(NumericVector(dimage->gediIO.nFiles), name);
  }
  for (int i=0;i<laiBins;i++) {
    sprintf(name, "gLAI%dt%d", i*(int)dimage->laiRes,(i+1)*(int)dimage->laiRes);
    df.push_back(NumericVector(dimage->gediIO.nFiles), name);
  }
  for (int i=0;i<laiBins;i++) {
    sprintf(name, "hgLAI%dt%d", i*(int)dimage->laiRes,(i+1)*(int)dimage->laiRes);
    df.push_back(NumericVector(dimage->gediIO.nFiles), name);
  }
  for (int i=0;i<laiBins;i++) {
    sprintf(name, "hiLAI%dt%d", i*(int)dimage->laiRes,(i+1)*(int)dimage->laiRes);
    df.push_back(NumericVector(dimage->gediIO.nFiles), name);
  }
  for (int i=0;i<laiBins;i++) {
    sprintf(name, "hmLAI%dt%d", i*(int)dimage->laiRes,(i+1)*(int)dimage->laiRes);
    df.push_back(NumericVector(dimage->gediIO.nFiles), name);
  }
  return (df);
}

void writeMetricsDataFrame(dataStruct *data,control *dimage,metStruct *metric,int numb,float *denoised,float *processed, DataFrame output) {
  int colIndex = 0;
  CharacterVector waveIDs = output[colIndex];
  waveIDs[numb] = data->waveID;
  output[colIndex++] = waveIDs;
  //true ground
  ((NumericVector)(output[colIndex++]))[numb] = data->gElev;
  //true top
  ((NumericVector)(output[colIndex++]))[numb] = data->tElev;
  //ground slope
  ((NumericVector)(output[colIndex++]))[numb] = data->slope;
  //ALS cover
  ((NumericVector)(output[colIndex++]))[numb] = data->cov;
  //gHeight
  ((NumericVector)(output[colIndex++]))[numb] = metric->gHeight;
  //maxGround
  ((NumericVector)(output[colIndex++]))[numb] = metric->maxGround;
  //inflGround
  ((NumericVector)(output[colIndex++]))[numb] = metric->inflGround;
  //signal top
  ((NumericVector)(output[colIndex++]))[numb] = metric->tElev;
  //signal bottom
  ((NumericVector)(output[colIndex++]))[numb] = metric->bElev;
  //cover
  ((NumericVector)(output[colIndex++]))[numb] = metric->cov;
  //leading edge ext
  ((NumericVector)(output[colIndex++]))[numb] = metric->leExt;
  //trailing edge extent
  ((NumericVector)(output[colIndex++]))[numb] = metric->teExt;
  //rhGauss ...
  for (int i=0; i < metric->nRH; i++)
    ((NumericVector)(output[colIndex++]))[numb] = metric->rh[i];
  //rhMax ...
  for (int i=0; i < metric->nRH; i++)
    ((NumericVector)(output[colIndex++]))[numb] = metric->rhMax[i];
  //rhInfl ...
  for (int i=0; i < metric->nRH; i++)
    ((NumericVector)(output[colIndex++]))[numb] = metric->rhInfl[i];
  //rhReal ...
  for (int i=0; i < metric->nRH; i++)
    ((NumericVector)(output[colIndex++]))[numb] = metric->rhReal[i];
  //gaussHalfCov
  ((NumericVector)(output[colIndex++]))[numb] = metric->covHalfG;
  //maxHalfCov
  ((NumericVector)(output[colIndex++]))[numb] = metric->covHalfM;
  //infHalfCov
  ((NumericVector)(output[colIndex++]))[numb] = metric->covHalfI;
  //bayHalfCov
  ((NumericVector)(output[colIndex++]))[numb] = metric->covHalfB;
  //pSigma
  ((NumericVector)(output[colIndex++]))[numb] = data->pSigma;
  //fSigma
  ((NumericVector)(output[colIndex++]))[numb] = data->fSigma;
  if(dimage->noise.linkNoise) {
    //linkM
    ((NumericVector)(output[colIndex++]))[numb] = dimage->noise.linkM;
    //linkCov
    ((NumericVector)(output[colIndex++]))[numb] = dimage->noise.linkCov;
  } else {
    ((NumericVector)(output[colIndex++]))[numb] = NA_REAL;
    //linkCov
    ((NumericVector)(output[colIndex++]))[numb] = NA_REAL;
  }
  //lon
  ((NumericVector)(output[colIndex++]))[numb] = data->lon;
  //lat
  ((NumericVector)(output[colIndex++]))[numb] = data->lat;
  //groundOverlap
  ((NumericVector)(output[colIndex++]))[numb] = data->gLap;
  //groundMin
  ((NumericVector)(output[colIndex++]))[numb] = data->gMinimum;
  //groundInfl
  ((NumericVector)(output[colIndex++]))[numb] = data->gInfl;
  //waveEnergy
  ((NumericVector)(output[colIndex++]))[numb] = metric->totE;
  //blairSense
  ((NumericVector)(output[colIndex++]))[numb] = metric->blairSense;
  //pointDense
  ((NumericVector)(output[colIndex++]))[numb] = data->pointDense;
  //beamDense
  ((NumericVector)(output[colIndex++]))[numb] = data->beamDense;
  //zenith
  ((NumericVector)(output[colIndex++]))[numb] = data->zen;
  //FHD
  ((NumericVector)(output[colIndex++]))[numb] = metric->FHD;
  //niM2
  ((NumericVector)(output[colIndex++]))[numb] = metric->niM2;
  //niM2.1
  ((NumericVector)(output[colIndex++]))[numb] = metric->niM21;
  //meanNoise
  ((NumericVector)(output[colIndex++]))[numb] = dimage->gediIO.den->meanN;
  //noiseStdev
  ((NumericVector)(output[colIndex++]))[numb] = (dimage->gediIO.den->thresh-dimage->gediIO.den->meanN)/dimage->gediIO.den->threshScale;
  //noiseThresh
  ((NumericVector)(output[colIndex++]))[numb] = dimage->gediIO.den->thresh;
  //FHDhist
  ((NumericVector)(output[colIndex++]))[numb] = metric->FHDhist;
  //FHDcan
  ((NumericVector)(output[colIndex++]))[numb] = metric->FHDcan;
  //FHDcanHist
  ((NumericVector)(output[colIndex++]))[numb] = metric->FHDcanH;
  //FHDcanGauss
  ((NumericVector)(output[colIndex++]))[numb] = metric->FHDcanGauss;
  //FHDcanGhist
  ((NumericVector)(output[colIndex++]))[numb] = metric->FHDcanGhist;
  //tLAI ...
  for (int i=0; i < metric->laiBins; i++)
    ((NumericVector)(output[colIndex++]))[numb] = metric->tLAI[i];
  //gLAI ...
  for (int i=0; i < metric->laiBins; i++)
    ((NumericVector)(output[colIndex++]))[numb] = metric->gLAI[i];
  //hgLAI ...
  for (int i=0; i < metric->laiBins; i++)
    ((NumericVector)(output[colIndex++]))[numb] = metric->hgLAI[i];
  //hiLAI ...
  for (int i=0; i < metric->laiBins; i++)
    ((NumericVector)(output[colIndex++]))[numb] = metric->hiLAI[i];
  //hmLAI ...
  for (int i=0; i < metric->laiBins; i++)
    ((NumericVector)(output[colIndex++]))[numb] = metric->hmLAI[i];
}
