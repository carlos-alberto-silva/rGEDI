#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "hdf5.h"
extern "C" {
#include "tools/tools.h"
#include "tools/msgHandling.h"
#include "libclidar/libLasRead.h"
#include "libclidar/libLidVoxel.h"
#include "libclidar/libLasProcess.h"
#include "libclidar/libLidarHDF.h"
#include "libclidar/libOctree.h"
#include "gedisimulator/gediIO.h"
#include "gedisimulator/gediRat.h"
/*#define USEPHOTON*/

  #ifdef USEPHOTON
  #include "photonCount.h"
  #endif

}

#include <Rcpp.h>
#include <gediRat.hpp>

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

rat_control* rat_makeControl(
  // Input output filenames and format
  CharacterVector input,
  CharacterVector output,
  bool            inList          = false,
  bool            ground          = false,
  bool            hdf             = false,
  CharacterVector waveID          = CharacterVector(0),

  // Single footprint, list of footprints, or grid of footprints
  NumericVector   coords          = NumericVector(0),
  CharacterVector listCoord       = CharacterVector(0),
  NumericVector   gridBound       = NumericVector(0),
  float           gridStep        = 30.0,

  // Lidar characteristics. Defaults are expected GEDI values.
  float           pSigma          = -1.0,
  float           pFWHM           = 15.0,
  CharacterVector readPulse       = CharacterVector(0),
  float           fSigma          = 5.5,
  CharacterVector wavefront       = CharacterVector(0),
  float           res             = 0.15,
  bool            LVIS            = false,
  bool            topHat          = false,
  bool            sideLobe        = false,
  float           lobeAng         = 0.0,


  // Input data quality filters
  bool            checkCover      = false,
  float           maxScanAng      = 1000000.0,
  float           decimate        = 1.0,

  // Computational speed options
  uint64_t        pBuff           = 0.2,
  int             maxBins         = 1024,
  bool            countOnly       = false,
  bool            pulseAfter      = false,
  bool            pulseBefore     = false,
  bool            noNorm          = false,

  // Octree
  bool            noOctree        = false,
  int             octLevels       = 0,
  int             nOctPix         = 40,

  // Using full-waveform input data (not tested)
  bool            decon           = false,
  bool            indDecon        = false,
  bool            readWave        = false,

  // Miscellaneous
  bool            listFiles       = false,
  bool            keepOld         = false,
  bool            useShadow       = false,
  bool            polyGround      = false,
  bool            nnGround        = false,
  NumericVector   seed            = NumericVector(0)
)
{
  rat_control *dimage=NULL;

  if(!(dimage=(rat_control *)calloc(1,sizeof(rat_control)))){
    errorf("error control allocation.\n");
    exit(1);
  }

  dimage->gediIO.nFiles=1;
  dimage->inList=chChalloc(dimage->gediIO.nFiles,"inList",0);
  dimage->inList[0]=challoc(200,"inList",0);
  strcpy(&(dimage->inList[0][0]),"/Users/dill/data/teast/maryland_play/sc_79_112_1.las");
  strcpy(dimage->outNamen,"teast.wave");
  dimage->gediIO.pRes=0.01;
  dimage->gediRat.coord[0]=624366.0;
  dimage->gediRat.coord[1]=3.69810*pow(10.0,6.0);
  dimage->gediRat.decon=NULL;

  /*switches*/
  dimage->gediRat.readWave=0;
  dimage->gediIO.ground=0;
  dimage->gediRat.sideLobe=0;   /*no side lobes*/
  dimage->gediRat.pulseAfter=1;  /*smooth after for speed*/
  dimage->listFiles=0;
  dimage->pBuffSize=(uint64_t)200000000;
  dimage->gediRat.checkCover=0;
  dimage->gediRat.normCover=1;
  dimage->gediRat.cleanOut=0;
  dimage->gediRat.topHat=0;
  dimage->useID=0;
  dimage->gediIO.readPulse=0;
  dimage->gediRat.useShadow=0;
  dimage->gediRat.vRes[0]=dimage->gediRat.vRes[1]=dimage->gediRat.vRes[2]=1.0;
  dimage->gediRat.beamRad=0.165;    /*33 cm*/
  dimage->polyGr=0;     /*don't fit a polynomial through the ground*/
  dimage->nnGr=0;       /*don't make a DEM from nearest neighbour*/
  dimage->overWrite=1;  /*over write any files with the same name if they exist*/
  dimage->gediRat.readALSonce=0;/*read each footprint separately*/
  dimage->writeHDF=0;   /*write output as ascii*/
  dimage->gediRat.defWfront=0;   /*Gaussian footprint*/
  dimage->gediRat.wavefront=NULL;

  /*beams*/
  dimage->gediIO.useCount=dimage->gediIO.useFrac=dimage->gediIO.useInt=1;

  /*octree*/
  dimage->gediRat.useOctree=1;
  dimage->gediRat.octree=NULL;

  /*gridding options*/
  dimage->gediRat.doGrid=0;           /*gridded switch*/
  dimage->gediRat.gNx=dimage->gediRat.gNy=0;
  /*batch*/
  dimage->gediRat.coords=NULL;       /*list of coordinates*/
  dimage->gediRat.waveIDlist=NULL;   /*list of waveform IDs*/
  dimage->gediIO.nMessages=200;
  strcpy(dimage->waveID,"gediWave");

  dimage->gediRat.iThresh=0.0006;
  dimage->gediRat.meanN=12.0;
  dimage->gediRat.doDecon=0;
  dimage->gediRat.indDecon=0;

  /*read the command line*/
  std::string inputStr = Rcpp::as<std::string>(input[0]);
  char* inputCStr = const_cast<char*>(inputStr.c_str());
  if (inList) {

    TTIDY((void **)dimage->inList,dimage->gediIO.nFiles);
    dimage->inList=readInList(&dimage->gediIO.nFiles, inputCStr);
  } else {
    TTIDY((void **)dimage->inList,dimage->gediIO.nFiles);
    dimage->gediIO.nFiles=1;
    dimage->gediIO.inList=chChalloc(dimage->gediIO.nFiles,"input name list",0);
    dimage->gediIO.inList[0]=challoc((uint64_t)strlen(inputCStr)+1,"input name list",0);
    strcpy(dimage->gediIO.inList[0],inputCStr);
  }

  std::string outputStr = Rcpp::as<std::string>(output[0]);
  char* outputCStr = const_cast<char*>(outputStr.c_str());
  strcpy(dimage->outNamen, outputCStr);

  if (coords.length() == 2) {
    dimage->gediRat.coord[0]=coords[0];
    dimage->gediRat.coord[1]=coords[1];
    dimage->gediRat.useOctree=0;    /*no point using octree for single*/
  }

  if (decon) {
    dimage->gediRat.doDecon=1;
  }
  if (indDecon) {
    dimage->gediRat.indDecon=1;
    dimage->gediRat.doDecon=1;
  }

  if (LVIS) {
    pSigma=0.6893;  /*two way trip*/
    fSigma=6.25;
  }

  dimage->gediIO.pSigma=pSigma;
  dimage->gediIO.pFWHM=pFWHM;
  dimage->gediIO.fSigma=fSigma;

  if (readWave) {
    dimage->gediRat.readWave=1;
  }
  if (ground) {
    dimage->gediIO.ground=1;
    dimage->gediRat.cleanOut=1;
  }

  if (sideLobe) {
    dimage->gediRat.sideLobe=1;
  }

  if (listFiles) {
    dimage->listFiles=1;
  }
  dimage->pBuffSize=(uint64_t)(pBuff*1000000000.0);

  if (noNorm) {
    dimage->gediRat.normCover=0;
  }
  if (checkCover){
    dimage->gediRat.checkCover=1;
  }
  if (topHat) {
    dimage->gediRat.topHat=1;
  }
  if (waveID.length() == 1) {
      dimage->useID=1;
      strcpy(dimage->waveID, waveID[0]);
  }
  if (readPulse.length() == 1) {
    dimage->gediIO.readPulse=1;
    strcpy(dimage->gediIO.pulseFile,readPulse[0]);
  }
  if(pulseAfter){
    dimage->gediRat.pulseAfter=1;
  }
  if(pulseBefore){
    dimage->gediRat.pulseAfter=0;
  }

  dimage->gediRat.maxScanAng=maxScanAng;

  if(useShadow){
    dimage->gediRat.useShadow=1;
  }

  if(polyGround){
    dimage->polyGr=1;
  }

  if(nnGround){
    dimage->nnGr=1;
  }

  dimage->gediIO.res=res;

  if(gridBound.length() == 4){
    dimage->gediRat.doGrid=1;
    dimage->useID=1;
    dimage->gediRat.gMinX=gridBound[0];
    dimage->gediRat.gMaxX=gridBound[1];
    dimage->gediRat.gMinY=gridBound[2];
    dimage->gediRat.gMaxY=gridBound[3];
  }
  dimage->gediRat.gRes=gridStep;
  if(keepOld){
    dimage->overWrite=0;
  }
  if(listCoord.length() > 0){
    dimage->gediRat.readALSonce=1;
    dimage->useID=1;
    strcpy(dimage->gediRat.coordList,listCoord[0]);
  }
  if(hdf){
    dimage->writeHDF=1;
  }
  dimage->maxBins=maxBins;
  if(wavefront.length() == 1){
    dimage->gediRat.defWfront=1;
    dimage->gediRat.wavefront=copyFrontFilename(wavefront[0]);
  }

  if(noOctree){
    dimage->gediRat.useOctree=0;
  }

  dimage->gediRat.octLevels=octLevels;

  dimage->gediRat.nOctTop=nOctPix;

  if(countOnly){
    dimage->gediIO.useCount=1;
    dimage->gediIO.useInt=0;
    dimage->gediIO.useFrac=0;
  }

  if(seed.length() == 1){
    srand(seed[0]);
  }

  dimage->gediRat.decimate=decimate;

  /*total number of beams*/
  dimage->gediIO.nTypeWaves=dimage->gediIO.useCount+dimage->gediIO.useFrac+dimage->gediIO.useInt;

  return(dimage);
}/*readCommands*/

#ifdef DEBUG
void gediRat(CharacterVector input,
  CharacterVector output,
  bool            inList          = false,
  bool            ground          = false,
  bool            hdf             = false,
  CharacterVector waveID          = CharacterVector(0),
  NumericVector   coords          = NumericVector(0),
  CharacterVector listCoord       = CharacterVector(0),
  NumericVector   gridBound       = NumericVector(0),
  float           gridStep        = 30.0,
  float           pSigma          = -1.0,
  float           pFWHM           = 15.0,
  CharacterVector readPulse       = CharacterVector(0),
  float           fSigma          = 5.5,
  CharacterVector wavefront       = CharacterVector(0),
  float           res             = 0.15,
  bool            LVIS            = false,
  bool            topHat          = false,
  bool            sideLobe        = false,
  float           lobeAng         = 0.0,
  bool            checkCover      = false,
  float           maxScanAng      = 1000000.0,
  float           decimate        = 1.0,
  uint64_t        pBuff           = 0.2,
  int             maxBins         = 1024,
  bool            countOnly       = false,
  bool            pulseAfter      = false,
  bool            pulseBefore     = false,
  bool            noNorm          = false,
  bool            noOctree        = false,
  int             octLevels       = 0,
  int             nOctPix         = 40,
  bool            decon           = false,
  bool            indDecon        = false,
  bool            readWave        = false,
  bool            listFiles       = false,
  bool            keepOld         = false,
  bool            useShadow       = false,
  bool            polyGround      = false,
  bool            nnGround        = false,
  NumericVector   seed            = NumericVector(0)      );
int main() {
  CharacterVector vec = CharacterVector(1);
  CharacterVector out = CharacterVector(1);
  vec[0] = "E:\\Documentos\\sample.las";
  out[0] = "E:\\Documentos\\saida.txt";

  gediRat(vec, out);
  return 0;
}
#endif

// [[Rcpp::export]]
void gediRat(
  CharacterVector input,
  CharacterVector output,
  bool            inList,
  bool            ground,
  bool            hdf,
  CharacterVector waveID,
  NumericVector   coords,
  CharacterVector listCoord,
  NumericVector   gridBound,
  float           gridStep,
  float           pSigma,
  float           pFWHM,
  CharacterVector readPulse,
  float           fSigma,
  CharacterVector wavefront,
  float           res,
  bool            LVIS,
  bool            topHat,
  bool            sideLobe,
  float           lobeAng,
  bool            checkCover,
  float           maxScanAng,
  float           decimate,
  uint64_t        pBuff,
  int             maxBins,
  bool            countOnly,
  bool            pulseAfter,
  bool            pulseBefore,
  bool            noNorm,
  bool            noOctree,
  int             octLevels,
  int             nOctPix,
  bool            decon,
  bool            indDecon,
  bool            readWave,
  bool            listFiles,
  bool            keepOld,
  bool            useShadow,
  bool            polyGround,
  bool            nnGround,
  NumericVector   seed)
{
  int i=0,j=0;
  rat_control *dimage=NULL;
  lasFile *las=NULL;
  pCloudStruct **data=NULL;
  waveStruct *waves=NULL;
  gediHDF *hdfData=NULL;
  void writeGEDIwave(rat_control *,waveStruct *,int);
  void tidySMoothPulse();
  void checkThisFile(lasFile *,rat_control *,int);
  void groundFromDEM(pCloudStruct **,rat_control *,waveStruct *);
  void checkWaveOverwrite(rat_control *,int);


  /*read command line*/
  dimage=rat_makeControl(
    input,
    output,
    inList,
    ground,
    hdf,
    waveID,
    coords,
    listCoord,
    gridBound,
    gridStep,
    pSigma,
    pFWHM,
    readPulse,
    fSigma,
    wavefront,
    res,
    LVIS,
    topHat,
    sideLobe,
    lobeAng,
    checkCover,
    maxScanAng,
    decimate,
    pBuff,
    maxBins,
    countOnly,
    pulseAfter,
    pulseBefore,
    noNorm,
    noOctree,
    octLevels,
    nOctPix,
    decon,
    indDecon,
    readWave,
    listFiles,
    keepOld,
    useShadow,
    polyGround,
    nnGround,
    seed
  );

  /*set up the pulse*/
  setGediPulse(&dimage->gediIO,&dimage->gediRat);

  /*set up grid or batch if needed*/
  setGediGrid(&dimage->gediIO,&dimage->gediRat);

  /*loop over las files and read*/
  if(!(data=(pCloudStruct **)calloc(dimage->gediIO.nFiles,sizeof(pCloudStruct *)))){
    errorf("error waveStruct allocation.\n");
    exit(1);
  }
  for(i=0;i<dimage->gediIO.nFiles;i++){
    /*report progress if reading all data here*/
    if(dimage->gediRat.doGrid||dimage->gediRat.readALSonce)msgf("File %d of %d",i+1,dimage->gediIO.nFiles);
    /*read lasFile*/
    las=readLasHead(dimage->inList[i],dimage->pBuffSize);

    /*read data or write filename if needed*/
    if(dimage->listFiles==0)data[i]=readALSdata(las,&dimage->gediRat,i);
    else                    checkThisFile(las,dimage,i);
    if(dimage->gediRat.doGrid||dimage->gediRat.readALSonce)msgf(" nPoints %u\n",data[i]->nPoints);

    /*tidy lasFIle*/
    las=tidyLasFile(las);
  }/*file loop*/

  /*set up HDF5 if needed*/
  if(dimage->writeHDF)hdfData=setUpHDF(&dimage->gediIO,&dimage->gediRat,dimage->useID,dimage->waveID,&dimage->hdfCount,dimage->maxBins);

  /*make waveforms*/
  if(dimage->listFiles==0){
    /*loop over waveforms*/
    for(i=0;i<dimage->gediRat.gNx;i++){
      for(j=0;j<dimage->gediRat.gNy;j++){
        /*update centre coord*/
        updateGediCoord(&dimage->gediRat,i,j);

        /*see if that file already exists*/
        checkWaveOverwrite(dimage,i);

        /*if it is not to be overwritten*/
        if(dimage->gediRat.useFootprint){
          /*set up footprint*/
          setGediFootprint(&dimage->gediRat,&dimage->gediIO);

          /*make waveforms*/
          waves=makeGediWaves(&dimage->gediRat,&dimage->gediIO,data);
        }

        /*if it is usable*/
        if(dimage->gediRat.useFootprint){
          /*progress report*/
          if(dimage->writeHDF){
            if((dimage->hdfCount%dimage->gediIO.nMessages)==0){
              msgf("Wave %d of %d\n",i*dimage->gediRat.gNy+j,dimage->gediRat.gNx*dimage->gediRat.gNy);
            }
          }
          /*find the ground if needed*/
          if(dimage->gediIO.ground&&(dimage->polyGr||dimage->nnGr))groundFromDEM(data,dimage,waves);

          /*output results*/
          if(dimage->writeHDF)packGEDIhdf(waves,hdfData,i+j*dimage->gediRat.gNx,&dimage->gediIO,&dimage->gediRat,&dimage->hdfCount,dimage->useID,dimage->waveID);
          else                writeGEDIwave(dimage,waves,i+j*dimage->gediRat.gNx);
        }

        /*tidy up*/
        TIDY(dimage->gediRat.nGrid);
        TIDY(dimage->gediRat.lobe);
        if(waves){
          TTIDY((void **)waves->wave,waves->nWaves);
          TTIDY((void **)waves->canopy,waves->nWaves);
          TTIDY((void **)waves->ground,waves->nWaves);
          TIDY(waves);
        }
      }/*grid y loop*/
    }/*grid x loop*/

    /*write HDF if needed and not blank*/
    if(dimage->writeHDF&&(dimage->hdfCount>0)){
      hdfData->nWaves=dimage->hdfCount;  /*account for unusable footprints*/
      writeGEDIhdf(hdfData,dimage->outNamen,&(dimage->gediIO));
    }
  }/*make and write a waveform if needed*/


  /*tidy up*/
  if(data){
    for(i=0;i<dimage->gediIO.nFiles;i++)data[i]=tidyPointCloud(data[i]);
    TIDY(data);
  }
  hdfData=tidyGediHDF(hdfData);
  if(dimage){
    if(dimage->gediIO.pulse){
      TIDY(dimage->gediIO.pulse->y);
      TIDY(dimage->gediIO.pulse->x);
      TIDY(dimage->gediIO.pulse);
    }
    TTIDY((void **)dimage->gediRat.coords,dimage->gediRat.gNx);
    TTIDY((void **)dimage->gediRat.geoCoords,dimage->gediRat.gNx);
    TTIDY((void **)dimage->gediRat.waveIDlist,dimage->gediRat.gNx);
    TTIDY((void **)dimage->inList,dimage->gediIO.nFiles);
    if(dimage->gediRat.wavefront){
      TTIDY((void **)dimage->gediRat.wavefront->front,dimage->gediRat.wavefront->nX);
      TIDY(dimage->gediRat.wavefront);
    }
    dimage->gediRat.octree=tidyOctree(dimage->gediRat.octree);
    TIDY(dimage);
  }
  tidySMoothPulse();
  return;
}

DataFrame rat_createMetricsDataFrame(rat_control* dimage) {
  // char name[22];
  DataFrame df = DataFrame::create();


  return(df);
}
