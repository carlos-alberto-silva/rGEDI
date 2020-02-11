#include <R.h>
#include <Rinternals.h>

#include <gediIO.h>
#include <gediRat.h>
#include <tools.h>

#define msgf Rprintf
#define errorf REprintf
#define exit(n) return(ScalarInteger(n))

rat_control* rat_makeControl(// Input output filenames and format
  const char*     input,
  const char*     output,
  const char*     inList,
  char            ground,
  char            hdf,
  const char*     waveID,

  // Single footprint, list of footprints, or grid of footprints
  double*         coords,
  const char*     listCoord,
  double*         gridBound,
  float           gridStep,

  // Lidar characteristics. Defaults are expected GEDI values.
  float           pSigma,
  float           pFWHM,
  const char*     readPulse,
  float           fSigma,
  const char*     wavefront,
  float           res,
  char            LVIS,
  char            topHat,
  char            sideLobe,
  float           lobeAng,


  // Input data quality filters
  char            checkCover,
  float           maxScanAng,
  float           decimate,

  // Computational speed options
  uint64_t        pBuff,
  int             maxBins,
  char            countOnly,
  char            pulseAfter,
  char            pulseBefore,
  char            noNorm,

  // Octree
  char            noOctree,
  int             octLevels,
  int             nOctPix,

  // Using full-waveform input data (not tested)
  char            decon,
  char            indDecon,
  char            readWave,

  // Miscellaneous
  char            listFiles,
  char            keepOld,
  char            useShadow,
  char            polyGround,
  char            nnGround,
  uint32_t*       seed);


void R_message(const char *txt) {
  SEXP r_msg = install("message");
  SEXP call_msg = PROTECT(lang2(r_msg, ScalarString(mkChar(txt))));
  eval(call_msg, R_BaseEnv);
  UNPROTECT(1);
}

SEXP C_getGediPulse(
    SEXP r_pSigma,
    SEXP r_pRes,
    SEXP r_pFWHM,
    SEXP r_iThresh,
    SEXP r_linkPsig,
    SEXP r_readPulse,
    SEXP r_pulseFile
) {
  const char* pulseFile = CHAR(asChar(r_pulseFile));
  char readPulse = (char)asLogical(r_readPulse);
  float pSigma = (float)asReal(r_pSigma);
  float pRes = (float)asReal(r_pRes);
  float pFWHM = (float)asReal(r_pFWHM);
  float iThresh = (float)asReal(r_iThresh);
  float linkPsig = 0.0f;
  pulseStruct *pulse = NULL;

  // Data frame setup
  SEXP result = PROTECT(allocVector(VECSXP, 5));

  setGediPulse(
      &pulse,
      &pSigma,
      &pRes,
      &linkPsig,
      readPulse,
      (char*)pulseFile,
      pFWHM,
      iThresh);

  SEXP x = PROTECT(allocVector(REALSXP, pulse->nBins));
  SEXP y = PROTECT(allocVector(REALSXP, pulse->nBins));
  for (int i=0; i<pulse->nBins; i++) {
    REAL(x)[i] = (double)pulse->x[i];
    REAL(y)[i] = (double)pulse->y[i];
  }
  free(pulse->x);
  free(pulse->y);
  free(pulse);


  SEXP res_pSigma = PROTECT(allocVector(REALSXP, 1));
  REAL(res_pSigma)[0] = (double)pSigma;

  SEXP res_pRes = PROTECT(allocVector(REALSXP, 1));
  REAL(res_pRes)[0] = (double)pRes;

  SEXP res_linkPsig = PROTECT(allocVector(REALSXP, 1));
  REAL(res_linkPsig)[0] = (double)linkPsig;

  SET_VECTOR_ELT(result, 0, x);
  SET_VECTOR_ELT(result, 1, y);
  SET_VECTOR_ELT(result, 2, res_pSigma);
  SET_VECTOR_ELT(result, 3, res_pRes);
  SET_VECTOR_ELT(result, 4, res_linkPsig);

  SEXP nam = PROTECT(allocVector(STRSXP, 5)); // names attribute (column names)
  SET_STRING_ELT(nam, 0, mkChar("x"));
  SET_STRING_ELT(nam, 1, mkChar("y"));
  SET_STRING_ELT(nam, 2, mkChar("pSigma"));
  SET_STRING_ELT(nam, 3, mkChar("pRes"));
  SET_STRING_ELT(nam, 4, mkChar("linkPsig"));
  namesgets(result, nam);

  UNPROTECT(7);

  return result;
}


SEXP C_gedisimulate(
  SEXP input,
  SEXP output,
  SEXP inList,
  SEXP ground,
  SEXP hdf,
  SEXP waveID,
  SEXP coords,
  SEXP listCoord,
  SEXP gridBound,
  SEXP gridStep,
  SEXP pSigma,
  SEXP pFWHM,
  SEXP readPulse,
  SEXP fSigma,
  SEXP wavefront,
  SEXP res,
  SEXP LVIS,
  SEXP topHat,
  SEXP sideLobe,
  SEXP lobeAng,
  SEXP checkCover,
  SEXP maxScanAng,
  SEXP decimate,
  SEXP pBuff,
  SEXP maxBins,
  SEXP countOnly,
  SEXP pulseAfter,
  SEXP pulseBefore,
  SEXP noNorm,
  SEXP noOctree,
  SEXP octLevels,
  SEXP nOctPix,
  SEXP decon,
  SEXP indDecon,
  SEXP readWave,
  SEXP listFiles,
  SEXP keepOld,
  SEXP useShadow,
  SEXP polyGround,
  SEXP nnGround,
  SEXP seed
) {
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
  
  double* c_coords = NULL;
  if (! isNull(coords)) {
    c_coords = malloc(2*sizeof(c_coords));
    c_coords[0] = REAL(coords)[0];
    c_coords[1] = REAL(coords)[1];
  }
  double* c_gridBound = NULL;
  if (! isNull(gridBound)) {
    c_gridBound = malloc(4*sizeof(c_gridBound));
    c_gridBound[0] = REAL(gridBound)[0];
    c_gridBound[1] = REAL(gridBound)[1];
    c_gridBound[2] = REAL(gridBound)[2];
    c_gridBound[3] = REAL(gridBound)[3];
  }
   
  /*read command line*/
  dimage=rat_makeControl(
    CHAR(asChar(input)),
    CHAR(asChar(output)),
    isNull(inList) ? NULL : CHAR(asChar(inList)),
    (char)asLogical(ground),
    (char)asLogical(hdf),
    isNull(waveID) ? NULL : CHAR(asChar(waveID)),
    c_coords,
    isNull(listCoord) ? NULL : CHAR(asChar(listCoord)),
    c_gridBound,
    (float)asReal(gridStep),

    (float)asReal(pSigma),
    (float)asReal(pFWHM),
    isNull(readPulse) ? NULL : CHAR(asChar(readPulse)),
    (float)asReal(fSigma),
    isNull(wavefront) ? NULL : CHAR(asChar(wavefront)),
    (float)asReal(res),
    (char)asLogical(LVIS),
    (char)asLogical(topHat),
    (char)asLogical(sideLobe),
    (float)asReal(lobeAng),

    (char)asLogical(checkCover),
    (float)asReal(maxScanAng),
    (float)asReal(decimate),

    (uint64_t)asInteger(pBuff),
    asInteger(maxBins),
    (char)asLogical(countOnly),
    (char)asLogical(pulseAfter),
    (char)asLogical(pulseBefore),
    (char)asLogical(noNorm),

    (char)asLogical(noOctree),
    asInteger(octLevels),
    asInteger(nOctPix),

    (char)asLogical(decon),
    (char)asLogical(indDecon),
    (char)asLogical(readWave),

    (char)asLogical(listFiles),
    (char)asLogical(keepOld),
    (char)asLogical(useShadow),
    (char)asLogical(polyGround),
    (char)asLogical(nnGround),
    isNull(seed) ? NULL : (uint32_t)asInteger(seed)
  );

  /*set up the pulse*/
  setGediPulse(&dimage->gediIO.pulse,
    &dimage->gediIO.pSigma,
    &dimage->gediIO.pRes,
    &dimage->gediIO.linkPsig,
    dimage->gediIO.readPulse,
    dimage->gediIO.pulseFile,
    dimage->gediIO.pFWHM,
    dimage->gediRat.iThresh);

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
              msgf("Wave %d of %d\r",i*dimage->gediRat.gNy+j,dimage->gediRat.gNx*dimage->gediRat.gNy);
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
  return ScalarInteger(0);
}/*main*/

rat_control* rat_makeControl(// Input output filenames and format
  const char*     input,
  const char*     output,
  const char*     inList,
  char            ground,
  char            hdf,
  const char*     waveID,

  // Single footprint, list of footprints, or grid of footprints
  double*         coords,
  const char*     listCoord,
  double*         gridBound,
  float           gridStep,

  // Lidar characteristics. Defaults are expected GEDI values.
  float           pSigma,
  float           pFWHM,
  const char*     readPulse,
  float           fSigma,
  const char*     wavefront,
  float           res,
  char            LVIS,
  char            topHat,
  char            sideLobe,
  float           lobeAng,


  // Input data quality filters
  char            checkCover,
  float           maxScanAng,
  float           decimate,

  // Computational speed options
  uint64_t        pBuff,
  int             maxBins,
  char            countOnly,
  char            pulseAfter,
  char            pulseBefore,
  char            noNorm,

  // Octree
  char            noOctree,
  int             octLevels,
  int             nOctPix,

  // Using full-waveform input data (not tested)
  char            decon,
  char            indDecon,
  char            readWave,

  // Miscellaneous
  char            listFiles,
  char            keepOld,
  char            useShadow,
  char            polyGround,
  char            nnGround,
  uint32_t*       seed)
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
  if (inList != NULL) {
    TTIDY((void **)dimage->inList,dimage->gediIO.nFiles);
    dimage->inList=readInList(&dimage->gediIO.nFiles, inList);
  } else {
    TTIDY((void **)dimage->inList,dimage->gediIO.nFiles);
    dimage->gediIO.nFiles=1;
    dimage->inList=chChalloc(dimage->gediIO.nFiles,"input name list",0);
    dimage->inList[0]=challoc((uint64_t)strlen(input)+1,"input name list",0);
    strcpy(dimage->inList[0],input);
  }
  strcpy(dimage->outNamen, output);

  if (coords != NULL) {
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
  if (waveID != NULL) {
      dimage->useID=1;
      strcpy(dimage->waveID, waveID[0]);
  }
  if (readPulse != NULL) {
    dimage->gediIO.readPulse=1;
    strcpy(dimage->gediIO.pulseFile,readPulse);
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

  if(gridBound != NULL){
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
  if(listCoord != NULL){
    dimage->gediRat.readALSonce=1;
    dimage->useID=1;
    strcpy(dimage->gediRat.coordList,listCoord);
  }
  if(hdf){
    dimage->writeHDF=1;
  }
  dimage->maxBins=maxBins;
  if(wavefront != NULL){
    dimage->gediRat.defWfront=1;
    dimage->gediRat.wavefront=copyFrontFilename(wavefront);
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

  if(seed != NULL){
    srand(seed[0]);
  }

  dimage->gediRat.decimate=decimate;

  /*total number of beams*/
  dimage->gediIO.nTypeWaves=dimage->gediIO.useCount+dimage->gediIO.useFrac+dimage->gediIO.useInt;

  return(dimage);
}

#undef exit
#undef msgf
#undef errorf
