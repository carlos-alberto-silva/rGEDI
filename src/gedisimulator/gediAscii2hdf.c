#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "hdf5.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "libLasRead.h"
#include "libLasProcess.h"
#include "libLidarHDF.h"
#include "libOctree.h"
#include "gediIO.h"

/*tolerances*/
#define TOL 0.00001

/*####################################*/
/*control structure*/

typedef struct{
  /*input/output*/
  gediIOstruct gediIO; /*input/output structure*/
  char outNamen[1000];
  int maxGauss;     /*maximum number of Gaussians for output*/

  /*switches*/
  char useBounds;    /*when we will process only a subset of bounds*/

  /*noise parameters*/
  float meanN;
  float nSig;
  float bThresh;   /*bounds threshold*/
  float hNoise;    /*hard threshold noise as a fraction of integral*/
  char missGround; /*force to miss ground to get RH errors*/
  float minGap;    /*minimum detectavle gap fraction*/
  char linkNoise;  /*use link noise or not*/
  float linkM;     /*link margin*/
  float linkCov;   /*cover at which link margin is defined*/
  float linkSig;   /*link noise sigma*/
  float linkFsig;  /*footprint sigma used for link margin*/
  float linkPsig;  /*pulse sigma used for link margin*/
  float trueSig;   /*true noise sigma in DN*/
  float deSig;     /*detector sigma*/
  char bitRate;   /*digitiser bit rate*/
  float maxDN;    /*maximum DN we need to digitise*/
  float offset;   /*waveform DN offset*/

  /*bounds for subsets*/
  double minX;
  double maxX;
  double minY;
  double maxY;
}control;


/*#####################################*/
/*main*/

int main(int argc,char **argv)
{
  control *dimage=NULL;
  control *readCommands(int,char **);
  dataStruct **data=NULL;
  dataStruct **readMultiData(control *);
  gediHDF *hdfData=NULL;


  /*read command line*/
  dimage=readCommands(argc,argv);

  /*read data*/
  data=readMultiData(dimage);

  /*copy into HDF structure*/
  hdfData=arrangeGEDIhdf(data,&(dimage->gediIO));

  /*tidy data array*/
  data=tidyAsciiStruct(data,dimage->gediIO.nFiles);

  /*write HDF5*/
  writeGEDIhdf(hdfData,dimage->outNamen);

  /*tidy arrays*/
  hdfData=tidyGediHDF(hdfData);
  if(dimage){
    TTIDY((void **)dimage->gediIO.inList,dimage->gediIO.nFiles);
    TIDY(dimage->gediIO.gFit);
    TIDY(dimage->gediIO.den);
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*####################################################*/
/*read multiple data files*/

dataStruct **readMultiData(control *dimage)
{
  int i=0;
  dataStruct **data=NULL;

  /*allocate space for all*/
  if(!(data=(dataStruct **)calloc(dimage->gediIO.nFiles,sizeof(dataStruct *)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }

  /*read data*/
  for(i=0;i<dimage->gediIO.nFiles;i++)data[i]=readASCIIdata(dimage->gediIO.inList[i],&(dimage->gediIO));

  return(data);
}/*readMultiData*/


/*####################################################*/
/*read command line*/

control *readCommands(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;
  void setDenoiseDefault(denPar *);
  void readPulse(denPar *);

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
  dimage->gediIO.inList=chChalloc(dimage->gediIO.nFiles,"inList",0);
  dimage->gediIO.inList[0]=challoc(200,"inList",0);
  dimage->gediIO.useInt=1;
  dimage->gediIO.useCount=1;
  dimage->gediIO.useFrac=1;
  strcpy(&(dimage->gediIO.inList[0][0]),"/Users/stevenhancock/data/gedi/simulations/USDA_CO/gedi.USDA_CO.4112.wave");
  strcpy(&(dimage->outNamen[0]),"teast.h5");
  dimage->gediIO.ground=1;


  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        TTIDY((void **)dimage->gediIO.inList,dimage->gediIO.nFiles);
        dimage->gediIO.inList=NULL;
        dimage->gediIO.nFiles=1;
        dimage->gediIO.inList=chChalloc(dimage->gediIO.nFiles,"input name list",0);
        dimage->gediIO.inList[0]=challoc((uint64_t)strlen(argv[++i])+1,"input name list",0);
        strcpy(dimage->gediIO.inList[0],argv[i]);
      }else if(!strncasecmp(argv[i],"-inList",7)){
        checkArguments(1,i,argc,"-inList");
        TTIDY((void **)dimage->gediIO.inList,dimage->gediIO.nFiles);
        dimage->gediIO.inList=readInList(&dimage->gediIO.nFiles,argv[++i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(dimage->outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n#####\nProgram to convert ASCII GEDI waveforms to HDF5\n#####\n\n-input name;     single input filename\n-output name;    output filename\n-inList list;    input file list for multiple files\n\n");
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
/*#####################################*/

