#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "hdf5.h"
#include "tools.h"
#include "tools.c"
#include "libLasRead.h"
#include "libLasProcess.h"
#include "libLidarHDF.h"
#include "libOctree.h"
#include "gediIO.h"


/*############################*/
/*# Converts LVIS lgw format #*/
/*# format into HDF5    2017 #*/
/*# svenhancock@gmail.com    #*/
/*############################*/

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




/*####################################*/
/*control structure*/

typedef struct{
  char inNamen[1000];   /*input filename*/
  char outNamen[1000];  /*output filename*/
}control;


/*###########################################################*/
/*main*/

int main(int argc,char **argv)
{
  control *dimage=NULL;
  control *readCommands(int,char **);
  lvisLGWstruct lgw;
  lvisHDF *hdf=NULL;
  lvisHDF *copyLVISdata(lvisLGWstruct *);
  void tidyLGW(lvisLGWstruct *);
  void writeHDFlvis(lvisHDF *,char *);

  /*read command Line*/
  dimage=readCommands(argc,argv);

  /*initialise*/
  lgw.ipoo=NULL;
  lgw.data=NULL;

  /*read LVIS data*/
  lgw.data=readLVISlgw(dimage->inNamen,&lgw);

  /*pack into HDF5 structure*/
  hdf=copyLVISdata(&lgw);
  tidyLGW(&lgw);

  /*write HDF5*/
  writeHDFlvis(hdf,dimage->outNamen);

  /*tidy up*/
  hdf->nWaves=1;  /*to trick it into the reorganised arrays for writing*/
  hdf=tidyLVISstruct(hdf);
  if(dimage){
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*###########################################################*/
/*write HDF data*/

void writeHDFlvis(lvisHDF *hdf,char *outNamen)
{
  hid_t file;         /* Handles */

  file=H5Fcreate(outNamen,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

  /*write data*/
  write1dDoubleHDF5(file,"LON0",hdf->lon0,hdf->nWaves);
  write1dDoubleHDF5(file,"LAT0",hdf->lat0,hdf->nWaves);
  write1dDoubleHDF5(file,"LON1023",hdf->lon1023,hdf->nWaves);
  write1dDoubleHDF5(file,"LAT1023",hdf->lat1023,hdf->nWaves);
  write1dDoubleHDF5(file,"TIME",hdf->time,hdf->nWaves);
  write1dFloatHDF5(file,"INCIDENTANGLE",hdf->zen,hdf->nWaves);
  write1dFloatHDF5(file,"Z0",hdf->z0,hdf->nWaves);
  write1dFloatHDF5(file,"SIGMEAN",hdf->sigmean,hdf->nWaves);
  write1dFloatHDF5(file,"Z1023",hdf->z1023,hdf->nWaves);
  write1dUint32HDF5(file,"LFID",hdf->lfid,hdf->nWaves);
  write1dUint32HDF5(file,"SHOTNUMBER",hdf->shotN,hdf->nWaves);
  write2dUint16HDF5(file,"RXWAVE",hdf->wave[0],hdf->nWaves,hdf->nBins);
  if(hdf->pulse)write2dUint16HDF5(file,"TXWAVE",hdf->pulse[0],hdf->nWaves,hdf->pBins);

  /*close file*/
  if(H5Fclose(file)){
    fprintf(stderr,"Issue closing file\n");
    exit(1);
  }
  fprintf(stdout,"Waveforms to %s\n",outNamen);
  return;
}/*writeHDFlvis*/


/*###########################################################*/
/*copy data from LGW to HDF5 structure*/

lvisHDF *copyLVISdata(lvisLGWstruct *lgw)
{
  int i=0,j=0;
  lvisHDF *hdf=NULL;
  float dz=0,scale=0;
  double dx=0,dy=0;
  char usePulse=0;

  /*is there a pulse?*/
  if(lgw->data[0].txwave)usePulse=1;
  else                   usePulse=0;

  /*allocate structure*/
  if(!(hdf=(lvisHDF *)calloc(1,sizeof(lvisHDF)))){
    fprintf(stderr,"error in HDF structure allocation.\n");
    exit(1);
  }
  hdf->nWaves=lgw->nWaves;
  hdf->nBins=1024; /*lgw->nBins; must be 1023 to match new LVIS*/
  hdf->pBins=lgw->nTxBins;
  fprintf(stdout,"Bins %d %d\n",lgw->nWaves,lgw->nBins);

  /*allocate arrays*/
  hdf->lon0=dalloc(hdf->nWaves,"lon0",0);
  hdf->lat0=dalloc(hdf->nWaves,"lat0",0);
  hdf->lon1023=dalloc(hdf->nWaves,"lon0",0);
  hdf->lat1023=dalloc(hdf->nWaves,"lat0",0);
  if(!(hdf->lfid=(uint32_t *)calloc(hdf->nWaves,sizeof(uint32_t)))){
    fprintf(stderr,"error in lfid allocation, allocating %ld.\n",hdf->nWaves*sizeof(uint32_t));
    exit(1);
  }
  if(!(hdf->shotN=(uint32_t *)calloc(hdf->nWaves,sizeof(uint32_t)))){
    fprintf(stderr,"error in shotN allocation, allocating %ld\n",hdf->nWaves*sizeof(uint32_t));
    exit(1);
  }
  if(!(hdf->wave=(uint16_t **)calloc(1,sizeof(uint16_t *)))){
    fprintf(stderr,"error in wave array allocation, allocating 1.\n");
    exit(1);
  }
  if(!(hdf->wave[0]=(uint16_t *)calloc((uint64_t)hdf->nWaves*(uint64_t)hdf->nBins,(uint64_t)sizeof(uint16_t)))){
    fprintf(stderr,"error in wave array allocation, allocating %ld.\n",(uint64_t)hdf->nWaves*(uint64_t)hdf->nBins*(uint64_t)sizeof(uint16_t));
    exit(1);
  }
  if(usePulse){
    if(!(hdf->pulse=(uint16_t **)calloc(1,sizeof(uint16_t *)))){
      fprintf(stderr,"error pulse array allocation, allocating 1.\n");
      exit(1);
    }
    if(!(hdf->pulse[0]=(uint16_t *)calloc((uint64_t)hdf->nWaves*(uint64_t)hdf->pBins,(uint64_t)sizeof(uint16_t)))){
      fprintf(stderr,"error pulse array allocation, allocating %ld.\n",(uint64_t)hdf->nWaves*(uint64_t)hdf->pBins*(uint64_t)sizeof(uint16_t));
      exit(1);
    }
  }


  hdf->zen=falloc((uint64_t)hdf->nWaves,"zen",0);
  hdf->z0=falloc((uint64_t)hdf->nWaves,"z0",0);
  hdf->z1023=falloc((uint64_t)hdf->nWaves,"z1023",0);
  hdf->sigmean=falloc((uint64_t)hdf->nWaves,"sigmean",0);
  hdf->time=dalloc(hdf->nWaves,"time",0);


  /*copy data*/
  for(i=0;i<lgw->nWaves;i++){
    hdf->lfid[i]=lgw->data[i].lfid;
    hdf->shotN[i]=lgw->data[i].shotN;
    hdf->zen[i]=lgw->data[i].zen;
    hdf->time[i]=lgw->data[i].lvistime;
    hdf->lon0[i]=lgw->data[i].lon0;
    hdf->lat0[i]=lgw->data[i].lat0;
    hdf->z0[i]=lgw->data[i].z0;
    hdf->sigmean[i]=lgw->data[i].sigmean;

    /*to match new LVIS bin numbers, extrapolate to 1024 bins*/
    dx=lgw->data[i].lon431-lgw->data[i].lon0;
    dy=lgw->data[i].lat431-lgw->data[i].lat0;
    dz=lgw->data[i].z431-lgw->data[i].z0;
    scale=(float)hdf->nBins/(float)lgw->nBins;

    hdf->lon1023[i]=lgw->data[i].lon0+(double)scale*dx;
    hdf->lat1023[i]=lgw->data[i].lat0+(double)scale*dy;
    hdf->z1023[i]=lgw->data[i].z0+scale*dz;

    /*copy and reecast waveform*/
    for(j=0;j<lgw->nBins;j++){
      if(lgw->verMin<4)hdf->wave[0][i*hdf->nBins+j]=(uint16_t)lgw->data[i].rxwave[j];
      else             hdf->wave[0][i*hdf->nBins+j]=lgw->data[i].rxwave4[j];
    }
    /*pad end of waveform to make 1024 bins*/
    for(j=lgw->nBins;j<hdf->nBins;j++)hdf->wave[0][i*hdf->nBins+j]=(uint16_t)(hdf->sigmean[i]-0.5);

    /*copy pulse if needed*/
    if(usePulse){
      for(j=0;j<hdf->pBins;j++){
        if(lgw->verMin<4)hdf->pulse[0][i*hdf->pBins+j]=(uint16_t)lgw->data[i].txwave[j];
        else             hdf->pulse[0][i*hdf->pBins+j]=lgw->data[i].txwave4[j];
      }
    }
  }/*wave loop*/

  return(hdf);
}/*copyLVISdata*/


/*###########################################################*/
/*tidy LGW structure*/

void tidyLGW(lvisLGWstruct *lgw)
{
  int i=0;

  if(lgw->ipoo){
    fclose(lgw->ipoo);
    lgw->ipoo=NULL;
  }
  if(lgw->data){
    for(i=0;i<lgw->nWaves;i++){
      TIDY(lgw->data[i].rxwave);
      TIDY(lgw->data[i].txwave);
      TIDY(lgw->data[i].rxwave4);
      TIDY(lgw->data[i].txwave4);
    }
    TIDY(lgw->data);
  }
  return;
}/*tidyLGW*/


/*###########################################################*/
/*read command Line*/

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


  /*defaults*/
  strcpy(dimage->outNamen,"teast.h5");
  strcpy(dimage->inNamen,"/Users/dill/data/teast/lgw/LVIS_US_NH_2009_VECT_20100328.subset.lgw");

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        strcpy(dimage->inNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(dimage->outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n#####\nProgram to convert LVIS lgw files to HDF5\n#####\n\n-input name;     input LGW filename\n-output name;   output HDF5 filename\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry lgw2hdf -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }

  return(dimage);
}/*readCommands*/


/*the end*/
/*###########################################################*/

