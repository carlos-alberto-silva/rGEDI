

/*#######################*/
/*# A library for       #*/
/*# handling LVIS files #*/
/*# S Hancock, 2017     #*/
/*#######################*/

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




/*#######################################*/
/*LVIS LGW v1.03 and below structure*/

typedef struct{
   uint32_t lfid;       /* LVIS file identifier*/
   uint32_t shotN; /* LVIS shotnumber*/
   float az;            /* true heading from the aircraft to the ground (degrees)*/
   float zen;      /*zenith angle (deghrees)*/
   float range;      /* range from the aircraft to the ground (meters)*/
   double lvistime;  /* LVIS recorded UTC time (seconds of the day) when the shot was acquired*/
   double lon0;      /* longitude of the highest sample of the waveform (degrees east)*/
   double lat0;      /* latitude of the highest sample of the waveform (degrees north)*/
   float  z0;        /* elevation of the highest sample of the waveform (m)*/
   double lon431;    /* longitude of the lowest sample of the waveform (degrees east)*/
   double lat431;    /* latitude of the lowest sample of the waveform (degrees north)*/
   float  z431;      /* elevation of the lowest sample of the waveform (m)*/
   float  sigmean;   /* signal mean noise level, calculated in-flight (counts)*/
   unsigned char *txwave; /* transmit waveform, recorded in-flight (counts)*/
   unsigned char *rxwave; /* return   waveform, recorded in-flight (counts)*/
   uint16_t *txwave4;     /* transmit waveform, recorded in-flight (counts) version on 1.4+*/
   uint16_t *rxwave4;     /* return   waveform, recorded in-flight (counts) version 1.4+*/
}lvisLGWdata;


/*#######################################*/
/*LVIS overall structure*/

typedef struct{
  int verMaj;      /*major version*/
  int verMin;      /*minor version*/
  int nWaves;      /*number of waveforms*/
  int nBins;       /*number of rx waveform bins*/
  int nTxBins;     /*number of tx waveform bins*/
  FILE *ipoo;      /*input file*/
  lvisLGWdata *data;  /*data pointer*/
  char byteord;    /*byte order of this computer*/
}lvisLGWstruct;


/*#####################################*/
/*LVIS HDF5 structure*/

typedef struct{
  int nWaves;   /*number of waveforms*/
  int nBins;    /*number of waveform bins*/
  int pBins;     /*number of pulse bins*/
  /*data per wave*/
  double *lon0;       /*LON0*/
  double *lat0;       /*LAT0*/
  double *lon1023;    /*LON1023*/
  double *lat1023;    /*LAT1023*/
  uint32_t *lfid;     /*LFID*/
  uint32_t *shotN;    /*SHOTNUMBER*/
  uint16_t **wave;    /*RXWAVE*/
  uint16_t **pulse;   /*TXWAVE*/
  float *zen;         /*INCIDENTANGLE*/
  float *z0;         /*Z0*/
  float *z1023;      /*Z1023*/
  float *sigmean;     /*SIGMEAN*/
  double *time;       /*TIME*/
}lvisHDF;


/*#######################################*/
/*lgw types fpr reading*/

/*#pragma pack(1)*/
struct lvis_lgw_v1_00{
   double lon0;
   double lat0;
   float  z0;
   double lon431;
   double lat431;
   float  z431;
   float  sigmean;
   unsigned char wave[432];
};
typedef struct lgw_v1_00 * ptr_lgw_v1_00;

/*#pragma pack(1)*/
struct lvis_lgw_v1_01{
   uint32_t lfid;
   uint32_t shotnumber;
   double lon0;
   double lat0;
   float  z0;
   double lon431;
   double lat431;
   float  z431;
   float  sigmean;
   unsigned char wave[432];
};
typedef struct lgw_v1_01 * ptr_lgw_v1_01;

/*#pragma pack(1)*/
struct lvis_lgw_v1_02{
   uint32_t lfid;
   uint32_t shotnumber;
   double lvistime;
   double lon0;
   double lat0;
   float  z0;
   double lon431;
   double lat431;
   float  z431;
   float  sigmean;
   unsigned char wave[432];
};
typedef struct lgw_v1_02 * ptr_lgw_v1_02;

/*#pragma pack(1)*/
struct lvis_lgw_v1_03{
   uint32_t lfid;
   uint32_t shotnumber;
   float azimuth;
   float incidentangle;
   float range;
   double lvistime;
   double lon0;
   double lat0;
   float  z0;
   double lon431;
   double lat431;
   float  z431;
   float  sigmean;
   unsigned char txwave[80];
   unsigned char rxwave[432];
};
typedef struct lgw_v1_03 * ptr_lgw_v1_03;

/*#pragma pack(1)*/
struct lvis_lgw_v1_04{
   uint32_t lfid;
   uint32_t shotnumber;
   float azimuth;
   float incidentangle;
   float range;
   double lvistime;
   double lon0;
   double lat0;
   float  z0;
   double lon527;
   double lat527;
   float  z527;
   float  sigmean;
   uint16_t txwave[120];
   uint16_t rxwave[528];
};
typedef struct lgw_v1_04 * ptr_lgw_v1_04;


/*#######################################*/
/*functions*/

lvisLGWdata *readLVISlgw(char *,lvisLGWstruct *);
void checkLVISsizes();
lvisHDF *tidyLVISstruct(lvisHDF *);
lvisHDF *readLVIShdf(char *);
void write1dDoubleHDF5(hid_t,char *,double *,int);
void write1dFloatHDF5(hid_t,char *,float *,int);
void write2dFloatHDF5(hid_t,char *,float *,int,int);
void write2dCharHDF5(hid_t,char *,char *,int,int);
void write1dIntHDF5(hid_t,char *,int *,int);
void write1dInt64HDF5(hid_t,char *,int64_t *,int);
void write1dUint32HDF5(hid_t,char *,uint32_t *,int);
void write1dUint16HDF5(hid_t,char *,uint16_t *,int);
void write2dUint16HDF5(hid_t,char *,uint16_t *,int,int);
void writeComp2dInt8HDF5(hid_t,char *,int8_t *,int,int);
void writeComp2dFloatHDF5(hid_t,char *,float *,int,int);
void writeComp1dFloatHDF5(hid_t,char *,float *,int);
void writeComp1dUint8HDF5(hid_t,char *,uint8_t *,int);
void writeComp1dUint16HDF5(hid_t,char *,uint16_t *,int);
void writeComp1dUint32HDF5(hid_t,char *,uint32_t *,int);
void writeComp1dUint64HDF5(hid_t,char *,uint64_t *,int);
void writeComp1dInt8HDF5(hid_t,char *,int8_t *,int);
void writeComp1dInt32HDF5(hid_t,char *,int32_t *,int);
void writeComp1dDoubleHDF5(hid_t,char *,double *,int);
float *read1dFloatHDF5(hid_t,char *,int *);
float **read2dFloatHDF5(hid_t,char *,int *,int *);
double *read1dDoubleHDF5(hid_t,char *,int *);
uint8_t *read1dUint8HDF5(hid_t,char *,int *);
uint16_t *read1dUint16HDF5(hid_t,char *,int *);
uint32_t *read1dUint32HDF5(hid_t,char *,int *);
uint64_t *read1dUint64HDF5(hid_t,char *,int *);
int *read1dIntHDF5(hid_t,char *,int *); 
float *read15dFloatHDF5(hid_t,char *,int *,int *);
char *read15dCharHDF5(hid_t,char *,int *,int *);
uint16_t **read2dUint16HDF5(hid_t,char *,int *,int *);

/*the end*/
/*####################################################*/

