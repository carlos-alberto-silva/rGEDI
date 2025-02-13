

/*######################*/
/*# A library for      #*/
/*# handling las files #*/
/*# S Hancock, 2015    #*/
/*######################*/

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




/*###########################################*/
/*lasfile point bit field, 1.3 and earlier*/

typedef struct{
  unsigned char retNumb: 3;
  unsigned char nRet: 3;
  unsigned char sDir: 1;
  unsigned char edge: 1;
}b_field;


/*###########################################*/
/*lasfile point bit field, 1.4*/

typedef struct{
  unsigned char retNumb: 4;
  unsigned char nRet: 4;
  unsigned char classF: 4;
  unsigned char sChan: 2;
  unsigned char sDir: 1;
  unsigned char edge: 1;
}nb_field;


/*###########################################*/
/*a las ARSF file, version 1.0 to 1.3*/

typedef struct{
  /*file pointer*/
  FILE *ipoo;
  char namen[300];   /*filename*/

  /*public header*/
  uint8_t vMajor;    /*version major*/
  uint8_t vMinor;    /*version minor*/
  int headLen;
  uint16_t doy;
  uint16_t year;
  uint16_t headSize;
  uint32_t offsetToP; /*offset to point data*/
  uint32_t nVarRec;   /*number of variable length records*/
  unsigned char pointFormat;   /*point data format ID*/
  uint16_t pRecLen;   /*point record length*/
  uint32_t nPoints;   /*number of point records*/
  uint32_t nPbyRet[7];/*number of point records by return*/
  double posScale[3];
  double posOffset[3];
  double minB[3];
  double maxB[3];
  uint64_t waveStart;

  /*geolocation*/
  uint16_t epsg;       /*projection code*/

  /*point*/
  int32_t x;           /*index to calculate coordinate*/
  int32_t y;           /*index to calculate coordinate*/
  int32_t z;           /*index to calculate coordinate*/
  uint16_t refl;       /*point intensity*/
  b_field field;         /*bit field containing all sorts*/
  nb_field newField;      /*las 1.4 bit field*/

  /*bit field parameters*/
  unsigned char retNumb;
  unsigned char nRet;
  unsigned char classF;
  unsigned char sChan;
  unsigned char sDir;
  unsigned char edge;

  unsigned char classif;  /*point classification*/
  int16_t scanAng;        /*scan angle*/
  unsigned char userData; /*user data*/
  uint16_t psID;          /*point source ID, used for Leica's AGC*/
  double gpsTime;         /*GPS time*/
  uint16_t RGB[3];        /*RGB image*/
  unsigned char packetDes;/*waveform packed description*/
  uint64_t waveMap;       /*pointer to waveform in file*/
  uint32_t waveLen;       /*length of waveform in bins*/
  float time;             /*time in picoseconds of this wave*/
  float grad[3];          /*waveform gradient*/
  /*uint64_t counter;*/    /*I do not know why this is here?*/

  /*buffer to read multiple points at a time*/
  char *pointBuff;       /*buffer with multiple points*/
  uint32_t buffStart;    /*buffer start in number of points*/
  uint32_t buffLen;      /*buffer length in number of points*/
  uint64_t buffByte;     /*buffer length in number of bytes*/
  uint32_t maxBuffLen;   /*max buffer length in number of points*/
  uint64_t maxBuffSize;  /*max buffer length in number of bytes*/

  /*waveform*/
  unsigned char *wave;

  /*GBIC for instruments with variable gain*/
  int gbLen;            /*array length*/
  float *gbic;          /*array of gain factors*/
  float flightBal;      /*balance factor between flights*/
}lasFile;


/*###########################################*/
/*list of las files, set up for Leica specific ALS information*/

typedef struct{
  int nFiles;       /*number of files*/
  char **nameList;  /*file names*/
  char **gbicNamen; /*GBIC input*/
  float *flightBal; /*intensity balance between flights*/

  /*only for use with TLS files*/
  double **scanCent;    /*scan centres*/
  double **align;       /*translation to align non-geolocated files*/
}listLas;


/*#####################################*/
/*structure to hold a point cloud*/

typedef struct{
  uint32_t nPoints;
  double *x;
  double *y;
  double *z;
  double bounds[6];         /*minX minY minZ maxX maxY maxZ*/
  unsigned char *class;     /*point classification*/
  int *refl;
  char *nRet;               /*number of discrete returns per beam*/
  char *retNumb;            /*this point's return number*/
  int16_t *scanAng;         /*scan angle*/
  unsigned char *packetDes; /*waveform or not*/
  float **grad;             /*Poynting vector*/
  float *time;              /*time in picoseconds of this wave*/
  uint64_t *waveMap;        /*pointer to waveform in file*/
  uint32_t *waveLen;        /*length of waveform in bins*/
  uint64_t waveStart;       /*offset to waveform data*/
  FILE *ipoo;               /*file pointer*/
  char hasWave;             /*waveform included*/
  float *gap;               /*gap fraction up to a point*/
  float *range;             /*range fram scan centre to this point*/
}pCloudStruct;


/*###########################################*/
/*function definitions*/

lasFile *readLasHead(char *,uint64_t);
void readLasGeo(lasFile *);
void readLasPoint(lasFile *,uint32_t);
unsigned char *readLasWave(uint64_t,int32_t,FILE *,uint64_t);
char **readInList(int *,char *);
lasFile *tidyLasFile(lasFile *);
void setCoords(double *,double *,double *,lasFile *);
void binPosition(double *,double *,double *,int,double,double,double,float,float *);
char checkOneWave(lasFile *);
char checkFileBounds(lasFile *,double,double,double,double);
void readGBIC(char,char,lasFile **,listLas *);
lasFile **lfalloc(int);
listLas *readLasList(char *);
void tidyListLas(listLas *);
double **readCoordList(int,char **,char *);
pCloudStruct *tidyPointCloud(pCloudStruct *);

/*the end*/
/*###########################################*/

