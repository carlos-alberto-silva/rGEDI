/*########################*/
/*# structures for voxel #*/
/*# lidar programs       #*/
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



/*###################################################*/
/*voxels*/

typedef struct{
  int nScans;        /*number of contributing files*/
  float **hits;      /*beams that hit voxel, per scan location*/
  float **miss;      /*beams blocked before voxel, per scan location*/
  float **inHit;     /*hits within voxel, per scan location*/
  float **inMiss;    /*missed within voxel, per scan location*/
  float **sampVol;   /*volume of voxel sampled, per scan location*/
  float **totVol;    /*volume of voxel passed through, per scan location*/
  float **sumRsq;    /*sum of radius squared, for point area*/
  float **meanRefl;  /*mean reflectances from all returns*/
  float **meanZen;   /*mean zenith angle of all beams intersecting*/
  float **meanRange; /*mean range of all beams intersecting*/
  float *rmse;       /*rmse of signal going in*/
  int *contN;        /*for normalising ALS*/
  int nVox;          /*total number of voxels*/
  int nX;            /*number of voxels*/
  int nY;            /*number of voxels*/
  int nZ;            /*number of voxels*/
  double res[3];     /*voxel resolution*/
  double bounds[6];  /*voxel bounds minX minY minZ maxX maxY maxZ*/
  float volume;      /*volume of the voxel*/
  /*derived values, not always used*/
  float **gap;       /*gap fraction within voxel*/
  float **gapTo;     /*gap fraction to a voxel*/
  float **PAIb;      /*PAI according to Beland*/
  float **PAIrad;    /*PAI from adding up spheres*/
  /*switches*/
  char useRMSE;      /*switch to save RAM in voxelate*/
  char savePts;      /*save points switch*/
  float maxZen;      /*maximum absolute zenith to allow. To filter tilt mount*/
  /*DTM for removing the ground*/
  demStruct *dem;    /*DEM structure*/
  float demTol;      /*tolerance to separate ground and non-ground points*/
}voxStruct;


/*###################################################*/
/*lidar radiometric parameters for voxels*/

typedef struct{
  float minRefl;    /*minimum refletance value to scale between 0 and 1*/
  float maxRefl;    /*maximum refletance value to scale between 0 and 1*/
  float appRefl;    /*scale between TLS reflectance and size*/
  float beamTanDiv; /*tan of tls beam divergence*/
  float beamRad;    /*TLS start radius*/
  float minGap;     /*minimum gap fraction correction to apply*/
  char correctR;    /*correct reflectance for range to get gap fractions*/
}lidVoxPar;


/*#####################################*/
/*structure to hold a range image*/

typedef struct{
  int nBins;        /*number of range bins*/
  int nX;           /*number of x bins*/
  int nY;           /*number of y bins*/
  float rRes;       /*range resolution*/
  float iRes;       /*image resolution*/
  double bounds[6]; /*minX minY minZ, maxX maxY maxZ*/
  double x0;        /*central coordinate*/
  double y0;        /*central coordinate*/
  double z0;        /*central coordinate*/
  char **image;     /*gap image at each range*/
  float grad[3];    /*vector along rage image*/
}rImageStruct;


/*####################################*/
/*voxel gap structure*/

typedef struct{
  /*the voxel*/
  voxStruct *vox;   /*TLS voxels. Gap within voxel*/
  voxStruct *toTLS; /*gap to TLS voxels*/
  float **meanGap;  /*mean minimum gap for voxels*/
  float **meanRefl; /*mean reflectance for voxels*/
  int **contN;      /*number of beams contributing to the above means*/

  /*TLS map*/
  int **mapFile;        /*file per voxel*/
  uint32_t **mapPoint;  /*point per voxel*/
  int *nIn;             /*number of points per voxel*/
}tlsVoxMap;

/*####################################*/
/*TLS beams, polar coords*/

typedef struct{
  float zen;     /*zenith*/
  float az;      /*azimuth*/
  float x;        /*beam origin*/
  float y;        /*beam origin*/
  float z;        /*beam origin*/
  uint32_t shotN; /*shot number within this scan*/
  uint8_t nHits;  /*number of hits of this beam*/
  float *r;       /*range*/
  float *refl;    /*reflectance*/
}tlsBeam;


/*##########################################*/
/*TLS point cloud*/

typedef struct{
  int bin;       /*bin number*/
  float x;       /*coordinate*/
  float y;       /*coordinate*/
  float z;       /*coordinate*/
  float gap;     /*voxel gap fraction*/
  float r;       /*range*/
  uint16_t refl; /*reflectance*/
  uint32_t hitN;  /*hit number*/
  uint8_t nHits;  /*number of hits of this beam*/
}tlsPoint;


/*##########################################*/
/*TLS scan*/

typedef struct{
  tlsBeam *beam;     /*array of beams*/
  tlsPoint *point;   /*array of points*/
  double xOff;       /*offset to allow coords to be floats*/
  double yOff;       /*offset to allow coords to be floats*/
  double zOff;       /*offset to allow coords to be floats*/
  uint32_t nBeams;   /*number of beams in this scan*/
  uint32_t nPoints;  /*number of points in this scan*/
  FILE *ipoo;        /*file pointer*/
  hid_t hdfFile;     /*HDF5 file pointer*/
  uint32_t pOffset;  /*current point position for buffering*/
  uint32_t nRead;    /*number of beams to read at once*/
  uint32_t maxRead;  /*maximum number of beams we could have in a region*/
  uint64_t totSize;  /*total file size*/
  uint64_t totRead;  /*amount of file read so far*/
  float **matrix;    /*matrix needed for ptx files*/
}tlsScan;


/*##########################################*/
/*function definitions*/

tlsScan *tidyTLScan(tlsScan *);
tlsScan *tidyTLScans(tlsScan *,int);
tlsScan *readTLSwithinVox(char **,int,voxStruct *,char,tlsVoxMap *);
tlsScan *readOneTLS(char *,voxStruct *,char,tlsVoxMap *,int,lidVoxPar *);
void readTLSpolarBinary(char *,uint32_t,tlsScan **);
void readPTXscan(char *,uint32_t,tlsScan **);
rImageStruct *allocateRangeImage(float,float,float,float *,double *,double *);
voxStruct *voxAllocate(int,float *,double *,char);
voxStruct *tidyVox(voxStruct *);
int *findVoxels(double *,double,double,double,double *,double *,int *,int,int,int,double **);
int *findAllVoxels(double *,double,double,double,voxStruct **,int,int **,double **,int *);
int *beamVoxels(float *,double,double,double,double *,double *,int,int,int,int *,double,double **,float);
float tlsPointSize(double,uint16_t,float,float,float,float,float,float);
double *findVoxelBounds(int *,int,voxStruct *,tlsVoxMap *,float *,double,double,double);
void tidyVoxelMap(tlsVoxMap *,int);
void silhouetteImage(int,pCloudStruct **,tlsScan *,rImageStruct *,lidVoxPar *,int *,int,tlsVoxMap *);
void waveFromImage(rImageStruct *,float **,char,float);
void fillInRimageGround(rImageStruct *);
void setWaveformRange(float *,double,float *,int,float);
void rotateX(double *,double);
void rotateZ(double *,double);
void readBoundsFromTLS(double *,char **,int);
void beamVoxelBounds(double *,float *,float,char,double *);
void makeBinImage(double *,float *,float *,char *,int,float,double,uint16_t,float,int,int,float,lidVoxPar *);
void readCanBounds(canBstruct *,char *,double *);
void setCanInd(int *,int *,double,double,double,lasFile *,canBstruct *);
void imposeCanBound(float *,int,int,int);
double setCanGround(double,double,canBstruct *);
void setTopVoxBlank(voxStruct *);
void voxelate(voxStruct *,float *,lasFile *,double,double,double,float);
void writeAsciiVox(voxStruct *,char *);

/*the end*/
/*###################################################*/

