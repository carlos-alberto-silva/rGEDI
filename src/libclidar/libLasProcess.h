
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


/*#############################################*/
/*strucure of denoising variables*/

typedef struct{
  /*denoising*/
  float thresh;      /*noise threshold*/
  float meanN;       /*mean noise level*/
  float tailThresh;  /*trailing edge threshold, meanN by default*/
  float sWidth;      /*smoothing length*/
  float msWidth;     /*middle smoothing width (after stats, before denoising*/
  float psWidth;     /*pre-smoothing width*/
  int minWidth;      /*min width above noise to accept*/
  int medLen;        /*median filter width*/
  char varNoise;     /*variable noise threshold switch*/
  char medStats;     /*use median rather than mean statistics (ARSF)*/
  float statsLen;    /*distance to calculate noise statistics over*/
  char noiseTrack;   /*noise tracking switch*/
  float threshScale; /*scale variable noise threshold*/
  char bitRate;      /*bit rate if using variable noise with floats*/
  float maxDN;       /*maxDN for digitisation*/
  char preMatchF;    /*matched filter before denoising*/
  char posMatchF;    /*matched filter after denoising*/
  char corrDrift;    /*correct detector drift switch*/
  char varDrift;     /*calculate drift parameter for every waveform*/
  float fixedDrift;  /*fixed drift correction factor if using*/

  /*deconvolution*/
  float pScale;      /*scale pulse length by*/
  int maxIter;       /*maximum number of iterations*/
  double deChang;    /*change between decon iterations to stop*/
  char deconMeth;    /*deconvolution method, -1 none, 0 Gold, 1 Richardson-Lucy*/
  char deconGauss;   /*is pulse a Gaussian or not*/
  char pNamen[400];  /*pulse filename*/
  float pSigma;      /*pulse width if a Gaussian is used*/
  float **pulse;     /*pulse to deconvolve by*/
  int pBins;         /*number of pulse bins*/
  float *matchPulse; /*matched filter pulse*/
  float res;         /*waveform resolution*/
  int maxPbin;       /*maximum pulse bin*/

  /*Gaussian fitting*/
  char fitGauss;     /*switch*/
  float gWidth;      /*smoothing width before picking features*/
  int nGauss;        /*number of Gaussians*/
  float *gPar;       /*Gaussian parameters*/
  char gaussPulse;   /*convert pulse to a Gaussian*/
  float minGsig;     /*minimum Gaussian width to fit*/

  /*Gaussians to identify hard targets*/
  char matchHard;    /*match check hard targets*/
  float hardThresh;  /*minimum match*/
  float *hardPulse; /*pulse to compare for hard targets*/
  /*OLD*/
  char gaussFilt;    /*switch*/
  float hardWidth;   /*maxWidth of hard feature*/
  float hardTol;     /*tolerance to scale width by*/
}denPar;


/*##################################################*/
/*global structure to save reallocating smoothing pulse*/

typedef struct{
  int nPulses;
  int *nBins;
  float *res;
  float *sWidth;
  float **pulse;
}smoothPulse;


/*###################################################*/
/*ground data structure*/

typedef struct{
  uint32_t nPoints;
  double *xUse;
  double *yUse;
  double *zUse;
  double *par;
  int nPoly[2];
}groundDstruct;


/*#####################################*/
/*canopy bounds*/

typedef struct{
  float *canMin;     /*minimum canopy height*/
  float *canMax;     /*maximum canopy height*/
  int cNx;           /*number of x elements*/
  int cNy;           /*number of y elements*/
  float cRes;        /*canopy bounds resolution*/
  double cUbound[4]; /*x, y bounds*/
}canBstruct;


/*#############################################*/
/*common function definitions*/

float *smooth(float,int,float *,float);

/*#############################################*/
/*functions*/

float *processWave(unsigned char *,int,denPar *,float);
float *processFloWave(float *,int,denPar *,float);
char checkWaveform(float *,uint32_t);
float *findRH(float *,double *,int,double,float,int *);
float foliageHeightDiversity(float *,int);
float foliageHeightDiversityHist(float *,int,float);
float *waveLmoments(float *,int,float,int);
float niMetric(float *,double *,int,float,double,float);
float determineGaussSep(float,float);
float *matchedFilter(float *,int,denPar *,float);
float *canProfile(float *,float *,int);
float *subtractGroundFromCan(float *,float *,int);
float *subtractGaussFromCan(float *,int,float,float,float,double *);

/*the end*/
/*#############################################*/

