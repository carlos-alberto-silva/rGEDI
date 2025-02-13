

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


/*###########################################################*/
/*hold photon coiunting data for binary output*/

typedef struct{
  uint64_t nPhots;   /*number of photons*/
  double *x;         /*x coord*/
  double *y;         /*y coord*/
  double *z;         /*z coord*/
}photonHDF;


/*###########################################################*/
/*poton counting structure*/

typedef struct{
  /*photon rates*/
  float designval;     /*mean number of photons per footprint*/
  float *prob;         /*probability of each number of photons occuring*/
  int pBins;           /*number of probability bins*/
  float rhoVrhoG;      /*ratio of canopy to ground reflectance for weighting*/
  float nPhotC;        /*mean number of canopy photons per footprint*/
  float nPhotG;        /*mean number of ground photons per footprint*/
  char reflDiff;       /*are we accounting for a difference in reflectance?*/
  /*noise*/
  float noise_mult;    /*noise scaling factor*/
  float H;             /*search window length, metres*/
  /*IO*/
  FILE *opoo;          /*output file*/
  char outNamen[1000];  /*output filename*/
  char writeHDF;       /*writeHDF switch*/
  /*HDF file*/
  photonHDF *hdf;      /*HDF file, if needed*/
}photonStruct;


/*###########################################################*/
/*global functions*/

void setPhotonRates(photonStruct *);
float *uncompressPhotons(float *,dataStruct *,photonStruct *,noisePar *,gediIOstruct *);
float **countPhotons(float *,dataStruct *,photonStruct *,int *,denPar *,noisePar *,char);
float *countWaveform(float *,dataStruct *,photonStruct *,denPar *,noisePar *);
float pickArrayElement(float,float *,int,char);
void removeAsymmetryPCL(float *,int);
void setPhotonProb(photonStruct *);

/*###########################################################*/

