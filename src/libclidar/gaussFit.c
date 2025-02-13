#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "mpfit.h"
#include "libLasProcess.h"




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


/*###############################################*/
/*structure for mpfit*/

typedef struct{  /* EXAMPLE: fitting y(x) */
  float *x;         /* x - independent variable of model */
  float *y;         /* y - measured "y" values */
  int nGauss;
}dataStruct;


/*###############################################*/
/*structure for multiple returns*/

typedef struct{
  int nFeat;     /*number of feautures*/
  float **temp;  /*feature waveform*/
  int *sBin;     /*start bins, to get range*/
  int *width;    /*number of bins for each feature*/
}multRet;


/*###############################################*/
/*structure for turning points*/

typedef struct{
  int nFeats;
  int *bound;
}turnStruct;


/*###############################################*/
/*global functions*/

float gauss(float,float,float);
int *markInt(int,int *,int);
int gaussBuffBins;


/*##################################################*/
/*fit multiple Gaussians*/

float *fitMultiGauss(float *x,float *decon,int nBins,float gSmooth,int *totGauss,float minGsig)
{
  int i=0,j=0,k=0,ret=0;
  int *nGauss=NULL;     /*number of Gaussians per feature*/
  float **params=NULL,*fitted=NULL;
  float *padX=NULL,gRes=0;
  float minErr=0;
  float *fitGauss(float *,float *,int,float,turnStruct *,float);
  turnStruct *turnings=NULL;
  turnStruct *findTurning(float *,int,float,float *);
  multRet *returns=NULL;
  multRet *filterData(float *,int,float);

  minErr=0.01;

  /*split up separate returns*/
  returns=filterData(decon,nBins,0.0001);

  /*set results blank*/
  fitted=falloc((uint64_t)nBins,"",0);
  for(i=0;i<nBins;i++)fitted[i]=0.0;
  params=fFalloc(returns->nFeat,"Gaussian parameters",0);
  nGauss=ialloc(returns->nFeat,"number of Gaussians",0);

  for(ret=0;ret<returns->nFeat;ret++){
    if(returns->sBin[ret]<0){  /*in case all bounds are of interest*/
      returns->sBin[ret]=0;
      returns->temp[ret]=&(returns->temp[ret][0]);
      returns->width[ret]=nBins;
    }

    /*determine number of Gaussians in this return*/
    turnings=findTurning(&(returns->temp[ret][0]),returns->width[ret],gSmooth,x);

    /*make a padded x arrays*/
    padX=falloc((uint64_t)returns->width[ret],"padded x",0);
    gRes=fabs(x[2]-x[1]);
    for(i=0;i<returns->width[ret];i++)padX[i]=x[returns->sBin[ret]]+(float)(i-gaussBuffBins)*gRes;

    /*do the fitting*/
    params[ret]=fitGauss(&(padX[0]),returns->temp[ret],returns->width[ret],minErr,turnings,minGsig);
    TIDY(padX);

    /*load up the final wave*/
    for(i=0;i<nBins;i++){
      for(j=0;j<turnings->nFeats-1;j++){
        if(params[ret][3*j+1]>0.0){
          fitted[i]+=params[ret][3*j+1]*gauss(x[i],params[ret][3*j+2],params[ret][3*j+0]);
        }
      }
    }/*load up final wave bin loop*/
    nGauss[ret]=turnings->nFeats-1;
    (*totGauss)+=nGauss[ret];
    if(turnings){
      TIDY(turnings->bound);
      TIDY(turnings);
    }
  }/*return loop*/


  /*copy Gaussian parameters to end of fitted array*/
  if((*totGauss)>0){
    if(!(fitted=(float *)realloc(fitted,(nBins+3*(*totGauss))*sizeof(float)))){
      fprintf(stderr,"Error in fitted Gaussian realoocation %lu\n",(nBins+3*(*totGauss))*sizeof(float));
      exit(1);
    }
    k=0;
    for(i=0;i<returns->nFeat;i++){
      for(j=0;j<3*nGauss[i];j++){
        fitted[k+nBins]=params[i][j];
        k++;
      }
    }
  }else TIDY(fitted);
  TIDY(nGauss);
  TTIDY((void **)params,returns->nFeat);

  /*tidy arrays*/
  if(returns){
    TIDY(returns->sBin);
    TIDY(returns->width);
    TTIDY((void **)returns->temp,returns->nFeat);
    TIDY(returns);
  }
  return(fitted);
}/*gaussFit*/


/*##################################################*/
/*fit multiple Gaussians*/

float *fitSingleGauss(float *x,float *decon,int nBins,float gSmooth,int *totGauss,float **gaussPar)
{
  int i=0;
  int nParams=0;
  int fitCheck=0;
  int mpfit(mp_func,int,int,double *, mp_par *,mp_config *,void *,mp_result *);
  int gaussErr(int,int,double *,double *,double **,void *);
  float *fitted=NULL;
  double *params=NULL;
  double *initialSingleGauss(float *,float *,int,float);
  mp_par *parStruct=NULL;
  mp_result *result=NULL;
  mp_config *config=NULL;
  dataStruct data;     /*data structure*/

  data.x=x;
  data.y=decon;
  data.nGauss=1;

  nParams=3;

  /*initial estimates*/
  params=initialSingleGauss(x,decon,nBins,0.0);

  /*allocate arrays*/
  if(!(config=(mp_config *)calloc(1,sizeof(mp_config)))){
    fprintf(stderr,"error in control structure.\n");
    exit(1);
  }
  config->nofinitecheck=1;
  config->maxiter=1000;
  if(!(result=(mp_result *)calloc(1,sizeof(mp_result)))){
    fprintf(stderr,"error in mpfit structure.\n");
    exit(1);
  }
  result->resid=dalloc(nBins,"",0);
  result->xerror=dalloc(nParams,"",0);
  result->covar=dalloc(nParams*nParams,"",0);

  /*parameter bounds*/
  if(!(parStruct=(mp_par *)calloc(nParams,sizeof(mp_par)))){
    fprintf(stderr,"error in bound structure.\n");
    exit(1);
  }
  /*mu*/
  parStruct[0].fixed=0;
  parStruct[0].side=3;
  parStruct[0].limited[0]=1;
  parStruct[0].limits[0]=params[0]-4.0*params[2];
  parStruct[0].limited[1]=1;
  parStruct[0].limits[1]=params[0]+4.0*params[2];
  /*A*/
  parStruct[1].fixed=0;
  parStruct[1].side=3;
  parStruct[1].limited[0]=1;
  parStruct[1].limits[0]=params[1]/10.0;
  parStruct[1].limited[1]=1;
  parStruct[1].limits[1]=params[1]*10.0;
  /*sig*/
  parStruct[2].fixed=0;
  parStruct[2].side=3;
  parStruct[2].limited[0]=1;
  parStruct[2].limits[0]=0.00001;  /*a small number*/
  parStruct[2].limited[1]=1;
  parStruct[2].limits[1]=params[2]*4.0;


  /*the fitting*/
  fitCheck=mpfit(gaussErr,nBins,nParams,params,parStruct,config,(void *)(&data),result);
  /*if(fitCheck<0){
    fprintf(stderr,"fitCheck %d\n",fitCheck);
  }*/

  /*tidy up*/
  if(result){
    TIDY(result->resid);
    TIDY(result->xerror);
    TIDY(result->covar);
    TIDY(result);
  }
  TIDY(config);
  TIDY(parStruct);

  (*totGauss)=data.nGauss;
  data.x=NULL;
  data.y=NULL;

  /*generate fitted wave*/
  fitted=falloc((uint64_t)nBins,"Gaussian fit",0);
  for(i=0;i<nBins;i++)fitted[i]=params[1]*gauss(x[i],params[2],params[0]);

  /*copy parameters*/
  (*gaussPar)=falloc((uint64_t)nParams,"Gaussian parameters",0);
  for(i=0;i<nParams;i++)(*gaussPar)[i]=params[i];
  TIDY(params);
  return(fitted);
}/*fitSingleGauss*/


/*###############################################*/
/*initial guess for single Gaussian*/

double *initialSingleGauss(float *x,float *y,int nBins,float gWidth)
{
  int i=0,centBin=0;
  float start=0,end=0;
  float thresh=0;
  float *smoothed=NULL;
  double *params=NULL;

  params=dalloc(3,"Gaussian parameters",0);
  params[1]=-1000.0;   /*A*/

  /*smooth to taste*/
  if(gWidth>0.0)smoothed=smooth(gWidth,nBins,y,0.15);
  else          smoothed=&(y[0]);

  /*find maximum*/
  for(i=0;i<nBins;i++){
    if(smoothed[i]>params[1]){
      params[0]=x[i];  /*mu*/
      params[1]=smoothed[i];  /*A*/
      centBin=i;
    }
  }

  /*trace fore and aft to get width*/
  thresh=params[1]*exp(-0.5);
  for(i=centBin;i>=0;i--){
    if(smoothed[i]<=thresh){
      start=x[i];
      break;
    }
  }
  for(i=centBin+1;i<nBins;i++){
    if(smoothed[i]<=thresh){
      end=x[i];
      break;
    }
  }
  params[2]=fabs(end-start)/2.0;  /*sigma*/

  if(gWidth>0.0){
    TIDY(smoothed);
  }else{
    smoothed=NULL;
  }

  return(params);
}/*initialSingleGauss*/


/*###############################################*/
/*fit Gaussians, from testDiffWave.c*/

float *fitGauss(float *x,float *y,int numb,float minErr,turnStruct *turnings,float minGsig)
{
  int i=0;
  int nGauss=0,nParams=0;
  float *params=NULL;
  float *initialGuess(float *,float *,int,int,int *,float);
  dataStruct *data=NULL;
  mp_par *parStruct=NULL;
  mp_par *setGaussBounds(float *,float *,int,double *,int,int,float);
  mp_result *result=NULL;
  mp_config *config=NULL;
  int mpfit(mp_func,int,int,double *, mp_par *,mp_config *,void *,mp_result *);
  int fitCheck=0;
  int gaussErr(int,int,double *,double *,double **,void *);
  void checkGaussBounds(int ,mp_par *,double *);
  void cleanWeakGauss(double *,int *);
  double *doParams=NULL;


  nGauss=turnings->nFeats-1;
  nParams=3*nGauss;
  params=initialGuess(x,y,numb,nGauss,turnings->bound,minGsig);
  doParams=dalloc(nParams,"",0);  /*convert to doubles and array base 0 for mpfit*/
  for(i=0;i<nParams;i++)doParams[i]=(double)params[i];
  minErr=0.0000001;

  /*load observations into a structure for passing about*/
  if(!(data=(dataStruct *)calloc(1,sizeof(dataStruct)))){
    fprintf(stderr,"error in control structure.\n");
    exit(1);
  }
  data->x=x;
  data->y=y;
  data->nGauss=nGauss;
  /*set bounds*/
  parStruct=setGaussBounds(x,y,numb,doParams,nParams,nGauss,minGsig);
  if(numb<2){  /*one point above noise, fix max and width*/
    params[3]=0.17;
    fprintf(stderr,"Insufficient parameters for a Gaussian fit\n");
    return(params);
  }

  /*check that the guesses are within the bounds*/
  checkGaussBounds(nGauss,parStruct,doParams);

  if(!(config=(mp_config *)calloc(1,sizeof(mp_config)))){
    fprintf(stderr,"error in control structure.\n");
    exit(1);
  }
  config->nofinitecheck=1;
  config->maxiter=2000;
  if(!(result=(mp_result *)calloc(1,sizeof(mp_result)))){
    fprintf(stderr,"error in mpfit structure.\n");
    exit(1);
  }
  result->resid=dalloc(numb,"",0);
  result->xerror=dalloc(nParams,"",0);
  result->covar=dalloc(nParams*nParams,"",0);
  fitCheck=mpfit(gaussErr,numb,nParams,doParams,parStruct,config,data,result);
  if(fitCheck<0)fprintf(stderr,"Fit check %d numb %d nGauss %d\n",fitCheck,numb,nGauss);


  if(fitCheck>=0){
    cleanWeakGauss(doParams,&nGauss); /*delete Gaussians with very low amplitude*/
    for(i=0;i<nParams;i++)params[i]=(float)doParams[i]; /*2 is a success flag for mpfit*/
  }else           for(i=0;i<nParams;i++)params[i]=-1.0;

  TIDY(parStruct);
  if(result){
    TIDY(result->resid);
    TIDY(result->xerror);
    TIDY(result->covar);
    TIDY(result);
  }
  TIDY(config);
  TIDY(doParams);
  TIDY(data);

  return(params);
}/*fitGauss*/


/*###############################################*/
/*clean weak Gaussians*/

void cleanWeakGauss(double *doParams,int *nGauss)
{
  int i=0;
  float maxA=0,thresh=0;

  /*find max Gaussian amplitude*/
  for(i=0;i<(*nGauss);i++)if(doParams[3*i+1]>maxA)maxA=doParams[3*i+1];

  thresh=0.0005*maxA;   /*0.05% of max amplitude*/
  for(i=0;i<(*nGauss);i++){
    if(doParams[3*i+1]<thresh){
      doParams[3*i]=doParams[3*i+1]=doParams[3*i+2]=0.0;
    }
  }

  return;
}/*cleanWeakGauss*/


/*###############################################*/
/*check that initial guesses are within bounds*/

void checkGaussBounds(int nGauss,mp_par *parStruct,double *doParams)
{
  int i=0;

  /*loop over Gaussians*/
  for(i=0;i<nGauss;i++){
    /*mu*/
    if(doParams[3*i]<=parStruct[3*i].limits[0])doParams[3*i]=parStruct[3*i].limits[0]+0.05;
    if(doParams[3*i]>=parStruct[3*i].limits[1])doParams[3*i]=parStruct[3*i].limits[1]-0.05;
    /*A*/
    if(doParams[3*i+1]<=parStruct[3*i+1].limits[0])doParams[3*i+1]=parStruct[3*i].limits[0]+0.0005;
    if(doParams[3*i+1]>=parStruct[3*i+1].limits[1])doParams[3*i+1]=parStruct[3*i].limits[1]-0.0005;
    /*sig*/
    if(doParams[3*i+2]<=parStruct[3*i+2].limits[0])doParams[3*i+2]=parStruct[3*i+2].limits[0]+0.05;
    if(doParams[3*i+2]>=parStruct[3*i+2].limits[1])doParams[3*i+2]=parStruct[3*i+2].limits[1]-0.05;
  }/*Gaussian loop*/

  return;
}/*checkBounds*/


/*###############################################*/
/*extract relevant data and determine width*/

multRet *filterData(float *y,int numb,float offset)
{
  int i=0,j=0,place=0;
  float mean=0; /*,stdev=0;*/
  float waveThresh=0;
  float tol=0,max=0;
  multRet *returns=NULL;
  char inFeat=0,markStart=0;
  char found=0;

  max=-1000.0;
  for(i=0;i<numb;i++)if(y[i]>max)max=y[i];
  tol=max/800.0; //0.000001; /*smallest floating point number*/

  if(!(returns=(multRet *)calloc(1,sizeof(multRet)))){
    fprintf(stderr,"error in multiple return structure.\n");
    exit(1);
  }
  returns->nFeat=0;
  returns->temp=NULL;
  returns->sBin=NULL;
  returns->width=NULL;

  /*get background signal*/
  /*mean=0.0;
  for(i=numb/2;i<numb;i++)mean+=y[i];
  mean/=(float)(numb-numb/2);
  stdev=0.0;
  for(i=numb/2;i<numb;i++)stdev+=(y[i]-mean)*(y[i]-mean);
  stdev=sqrt(stdev/(float)(numb-numb/2));
  if(stdev>tol)waveThresh=mean+4.0*stdev;
  else         waveThresh=mean+tol;
  if(waveThresh<tol)waveThresh=tol;*/

  waveThresh=tol;
  mean=0.0;

  inFeat=0;
  for(i=0;i<numb;i++){
    if(y[i]>waveThresh){  /*in a feature*/
      if(inFeat==0)markStart=1;
      else         markStart=0;
      inFeat=1;
      if(markStart){
        found=0;
        for(j=i;j>=0;j--){  /*track to feature start*/
          if(y[j]<=mean){
            returns->sBin=markInt(returns->nFeat,returns->sBin,j);
            found=1;
            break;
          }
        }/*track to start*/
        if(found==0)returns->sBin=markInt(returns->nFeat,returns->sBin,1);
      }
    }else if((y[i]<=waveThresh)&&(inFeat==1)){
      inFeat=0;
      found=0;
      for(;i<numb;i++){ /*track to end*/
        if(y[i]<=mean){
          returns->width=markInt(returns->nFeat,returns->width,i-returns->sBin[returns->nFeat]+1);
          returns->nFeat++;
          found=1;
          break;
        }
      }/*track to end*/
      if(found==0){
        returns->width=markInt(returns->nFeat,returns->width,(numb/2)-returns->sBin[returns->nFeat]+1);
        returns->nFeat++;
      }
    }/*left feature*/
  }/*feature counting loop*/
  if(inFeat){  /*for low tails*/
    returns->width=markInt(returns->nFeat,returns->width,i-returns->sBin[returns->nFeat]+1);
    returns->nFeat++;
  }
  /*copy intensity data*/
  gaussBuffBins=130;
  returns->temp=fFalloc(returns->nFeat,"temp feature",0);
  for(i=0;i<returns->nFeat;i++){
    returns->temp[i]=falloc((uint64_t)returns->width[i]+2*(uint64_t)gaussBuffBins,"temp feature",i+1);
    for(j=0;j<gaussBuffBins;j++)returns->temp[i][j]=0.0;
    for(j=0;j<returns->width[i];j++){
      place=j+returns->sBin[i];
      if((place>=0)&&(place<numb)){
        returns->temp[i][j+gaussBuffBins]=y[place];
      }else fprintf(stderr,"Bawbag %d of %d\n",place,numb);
    }
    for(j=gaussBuffBins+returns->width[i];j<2*gaussBuffBins+returns->width[i];j++)returns->temp[i][j]=0.0;
    returns->width[i]+=2*gaussBuffBins;
  }/*intensity data copying*/
  return(returns);
}/*filterData*/


/*###############################################*/
/*Gaussian error function*/

int gaussErr(int numb, int npar, double *p, double *deviates,double **derivs, void *private)
{
  int i=0,j=0;
  float y=0,arg=0,expo=0;
  float A=0,mu=0,sig=0;
  dataStruct *data=NULL;

  data=(dataStruct *)private;

  for(j=0;j<data->nGauss;j++){
    mu=(float)p[3*j+0];
    A=(float)p[3*j+1];
    sig=(float)p[3*j+2];
  }

  for(i=0;i<numb;i++){
    y=0.0;
    for(j=0;j<data->nGauss;j++){
      mu=(float)p[3*j+0];
      A=(float)p[3*j+1];
      sig=(float)p[3*j+2];
      y+=A*gauss(data->x[i],sig,mu);

      /*calculate gradients*/
      if(derivs){
        arg=-1.0*(data->x[i]-mu)*(data->x[i]-mu)/(2.0*sig*sig);
        if(arg>=-730.0)expo=exp(arg);
        else           expo=0.0;
        derivs[3*j+0][i]=(double)(A*(data->x[i]-mu)/(sig*sig)*expo);
        derivs[3*j+1][i]=(double)expo;
        derivs[3*j+2][i]=(double)(A*(data->x[i]-mu)*(data->x[i]-mu)/(sig*sig*sig)*expo);
      }
    }/*Gaussian loop*/
    deviates[i]=(double)(y-data->y[i]);
  }/*time loop*/


  return(0);
}/*gaussErr*/


/*###############################################*/
/*initial Gaussian parameter guess*/

float *initialGuess(float *x,float *y,int numb,int nGauss,int *bound,float minGsig)
{
  int i=0,j=0,maxBin=0;
  float *params=NULL;
  float contN=0,thisThresh=0;

  /*variablaes are 3*j+1 mu, 3*j+2 A, 3*j+3 sig*/
  params=falloc(3*(uint64_t)nGauss,"parameters",0);

  for(j=0;j<nGauss;j++){
    if((bound[j+1]-bound[j])>=3){
      contN=0.0;
      params[3*j+0]=0.0;
      params[3*j+1]=-10.0;
      for(i=bound[j];i<=bound[j+1];i++){
        params[3*j+0]+=x[i]*y[i];  /*position is weighted by intensity*/
        contN+=y[i];
        if(y[i]>params[3*j+1]){
          params[3*j+1]=y[i];
          maxBin=i;           /*to calculate the width from */
        }
      }
      if(contN>0.0)params[3*j+0]/=contN;       /*normalise mean range*/
      else         params[3*j+0]=(x[bound[j+1]]-x[bound[j]])/2.0;

      thisThresh=params[3*j+1]*exp(-0.5);  /*1/e^2 of maximum*/
      params[3*j+2]=-1.0;              /*nonsense value*/
      for(i=maxBin;i<=numb;i++){
        if(y[i]<=thisThresh){
          params[3*j+2]=(x[i]-params[3*j+0])/2.0;
          break;
        }
      }
    }else{
      params[3*j+0]=(x[bound[j+1]]-x[bound[j]])/2.0;
      params[3*j+1]=(y[bound[j+1]]-y[bound[j]])/2.0;
      params[3*j+2]=minGsig+0.02;
    }
    /*to prevent daft values*/
    if(params[3*j+2]<=minGsig)params[3*j+2]=minGsig+0.02;
    if(params[3*j+1]<=0.0)params[3*j+1]=0.001;

  }/*sub feature loop*/
  return(params);
}/*initialGuess*/



/*###############################################*/
/*find turning points for other methods*/

turnStruct *findTurning(float *y,int width,float preSmooth,float *x)
{
  int i=0;
  float *smoothed=NULL;
  float *d2x=NULL,*temp=NULL;
  turnStruct *turnings=NULL;

  /*smooth to taste*/
  if(preSmooth>0.0)smoothed=smooth(preSmooth,width,y,0.15);
  else             smoothed=&(y[0]);


  /*determine second derivative*/
  d2x=falloc((uint64_t)width,"d2x",0);
  for(i=1;i<width-1;i++)d2x[i]=2.0*smoothed[i]-(smoothed[i+1]+smoothed[i-1]);

  temp=smooth(preSmooth,width,d2x,0.15);
  TIDY(d2x);
  d2x=temp;
  temp=NULL;

  if(!(turnings=(turnStruct *)calloc(1,sizeof(turnStruct)))){
    fprintf(stderr,"error in multiple return structure.\n");
    exit(1);
  }
  turnings->nFeats=0;
  turnings->bound=NULL;

  /*mark the first bin*/
  i=0;
  turnings->bound=markInt(turnings->nFeats,turnings->bound,i);
  turnings->nFeats++;
  for(i=2;i<width-2;i++){
    if(smoothed[i]>0.0){
      if((d2x[i]<=d2x[i-1])&&(d2x[i]<d2x[i+1])){  /*minimum of the second derivative*/
        turnings->bound=markInt(turnings->nFeats,turnings->bound,i-1);
        turnings->nFeats++;
      }
    }
  }
  /*and then mark the last bin*/
  turnings->bound=markInt(turnings->nFeats,turnings->bound,width-1);
  turnings->nFeats++;
  if(preSmooth>0.0){
    TIDY(smoothed);
  }else{
    smoothed=NULL;
  }
  TIDY(d2x);
  return(turnings);
}/*findTurning*/


/*###############################################*/
/*set bounds for Gaussian parameters*/

mp_par *setGaussBounds(float *x,float *y,int numb,double *doParams,int nParams,int nGauss,float minGsig)
{
  int i=0;
  mp_par *parStruct=NULL;

  if(!(parStruct=(mp_par *)calloc(nParams,sizeof(mp_par)))){
    fprintf(stderr,"error in bound structure.\n");
    exit(1);
  }

  for(i=0;i<nGauss;i++){
    /*mu*/
    parStruct[3*i].fixed=0;
    parStruct[3*i].side=3;
    parStruct[3*i].limited[0]=1;
    parStruct[3*i].limits[0]=x[0];
    parStruct[3*i].limited[1]=1;
    parStruct[3*i].limits[1]=x[numb-1];
    /*A*/
    parStruct[3*i+1].fixed=0;
    parStruct[3*i+1].side=3;
    parStruct[3*i+1].limited[0]=1;
    parStruct[3*i+1].limits[0]=0.0; //doParams[3*i+1]/20.0;
    parStruct[3*i+1].limited[1]=1;
    parStruct[3*i+1].limits[1]=doParams[3*i+1]*6.6;
    /*sig*/
    parStruct[3*i+2].fixed=0;
    parStruct[3*i+2].side=3;
    parStruct[3*i+2].limited[0]=1;
    parStruct[3*i+2].limits[0]=minGsig;  /*system pulse*/
    parStruct[3*i+2].limited[1]=0;
    parStruct[3*i+2].limits[1]=minGsig+50.0; //(x[numb-1]-x[0]);
  }/*Gaussian loop*/

  return(parStruct);
}/*setGaussBounds*/


/*###############################################*/
/*non normalised Gaussian*/

float gauss(float x,float sigma,float offset)
{
  float y=0,arg=0;

  arg=-1.0*(x-offset)*(x-offset)/(2.0*sigma*sigma);
  if(arg>=-700.0)y=exp(-1.0*(x-offset)*(x-offset)/(2.0*sigma*sigma));
  else           y=0.0;

  return(y);
}/*gauss*/

/*the end*/
/*####################################################################################*/

