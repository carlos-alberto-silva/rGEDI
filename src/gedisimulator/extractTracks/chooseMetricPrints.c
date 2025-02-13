#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"


/*###########################################*/
/*control structure*/

typedef struct{
  char metricNamen[1000];   /*metric file input*/
  char trackNamen[1000];    /*GEDI tracks to use*/
  char outNamen[1000];      /*output filename*/
  double minSep;           /*minimum separation to accept*/
  float gridRes;           /*grid resolution*/
}control;


/*###########################################*/
/*track structure*/

typedef struct{
  double *x;
  double *y;
  double minX;   /*x corner*/
  double minY;   /*y corner*/
  double maxX;   /*x corner*/
  double maxY;   /*y corner*/
  int nTracks;
}trackStruct;

/*###########################################*/
/*grid structure for mapping tacks*/

typedef struct{
  int nX;        /*number of x grid cells*/
  int nY;        /*number of y grid cells*/
  double minX;   /*x corner*/
  double minY;   /*y corner*/
  double maxX;   /*x corner*/
  double maxY;   /*y corner*/
  double res;    /*grid resolution*/
  int *nIn;      /*number of footprints per cell*/
  int **map;     /*list of track indexes per grid cell*/
}gridStruct;


/*###########################################*/
/*main*/

int main(int argc,char **argv)
{
  control *dimage=NULL;
  control *readCommands(int,char **);
  trackStruct *tracks=NULL;
  trackStruct *readTracks(char *);
  gridStruct *grid=NULL;
  gridStruct *setMapGrid(control *,trackStruct *);
  void selectMetrics(trackStruct *,char *,char *,double,gridStruct *);


  /*read commands*/
  dimage=readCommands(argc,argv);

  /*read track locations*/
  tracks=readTracks(dimage->trackNamen);

  /*set up mapping grid*/
  grid=setMapGrid(dimage,tracks);

  /*read metris and output needed*/
  selectMetrics(tracks,dimage->metricNamen,dimage->outNamen,dimage->minSep,grid);

  /*tidy up*/
  if(grid){
    TIDY(grid->nIn);
    TTIDY((void **)grid->map,grid->nX*grid->nY);
    TIDY(grid);
  }
  if(tracks){
    TIDY(tracks->x);
    TIDY(tracks->y);
    TIDY(tracks);
  }
  TIDY(dimage);
  return(0);
}/*main*/


/*###########################################*/
/*set up map grid structure*/

gridStruct *setMapGrid(control *dimage,trackStruct *tracks)
{
  int i=0,xBin=0,yBin=0,place=0;
  gridStruct *grid=NULL;

  /*allocate space*/
  if(!(grid=(gridStruct *)calloc(1,sizeof(gridStruct)))){
    fprintf(stderr,"error grid map allocation.\n");
    exit(1);
  }

  /*set values*/
  grid->minX=tracks->minX;
  grid->maxX=tracks->maxX;   
  grid->minY=tracks->minY;   
  grid->maxY=tracks->maxY;   
  grid->res=(double)dimage->gridRes;
  grid->nX=(int)((grid->maxX-grid->minX)/grid->res+1.0);
  grid->nY=(int)((grid->maxY-grid->minY)/grid->res+1.0);

  grid->nIn=ialloc(grid->nX*grid->nY,"nIn grid",0);
  grid->map=iIalloc(grid->nX*grid->nY,"grid map",0);
  for(i=grid->nX*grid->nY-1;i>=0;i--){
    grid->nIn[i]=0;
    grid->map[i]=NULL;
  }

  /*map footprints*/
  for(i=0;i<tracks->nTracks;i++){
    xBin=(int)((tracks->x[i]-grid->minX)/grid->res+0.5);
    yBin=(int)((tracks->y[i]-grid->minY)/grid->res+0.5);

    if((xBin>=0)&&(xBin<grid->nX)&&(yBin>=0)&&(yBin<grid->nY)){
      place=yBin*grid->nX+xBin;
      grid->map[place]=markInt(grid->nIn[place],grid->map[place],i);
      grid->nIn[place]++;
    }
  }/*grid mapping loop*/

  return(grid);
}/*setMapGrid*/


/*###########################################*/
/*select and write metrics*/

void selectMetrics(trackStruct *tracks,char *metricNamen,char *outNamen,double minSep,gridStruct *grid)
{
  int j=0,xCol=0,yCol=0,nUse=0;
  int xCent=0,yCent=0,place=0;     /*to choose grid locations*/
  int xBin=0,yBin=0,k=0;           /*to choose grid locations*/
  int xS=0,xE=0,yS=0,yE=0;         /*to choose grid locations*/
  int64_t i=0,*useList=NULL;
  double sepSq=0,*minSepSq=NULL;
  double x=0,y=0,dx=0,dy=0;
  double thresh=0;
  char line[20000],*token=NULL;
  char temp[20000];
  char lastTok[100];
  char writtenHead=0;
  char **useLines=NULL;
  FILE *ipoo=NULL,*opoo=NULL;

  /*open metrics*/
  if((ipoo=fopen(metricNamen,"r"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",metricNamen);
    exit(1);
  }
  thresh=minSep*minSep;
  writtenHead=0;

  /*open output*/
  if((opoo=fopen(outNamen,"w"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",outNamen);
    exit(1);
  }

  /*choose metrics*/
  if(!(useList=(int64_t *)calloc(tracks->nTracks,sizeof(int64_t)))){
    fprintf(stderr,"error useList allocation.\n");
    exit(1);
  } 

  /*allocate distance array and set to blank*/
  minSepSq=dalloc(tracks->nTracks,"useList",0);
  useLines=chChalloc(tracks->nTracks,"useLines",0);
  for(i=0;i<tracks->nTracks;i++){
    minSepSq[i]=10000000.0;
    useList[i]=-1;
    useLines[i]=NULL;
  }

  /*search for closest footprints to tracks*/
  i=0;
  while(fgets(line,10000,ipoo)!=NULL){
    strcpy(temp,line);
    if(strncasecmp(line,"#",1)){  /*read data*/
      x=y=1000000000000.0;  /*silly values*/
      /*read coords*/
      j=1;
      token=strtok(line," ");
      while(token){
        if(j==yCol)y=atof(token);
        else if(j==xCol)x=atof(token);
        else if(j>yCol)break;
        token=strtok(NULL," ");
        j++;
      }

      /*see which grid cells it overlaps*/
      xCent=(int)((x-grid->minX)/grid->res+0.5);
      yCent=(int)((y-grid->minY)/grid->res+0.5);
      if((x-(double)xCent*grid->res+grid->minX)<minSep)xS=xCent-1;
      else                                             xS=xCent;
      if(((double)(xCent+1)*grid->res+grid->minX-x)<minSep)xE=xCent+2;
      else                                                 xE=xCent+1;
      if((y-(double)yCent*grid->res+grid->minY)<minSep)yS=yCent-1;
      else                                             yS=yCent;
      if(((double)(yCent+1)*grid->res+grid->minY-y)<minSep)yE=yCent+2;
      else                                                 yE=yCent+1;
      if(xS<0)xS=0;
      if(yS<0)yS=0;
      if(xE>grid->nX)xE=grid->nX;
      if(yE>grid->nY)yE=grid->nY;

      /*loop through intersected grid cells*/
      for(xBin=xS;xBin<xE;xBin++){
        for(yBin=yS;yBin<yE;yBin++){
          place=yBin*grid->nX+xBin;
          for(k=0;k<grid->nIn[place];k++){
            /*look at map*/
            j=grid->map[place][k];
            /*determine separation*/
            dx=x-tracks->x[j];
            dy=y-tracks->y[j];
            sepSq=dx*dx+dy*dy;
            /*record closes suitable*/
            if((sepSq<minSepSq[j])&&(sepSq<=thresh)){
              if(useList[j]<=0)nUse++;
              useList[j]=i;
              minSepSq[j]=sepSq;
              /*copy line*/
              TIDY(useLines[j]);
              useLines[j]=challoc(strlen(temp)+1,"useLines",j+1);
              strcpy(useLines[j],temp);
            }
          }/*point within grid loop*/
        }/*intersected track grid loop*/
      }/*intersected track grid loop*/

      i++;
    }else{  /*read header*/
      j=1;
      token=strtok(line," ");
      while(token){
        if(!strncasecmp(token,"lon,",4))xCol=atoi(lastTok);
        else if(!strncasecmp(token,"lat,",4))yCol=atoi(lastTok);
        strcpy(lastTok,token);
        token=strtok(NULL," ");
        j++;
      }
      fprintf(stdout,"Metric coord cols %d %d\n",xCol,yCol);

      /*write header to output*/
      if(writtenHead==0){
        fprintf(opoo,"%s",temp);
        writtenHead=1;
      }
    }
  }
  TIDY(minSepSq);
  fprintf(stdout,"Read %lld metrics and selected %d\n",(long long int)i,nUse);
  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }

  /*write out output*/
  for(j=0;j<tracks->nTracks;j++){
    if(useList[j]>=0)fprintf(opoo,"%s",useLines[j]);
  }

  /*tidy up*/
  TTIDY((void **)useLines,tracks->nTracks);
  TIDY(useList);
  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Written to %s\n",outNamen);
  return;
}/*selectMetrics*/


/*###########################################*/
/*read track locations*/

trackStruct *readTracks(char *namen)
{
  int i=0;
  trackStruct *tracks=NULL;
  char line[400],temp1[200],temp2[200];
  FILE *ipoo=NULL;

  /*allocate*/
  if(!(tracks=(trackStruct *)calloc(1,sizeof(trackStruct)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  tracks->x=NULL;
  tracks->y=NULL;
  tracks->nTracks=0;
  tracks->minX=tracks->minY=1000000000000.0;
  tracks->maxX=tracks->maxY=-1000000000000.0;

  /*open data*/
  if((ipoo=fopen(namen,"r"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",namen);
    exit(1);
  }

  /*count number of lines*/
  while(fgets(line,400,ipoo)!=NULL)if(strncasecmp(line,"#",1))tracks->nTracks++;
  tracks->x=dalloc(tracks->nTracks,"x tracks",0);
  tracks->y=dalloc(tracks->nTracks,"y tracks",0);

  /*rewind to start of file*/
  if(fseek(ipoo,(long)0,SEEK_SET)){
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  /*read data*/
  i=0;
  while(fgets(line,400,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      if(sscanf(line,"%s %s",temp1,temp2)==2){
        tracks->x[i]=atof(temp1);
        tracks->y[i]=atof(temp2);
        /*determine area bounds*/
        if(tracks->x[i]<tracks->minX)tracks->minX=tracks->x[i];
        if(tracks->x[i]>tracks->maxX)tracks->maxX=tracks->x[i];
        if(tracks->y[i]<tracks->minY)tracks->minY=tracks->y[i];
        if(tracks->y[i]>tracks->maxY)tracks->maxY=tracks->y[i];
        i++;
      }else{
        fprintf(stderr,"Error reading %s\n",namen);
        exit(1);
      }
    }
  }/*data reading loop*/

  /*close file*/
  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  fprintf(stdout,"There will be %d footprints\n",tracks->nTracks);
  return(tracks);
}/*readTracks*/


/*###########################################*/
/*read commands*/

control *readCommands(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;

  /*allocate*/
  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }

  /*defaults*/
  dimage->gridRes=1000.0;

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-metric",7)){
        checkArguments(1,i,argc,"-metric");
        strcpy(dimage->metricNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-tracks",7)){
        checkArguments(1,i,argc,"-tracks");
        strcpy(dimage->trackNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-minSep",7)){
        checkArguments(1,i,argc,"-minSep");
        dimage->minSep=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(dimage->outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-gridRes",8)){
        checkArguments(1,i,argc,"-gridRes");
        dimage->gridRes=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry chooseMetricPrints -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }

  return(dimage);
}/*readCommands*/

/*the end*/
/*###########################################*/
