#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "geotiffio.h"
#include "xtiffio.h"
#include "tools.h"
#include "tiffRead.h"


/*############################################*/
/*read a geotiff*/

void readGeotiff(geot *geotiff,char *namen,char readData)
{
  int i=0,j=0,x=0,y=0,place=0;
  uint32_t tileWidth=0,tileLength=0;
  unsigned int type=0;
  short tiepointsize=0,scalesize=0;
  double *tiepoints=NULL,*scale=NULL;
  float *buff=NULL;
  //GTIF *tiffStruct=(GTIF*)0;
  TIFF *tiffIn=(TIFF*)0;

  if((tiffIn=XTIFFOpen(namen,"r"))==NULL){
    fprintf(stderr,"Damn, no %s\n",namen);
    exit(1);
  }
  geotiff->tiepoints=dalloc(6,"",0);
  geotiff->scale=dalloc(3,"",0);
  geotiff->fImage=NULL;
  geotiff->dImage=NULL;
  geotiff->image=NULL;

  /*for now we have assumed OSNG, 27700, projection*/
  //tiffStruct=GTIFNew(tiffIn);
  TIFFGetField(tiffIn,TIFFTAG_IMAGEWIDTH,&geotiff->nX);
  TIFFGetField(tiffIn,TIFFTAG_IMAGELENGTH,&geotiff->nY);
  TIFFGetField(tiffIn,TIFFTAG_DATATYPE,&type);


  /*if the data is tiled*/
  TIFFGetField(tiffIn,TIFFTAG_TILEWIDTH,&tileWidth);
  if(tileWidth>0){
    TIFFGetField(tiffIn,TIFFTAG_TILELENGTH,&tileLength);
    buff=falloc((uint64_t)tileWidth*(uint64_t)tileLength,"buffer",0);
  }


  //GTIFPrint(tiffStruct,0,0);

  TIFFGetField(tiffIn,TIFFTAG_GEOPIXELSCALE,&scalesize,&scale);
  TIFFGetField(tiffIn,TIFFTAG_GEOTIEPOINTS,&tiepointsize,&tiepoints);
  for(i=0;i<3;i++) geotiff->scale[i]=scale[i];
  for(i=0;i<6;i++) geotiff->tiepoints[i]=tiepoints[i];

  if(readData){
    if(type==0){  /*unsigned char*/
      geotiff->image=uchalloc((uint64_t)geotiff->nX*(uint64_t)geotiff->nY,namen,0);
      for(i=0;i<geotiff->nY;i++){                  /*looping along the lattitude*/
        if(TIFFReadScanline(tiffIn,&(geotiff->image[i*geotiff->nX]),i,1)!=1){
          fprintf(stderr,"Error reading scan line %d from tiff image\n",i);
          exit(1);
        }
      }
    }else if(type==3){ /*float*/
      if(((int)TIFFScanlineSize(tiffIn)/geotiff->nX)==4){
        geotiff->fImage=falloc((uint64_t)geotiff->nX*(uint64_t)geotiff->nY,namen,0);

        if(tileWidth==0){   /*read scan lines*/
          for(i=0;i<geotiff->nY;i++){                  /*looping along the lattitude*/
            if(TIFFReadScanline(tiffIn,&(geotiff->fImage[i*geotiff->nX]),i,1)!=1){
              fprintf(stderr,"Error reading scan line %d from tiff image\n",i);
              exit(1);
            }
          }
        }else{ /*read tiled data*/
          for(y=0;y<geotiff->nY;y+=tileLength){                  /*looping along the lattitude*/
            for(x=0;x<geotiff->nX;x+=tileWidth){
              TIFFReadTile(tiffIn,buff,x,y,0,0);

              /*pack into results*/
              for(i=0;i<tileWidth;i++){
                if((i+x)>=geotiff->nX)continue;
                for(j=0;j<tileLength;j++){
                  if((j+y)>=geotiff->nY)continue;
                  place=x+i+(j+y)*geotiff->nX;
                  geotiff->fImage[place]=buff[j*tileWidth+i];
                }
              }
            }
          }
        }
      }else if(((int)TIFFScanlineSize(tiffIn)/geotiff->nX)==8){
        geotiff->dImage=dalloc(geotiff->nX*geotiff->nY,namen,0);
        for(i=0;i<geotiff->nY;i++){                  /*looping along the lattitude*/
          if(TIFFReadScanline(tiffIn,&(geotiff->dImage[i*geotiff->nX]),i,1)!=1){
            fprintf(stderr,"Error reading scan line %d from tiff image\n",i);
            exit(1);
          }
        }
      }else{
        fprintf(stderr,"What do you think you're doing!?!\n");
        exit(1);
      }
    }else{
      fprintf(stderr,"Cannot handle type %d\n",type);
      exit(1);
    }
  }/*read data question*/

  TIDY(buff);
  XTIFFClose(tiffIn);
  return;
}/*readGeotiff*/


/*############################################*/
/*tidy tiff file structure*/

geot *tidyTiff(geot *tiff)
{
  if(tiff){
    TIDY(tiff->image);
    TIDY(tiff->fImage);
    TIDY(tiff->dImage);
    TIDY(tiff->tiepoints);
    TIDY(tiff->scale);
    TIDY(tiff);
  }

  return(tiff);
}/*tidyTiff*/
/*the end*/
/*############################################*/

