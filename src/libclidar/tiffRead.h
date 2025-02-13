


/*############################################*/
/*a geotiff*/

typedef struct{
  unsigned char *image; /*image*/
  float *fImage;        /*float image*/
  double *dImage;       /*double image*/
  int nX;          /*image size*/
  int nY;          /*image size*/
  double *tiepoints;  /*geolocation*/
  double *scale;      /*scaling values*/
}geot;  /*a geotiff*/

/*############################################*/

/*global functions*/
geot *tidyTiff(geot *);
void readGeotiff(geot *,char *,char);


/*############################################*/

