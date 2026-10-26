/*
 * Rebuild tests/testthat/fixtures/hancock-reference.csv.
 *
 * This small validation program is a direct extraction of the supported
 * Gaussian/count/pulse-after path in Steven Hancock's gedisimulator gediIO.c,
 * commit 49ad4f26b03083b61ee0bd424e8a7e4a8f2c301d.  It intentionally uses C
 * float arrays where the upstream source does.  Compile with:
 *   gcc -O2 tools/hancock_reference_fixture.c -o hancock_fixture.exe -lm
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct { double x, y, z; float refl; int cls, ret, nret; } Point;

static const Point pts[] = {
  {-3.0,-2.0, 0.8,10,2,1,1}, {-1.2, 0.4, 1.1,12,2,1,1},
  { 0.2,-0.3, 1.0,14,2,1,1}, { 1.1, 1.0, 1.3,11,2,1,1},
  { 2.7,-1.4, 0.9, 9,2,1,1}, {-2.1, 1.7, 8.0,18,5,1,2},
  {-2.1, 1.7,15.0,22,5,2,2}, {-0.4,-1.1,12.5,20,5,1,1},
  { 0.8, 0.5,18.0,25,5,1,2}, { 0.8, 0.5,22.0,27,5,2,2},
  { 2.2, 1.9,16.0,19,5,1,1}, { 3.4, 0.2,20.0,24,5,1,1}
};

static void simulate(int method, int normalize, const char *scenario) {
  const int np = (int)(sizeof(pts) / sizeof(pts[0]));
  const float fsigma=5.5f, res=0.15f, pres=0.15f, thresh=0.0006f;
  float psigma=(15.0f*0.2998f/2.0f)/2.35482f;
  int pn=0, j, i, centre, nb, gxN, gyN;
  float y, x=0.0f, total=0.0f;
  double minz=1e100,maxz=-1e100,minx=1e100,miny=1e100,maxx=-1e100,maxy=-1e100;
  do { y=(float)exp(-0.5*(x/psigma)*(x/psigma)); x+=pres; pn+=2; } while(y>=thresh);
  float *py=(float *)calloc((size_t)pn,sizeof(float));
  centre=pn/2; x=-(float)centre*pres;
  for(i=0;i<pn;i++){ py[i]=(float)exp(-0.5*(x/psigma)*(x/psigma)); total+=py[i]; x+=pres; }
  for(i=0;i<pn;i++)py[i]/=total;
  for(i=0;i<np;i++){
    if(pts[i].z<minz)minz=pts[i].z; if(pts[i].z>maxz)maxz=pts[i].z;
    if(pts[i].x<minx)minx=pts[i].x; if(pts[i].x>maxx)maxx=pts[i].x;
    if(pts[i].y<miny)miny=pts[i].y; if(pts[i].y>maxy)maxy=pts[i].y;
  }
  double buff=35.0+(double)pn*(double)pres/2.0;
  minz-=buff; maxz+=buff; nb=(int)((maxz-minz)/(double)res);
  float *wave=(float *)calloc((size_t)nb,sizeof(float));
  float *ground=(float *)calloc((size_t)nb,sizeof(float));
  float *sm=(float *)calloc((size_t)nb,sizeof(float));
  float *smg=(float *)calloc((size_t)nb,sizeof(float));
  const float grid=1.5f; double ox=minx+grid, oy=miny+grid;
  gxN=(int)((maxx-minx)/grid)+2; gyN=(int)((maxy-miny)/grid)+2;
  int *density=(int *)calloc((size_t)(gxN*gyN),sizeof(int));
  if(normalize){
    for(i=0;i<np;i++) if(pts[i].ret==pts[i].nret){
      int gx=(int)((pts[i].x-ox)/grid), gy=(int)((pts[i].y-oy)/grid);
      if(gx>=0&&gx<gxN&&gy>=0&&gy<gyN)density[gy*gxN+gx]++;
    }
  }
  for(i=0;i<np;i++){
    double sep=sqrt(pts[i].x*pts[i].x+pts[i].y*pts[i].y);
    float scale=(float)exp(-0.5*(sep/fsigma)*(sep/fsigma));
    if(scale>thresh){
      if(normalize){
        int gx=(int)((pts[i].x-ox)/grid), gy=(int)((pts[i].y-oy)/grid);
        if(gx>=0&&gx<gxN&&gy>=0&&gy<gyN&&density[gy*gxN+gx]>0)
          scale/=(float)density[gy*gxN+gx];
      }
      float value=scale;
      if(method==1)value=(float)(pts[i].refl*scale);
      else if(method==2)value=(float)(scale*(1.0f/(float)pts[i].nret));
      int bin=(int)floor((maxz-pts[i].z)/(double)res);
      if(bin>=0&&bin<nb){ wave[bin]+=value; if(pts[i].cls==2)ground[bin]+=value; }
    }
  }
  for(i=0;i<nb;i++){
    float n=0.0f;
    for(j=0;j<pn;j++){
      int bin=i+(int)floor((float)(centre-j)*pres/res);
      if(bin>=0&&bin<nb){ sm[i]+=wave[bin]*py[j]; smg[i]+=ground[bin]*py[j]; n+=1.0f; }
    }
    if(n>0.0f){sm[i]/=n;smg[i]/=n;}
  }
  total=0.0f; for(i=0;i<nb;i++)total+=sm[i]*res;
  for(i=0;i<nb;i++){sm[i]/=total;smg[i]/=total; printf("%s,%d,%.10f,%.10g,%.10g\n",scenario,i+1,maxz-(double)i*res,sm[i],smg[i]);}
  free(py);free(wave);free(ground);free(sm);free(smg);free(density);
}

int main(void){
  puts("scenario,bin,elevation,wave,ground");
  simulate(0,0,"count_raw"); simulate(0,1,"count_normalized");
  simulate(1,0,"intensity_raw"); simulate(1,1,"intensity_normalized");
  simulate(2,0,"fraction_raw"); simulate(2,1,"fraction_normalized");
  return 0;
}
