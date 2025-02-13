import h5py
import argparse
import string
import numpy as np
from osgeo import gdal
import random
from pyproj import Proj, transform
from scipy import pi
from math import sqrt


# read MODIS
def readMODIS(modNamen):
  driver = gdal.GetDriverByName('GTiff')
  dataset_MOD = gdal.Open(modNamen)
  b_landcover = dataset_MOD.GetRasterBand(1)
  b_vcf = dataset_MOD.GetRasterBand(2)
  b_startgrowing = dataset_MOD.GetRasterBand(3)
  b_endgrowing = dataset_MOD.GetRasterBand(4)
  cols = dataset_MOD.RasterXSize
  rows = dataset_MOD.RasterYSize
  transMod = dataset_MOD.GetGeoTransform()
  xOrigin = transMod[0]
  yOrigin = transMod[3]
  pixelWidth = transMod[1]
  pixelHeight = transMod[5]
  #prepare raster data set
  data_lc = b_landcover.ReadAsArray(0, 0, cols, rows)
  data_vcf = b_vcf.ReadAsArray(0, 0, cols, rows)
  data_leafon = b_startgrowing.ReadAsArray(0, 0, cols, rows)
  data_leafoff = b_endgrowing.ReadAsArray(0, 0, cols, rows)
  return(data_lc,data_vcf,data_leafon,data_leafoff,xOrigin,yOrigin,pixelWidth,pixelHeight,cols,rows)


# intersect MODIS pixels
def intersectLines(x1,x0,y1,y0,xOrigin,yOrigin,pixelWidth,pixelHeight):
  # set arrays with oprigin as first element
  xi=[]
  yi=[]
  xi=np.append(xi,x0)
  yi=np.append(yi,y0)
  tx=x0
  ty=y0
  # vector of track
  dx=x1-x0
  dy=y1-y0
  vectLen=sqrt(dx**2+dy**2)
  dx=dx/vectLen
  dy=dy/vectLen
  xStep=1 if(dx>0) else -1
  yStep=1 if(dy>0) else -1
  # distance to end of track
  distToX=(x1-tx)*xStep
  distToY=(y1-ty)*yStep
  # loop along track
  while((distToX>0)&(distToY>0)):
    # coordinates of next intersection along x or y
    xT=(int((tx-xOrigin)/pixelWidth+0.00001)+xStep)*pixelWidth+xOrigin
    yT=(int((ty-yOrigin)/pixelHeight+0.00001)-yStep)*pixelHeight+yOrigin  # y pixels measured from top
    # distance along vector of intersection
    if(dx!=0.0):
      xDist=(xT-tx)/dx
    else:
      xDist=100000.0
    if(dy!=0.0):
      yDist=(yT-ty)/dy
    else:
      yDist=100000.0;
    # which intersection is closer?
    if(xDist<=yDist):
      tx=xT
      ty=ty+dy*xDist
    else:
      tx=tx+dx*yDist
      ty=yT
    # update distance to end
    distToX=(x1-tx)*xStep
    distToY=(y1-ty)*yStep
    if((distToX<0)|(distToY<0)):
      break
    # record closest intersection
    xi=np.append(xi,tx)
    yi=np.append(yi,ty)
  return(xi,yi)


# read command line
def getCmdArgs():
    '''
    Get commdnaline arguments
    GEDI grid resolution
    latitude band
    '''
    p = argparse.ArgumentParser(description=("Print out GEDI footprint locations from orbital traks\n"))
    p.add_argument("--modis", dest ="modNamen", type=str, default='/gpfs/data1/vclgp/htang/GEDI/ISS/MODIS_ancillary_sets.tif', help=("MODIS ancillary filename."))
    p.add_argument("--output",type=str,dest="outNamen",help="output filename",default='teast.coords')
    p.add_argument("--cloud",type=float,dest="cFrac",help="cloud fraction",default=0.0)
    p.add_argument("--seed",type=int,dest="seed",help="random number seed",default=0)
    p.add_argument("--bounds",type=float,nargs='*',dest="bounds",help="Bounds of interest in input EPSG",default=[-180,180,-1,1])
    p.add_argument("--iEPSG",type=int,dest="iEPSG",help="Input EPSG",default=4326)
    p.add_argument("--oEPSG",type=int,dest="oEPSG",help="Output EPSG",default=4326)
    p.add_argument("--mission_length",dest = "miss_len", type=int, default=395, help=("Total mission length in days since 01 Nov 2018 00:00:00 UTC.\nDefault for the full nominal 2-yr mission and a valid input should be less than 365*2."))
    p.add_argument("--orbit_sim_dir",dest = "orbit_sim_dir", type=str, default='/gpfs/data1/vclgp/htang/GEDI/ISS/BaselineSIM', help=("Directory containing orbital simulations from cott Luthcke.\nDefault /gpfs/data1/vclgp/htang/GEDI/ISS/BaselineSIM."))#added on 2nd Apr 2018
    p.add_argument("-s","--skip_track",dest = "skip_track_list", type=int, nargs='*', default=[-1], help=("Track numbering to be dropped in order to assess its impact on mission level 1 requirement"))#added on Feb 23-2018
    p.add_argument("-p","--power_beams",dest = "pow_beam_list", type=int, nargs='*', default=[5,6], help=("Track numbers of power beams"))
    p.add_argument("--usePhen",dest = "usePhen", type=int,default=0, help=("Use phenology sqitch"))
    cmdargs = p.parse_args()

    return cmdargs


## Main ##
if __name__ == '__main__':
  # read command line
  options=getCmdArgs()
  modNamen=options.modNamen
  orbit_sim_dir=options.orbit_sim_dir
  bounds=options.bounds
  iEPSG=options.iEPSG
  oEPSG=options.oEPSG
  pow_beam_list=options.pow_beam_list
  cFrac=options.cFrac
  random.seed(options.seed)
  miss_len=options.miss_len
  outNamen=options.outNamen
  skip_track_list=options.skip_track_list
  usePhen=options.usePhen

  # mission 0 time in days
  t0=11*30.4375+1

  # tanslate bounds from ALS to orbit projection
  epsg="epsg:"+str(iEPSG)
  inProj=Proj(init=epsg)
  epsg="epsg:"+str(oEPSG)
  outProj=Proj(init=epsg)
  bX,bY=transform(outProj,inProj,bounds[0:2],bounds[2:4])

  # resolution in metres
  #sRes=transform(outProj,inProj,0,60)[1]
  sRes=60.0 #(60/6371007.181)*180/pi

  # read MODIS
  data_lc,data_vcf,data_leafon,data_leafoff,xOrigin,yOrigin,pixelWidth,pixelHeight,cols,rows=readMODIS(modNamen)

  # open output
  fOut=open(outNamen, 'w')

  # loop over GEDI tracks
  for beamid in range(1,11):
    # do we need this beam?
    if beamid in skip_track_list:
      print("Skipping track",beamid)
      continue

    # read orbital data
    inNamen="%s/GEDI_2yr_BaselineSIM_beam%d_v01.h5" % (orbit_sim_dir,beamid)
    f=h5py.File(inNamen,'r')
    y=np.array( f['/GT/latitude'])
    x=np.array(f['/GT/longitude'])
    sim_sun_el=np.array(f['/GT/sun_el'])
    delta_time=f['/delta_time'].value
    f.close()
 
    # which coords to use
    buff_res=(0.066667 * np.sin(52 * np.pi /180.0))*4
    useInd=np.sort(np.where((y>=(bY[0]-buff_res))&(y<(bY[1]+buff_res))&(x>=(bX[0]-buff_res))&(x<(bX[1]+buff_res)))[0])

    # loop over usable tracks
    for i in range(1,len(useInd)):
      ind=useInd[i]

      if(useInd[i-1]!=(useInd[i]-1)):
        continue

      crossing_lon,crossing_lat=intersectLines(x[ind],x[ind-1],y[ind],y[ind-1],xOrigin,yOrigin,pixelWidth,pixelHeight)
      tX,tY=transform(inProj,outProj,crossing_lon,crossing_lat)
      print("nIn",len(crossing_lon))

      for j in range(0,len(crossing_lon)):
        line=str(tX[j])+" "+str(tY[j])+"\n"
        fOut.write(line)


  fOut.close()

  print("Written to ",outNamen)

