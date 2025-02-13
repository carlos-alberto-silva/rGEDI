import h5py
import argparse
import string
import numpy as np
from osgeo import gdal
from pyproj import Proj, transform
from scipy import pi
from math import sqrt


####################################################################
# pointing uncertainty

def pointingUncertainty(errSigma,dX,dY):
  shift=np.random.normal(loc=0.0,scale=errSigma)
  uncX=-1*dY*shift  # note that we rotate 90 degrees to be normal to track
  uncY=dX*shift     # note that we rotate 90 degrees to be normal to track
  return(uncX,uncY)


####################################################################
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


####################################################################
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


####################################################################
# read command line

def getCmdArgs():
    '''
    Get commdnaline arguments
    GEDI grid resolution
    latitude band
    '''
    p = argparse.ArgumentParser(description=("Print out GEDI footprint locations from orbital traks\n"))
    p.add_argument("--modis", dest ="modNamen", type=str, default='/gpfs/data1/vclgp/data/phenology/leaf_on_off/forestonly/MODIS_land_ancillary_sets6.tif', help=("MODIS ancillary filename."))
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
    p.add_argument("--useWeak",dest = "useWeak", type=int,default=0, help=("Use weak beam by day"))
    p.add_argument("--pointErr",type=float,dest="errSigma",help="pointing uncetainty in metres, 1 sigma",default=0)
    cmdargs = p.parse_args()
    return cmdargs


####################################################################
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
  np.random.seed(options.seed)
  miss_len=options.miss_len
  outNamen=options.outNamen
  skip_track_list=options.skip_track_list
  usePhen=options.usePhen
  useWeak=options.useWeak
  errSigma=options.errSigma

  # mission 0 time in days
  t0=11*30.4375+1

  # tanslate bounds from ALS to orbit projection
  epsg="epsg:"+str(iEPSG)
  inProj=Proj(init=epsg)
  epsg="epsg:"+str(oEPSG)
  outProj=Proj(init=epsg)
  bX,bY=transform(outProj,inProj,bounds[0:2],bounds[2:4])

  # resolution in metres
  sRes=60.0

  # read MODIS
  data_lc,data_vcf,data_leafon,data_leafoff,xOrigin,yOrigin,pixelWidth,pixelHeight,cols,rows=readMODIS(modNamen)

  # open output
  fOut=open(outNamen, 'w')

  # position shifts
  if(errSigma>0.0):
    errX=np.empty(shape=(0),dtype=float)
    errY=np.empty(shape=(0),dtype=float)
    dayList=np.empty(shape=(0),dtype=int)
  else:
    errX=[0.0]
    errY=[0.0]

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

    # blank arrays
    uX=np.empty(shape=(0))
    uY=np.empty(shape=(0))
    waveID=[]
    numb=0
 
    # loop over usable tracks
    for i in range(1,len(useInd)):
      ind=useInd[i]

      # only look at sections of lines. Reset counters as a line passes
      if(useInd[i]!=(useInd[i-1]+1)):
        numb=0
        continue

      # if over mission length, skip
      if((delta_time[ind]/(24*60*60))>miss_len):
        break
 
      # Which MODIS pixels are intersected
      crossing_lon,crossing_lat=intersectLines(x[ind],x[ind-1],y[ind],y[ind-1],xOrigin,yOrigin,pixelWidth,pixelHeight)

      # loop over intersected MODIS pixels
      for j in range(0,len(crossing_lon)-1):
 
        # cloud cover
        if(np.random.random()<cFrac):
          #print("Cloudy")
          continue

        # determine modis properties
        xInd=int((crossing_lon[j]-xOrigin)/pixelWidth)
        yInd=int((crossing_lat[j]-yOrigin)/pixelHeight)
        if((xInd<0)|(xInd>=cols)|(yInd<0)|(yInd>=rows)):
          print("No modis")
          continue

        vcf=data_vcf[yInd][xInd]
        lc=data_lc[yInd][xInd]
        onDat=data_leafon[yInd][xInd]
        offDat=data_leafoff[yInd][xInd]
        power =1 if(beamid in pow_beam_list) else 0
 

        # check is leaf-on
        c_doy=(delta_time[ind]/(24*60*60))+t0
        while( c_doy > 365.25 ):
          c_doy-=365.25
        if(lc==1)|(lc==2)|(lc>5)|(lc==0):  # evergreen or unreliable phenology
          condition_doy=True
        elif (onDat==0)&(offDat==0):  # no data
          condition_doy=True
        elif onDat < offDat: #mostly nothern hemi
          condition_doy=(c_doy>onDat)&(c_doy<offDat)
        else: #growing period extending from year1 to the next year
          condition_doy=(c_doy>offDat)&(c_doy<onDat)
 
        # is the beam usable?
        if(((usePhen==0)|(condition_doy==True))&((useWeak==1)|((power==1)|(sim_sun_el[ind]>=83)|(vcf<70)))):
          if(condition_doy==True):
            phenLab=1
          else:
            phenLab=0

          # footprint labels
          if(sim_sun_el[ind]>=83):
            tim="night"
          else:
            tim="day"

          # date
          ye=int(((delta_time[ind]/(24*60*60))+t0)/365.25)+18
          doy=int(((delta_time[ind]/(24*60*60))+t0)%365.25)
          day=int(delta_time[ind]/(24*60*60))
          minute=(delta_time[ind]/60.0)-(day*24.0*60.0)
 
          # how many footprints within this MODIS cell?
          xS,yS=transform(inProj,outProj,crossing_lon[j],crossing_lat[j])
          xE,yE=transform(inProj,outProj,crossing_lon[j+1],crossing_lat[j+1])
          nIn=int(sqrt((xE-xS)**2+(yE-yS)**2)/sRes)

          # direction of track for footprints
          dX=xE-xS
          dY=yE-yS
          if((np.fabs(dX)>2000)|(np.fabs(dY)>2000)): # this step is longer than one MODIS pixel
            print("A step too far",dX,dY,crossing_lon[j],crossing_lat[j],crossing_lon[j+1],crossing_lat[j+1])
          totLen=sqrt(dX**2+dY**2)
          dX=dX/totLen
          dY=dY/totLen

          # minimum 60 m spacing
          if(numb>0):
            diff=sqrt((xS-uX[len(uX)-1])**2+(yS-uY[len(uY)-1])**2)
          else:
            diff=sRes*4.0
          if(diff<sRes):
            x00=uX[len(uX)-1]+dX*sRes
            y00=uY[len(uY)-1]+dY*sRes
          else:
            x00,y00=transform(inProj,outProj,crossing_lon[j],crossing_lat[j])

          # position uncertainty
          if(errSigma>0.0):
            dayStart=int(delta_time[ind]/(24*60*60))
            if dayStart in dayList:  # exists
              doyInd=np.where(dayStart==dayList)[0][0]
            else: # new day, new uncertainty
              dayList=np.append(dayList,dayStart)
              uncX,uncY=pointingUncertainty(errSigma,dX,dY)
              errX=np.append(errX,uncX)
              errY=np.append(errY,uncY)
              doyInd=len(dayList)-1
          else:  # no uncertainty
            doyInd=0

          # loop over footprints within cell
          for p in range(0,nIn+1):
            uX=np.append(uX,x00+dX*p*sRes+errX[doyInd])
            uY=np.append(uY,y00+dY*p*sRes+errY[doyInd])
            waveID.append("%d.%d.%d.%s.zen.%d.y.%d.doy.%d.pow.%d.phen.%d.vcf.%d.day.%d.min.%d"%(beamid,j,numb,tim,int(sim_sun_el[ind]),ye,doy,power,phenLab,vcf,day,minute))
            numb=numb+1  # increment footprint counter along
 
    # write data
    useCoord=np.where((uX>=bounds[0])&(uX<=bounds[1])&(uY>=bounds[2])&(uY<=bounds[3]))[0]
    for p in useCoord:
      line=str(uX[p])+" "+str(uY[p])+" "+waveID[p]+"\n"
      fOut.write(line)
    # wipe arrays
    uX=np.empty(shape=(0))
    uY=np.empty(shape=(0))
    waveID=[]
    numb=0

  fOut.close()

  print("Written to ",outNamen)

