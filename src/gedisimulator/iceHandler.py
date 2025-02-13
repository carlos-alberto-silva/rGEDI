
'''
Handles ICESat-2 simulations
'''


#################
# Packages

import numpy as np
import h5py
from math import sqrt
from pyproj import Proj, transform
if __name__ == '__main__':
  import argparse



########################################

class ice2(object):
  '''
  Handles real ICESat-2 data
  '''

  #################################

  def __init__(self,namen,epsg=4326,minX=-100000000,maxX=100000000,minY=-1000000000,maxY=100000000):
    '''Class initialiser'''
    self.readPhotons(namen,epsg=epsg,minX=minX,maxX=maxX,minY=minY,maxY=maxY)

  #################################

  def readPhotons(self,namen,epsg=4326,minX=-100000000,maxX=100000000,minY=-1000000000,maxY=100000000):
    '''Read ICESat-2 HDF5 file'''
    f=h5py.File(namen,'r')
    lon=np.array(f['gt1l']['heights']['lon_ph'])
    lat=np.array(f['gt1l']['heights']['lat_ph'])
    z=np.array(f['gt1l']['heights']['h_ph'])
    # truth if this is a simulation
    if(('gt1l/veg_truth/refDEM' in f)&('gt1l/heights/signal_conf_photon' in f)&('gt1l/veg_truth/isground' in f)):
      self.sim=1
      s=np.array(f['gt1l']['heights']['signal_conf_photon'])
      g=np.array(f['gt1l']['veg_truth']['refDEM'])
      isGr=np.array(f['gt1l']['veg_truth']['isground'])
    else:
      self.sim=0
    # reproject
    if(epsg!=4326):
      inProj=Proj(init="epsg:4326")
      outProj=Proj(init="epsg:"+str(epsg))
      x,y=transform(inProj, outProj, lon, lat)
    else:
      x=lon
      y=lat
    # filter if needed
    useInds=np.where((x>=minX)&(x<=maxX)&(y>=minY)&(y<=maxY))
    if(len(useInds)>0):
      useInds=useInds[0]
      self.x=x[useInds]
      self.y=y[useInds]
      self.z=z[useInds]
      if(self.sim==1):
        self.s=s[useInds]
        self.g=g[useInds]
        self.isGr=isGr[useInds]


  #################################

  def writeCoords(self,outNamen="test.pts"):
    '''Write out coordinates and photons'''
    f=open(outNamen,'w')
    for i in range(0,len(self.x)):
      line=str(self.x[i])+" "+str(self.y[i])+" "+str(self.z[i])+"\n"
      f.write(line)
    f.close()
    print("Written to",outNamen)


########################################

class iceSim(object):
  '''
  Reads and acts upon
  an ICEsat-2 simulation
  '''

  #################################

  def __init__(self,namen,epsg):
    '''Class initialiser'''
    temp=np.loadtxt(namen,unpack=True, dtype=float,comments='#',delimiter=' ')
    self.x=temp[0]
    self.y=temp[1]
    self.z=temp[2]
    self.minht=temp[3]
    self.WFGroundZ=temp[4]
    self.RH50=temp[5]
    self.RH60=temp[6]
    self.RH75=temp[7]
    self.RH90=temp[8]
    self.RH95=temp[9]
    self.CanopyZ=temp[10]
    self.canopycover=temp[11]
    self.shotN=temp[12]
    self.photonN=temp[13]
    self.iterationN=temp[14]
    self.refdem=temp[15]
    self.noiseInt=temp[16]
    self.signal=np.array(temp[17],dtype=np.int16)
    self.ground=temp[18]
    self.epsg=epsg
    return


  #################################

  def appendFile(self,namen):
    '''Append a file to the data'''
    temp=np.loadtxt(namen,unpack=True, dtype=float,comments='#',delimiter=' ')
    self.x=np.append(self.x,temp[0])
    self.y=np.append(self.y,temp[1])
    self.z=np.append(self.z,temp[2])
    self.minht=np.append(self.minht,temp[3])
    self.WFGroundZ=np.append(self.WFGroundZ,temp[4])
    self.RH50=np.append(self.RH50,temp[5])
    self.RH60=np.append(self.RH60,temp[6])
    self.RH75=np.append(self.RH75,temp[7])
    self.RH90=np.append(self.RH90,temp[8])
    self.RH95=np.append(self.RH95,temp[9])
    self.CanopyZ=np.append(self.CanopyZ,temp[10])
    self.canopycover=np.append(self.canopycover,temp[11])
    self.shotN=np.append(self.shotN,temp[12])
    self.photonN=np.append(self.photonN,temp[13])
    self.iterationN=np.append(self.iterationN,temp[14])
    self.refdem=np.append(self.refdem,temp[15])
    self.noiseInt=np.append(self.noiseInt,temp[16])
    self.signal=np.append(self.signal,temp[17])
    self.ground=np.append(self.signal,temp[18])

  #################################

  def writeHDF(self,outNamen):
    '''Write the output HDF5 file'''
    # convert some parameters
    numb=self.x.shape[0]
    self.setDists()
    delta_time=self.dists*0.7/7599.68
    self.setPhtCount()
    inProj=Proj(init="epsg:"+str(self.epsg))
    outProj=Proj(init="epsg:"+str(4326))
    lon,lat=transform(inProj,outProj,self.x,self.y)
    segment_dist_x=np.remainder(self.dists,20)
    segment_id=np.around(np.array(self.dists/20))
    # open output
    f=h5py.File(outNamen,'w')
    # make top level groups
    f.create_group('#ref#')
    f.create_group('gt1l')
    f['gt1l'].create_group('bckgrd_atlas')
    f['gt1l'].create_group('geolocation')
    f['gt1l'].create_group('heights')
    f['gt1l'].create_group('orbit_info')
    f['gt1l'].create_group('veg_truth')
    # populate data
    f['#ref#'].create_dataset('a',data=[0,0],compression='gzip')
    f['gt1l']['bckgrd_atlas'].create_dataset('bckgrd_rate',data=np.full(numb,1),compression='gzip')
    f['gt1l']['geolocation'].create_dataset('segment_dist_x',data=segment_dist_x,compression='gzip')
    f['gt1l']['geolocation'].create_dataset('segment_id',data=segment_id,compression='gzip')
    f['gt1l']['geolocation'].create_dataset('sigma_h',data=np.full(numb,0.4),compression='gzip')
    f['gt1l']['geolocation'].create_dataset('surf_type',data=np.full(numb,1),compression='gzip')
    f['gt1l']['geolocation'].create_dataset('segment_length',data=np.full(numb,20),compression='gzip')
    f['gt1l']['geolocation'].create_dataset('ph_index_beg',data=np.full(numb,1),compression='gzip')
    f['gt1l']['geolocation'].create_dataset('segment_ph_cnt',data=self.seg_phtcount,compression='gzip')
    f['gt1l']['geolocation'].create_dataset('solar_elevation',data=np.full(numb,20),compression='gzip')
    f['gt1l']['geolocation'].create_dataset('delta_time',data=delta_time,compression='gzip')
    f['gt1l']['heights'].create_dataset('delta_time',data=delta_time,compression='gzip')
    f['gt1l']['heights'].create_dataset('h_ph',data=self.z,compression='gzip')
    f['gt1l']['heights'].create_dataset('lon_ph',data=lon,compression='gzip')
    f['gt1l']['heights'].create_dataset('lat_ph',data=lat,compression='gzip')
    f['gt1l']['heights'].create_dataset('signal_conf_photon',data=self.signal,compression='gzip')
    f['gt1l']['heights'].create_dataset('dist_ph_across',data=np.full(numb,5),compression='gzip')
    f['gt1l']['heights'].create_dataset('dist_ph_along',data=np.full(numb,10),compression='gzip')
    f['gt1l']['heights'].create_dataset('pce_mframe_cnt',data=np.full(numb,1),compression='gzip')
    f['gt1l']['orbit_info'].create_dataset('rgt',data=(1))
    f['gt1l']['orbit_info'].create_dataset('cycle_number',data=np.full(numb,1),compression='gzip')
    f['gt1l']['veg_truth'].create_dataset('x_utm',data=self.x,compression='gzip')
    f['gt1l']['veg_truth'].create_dataset('y_utm',data=self.y,compression='gzip')
    f['gt1l']['veg_truth'].create_dataset('refDEM',data=self.refdem,compression='gzip')
    f['gt1l']['veg_truth'].create_dataset('rh50',data=self.RH50,compression='gzip')
    f['gt1l']['veg_truth'].create_dataset('rh60',data=self.RH60,compression='gzip')
    f['gt1l']['veg_truth'].create_dataset('rh75',data=self.RH75,compression='gzip')
    f['gt1l']['veg_truth'].create_dataset('rh90',data=self.RH90,compression='gzip')
    f['gt1l']['veg_truth'].create_dataset('rh95',data=self.RH95,compression='gzip')
    f['gt1l']['veg_truth'].create_dataset('canopyz',data=self.CanopyZ,compression='gzip')
    f['gt1l']['veg_truth'].create_dataset('canopy_cover',data=self.canopycover,compression='gzip')
    f['gt1l']['veg_truth'].create_dataset('isground',data=self.ground,compression='gzip')

    # close up
    f.close()
    print("Written to",outNamen)
    return


  #################################

  def setPhtCount(self):
    '''Set photon count rates per 20 m'''
    minD=np.min(self.dists)
    maxD=np.max(self.dists)
    self.seg_phtcount=np.histogram(self.dists,int((maxD-minD)/20.0+1))[0]

  #################################

  def setDists(self):
    '''Set distances'''
    self.minX=np.min(self.x)
    self.minY=np.min(self.y)
    self.dists=np.sqrt((self.x-self.minX)**2+(self.y-self.minY)**2)

  #################################

  def padData(self,minLen):
    '''Pad data to make a minimum transect length'''

    # footprint spacing
    step=0.7

    # determine track length
    dx=self.x[-1]-self.x[0]
    dy=self.y[-1]-self.y[0]
    trackLen=sqrt(dx**2+dy**2)
    dx=dx/self.x.shape[0]
    dy=dy/self.y.shape[0]

    # see if 
    if(trackLen<minLen):  # then we need to pad
      vectX=dx/trackLen
      vectY=dy/trackLen

      nPer=self.x.shape[0]
      nExtras=int(minLen/trackLen+1)
      print('Padding',nExtras,'times')

      # make copies of originals
      x=np.copy(self.x)
      y=np.copy(self.y)
      z=np.copy(self.z)
      minht=np.copy(self.minht)
      WFGroundZ=np.copy(self.WFGroundZ)
      RH50=np.copy(self.RH50)
      RH60=np.copy(self.RH60)
      RH75=np.copy(self.RH75)
      RH90=np.copy(self.RH90)
      RH95=np.copy(self.RH95)
      CanopyZ=np.copy(self.CanopyZ)
      canopycover=np.copy(self.canopycover)
      shotN=np.copy(self.shotN)
      photonN=np.copy(self.photonN)
      iterationN=np.copy(self.iterationN)
      refdem=np.copy(self.refdem)
      noiseInt=np.copy(self.noiseInt)
      signal=np.copy(self.signal)
      ground=np.copy(self.ground)

      for i in range(1,nExtras):
        # are we reversing?
        isOdd=i%2

        # make up new coordinates
        self.x=np.append(self.x,np.linspace(start=self.x[-1],stop=self.x[-1]+dx*nPer,num=nPer))
        self.y=np.append(self.y,np.linspace(start=self.y[-1],stop=self.y[-1]+dy*nPer,num=nPer))

        # alternately reverse data
        if(isOdd==0):   # even, keep in this order
          self.z=np.append(self.z,z)
          self.minht=np.append(self.minht,minht)
          self.WFGroundZ=np.append(self.WFGroundZ,WFGroundZ)
          self.RH50=np.append(self.RH50,RH50)
          self.RH60=np.append(self.RH60,RH60)
          self.RH75=np.append(self.RH75,RH75)
          self.RH90=np.append(self.RH90,RH90)
          self.RH95=np.append(self.RH95,RH95)
          self.CanopyZ=np.append(self.CanopyZ,CanopyZ)
          self.canopycover=np.append(self.canopycover,canopycover)
          self.shotN=np.append(self.shotN,shotN)
          self.photonN=np.append(self.photonN,photonN)
          self.iterationN=np.append(self.iterationN,iterationN)
          self.refdem=np.append(self.refdem,refdem)
          self.noiseInt=np.append(self.noiseInt,noiseInt)
          self.signal=np.append(self.signal,signal)
          self.ground=np.append(self.ground,ground)
        else:          # odd, reverse
          self.z=np.append(self.z,np.flip(z,0))
          self.minht=np.append(self.minht,np.flip(minht,0))
          self.WFGroundZ=np.append(self.WFGroundZ,np.flip(WFGroundZ,0))
          self.RH50=np.append(self.RH50,np.flip(RH50,0))
          self.RH60=np.append(self.RH60,np.flip(RH60,0))
          self.RH75=np.append(self.RH75,np.flip(RH75,0))
          self.RH90=np.append(self.RH90,np.flip(RH90,0))
          self.RH95=np.append(self.RH95,np.flip(RH95,0))
          self.CanopyZ=np.append(self.CanopyZ,np.flip(CanopyZ,0))
          self.canopycover=np.append(self.canopycover,np.flip(canopycover,0))
          self.shotN=np.append(self.shotN,np.flip(shotN,0))
          self.photonN=np.append(self.photonN,np.flip(photonN,0))
          self.iterationN=np.append(self.iterationN,np.flip(iterationN,0))
          self.refdem=np.append(self.refdem,np.flip(refdem,0))
          self.noiseInt=np.append(self.noiseInt,np.flip(noiseInt,0))
          self.signal=np.append(self.signal,np.flip(signal,0))
          self.ground=np.append(self.ground,np.flip(ground,0))
    return


# iceSim() class end
########################################


def readCommands():
  '''Read the command line'''
  p = argparse.ArgumentParser(description=("Convert ICESat-2 .pts sims to HDF5"))
  p.add_argument("--input",dest="inNamen",type=str,help=("Input filename"))
  p.add_argument("--inList",dest="inList",default="none",type=str,help=("Input filename list"))
  p.add_argument("--output",dest="outNamen",type=str,default='ice2.h5',help=("Output filename"))
  p.add_argument("--epsg",dest="epsg",type=int,default=32632,help=("Input EPSG"))
  p.add_argument("--minLen",dest="minLen",type=float,default=0,help=("Minimum acceptable length"))
  cmdargs = p.parse_args()
  return cmdargs


########################################
# read mulriple ICESat-2 files

def readMultiSim(inList,epsg):
  '''Read multiple ICEsat-2 files'''
  # read list of files
  files=np.genfromtxt(inList,unpack=True,usecols=(0,),dtype=str,comments='#',delimiter=" ")

  # determine file order
  fileInds=np.empty(files.shape,dtype=int)
  for i in range(0,files.shape[0]):
    fileInds[i]=int(files[i].split('.')[-2])
  
  # determine index order
  indOrder=np.argsort(fileInds)

  # read first file
  data=iceSim(files[indOrder[0]],epsg)

  # loop and append rest
  for i in indOrder[1:]:
    data.appendFile(files[i])

  return(data)


########################################
# Main block

if __name__ == '__main__':
  # read commands
  cmdargs=readCommands()
  # read data
  if(cmdargs.inList=="none"):
    data=iceSim(cmdargs.inNamen,cmdargs.epsg)
  else:
    data=readMultiSim(cmdargs.inList,cmdargs.epsg)
  # pad data if needed
  data.padData(cmdargs.minLen)
  # write data
  data.writeHDF(cmdargs.outNamen)

# The end
########################################

