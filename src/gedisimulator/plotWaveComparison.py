import h5py
import numpy as np
import string
import matplotlib.pyplot as plt
import matplotlib
import argparse

# import a GEDI data handler from this repository
import sys
import os
sys.path.append(os.environ['GEDIRAT_ROOT'])
from gediHandler import gediData


###########################################

if __name__ == '__main__':
  def gediCommands():
    '''
    Read commandline arguments
    '''
    p = argparse.ArgumentParser(description=("Draws comparisons of GEDI and simulated waveforms. NOTE that the waves should be aligned with collocateWaves first."))
    p.add_argument("--simWave",dest="simWave",type=str,help=("Input simulated HDF5 filename"))
    p.add_argument("--gediWave",dest="gediWave",type=str,help=("Input GEDI HDF5 filename"))
    p.add_argument("--outRoot",dest="outRoot",type=str,default='test',help=("Output graph filename root"))
    p.add_argument("--offset",dest="offset",type=float,default=0,help=("Datum difference, if not already applied"))
    cmdargs = p.parse_args()
    return cmdargs


###########################################
# main block

if __name__ == '__main__':
  # read the command line
  cmd=cmdargs=gediCommands()

  # read sim HDF5
  sim=gediData(filename=cmd.simWave)

  # read GEDI HDF5
  gedi=gediData(filename=cmd.gediWave)
  matplotlib.rcParams.update({'font.size': 16})

  # loop through waves
  for i in range(0,sim.nWaves):

    # match IDs
    thisShotN=int(sim.waveID[i].split('.')[2])
    thisBeam=sim.waveID[i].split('.')[1]
    realInd=np.where((thisShotN==gedi.waveID)&(thisBeam==gedi.beamID))


    if(len(realInd)>0):
      if(len(realInd[0])>0):
        realInd=realInd[0][0]
      else:
        continue
    else:
      continue

    # set z arrays
    gedi.setOneZ(realInd)
    sim.setOneZ(i)


    # mean noise level
    noiseBins=int(10/0.15)
    meanN=np.mean(gedi.wave[realInd,:noiseBins])
    stdev=np.std(gedi.wave[realInd,:noiseBins])
    thresh=meanN+3.5*stdev
    totE=np.sum(gedi.wave[realInd,gedi.wave[realInd]>thresh]-meanN)


    plt.ioff()
    plt.plot(sim.wave[i]*totE/np.sum(sim.wave[i]),sim.z+cmd.offset,label="Simulation")
    plt.plot(gedi.wave[realInd,0:gedi.nBins]-meanN,gedi.z,label="GEDI")

    #plt.plot(alignG,z,label="ALS ground")
    plt.legend()
    plt.xlim(left=0)
    #plt.ylim((bot,top))
    plt.xlabel('DN')
    plt.ylabel('Height (m)')

    namen="%s.%s.%d.png" % (cmd.outRoot,thisBeam,thisShotN)

    plt.savefig(namen)
    plt.close()
    plt.clf()
    print("Written to",namen)

