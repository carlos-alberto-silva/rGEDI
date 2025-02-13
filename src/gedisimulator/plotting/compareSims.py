
'''
Compares simulated to real
GEDI waveforms

S. Hancock, 2019
'''


############################################

import numpy as np
import argparse
import matplotlib.pyplot as plt
from gediHandler import gediData


############################################

def readCommands():
  '''
  Read commandline arguments
  '''
  p = argparse.ArgumentParser(description=("Writes out properties of GEDI waveform files"))
  p.add_argument("--real",dest="realName",type=str,help=("Input real GEDI HDF5 filename"))
  p.add_argument("--sim",dest="simName",type=str,help=("Input simulated GEDI HDF5 filename"))
  p.add_argument("--bounds", dest ="bounds", type=float,nargs=4,default=[-100000000,-100000000,100000000000,10000000000], help=("Bounds to plot between. minX minY maxX maxY"))
  p.add_argument("--outRoot",dest="outRoot",type=str,default='test',help=("Output graph filename root"))
  cmdargs = p.parse_args()
  return cmdargs


############################################

if __name__ == '__main__':
  # read the command line
  cmdargs=readCommands()
  # read sim
  sim=gediData()
  # read real
