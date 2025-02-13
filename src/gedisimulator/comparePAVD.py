
'''
Compare PVD profiles from simulated
and real GEDI data. Aimed to solve the
scaling issue noticed by Pat Burns
'''

import numpy as np
import h5py
import matplotlib.pyplot as plt
import argparse
from sys import exit


########################################################

def readCommands():
  '''
  Read commandline arguments
  '''
  p = argparse.ArgumentParser(description=("Writes out properties of GEDI waveform files"))
  p.add_argument("--l2b",dest="l2bName",type=str,help=("Input GEDI L2B HDF5 filename"))
  p.add_argument("--real",dest="realName",type=str,help=("Input real data gediMetric metrics"))
  p.add_argument("--sim",dest="simName",type=str,help=("Input simulated data gediMetric metrics"))
  p.add_argument("--outRoot",dest="outRoot",type=str,default='test',help=("Output graph filename root"))
  cmdargs = p.parse_args()
  return cmdargs

########################################################

class l2bMetrics():
  '''Class t hold L2B metrics'''

  ######################

  def __init__(self,filename,mode=0):
    '''Class initialiser'''

    if(mode==1):
      self.readSim(filename)
    elif(mode==0):
      self.readReal(filename)
    elif(mode==2):
      self.readSim(filename,gCol=1,root=" tLAI")
    else:
      print("Mode",mode,"not recognised")
      exit(1)

    return


  ######################

  def readSim(self,filename,gCol=5,root=" gLAI"):
    '''Read simulated metrics'''

    f=open(filename, 'r')
    lines=f.readlines()
    self.nWaves=len(lines)-1

    # read header
    bits=lines[0].replace("#","").split(",")
    f.close()

    # find columns needed
    idCol=0

    colList=[]
    i=0
    self.nLAI=0
    for bit in bits:
      if( root in bit ):
        colList.append(i)
        self.nLAI+=1
      i+=1

    # read whole file
    self.pavd=np.full((self.nWaves,self.nLAI),-999.0)
    self.zG=np.full((self.nWaves),-999.0)
    self.waveID=np.genfromtxt(filename,usecols=(idCol),unpack=True,comments="#",dtype=str)

    f=open(filename, 'r')
    lines=f.readlines()
    i=0
    for line in lines:
      if "#" not in line:
        bits=line.split(" ")
        self.zG[i]=float(bits[gCol])
        for j in range(0,self.nLAI):
          self.pavd[i,j]=float(bits[colList[j]])
        i+=1

    f.close()
    return


  ######################

  def readReal(self,filename):
    '''Read an L2B file'''

    beamList=['BEAM0000', 'BEAM0001', 'BEAM0010', 'BEAM0011', 'BEAM0101', 'BEAM0110', 'BEAM1000', 'BEAM1011']


    self.waveID=np.array((),dtype=np.uint)
    self.beamID=np.array((),dtype=str)
    #self.zG=np.array((),dtype=float)   # need to read this from L2A file
    laiStarted=False

    f=h5py.File(filename,'r')

    # loop over beams
    for beam in beamList:
      # read shots in this beam

      tempPAVD=np.array(f[beam]['pavd_z'])
      tempID=np.array(f[beam]['shot_number'])
      tempBeam=np.repeat(beam,tempID.shape[0])
      #tempZG=np.array(f[beam]['shot_number'])

      self.waveID=np.append(self.waveID,tempID)
      self.beamID=np.append(self.beamID,tempBeam)
      if(laiStarted):
        self.pavd=np.append(self.pavd,tempPAVD,axis=0)
      else:
        self.pavd=tempPAVD

      laiStarted=True

    self.nLAI=tempPAVD.shape[1]
    self.nWaves=self.waveID.shape[0]
    f.close()

    return


########################################################

def comparePAVD(sim,real,l2b,outRoot):
  '''Compare PVD estimates'''

  nDP=3
  res=5
  tol=0.0000001

  # loop over one of the files
  for i in range(0,sim.nWaves):
    # find associated shots
    thisID=sim.waveID[i]
    beam=thisID.split(".")[1]
    shotN=np.uint(thisID.split(".")[2])
    realID="gedi."+str(beam)+"."+str(shotN)


    outname=outRoot+"."+str(beam)+"."+str(shotN)+".png"


    l2bPAVD=l2b.pavd[(l2b.beamID==beam)&(l2b.waveID==shotN)][0]
    simPAVD=sim.pavd[i]

    # total leaf area
    simA=np.sum(simPAVD)*res
    l2A=np.sum(l2bPAVD)*res
    simH=simPAVD[simPAVD>tol].argmax()*res
    l2H=l2bPAVD[l2bPAVD>tol].argmax()*res
    print("Areas",round(simA,3),round(l2A,3),round(simA/l2A,3),round(l2A/simA,3),"height",simH,l2H)

    realPAVD=real.pavd[real.waveID==realID][0]
    plt.plot(l2bPAVD,np.arange(0,l2bPAVD.shape[0]*res,res),label='l2b')
    plt.plot(simPAVD,np.arange(0,simPAVD.shape[0]*res,res),label='sim')
    #plt.plot(realPAVD,np.arange(0,realPAVD.shape[0]*res,res),label='real')
    plt.legend()
    plt.xlabel("PAVD (m^2/m^3)")
    plt.ylabel("Height (m)")
    plt.savefig(outname)
    plt.cla()
    plt.clf()
    #print("Drawn to",outname)

    #for j in range(0,sim.nLAI):
    #  print(i,j,'l2b',round(l2bPAVD[j],nDP),'sim',round(simPAVD[j],nDP),'real',round(realPAVD[j],nDP))

  return


########################################################

if __name__ == '__main__':
  cmd=readCommands()

  # read each dataset
  sim=l2bMetrics(cmd.simName,mode=2)
  real=l2bMetrics(cmd.realName,mode=1)
  l2b=l2bMetrics(cmd.l2bName,mode=0)

  # Align datasets by shot ID
  comparePAVD(sim,real,l2b,cmd.outRoot)


########################################################

