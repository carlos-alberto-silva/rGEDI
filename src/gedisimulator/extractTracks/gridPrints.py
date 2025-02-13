import numpy as np
import string
import optparse

# read command line
def clParser():
    desc = """Script extracts the MODIS binary layers from the HDFs and dumps and envi (*hopefully*) hdr file"""
    clp = optparse.OptionParser("Usage: %prog [options] filename", description = desc)
    clp.add_option("--input",type="string",dest="inNamen",help="input file name")
    clp.add_option("--output",type="string",dest="outNamen",help="output file name",default="teast")
    clp.add_option("--res",type="float",dest="res",help="grid resolution",default=1000)
    return clp.parse_args()
(options, args) = clParser()
res=options.res

# read data
x,y=np.loadtxt(options.inNamen, usecols=(0,1), unpack=True, dtype=float,comments='#')
waveID=np.genfromtxt(options.inNamen,usecols=(2,), unpack=True, dtype=str,comments='#')


# bounds
minX=min(x)
maxX=max(x)
minY=min(y)
maxY=max(y)
nX=int((maxX-minX)/res+1)
nY=int((maxY-minY)/res+1)

# open output
f=open(options.outNamen,"w")

# loop over tiles
for i in range(0,nX):
  x0=i*res+minX
  x1=x0+res
  for j in range(0,nY):
    y0=j*res+minY
    y1=y0+res
    useInd=np.where((x>=x0)&(x<x1)&(y>=y0)&(y<y1))
    if(len(useInd)>0):
      useInd=useInd[0]
      for k in useInd:
        line=str(x[k])+" "+str(y[k])+" "+waveID[k]+"\n"
        f.write(line)

f.close()
print("Written to",options.outNamen)

