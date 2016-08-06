from numpy import *

sDir = "/home/utsumi/test/AccEst"
#dattype = "IMERG"
#ny   = 522
#nx   = 662

dattype = "IMERG.MW"
ny   = 522
nx   = 662

#dattype = "RA"
#ny   = 2800
#nx   = 3200

#dattype = "GSMaP.MW"
#ny   = 522
#nx   = 662

#dattype = "GSMaP.IR"
#ny   = 522
#nx   = 662





Year = 2014
lMon = range(4,11+1)

a2out = zeros([ny,nx],float32)
for Mon in lMon:
  sPath= sDir + "/accEst.%s.%04d.%02d.%dx%d"%(dattype, Year,Mon,ny,nx)
  a    = fromfile(sPath, float32).reshape(ny,nx)
  a    = ma.masked_less(a,0.0)
  a2out= a2out + a
  print dattype,Mon, a.mean()

a2out = a2out / len(lMon)
print dattype,"Mean",a2out.mean()
