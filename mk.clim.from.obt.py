from cf2.io.GPM import GPM
from numpy import *
import datetime
import calendar
import h5py
import sys, os
Y  = 2014
lM = [6]

#prjName  = "GPM.KuPR"
#prdLv    = "L2"
#prdVer   = "02"
#varname  = "NS/SLV/precipRateESurface"

prjName  = "GPM.DPR"
prdLv    = "L2"
prdVer   = "03"
varname  = "NS/SLV/precipRateESurface"


gpm = GPM(prjName, prdLv, prdVer)
region = "JPN"

#---- BBox ---------
if region == "JPN":
  BBox= [[20.0, 118.0], [48.0, 150.0]]
else:
  print "check region !!"
  sys.exit()
#------------------
ny  = 280  # res=0.1
nx  = 320  # res=0.1
a2one = ones([ny,nx], float32)
a2s   = zeros([ny,nx], float32)
a2n   = zeros([ny,nx], float32)
miss  = -9999.0

for M in lM:
  sD  = 1
  eD = calendar.monthrange(Y,M)[1]
  #eD  = sD    # test
  for D in range(sD,eD+1):
    for H in range(0, 23+1):
      sDTime = datetime.datetime(Y,M,D,H)
      eDTime = sDTime + datetime.timedelta(hours=1)
      try:
        gpmmap = gpm(varname, sDTime, eDTime, BBox, res=0.1)
        a2dat  = array(gpmmap.griddata[0], float32) 
        a2dat  = ma.masked_less(gpmmap.griddata[0], 0.0)
        a2s    = a2s + a2dat.filled(0.0)
        a2n    = a2n + ma.masked_array( a2one, a2dat.mask).filled(0.0)
      except ValueError:
        print ""
        print "SKIP!! ValueError", sDTime
        print ""
      except IndexError :
        print ""
        print "SKIP!! IndexError", sDTime
  #--- convert unit ---------
  a2s      = a2s / 60./60.   # mm/hour --> mm/sec
  #--- mean precipitaion rate ----
  a2r      = (ma.masked_where(a2n ==0.0, a2s)/a2n).filled(miss)

  #--- write date to files --
  odir_base= "/tank/utsumi/PMM/CLIM.OBT/%s.%s.%s.%s"%(prjName, prdLv, prdVer, region)
  odir     = odir_base + "/%04d"%(Y)
  try:
    os.makedirs(odir)
  except OSError:
    pass
  print odir
  
  sumname  = odir + "/sum.%04d.%02d.%dx%d"%(Y,M,ny,nx)
  numname  = odir + "/num.%04d.%02d.%dx%d"%(Y,M,ny,nx)
  ratname  = odir + "/rat.%04d.%02d.%dx%d"%(Y,M,ny,nx)

  a2s      = array(a2s, float32)
  a2n      = array(a2n, float32)
  a2r      = array(a2r, float32)

  a2s.tofile(sumname)
  a2n.tofile(numname)
  a2r.tofile(ratname)
  print sumname


