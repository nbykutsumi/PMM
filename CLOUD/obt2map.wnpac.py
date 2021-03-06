from numpy import *
from PMM.pmm_fsub import *
from myfunc.IO import GPyM
from datetime import datetime, timedelta
import myfunc.util as util
import os

#iDTime  = datetime(2014,3,31,0)
#eDTime  = datetime(2016,1,2,0)

iDTime  = datetime(2014,5,27,0)
eDTime  = datetime(2014,6,1,0)

#iDTime  = datetime(2015,6,1,0)
#eDTime  = datetime(2015,7,1,0)




prj     = 'GPM.KuPR'
prdLv   = 'L2'
#prdVer  = '03'
prdVer  = '05'
var     = 'NS/SLV/precipRateESurface'

#prj     = 'GPM.KaPR'
#prdLv   = 'L2'
#prdVer  = '03'
##var     = 'HS/SLV/precipRateESurface'
#var     = 'MS/SLV/precipRateESurface'

#prj     = 'GPM.GMI'
#prdLv   = 'L2'
#prdVer  = '05'
#var     = 'S1/surfacePrecipitation'

#prj     = "TRMM.PR"
#prdLv   = "L2A25"
#prdVer  = "07"
#var     = "e_SurfRain"

#prj     = "TRMM.TMI"
#prdLv   = "L2A12"
#prdVer  = "07"
#var     = "surfacePrecipitation"

BBox    = [[-0.1, 113.875],[52.1, 180.125]]
ny,nx   = 261, 265
miss_in = -9999.9
miss_out= -9999.

dDTime  = timedelta(hours=0.5)
lDTime  = util.ret_lDTime(iDTime, eDTime, dDTime)
gpm     = GPyM.GPM(prj, prdLv, prdVer)


for DTime in lDTime:
  try:
    gpmobt  = gpm(var, DTime, DTime+dDTime, BBox=BBox)
  except IndexError:
    print "---------------"
    print "No data in BBox",DTime
    print "---------------"
    continue
  except ValueError:
    print "---------------"
    print "Value Error",DTime
    print "---------------"
    continue


  nl, nw        = gpmobt.data.shape
  a2sum, a2num  = pmm_fsub.obt2wnpac261x265(gpmobt.data.T, gpmobt.lon.T, gpmobt.lat.T, nw, nl)

  a2sum   = a2sum.T
  a2num   = a2num.T
  a2pr    = ma.masked_invalid(a2sum/a2num)
  #a2pr    = a2pr / (60.*60.)   # mm/h --> mm/s
  a2pr    = a2pr.filled(miss_out)
  #--- save ---------
  Year    = DTime.year
  Mon     = DTime.month
  Day     = DTime.day
  Hour    = DTime.hour
  Minute  = DTime.minute
  baseDir = "/tank/utsumi/PMM/WNP.261x265/%s/%s/%s"%(prj,prdLv,prdVer)
  dataDir = baseDir + "/%04d/%02d"%(Year,Mon)
  dataPath= dataDir + "/pr.%02d.%02d.%02d.%dx%d"%(Day,Hour,Minute,ny,nx)

  util.mk_dir(dataDir)
  a2pr.astype(float32).tofile(dataPath)
  print "OUTPUT"
  print dataPath
    
  



