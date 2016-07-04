from numpy import *
from cf2.io import GPM
from datetime import datetime, timedelta
import myfunc.util as util
import os

prj     = 'GPM.KuPR'
prdLv   = 'L2'
prdVer  = '03'
var     = 'NS/SLV/precipRateESurface'

#prj     = 'GPM.KaPR'
#prdLv   = 'L2'
#prdVer  = '03'
##var     = 'HS/SLV/precipRateESurface'
#var     = 'MS/SLV/precipRateESurface'

#prj     = 'GPM.GMI'
#prdLv   = 'L2'
#prdVer  = '02'
#var     = 'S1/surfacePrecipitation'

#prj     = "TRMM.PR"
#prdLv   = "L2A25"
#prdVer  = "07"
#var     = "e_SurfRain"

#prj     = "TRMM.TMI"
#prdLv   = "L2A12"
#prdVer  = "07"
#var     = "surfacePrecipitation"

BBox    = [[20.,118.],[48.,150.]]  # RadarAMeDAS
res     = 0.1
miss_in = -9999.9
miss_out= -9999.

if res == 0.1:
  ny,nx = 280,320

#iDTime  = datetime(2014,12,1,0)
#eDTime  = datetime(2015,1,1,0)

iDTime  = datetime(2014,12,1,0)
eDTime  = datetime(2015,1,1,0)
dDTime  = timedelta(hours=1)
lDTime  = util.ret_lDTime(iDTime, eDTime, dDTime)
gpm     = GPM.GPM(prj, prdLv, prdVer)

baseDir = "/tank/utsumi/PMM/RAdom/%s/%s/%s"%(prj,prdLv,prdVer)

for DTime in lDTime:
  try:
    gpmdom  = gpm(var, DTime, DTime+dDTime, BBox, res)
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
  #
  a2dat = ma.masked_equal(gpmdom.griddata, miss_in).mean(axis=0).filled(miss_out)

  Year, Mon, Day, Hour, Min =\
    DTime.year, DTime.month, DTime.day, DTime.hour, DTime.minute 
  oDir  = os.path.join(baseDir,"%04d%02d"%(Year,Mon))
  oPath = os.path.join(oDir, "mmh.%04d%02d%02d%02d.%dx%d"%(Year,Mon,Day,Hour,ny,nx))
  util.mk_dir(oDir)
  a2dat.astype(float32).tofile(oPath)
  print oPath
    
  



