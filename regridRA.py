from numpy import *
from myfunc.IO import RadarAMeDAS
from datetime import datetime, timedelta
import myfunc.util as util
import myfunc.regrid.Regrid as Regrid
import os

iDTime = datetime(2014,1,1,0)
#eDTime = datetime(2014,1,1,3)
eDTime = datetime(2014,12,31,23)
dDTime = timedelta(hours=1)
lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)
miss_out = -999.
resUp  = 0.1
ra = RadarAMeDAS.RadarAMeDAS(prj="ra_0.01")
reg= Regrid.UpScale()


BBox   = ra.BBox
LatOrg = ra.Lat
LonOrg = ra.Lon
LatUp  = arange(BBox[0][0]+resUp/2,
                BBox[1][0]+resUp/10,
                resUp, dtype='float64')

LonUp  = arange(BBox[0][1]+resUp/2,
                BBox[1][1]+resUp/10,
                resUp, dtype='float64')

nyUp, nxUp = len(LatUp), len(LonUp)
baseDir  = "/tank/utsumi/data/RadarAMeDAS/ra_%s"%(resUp)
#--- save meta data -----------
util.mk_dir(baseDir)
sLat = util.array2csv(LatUp)
sname= baseDir+"/lat.csv"
f=open(sname,"w");f.write(sLat); f.close()

sLon = util.array2csv(LonUp)
sname= baseDir+"/lon.csv"
f=open(sname,"w");f.write(sLon); f.close()
print sname

#------------------------------
def ret_oPath(DTime, resUp):
  Year, Mon, Day, Hour, Min =\
    DTime.year, DTime.month, DTime.day, DTime.hour, DTime.minute
  oDir     = os.path.join(baseDir, "%04d%02d"%(Year,Mon))
  oPath    = os.path.join(oDir, "RadarAmedas.%04d%02d%02d%02d%02d00.%dx%d"%(Year,Mon,Day,Hour,Min,nyUp, nxUp))
  return oDir, oPath

reg(LatOrg, LonOrg, LatUp, LonUp, globflag=False)
for DTime in lDTime:
  a2org = ra.loadForward_mmh(DTime)
  a2up  = reg.upscale(a2org, miss_in=ra.miss, miss_out = miss_out)
  oDir, oPath = ret_oPath(DTime, resUp)
  util.mk_dir(oDir)
  a2up.astype(float32).tofile(oPath)
  print oPath
