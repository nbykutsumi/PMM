from numpy import *
from datetime import datetime, timedelta
from myfunc.regrid import Regrid
import detect.Tag as Tag
import myfunc.IO.RadarAMeDAS as RadarAMeDAS
import myfunc.util as util
import os


iDTime = datetime(2014,4,1,0)
#eDTime = datetime(2014,12,31,18)
eDTime = datetime(2014,12,31,18)
#iDTime = datetime(2014,1,1,0)
#eDTime = datetime(2014,4,1,0)

dDTime = timedelta(hours=6)
lDTime = util.ret_lDTime(iDTime,eDTime,dDTime)
ra     = RadarAMeDAS.RadarAMeDAS(prj="ra_0.1")
BBox   = ra.BBox   # [[20.0, 118.0],[48.0, 150.0]]
res    = 0.1
miss_tag=0.0
model  = "JRA55"
TAG    = Tag.Tag(model=model, res="bn", miss=miss_tag)
ltag   = ["tc","fbc","c","cf","nbc"]
#ltag   = ["cf"]
LatIn  = TAG.Lat
LonIn  = TAG.Lon
LatOut = arange(BBox[0][0]+res/2,
                BBox[1][0]+res/10,
                res, dtype='float64')

LonOut = arange(BBox[0][1]+res/2,
                BBox[1][1]+res/10,
                res, dtype='float64')
nyOut  = len(LatOut)
nxOut  = len(LonOut)

baseDir = "/tank/utsumi/PMM/RAdom/tag"
#-------------------------
def biIntpMask(a2Mask, LatIn, LonIn, LatOut, LonOut, miss=miss_tag):
  aMaskFrac     = Regrid.biIntp(\
                    LatIn, LonIn   \
                  , a2Mask \
                  , LatOut, LonOut \
                  )[0]

  return          ma.masked_greater_equal(\
                      ma.masked_less(aMaskFrac, 0.5).filled(miss_tag) \
                    , 0.5 \
                   ).filled(1.0)  
#-------------------------
DTime = iDTime
MonPre = -9999
for DTime in lDTime:
  Year, Mon, Day, Hour = DTime.year, DTime.month, DTime.day, DTime.hour
  if MonPre != Mon:
    MonPre = Mon
    lYM  = [Year, Mon]
    TAG.init_cyclone(lYM, lYM, model=model, tctype="bst")
  for tag in ltag:
    a2MaskIn = TAG.mkMask(tag, DTime)
    a2Mask   = biIntpMask(a2MaskIn, LatIn, LonIn, LatOut, LonOut, miss=miss_tag)
    oDir     = os.path.join(baseDir, tag,"%04d%02d"%(DTime.year, DTime.month))
    oPath    = os.path.join(oDir, "%s.%04d%02d%02d%02d.%dx%d"%(tag,Year,Mon,Day,Hour,nyOut,nxOut))
    util.mk_dir(oDir)
    a2Mask.astype(float32).tofile(oPath)
    print oPath

