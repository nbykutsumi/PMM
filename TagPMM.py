from numpy import *
import detect.Tag as Tag
from datetime import datetime, timedelta
import myfunc.IO.RadarAMeDAS as RadarAMeDAS
import os
ra  = RadarAMeDAS.RadarAMeDAS(prj="ra_0.1")
TAG = Tag.Tag()

class TagPMM(object):
  def __init__(self):
    self.baseDir = "/tank/utsumi/PMM/RAdom/tag"
    self.res     = 0.1
    self.BBox    = ra.BBox   # [[20.0, 118.0],[48.0, 150.0]]
    self.Lat     = arange(self.BBox[0][0]+self.res/2,
                    self.BBox[1][0]+self.res/10,
                    self.res, dtype='float64')
    
    self.Lon     = arange(self.BBox[0][1]+self.res/2,
                    self.BBox[1][1]+self.res/10,
                    self.res, dtype='float64')
    self.ny  = len(self.Lat)
    self.nx  = len(self.Lon)
 
  def loadTag(self, tag, DTime):
    Year, Mon, Day, Hour = DTime.year, DTime.month, DTime.day, DTime.hour
    sDir  = os.path.join(self.baseDir,tag, "%04d%02d"%(Year,Mon))
    sPath = os.path.join(sDir, "%s.%04d%02d%02d%02d.%dx%d"%(tag,Year,Mon,Day,Hour,self.ny,self.nx))
    return fromfile(sPath, float32).reshape(self.ny,self.nx)

  def dictTagFrac(self, ltag, DTime, miss=0):
    dictTag = {}
    for tag in ltag:
      dictTag[tag] = self.loadTag(tag, DTime)
    return TAG.mkMaskFracCore(dictTag, miss=miss)

    
