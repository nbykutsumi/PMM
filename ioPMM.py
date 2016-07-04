from numpy import *
from datetime import datetime, timedelta
from TagPMM import TagPMM
import myfunc.IO.GSMaP as GSMaP
import os, sys

class ioRAdom(TagPMM):
  def __init__(self):
    TagPMM.__init__(self) 

  def init_GSMaP(self, prj="standard", ver="v6"):
    self.gsmap=GSMaP.GSMaP(prj=prj, ver=ver, BBox=self.BBox)

  def ret_pathMonRAmmh(self, Year, Mon):
    res  = 0.1
    baseDir = "/tank/utsumi/data/RadarAMeDAS/ra_%s"%(res)
    sDir    = os.path.join(baseDir, "%04d%02d"%(Year,Mon))
    sumPath = os.path.join(sDir, "RadarAmedas.Sum.%04d%02d.%dx%d"%(Year,Mon,self.ny, self.nx))
    numPath = os.path.join(sDir, "RadarAmedas.Num.%04d%02d.%dx%d"%(Year,Mon,self.ny, self.nx))
    avePath = os.path.join(sDir, "RadarAmedas.%04d%02d.%dx%d"%(Year,Mon,self.ny, self.nx))
    return baseDir, sDir, sumPath, numPath, avePath

  def loadMonRAmmh(self, Year, Mon, var="ave",maskflag=True):
    """
    var = "sum", "num", "ave"
    """
    miss = -999.
    baseDir, sDir, sumPath, numPath, avePath = self.ret_pathMonRAmmh(Year,Mon)
    if   var=="sum": sPath=sumPath
    elif var=="num": sPath=numPath
    elif var=="ave": sPath=avePath

    if maskflag ==True:
      return ma.masked_equal(fromfile(sPath, float32).reshape(self.ny, self.nx), miss)
    else:
      return fromfile(sPath, float32).reshape(self.ny, self.nx)

  def ret_pathRAmmh(self, DTime):
    res  = 0.1
    Year, Mon, Day, Hour, Min =\
       DTime.year, DTime.month, DTime.day, DTime.hour, DTime.minute
    baseDir = "/tank/utsumi/data/RadarAMeDAS/ra_%s"%(res)
    sDir    = os.path.join(baseDir, "%04d%02d"%(Year,Mon))
    sPath   = os.path.join(sDir, "RadarAmedas.%04d%02d%02d%02d%02d00.%dx%d"%(Year,Mon,Day,Hour,Min, self.ny, self.nx))
    return baseDir, sDir, sPath

  def loadRAmmh(self, DTime, maskflag=True):
    miss = -999.
    baseDir, sDir, sPath = self.ret_pathRAmmh(DTime)
    if maskflag ==True:
      return ma.masked_equal(fromfile(sPath, float32).reshape(self.ny, self.nx), miss)
    else:
      return fromfile(sPath, float32).reshape(self.ny, self.nx)

  def loadGPMmmh(self, DTime, prj="GPM.KuPR",maskflag=True, return_dummy=False):
    """ 
    return_dummy=False: Return TabError
    return_dummy=True:  Return dummy array if file does not exist
    """ 
    if prj   == 'GPM.KuPR':
      prdLv   = 'L2'
      prdVer  = '03'
      #var     = 'NS/SLV/precipRateESurface'

    elif prj == "GPM.KaPR":
      prdLv   = 'L2'
      prdVer  = '03'
      #var     = 'HS/SLV/precipRateESurface'
      #var     = 'MS/SLV/precipRateESurface'

    elif prj == 'GPM.GMI':
      prdLv   = 'L2'
      prdVer  = '02'
      #var     = 'S1/surfacePrecipitation'     
    else:
      print "in ioPMM.loadGPMmmh"
      print "invalid prj", prj
      sys.exit()

    miss = -9999.
    baseDir = "/tank/utsumi/PMM/RAdom/%s/%s/%s"%(prj,prdLv,prdVer)
    Year, Mon, Day, Hour, Min =\
      DTime.year, DTime.month, DTime.day, DTime.hour, DTime.minute
    sDir  = os.path.join(baseDir,"%04d%02d"%(Year,Mon))
    sPath = os.path.join(sDir, "mmh.%04d%02d%02d%02d.%dx%d"%(Year,Mon,Day,Hour,self.ny,self.nx))
    if os.path.exists(sPath):
      if maskflag ==True:
        return ma.masked_equal(fromfile(sPath, float32).reshape(self.ny, self.nx), miss)
      else:
        return fromfile(sPath, float32).reshape(self.ny, self.nx)
    else:
      if return_dummy != True:
        #raise Exception("NoObsError")
        raise TabError

      elif (os.path.exists(sDir) and (len(os.listdir(sDir))!=0)):
        print "no track in the domain: return miss"
        if maskflag ==True:
          return ma.masked_equal(ones([self.ny, self.nx], float32), 1.0)
        else:
          return ones([self.ny, self.nx], float32)*miss
      else:
        print "no file and no directory",sDir
        sys.exit()

  def loadGSMaPmmh(self, DTime, maskflag=True):
    if maskflag==True:
      if maskflag==True:
        return ma.masked_less(self.gsmap.load_mmh(DTime),0.0)
      else:
        return self.gsmap.load_mmh(DTime)

  
class matchRAdom(TagPMM):
  def __init__(self):
    TagPMM.__init__(self)

  def __call__(self, prtype1="RA", prtype2="GPM.KuPR", Trc="ALL", nclass=4):
    self.baseDir = "/tank/utsumi/PMM/RAdom/%s.vs.%s.on.%s.%sclass"%(prtype1, prtype2, Trc, nclass)
 
  def ret_monSumPath(self, prtype, tag, Year,Mon):
    ny,nx    = self.ny, self.nx
    sDir     = os.path.join(self.baseDir, "%04d"%(Year))
    sumPath  = os.path.join(sDir, "sum.%s.%s.%04d%02d.%dx%d"%(prtype,tag,Year,Mon,ny,nx))
    return sDir, sumPath
  
  def ret_monNumPath(self, tag, Year,Mon):
    ny,nx    = self.ny, self.nx
    sDir     = os.path.join(self.baseDir, "%04d"%(Year))
    numPath  = os.path.join(sDir, "num.%s.%04d%02d.%dx%d"%(tag,Year,Mon,ny,nx))
    return sDir, numPath
  
  def ret_monSum(self, prtype, tag, Year,Mon):
    ny,nx    = self.ny, self.nx
    sDir     = self.ret_monSumPath(prtype,tag,Year,Mon)[0]
    sumPath  = self.ret_monSumPath(prtype,tag,Year,Mon)[1]
    return fromfile(sumPath, float32).reshape(ny,nx)
  
  def ret_monNum(self, tag, Year,Mon):
    ny,nx    = self.ny, self.nx
    sDir     = self.ret_monNumPath(tag,Year,Mon)[0]
    numPath  = self.ret_monNumPath(tag,Year,Mon)[1]
  
    return fromfile(numPath, float32).reshape(ny,nx)
  
  def ret_monRte(self, prtype, tag, Year,Mon):
    ny,nx    = self.ny, self.nx
    sDir     = self.ret_monSumPath(prtype,tag,Year,Mon)[0]
    sumPath  = self.ret_monSumPath(prtype,tag,Year,Mon)[1]
    numPath  = self.ret_monNumPath(tag,Year,Mon)[1]
  
    a2sum    = fromfile(sumPath, float32).reshape(ny,nx)
    a2num    = fromfile(numPath, float32).reshape(ny,nx)
    return ma.masked_where(a2num==0.0, a2sum)/a2num

  def ret_monPrc(self, prtype, tag, Year,Mon):
    ny,nx    = self.ny, self.nx
    sDir     = self.ret_monSumPath(prtype,tag,Year,Mon)[0]
    sumPath  = self.ret_monSumPath(prtype,tag,Year,Mon)[1]
    numPath  = self.ret_monNumPath("plain",Year,Mon)[1]
  
    a2sum    = fromfile(sumPath, float32).reshape(ny,nx)
    a2num    = fromfile(numPath, float32).reshape(ny,nx)
    return ma.masked_where(a2num==0.0, a2sum)/a2num
  
