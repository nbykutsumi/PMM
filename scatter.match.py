from numpy import *
import matplotlib.pyplot as plt
import ioPMM
import myfunc.util as util
import os

iYear = 2014
eYear = 2014
lYear = range(iYear,eYear+1)
lMon  = [4,5,6]

ltagAll  = ["c","fbc","tc","ot"]
prtype1 = "RA"
prtype2 = "GPM.KuPR"
Trc     = "ALL"
IO    = ioPMM.ioRAdom()
ny,nx = IO.ny, IO.nx


baseDir = "/tank/utsumi/PMM/RAdom/%s.vs.%s.on.%s.%sclass"%(prtype1, prtype2, Trc, len(ltagAll))
def ret_monSumPath(prtype, tag, Year,Mon):
  sDir     = os.path.join(baseDir, "%04d"%(Year))
  sumPath  = os.path.join(sDir, "sum.%s.%s.%04d%02d.%dx%d"%(prtype,tag,Year,Mon,ny,nx))
  return sDir, sumPath

def ret_monNumPath(tag, Year,Mon):
  sDir     = os.path.join(baseDir, "%04d"%(Year))
  numPath  = os.path.join(sDir, "num.%s.%04d%02d.%dx%d"%(tag,Year,Mon,ny,nx))
  return sDir, numPath

def ret_monSum(prtype, tag, Year,Mon):
  sDir     = ret_monSumPath(prtype,tag,Year,Mon)[0]
  sumPath  = ret_monSumPath(prtype,tag,Year,Mon)[1]
  return fromfile(sumPath, float32).reshape(ny,nx)

def ret_monNum(tag, Year,Mon):
  sDir     = ret_monNumPath(tag,Year,Mon)[0]
  numPath  = ret_monNumPath(tag,Year,Mon)[1] 

  return fromfile(numPath, float32).reshape(ny,nx)

def ret_monAve(prtype, tag, Year,Mon):
  sDir     = ret_monSumPath(prtype,tag,Year,Mon)[0]
  sumPath  = ret_monSumPath(prtype,tag,Year,Mon)[1]
  numPath  = ret_monNumPath(prtype,tag,Year,Mon)[1] 

  a2sum    = fromfile(sumPath, float32).reshape(ny,nx)
  a2num    = fromfile(numPath, float32).reshape(ny,nx)
  return ma.masked_where(a2num==0.0, a2sum)/a2num

#********* main ********************************
for tag in ltagAll+["plain"]:
  a2sum1 = zeros([ny,nx],float32)
  a2sum2 = zeros([ny,nx],float32)
  a2num  = zeros([ny,nx],float32)
  for Year,Mon in [[Year,Mon] for Year in lYear for Mon in lMon]:
    a2sum1 = a2sum1 + ret_monSum(prtype1, tag ,Year,Mon)
    a2sum2 = a2sum2 + ret_monSum(prtype2, tag ,Year,Mon)
    a2num  = a2num  + ret_monNum(tag, Year,Mon)
    print a2sum1.sum(), a2sum2.sum(), a2num.sum()

  a2ave1 = ma.masked_where(a2num==0.0, a2sum1)/a2num
  a2ave2 = ma.masked_where(a2num==0.0, a2sum2)/a2num

  #****** Figure ***************
  figplot  = plt.figure(figsize=(5,5))
  axplot   = figplot.add_axes([0.13, 0.13, 0.77, 0.77])
  axplot.plot( a2ave1, a2ave2, "o", color="k", markersize=2)

  #-- axis -----
  axmax = 15.
  axplot.set_xlim(0.0, axmax)
  axplot.set_ylim(0.0, axmax)
 
  #-- axis label --
  axplot.set_xlabel("%s [mm/hour]"%(prtype1), fontsize=18)
  axplot.set_ylabel("%s [mm/hour]"%(prtype2), fontsize=18)

  #-- title -------
  Corr   = corrcoef(a2ave1.flatten(), a2ave2.flatten())[0][1]
  stitle = "X=%s Y=%s tag:%s"%(prtype1, prtype2,tag)
  stitle = stitle + "\n%04d-%04d M%02d-%02d R=%5.2f"%(iYear,eYear, lMon[0], lMon[-1], Corr)
  axplot.set_title(stitle, fontsize=15)

  #-- save -------
  figDir  = baseDir + "/pict"
  figPath = figDir  + "/scatter.%04d-%04d.M%02d-%02d.%s.png"%(iYear,eYear,lMon[0],lMon[-1],tag)
  util.mk_dir(figDir)
  plt.savefig(figPath)
  print figPath

