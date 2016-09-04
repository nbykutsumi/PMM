from numpy import *
from datetime import datetime, timedelta
from bisect import bisect, bisect_left, bisect_right
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import myfunc.util as util
import matplotlib.pyplot as plt
import calendar
import sys

clVer  = "JMA1"
clVer  = "MyWNP1"

iYM    = [2014,4]
eYM    = [2015,6]
lYM    = util.ret_lYM(iYM, eYM)
lYM    = [YM for YM in lYM if YM[1] not in [11,12,1,2,3]]
print lYM

#ldattype = ["RA","KuPR","GMI","IMERG","IMERG.IR","IMERG.MW","GSMaP","GSMaP.IR","GSMaP.MW"]
ldattype = ["RA"]

#
#rootDir = "/tank/utsumi"
rootDir = "/home/utsumi/mnt/well.share"
if clVer   == "JMA1":
  cl         = CLOUDTYPE.CloudWNP()
  ncltype = 8
  lcltype = range(ncltype)
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.JMA"

elif clVer == "MyWNP1":
  cl         = CLOUDTYPE.MyCloudWNP(ver=1)
  ncltype = 5
  lcltype = range(ncltype)
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.My1"

#
ny,nx = cl.ny, cl.nx   #261, 265
Lat   = cl.Lat
Lon   = cl.Lon

dommask  = "RA"   # Keep!
#dommask  = None

#BBox    = [[-0.1, 113.875],[52.1, 180.125]]  # WN.Pac
BBox    = [[20., 118.],[48., 150.]]    # RadarAMeDAS
miss  = -9999.

dclName = cl.dclName
dclShortName = cl.dclShortName

#----------------------
iX      = bisect_left (Lon,BBox[0][1])
eX      = bisect_right(Lon,BBox[1][1])
iY      = bisect_left (Lat,BBox[0][0])
eY      = bisect_left (Lat,BBox[1][0])
a2BBox  = ones([ny,nx],float32)*miss
a2BBox[iY:eY+1, iX:eX+1]= 1.0

if dommask =="RA":
  maskPath  = "/tank/utsumi/data/RadarAMeDAS/mask/RAmask.kubota.0.20x0.25WNP.261x265"
  a2dommask = fromfile(maskPath,float32).reshape(ny,nx)
  a2dommask = ma.masked_where(a2BBox ==miss, a2dommask)
  a3dommask = array([a2dommask for i in range(len(lcltype))])
else:
  a2dommask = a2BBox
  a3dommask = array([a2dommask for i in range(len(lcltype))])

#******************************************
def loadSum(dattype,Year,Mon,binPr):
  srcDir    = ibaseDir + "/CL.Pr.%s"%(dattype)
  iPath     = srcDir + "/sum.P%05.1f.%04d.%02d.%dx%dx%d"%(binPr, Year,Mon, ncltype,ny,nx)
  a3out = fromfile(iPath, float32).reshape(ncltype, ny, nx)

  return ma.masked_where(a3dommask==-9999., a3out)

def loadNum(dattype,Year,Mon,binPr):
  srcDir    = ibaseDir + "/CL.Pr.%s"%(dattype)
  iPath     = srcDir + "/num.P%05.1f.%04d.%02d.%dx%dx%d"%(binPr, Year,Mon, ncltype,ny,nx)
  a3out = fromfile(iPath, int32).reshape(ncltype, ny, nx)
  return ma.masked_where(a3dommask==-9999., a3out)


def accDat(dattype,lYM, binPr, sumnum="num"):
  dfunc = {"sum":loadSum
          ,"num":loadNum
          }
  ddtype= {"sum":"float32"
          ,"num":"int32"
          }

  accDat = zeros([ncltype, ny, nx], ddtype[sumnum])

  for (Year,Mon) in lYM:
    accDat = accDat + dfunc[sumnum](dattype,Year,Mon,binPr)
  return accDat

#******************************************
iDTime   = datetime(iYM[0],iYM[1],1,0)
eDTime   = datetime(eYM[0],eYM[1], calendar.monthrange(eYM[0],eYM[1])[1],23)
dDTime   = timedelta(hours=1)

a3numCL  = array([cl.loadNumAcc(iYM, eYM, cltype) for cltype in lcltype])

for dattype in ldattype:
  a3sum = accDat(dattype,lYM,999,"sum")
  a3num = accDat(dattype,lYM,999,"num")
  a3rat = ma.masked_invalid(a3sum / a3num)  # average rate (mm/h)

  a3acc = a3rat * a3numCL / len(lYM)  # mm/month
  a3acc = ma.masked_where(a3dommask==miss, a3acc)

  print dattype, a3acc.mean()

  #**** Figure: precipitation rate *****
  figplot = plt.figure(figsize=(4.1,3.2))
  axplot  = figplot.add_axes([0.11,0.10, 0.88, 0.8])

  # Each cloud type
  X  = [icl+1 for (icl,cltype) in enumerate(lcltype)]

#  """
#  CAUTION!: precip rate for Cb is multipled by 0.1
#  """
#  Y  = [a3rat[icl].mean()*0.1 if icl==1 else a3rat[icl].mean()
#              for (icl,cltype) in enumerate(lcltype)]

  Y  = [a3rat[icl].mean() 
              for (icl,cltype) in enumerate(lcltype)]

  # All
  X = [0] + X
  Y = [ma.masked_invalid(a3sum.sum(axis=0)/a3num.sum(axis=0)).mean()] + Y

  # Title
  plt.bar(X, Y, color=["k"]+["grey"]*(len(X)-1))
  plt.title(dattype, fontsize=20)

  # X-ticks
  plt.xticks(array(X)+0.5, [dclShortName[icl] for icl in [99]+lcltype]
             ,fontsize=18)

  # Y-ticks
  for tick in axplot.yaxis.get_major_ticks():
    tick.label.set_fontsize(15) 

  # Y limit
  #axplot.set_ylim(0,1.4)

  # zero line
  plt.plot([X[0]-0.1,X[-1]+1.1],[0,0],"-",color="k")

  # Names
  #sDir  = "/tank/utsumi/PMM/WNP.261x265/pict"
  #sDir  = "/home/utsumi/mnt/well.share/PMM/WNP.261x265/pict"
  sDir  = ibaseDir + "/pict"
  if dommask==False:  
    sPath = sDir + "/bar.prate.%s.png"%(dattype)
  else:
    sPath = sDir + "/%sdom.bar.prate.%s.png"%(dommask, dattype)
  figplot.savefig(sPath)
  print "all=",Y[0]
  print sPath 


  #**** Figure: Accumulated precipitation *****
  figplot = plt.figure(figsize=(4.1,3.2))
  axplot  = figplot.add_axes([0.11,0.10, 0.88, 0.8])

  # Each cloud type
  X  = [icl+1 for (icl,cltype) in enumerate(lcltype)]
  Y  = [a3acc[icl].mean() 
              for (icl,cltype) in enumerate(lcltype)]

  # All
  X = [0] + X
  Y = [a3acc.sum(axis=0).mean()] + Y

  # Title
  plt.bar(X, Y, color=["k"]+["grey"]*(len(X)-1))
  plt.title(dattype, fontsize=20)

  # X-ticks
  plt.xticks(array(X)+0.5, [dclShortName[icl] for icl in [99]+lcltype]
             ,fontsize=18)

  # Y-ticks
  for tick in axplot.yaxis.get_major_ticks():
    tick.label.set_fontsize(15) 

  # Y limit
  axplot.set_ylim(0,200)

  # zero line
  plt.plot([X[0]-0.1,X[-1]+1.1],[0,0],"-",color="k")

  # Names
  #sDir  = "/tank/utsumi/PMM/WNP.261x265/pict"
  #sDir  = "/home/utsumi/mnt/well.share/PMM/WNP.261x265/pict"
  sDir  = ibaseDir + "/pict"
  if dommask==False:  
    sPath = sDir + "/bar.acc.%s.png"%(dattype)
  else:
    sPath = sDir + "/%sdom.bar.acc.%s.png"%(dommask, dattype)
  figplot.savefig(sPath)
  print "all=",Y[0]
  print sPath 

