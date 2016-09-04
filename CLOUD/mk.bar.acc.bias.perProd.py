from numpy import *
from datetime import datetime, timedelta
from bisect import bisect, bisect_left, bisect_right
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import myfunc.util as util
import matplotlib.pyplot as plt
import calendar
import sys

iYM    = [2014,4]
eYM    = [2014,11]
lYM    = util.ret_lYM(iYM, eYM)
ldattype = ["KuPR","GMI","IMERG","IMERG.IR","IMERG.MW","GSMaP","GSMaP.IR","GSMaP.MW"]
#ldattype = ["GSMaP.IR"]

cl    = CLOUDTYPE.CloudWNP()
ny,nx = cl.ny, cl.nx   #261, 265
Lat   = cl.Lat
Lon   = cl.Lon

dommask  = "RA"  # Do NOT Change!!

#BBox    = [[-0.1, 113.875],[52.1, 180.125]]  # WN.Pac
BBox    = [[20., 118.],[48., 150.]]    # RadarAMeDAS
miss  = -9999.

lcltype = range(0,7+1)
ncltype = len(lcltype)
dclName ={0:"Clear Sky",   1:"Cumulonimbus(Cb)",  2:"High Cloud",3:"Mid Cloud"
         ,4:"Cumulus(Cu)", 5:"Stratocumulus(Sc)", 6:"Fog/St"    ,7:"Cloudy", 99:"All"}

dclShortName={0:"no", 1:"Cb",  2:"hi",3:"md"
             ,4:"Cu", 5:"Sc",  6:"St",7:"cw", 99:"All"}

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
  baseDir = "/tank/utsumi/PMM/WNP.261x265"
  srcDir    = baseDir + "/CL.Pr.%s"%(dattype)
  iPath     = srcDir + "/sum.P%05.1f.%04d.%02d.%dx%dx%d"%(binPr, Year,Mon, ncltype,ny,nx)
  a3out = fromfile(iPath, float32).reshape(ncltype, ny, nx)

  return ma.masked_where(a3dommask==-9999., a3out)

def loadNum(dattype,Year,Mon,binPr):
  baseDir = "/tank/utsumi/PMM/WNP.261x265"
  srcDir    = baseDir + "/CL.Pr.%s"%(dattype)
  iPath     = srcDir + "/num.P%05.1f.%04d.%02d.%dx%dx%d"%(binPr, Year,Mon, ncltype,ny,nx)
  a3out = fromfile(iPath, int32).reshape(ncltype, ny, nx)
  return ma.masked_where(a3dommask==-9999., a3out)


def accDat(dattype,iYM, eYM, binPr, sumnum="num"):
  dfunc = {"sum":loadSum
          ,"num":loadNum
          }
  ddtype= {"sum":"float32"
          ,"num":"int32"
          }

  accDat = zeros([ncltype, ny, nx], ddtype[sumnum])

  lYM = util.ret_lYM(iYM,eYM)
  for (Year,Mon) in lYM:
    accDat = accDat + dfunc[sumnum](dattype,Year,Mon,binPr)
  return accDat

#******************************************
iDTime   = datetime(iYM[0],iYM[1],1,0)
eDTime   = datetime(eYM[0],eYM[1], calendar.monthrange(eYM[0],eYM[1])[1],23)
dDTime   = timedelta(hours=1)

a3numCL  = array([cl.loadNumAcc(iYM, eYM, cltype) for cltype in lcltype])

a3sumRA = accDat("RA",iYM,eYM,999,"sum")
a3numRA = accDat("RA",iYM,eYM,999,"num")
a3ratRA = ma.masked_invalid(a3sumRA / a3numRA) # average rate (mm/h) including dry period
a3accRA = a3ratRA * a3numCL / len(lYM) # mm/month
a3accRA = ma.masked_where(a3dommask==-9999., a3accRA)

for dattype in ldattype:
  a3sum = accDat(dattype,iYM,eYM,999,"sum")
  a3num = accDat(dattype,iYM,eYM,999,"num")
  a3rat = ma.masked_invalid(a3sum / a3num)  # average rate (mm/h)
  a3acc = a3rat * a3numCL / len(lYM)  # mm/month
  a3acc = ma.masked_where(a3dommask==-9999., a3acc)

  figplot = plt.figure(figsize=(4.1,3.2))
  axplot  = figplot.add_axes([0.11,0.10, 0.88, 0.8])

  # Each cloud type
  X  = [icl+1 for (icl,cltype) in enumerate(lcltype)]
  Y  = [a3acc[icl].mean() - a3accRA[icl].mean()
              for (icl,cltype) in enumerate(lcltype)]

  # All
  X = [0] + X
  Y = [a3acc.sum(axis=0).mean() - a3accRA.sum(axis=0).mean()] + Y

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
  axplot.set_ylim(-60,60)

  # zero line
  plt.plot([X[0]-0.1,X[-1]+1.1],[0,0],"-",color="k")

  # Names
  #sDir  = "/tank/utsumi/PMM/WNP.261x265/pict"
  sDir  = "/home/utsumi/mnt/well.share/PMM/WNP.261x265/pict"
  if dommask==False:
    sPath = sDir + "/bar.bias.%s.png"%(dattype)
  else:
    sPath = sDir + "/%sdom.bar.bias.%s.png"%(dommask, dattype)
  figplot.savefig(sPath)
  print sPath 

