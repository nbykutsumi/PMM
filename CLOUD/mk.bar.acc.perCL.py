from numpy import *
from datetime import datetime, timedelta
from bisect import bisect, bisect_left, bisect_right
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import myfunc.util as util
import matplotlib.pyplot as plt
import calendar
import sys

#clVer  = "JMA1"
#clVer  = "MyWNP1"
clVer  = "MyWNP2"

iYM    = [2014,4]
eYM    = [2015,6]
#eYM    = [2014,5]
lYM    = util.ret_lYM(iYM, eYM)
lYM    = [YM for YM in lYM if YM[1] not in [11,12,1,2,3]]
print lYM

ldattype = ["RA","KuPR","GMI","IMERG","IMERG.IR","IMERG.MW","GSMaP","GSMaP.IR","GSMaP.MW"]

#
#rootDir = "/tank/utsumi"
rootDir = "/home/utsumi/mnt/well.share"
if clVer   == "JMA1":
  cl         = CLOUDTYPE.CloudWNP()
  ncltype = 8
  lcltype = range(ncltype)
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.JMA"

elif clVer[:5] == "MyWNP":
  ver        = int(clVer[5:])
  cl         = CLOUDTYPE.MyCloudWNP(ver=ver)
  ncltype = cl.ncl
  lcltype = cl.licl
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.My%d"%(ver)
  ibaseDirCL = rootDir + "/CLOUDTYPE/MyWNP%d"%(ver)

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

drat = {}
dacc = {}
for dattype in ldattype:
  a3sum = accDat(dattype,lYM,999,"sum")
  a3num = accDat(dattype,lYM,999,"num")
  a3rat = ma.masked_invalid(a3sum / a3num)  # average rate (mm/h)

  a3acc = a3rat * a3numCL / len(lYM)  # mm/month
  a3acc = ma.masked_where(a3dommask==miss, a3acc)

  for icl in range(ncltype): drat[dattype,icl]  = a3rat[icl].mean()
  for icl in range(ncltype): dacc[dattype,icl]  = a3acc[icl].mean()

  #drat[dattype,99]  = ma.masked_invalid(a3sum.sum(axis=0)/a3num.sum(axis=0)).mean()
  drat[dattype,99]  = ma.masked_invalid(a3acc.sum(axis=0) / (a3numCL / len(lYM)).sum(axis=0)).mean()

  dacc[dattype,99]  = a3acc.sum(axis=0).mean()

#**** Figure: precipitation rate *****
for icl in range(ncltype) + [99]:
  figplot = plt.figure(figsize=(4.1,1.0))
  axplot  = figplot.add_axes([0.11,0.2, 0.84, 0.58])

  X  = arange(len(ldattype))+0.5
  Y  = [drat[dattype, icl]
              for dattype in ldattype]

  bars = plt.bar(X, Y, color=["k"]+["grey"]*(len(X)-1))

  # Title
  plt.title(dclName[icl], fontsize=12)

  # Y limit
  if    clVer == "MyWNP1":
    if icl == 1:
      axplot.set_ylim(0,6)
    else:
      axplot.set_ylim(0,0.5)

  elif  clVer == "MyWNP2":
    if icl in [1,2]:
      axplot.set_ylim(0,10)
    else:
      axplot.set_ylim(0,0.6)

  # X-ticks
  plt.xticks(X+0.5, map(int, X-0.5), fontsize=10)

  # Y-ticks
  for itick, tick in enumerate(axplot.yaxis.get_major_ticks()):
    if itick%2 ==0:
      tick.label.set_fontsize(10) 
    else:
      tick.label.set_fontsize(0) 


  # zero line
  plt.plot([X[0]-0.1,X[-1]+1.1],[0,0],"-",color="k")

  # Names
  sDir  = ibaseDir + "/pict"
  if dommask==False:  
    sPath = sDir + "/bar.prate.perCL.%s.png"%(dclShortName[icl])
  else:
    sPath = sDir + "/%sdom.bar.prate.perCL.%s.png"%(dommask, dclShortName[icl])
  figplot.savefig(sPath)
  print sPath 



#**** Figure: Accumulated rate *****
for icl in range(ncltype) + [99]:
  figplot = plt.figure(figsize=(4.1,1.0))
  axplot  = figplot.add_axes([0.11,0.2, 0.84, 0.56])

  X  = arange(len(ldattype))+0.5
  Y  = [dacc[dattype, icl]
              for dattype in ldattype]

  plt.bar(X, Y, color=["k"]+["grey"]*(len(X)-1))

  # Title
  plt.title(dclName[icl], fontsize=12)

  # Y limit
  if    clVer[:5] == "MyWNP":
    if icl == 99:
      axplot.set_ylim(0,300)
    else:
      axplot.set_ylim(0,180)

  # X-ticks
  plt.xticks(X+0.5, map(int, X-0.5), fontsize=10)

  # Y-ticks
  for itick, tick in enumerate(axplot.yaxis.get_major_ticks()):
    if itick%2 ==0:
      tick.label.set_fontsize(10) 
    else:
      tick.label.set_fontsize(0) 


  # Names
  sDir  = ibaseDir + "/pict"
  if dommask==False:  
    sPath = sDir + "/bar.acc.perCL.%s.png"%(dclShortName[icl])
  else:
    sPath = sDir + "/%sdom.bar.acc.perCL.%s.png"%(dommask, dclShortName[icl])
  figplot.savefig(sPath)
  print sPath 

#**** Legend file ****
legPath = sDir + "/legend.bar.perCL.png"
sout    = "\n".join(["%d:  %s"%(i, ldattype[i]) for i in range(len(ldattype))])
figleg  = plt.figure(figsize=(2,3))
figleg.text(0.1,0.1,sout, fontsize=16)
figleg.savefig(legPath)



