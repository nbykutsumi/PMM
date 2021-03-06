from numpy import *
from datetime import datetime, timedelta
from bisect import bisect, bisect_left, bisect_right
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import myfunc.util as util
import matplotlib.pyplot as plt
import calendar
import sys

#TrackFlag = True
TrackFlag = False

#clVer  = "JMA1"
#clVer  = "MyWNP1"
#clVer  = "MyWNP2"
clVer  = "MyWNP.M.3"

iYM    = [2014,4]
eYM    = [2015,6]
#eYM    = [2014,4]
lYM    = util.ret_lYM(iYM, eYM)
lYM    = [YM for YM in lYM if YM[1] not in [11,12,1,2,3]]
print lYM

ldattype = ["RA","KuPR","GMI","IMERG","IMERG.IR","IMERG.MW","GSMaP","GSMaP.IR","GSMaP.MW"]
#ldattype = ["GSMaP.IR"]

lbin = [0,3,6,9,12,16,20,40]
lBin = map(list,zip(lbin[:-1],lbin[1:])) + [[lbin[-1],999]]

rootDir = "/home/utsumi/mnt/well.share"
if clVer   == "JMA1":
  cl         = CLOUDTYPE.CloudWNP()
  ncltype = 8
  lcltype = range(ncltype)
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.JMA"

elif clVer[:5] == "MyWNP":
  ver        = clVer[5:]
  cl         = CLOUDTYPE.MyCloudWNP(ver=ver)
  ncltype = cl.ncl
  lcltype = cl.licl
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.My%s"%(ver)
  ibaseDirCL = rootDir + "/CLOUDTYPE/MyWNP%s"%(ver)

#
ny,nx = cl.ny, cl.nx   #261, 265
Lat   = cl.Lat
Lon   = cl.Lon

dommask  = "RA"  # Do NOT Change!!

#BBox    = [[-0.1, 113.875],[52.1, 180.125]]  # WN.Pac
BBox    = [[20., 118.],[48., 150.]]    # RadarAMeDAS
miss  = -9999.

dclName      = cl.dclName
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


## land sea mask
#landfracDir = "/home/utsumi/mnt/well.share/PMM/WNP.261x265/MASK"
#landfracPath= landfracDir + "/landfrac.%dx%d"%(cl.ny, cl.nx)
#a2lndfrc    = fromfile(landfracPath, float32).reshape(261,265)
#
#a2lndmask   = ma.masked_less(a2lndfrc,1.0).mask  # mask sea&coast
#a2seamask   = ma.masked_greater(a2lndfrc,0.0).mask # mask land&coast
#a2cstmask   = ma.masked_equal(a2lndfrc, 1.0).mask \
#             +ma.masked_equal(a2lndfrc, 0.0).mask  # mask land&sea
#a2anymask   = ma.masked_equal(zeros([ny, nx]), 1.0).mask
#da2lsmask   = {"lnd":a2lndmask, "sea":a2seamask, "cst":a2cstmask,"any":a2anymask}
#
#llndsea     = ["lnd","sea","cst","any"]
#******************************************
def loadSum(dattype,Year,Mon,binPr):
  srcDir    = ibaseDir + "/CL.Pr.%s/%04d"%(dattype,Year)
  iPath     = srcDir + "/sum.P%05.1f.%04d.%02d.%dx%dx%d"%(binPr, Year,Mon, ncltype,ny,nx)
  a3out = fromfile(iPath, float32).reshape(ncltype, ny, nx)

  return ma.masked_where(a3dommask==-9999., a3out)

def loadNum(dattype,Year,Mon,binPr):
  srcDir    = ibaseDir + "/CL.Pr.%s/%04d"%(dattype,Year)
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
for dattype in ldattype:
  Val = empty([ncltype, len(lbin)], float32)
  for iBin, Bin in enumerate(lBin):
    binMin, binMax = Bin[0], Bin[1]
    a3numMin = accDat(dattype,lYM,binMin,"num")
    a3numMax = accDat(dattype,lYM,binMax,"num")
 
    a3num    = a3numMax - a3numMin 
    vTot     = a3num.sum() 

    for icl in range(ncltype):
      Val[icl,iBin] = a3num[icl].sum() / float(vTot)

  Bot  = empty([ncltype, len(lbin)], float32)
  Bot[0] = 0.0
  Bot[1:]= Val[:-1].cumsum(axis=0)
  #****** Figure ******************
  figplot = plt.figure(figsize=(3.0,1.0))
  axplot  = figplot.add_axes([0.11,0.20, 0.84, 0.55])

  Col= plt.matplotlib.cm.jet(linspace(0,1, ncltype))  
  X  = range(len(lbin))
  for icl in range(ncltype):
    axplot.bar(X, Val[icl], bottom=Bot[icl], color=Col[icl])

  print Val[1]

  # Title
  plt.title(dclName[icl], fontsize=12)

  # Add title
  stitle = "%04d/%02d-%04d/%02d %s"%(iYM[0],iYM[1],eYM[0],eYM[1],dattype)
  plt.title(stitle, fontsize=12)

  # X-ticks
  plt.xticks(X, map(str, lbin[:-1]) + ["%d<"%(lbin[-1])])


  # Y-ticks
  for itick, tick in enumerate(axplot.yaxis.get_major_ticks()):
    if itick%2 == 0: 
      tick.label.set_fontsize(10) 
    else:
      if icl == 99:
        tick.label.set_fontsize(10)
      else:
        tick.label.set_fontsize(0)

  # Y-limit
  axplot.set_ylim(0,1.0)

  # Names
  sDir  = ibaseDir + "/pict"
  if dommask==False:
    sPath = sDir + "/bar.ratioOfCL.%s.png"%(dattype)
  else:
    sPath = sDir + "/%sdom.bar.ratioOfCL.%s.png"%(dommask, dattype)
  figplot.savefig(sPath)
  print sPath 
  plt.close()

# Figure Legend
figplot = plt.figure(figsize=(3,2))
axplot  = figplot.add_axes([0,0,1,1])
for icl in range(ncltype):
  axplot.plot(0.5,     1+icl, "s", markerfacecolor=Col[icl], markersize=15)
  axplot.text(0.5+0.2, 1+icl,  dclName[icl], va="center",fontsize=15)

axplot.set_xlim(0,5)
axplot.set_ylim(0,ncltype+2)
legendPath = sDir + "/legend.ratioOfCL.png"
plt.savefig(legendPath)
print legendPath
