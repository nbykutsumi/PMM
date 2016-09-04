from numpy import *
from bisect import bisect, bisect_left, bisect_right
import myfunc.util as util
import myfunc.IO.CLOUDTYPE as CLOUDTYPE

clVer = "JMA1"
clVer = "MyWNP"

iYM   = [2014,4]
eYM   = [2014,4]
lYM   = util.ret_lYM(iYM, eYM)

ldattype= ["RA"]

if clVer   == "JMA1":
  cl         = CLOUDTYPE.CloudWNP()
  ncltype = 8
  lcltype = range(ncltype)
  ibaseDir   = "/tank/utsumi/PMM/WNP.261x265"
  ibaseDirCL = "/tank/utsumi/CLOUDTYPE/WNPAC"
  figDir     = "/home/utsumi/mnt/well.share/PMM/WNP.261x265/pict"
  def ret_CLPrDir(dattype):
    return ibaseDir + "/CL.Pr.%s"%(dattype)

elif clVer == "MyWNP":
  cl         = CLOUDTYPE.MyCloudWNP()
  ncltype = 6
  lcltype = range(ncltype)
  ibaseDir   = "/home/utsumi/mnt/well.share/PMM/WNP.261x265"
  ibaseDirCL = "/home/utsumi/mnt/well.share/PMM/WNP.261x265/MyCLTYPE"
  figDir     = "/home/utsumi/mnt/well.share/PMM/WNP.261x265/pict.MyCL"
  def ret_CLPrDir(dattype):
    return ibaseDir + "/MyCL.Pr.%s"%(dattype)

dclName = cl.dclName
dclShortName = cl.dclShortName

Lat   = cl.Lat
Lon   = cl.Lon
ny    = cl.ny
nx    = cl.nx

#dommask  = "RA"
dommask  = None

#BBox    = [[-0.1, 113.875],[52.1, 180.125]]  # WN.Pac
BBox    = [[20., 118.],[48., 150.]]    # RadarAMeDAS
miss  = -9999.
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

#----------------------
def loadNum(Year,Mon,binPr):
  srcDir    = ret_CLPrDir(dattype)
  iPath     = srcDir + "/num.P%05.1f.%04d.%02d.%dx%dx%d"%(binPr, Year,Mon, ncltype,ny,nx)
  if binPr == 0.0:
    a3out = zeros([ncltype,ny,nx],int32)
  else:
    a3out = fromfile(iPath, int32).reshape(ncltype, ny, nx)

  return ma.masked_where(a3dommask==-9999., a3out)

def loadSum(Year,Mon,binPr):
  srcDir    = ret_CLPrDir(dattype)
  srcDir    = baseDir + "/CL.Pr.%s"%(dattype)
  iPath     = srcDir + "/sum.P%05.1f.%04d.%02d.%dx%dx%d"%(binPr, Year,Mon, ncltype,ny,nx)
  if binPr == 0.0:
    a3out = zeros([ncltype,ny,nx],float32)
  else:
    a3out = fromfile(iPath, float32).reshape(ncltype, ny, nx)

  return ma.masked_where(a3dommask==-9999., a3out)

def accDat(iYM, eYM, binPr, sumnum="num"):
  dfunc = {"sum":loadSum
          ,"num":loadNum
          }
  ddtype= {"sum":"float32"
          ,"num":"int32"
          }

  accDat = zeros([ncltype, ny, nx], ddtype[sumnum])

  lYM = util.ret_lYM(iYM,eYM)
  for (Year,Mon) in lYM:
    accDat = accDat + dfunc[sumnum](Year,Mon,binPr)
  return accDat
#----------------------
da2num = {}
for dattype in ldattype:
  a3num = accDat(iYM, eYM, 999, sumnum="num")
  da2num[dattype] = a3num.sum(axis=0)
  print dattype, da2num[dattype].sum()




