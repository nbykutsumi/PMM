from numpy import *
from bisect import bisect, bisect_left, bisect_right
import myfunc.util as util
import myfunc.IO.CLOUDTYPE as CLOUDTYPE

iYM   = [2014,4]
eYM   = [2014,4]
lYM   = util.ret_lYM(iYM, eYM)

ldattype= ["GSMaP","GSMaP.IR","GSMaP.MW","IMERG","IMERG.IR","IMERG.MW"]

lcltype = range(0,7+1)
ncltype = len(lcltype)
dclid   = {0:0, 1:1, 2:201, 3:202, 4:4, 5:3, 6:204, 7:200}

dclName ={0:"Clear Sky",   1:"Cumulonimbus(Cb)",  2:"High Cloud",3:"Mid Cloud"
         ,4:"Cumulus(Cu)", 5:"Stratocumulus(Sc)", 6:"Fog/St"    ,7:"Cloudy"}

dclShortName={0:"no", 1:"Cb",  2:"hi",3:"md"
             ,4:"Cu", 5:"Sc",  6:"St",7:"cw"}

cl    = CLOUDTYPE.CloudWNP()
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
  baseDir = "/tank/utsumi/PMM/WNP.261x265"
  srcDir    = baseDir + "/CL.Pr.%s"%(dattype)
  iPath     = srcDir + "/num.P%05.1f.%04d.%02d.%dx%dx%d"%(binPr, Year,Mon, ncltype,ny,nx)
  if binPr == 0.0:
    a3out = zeros([ncltype,ny,nx],int32)
  else:
    a3out = fromfile(iPath, int32).reshape(ncltype, ny, nx)

  return ma.masked_where(a3dommask==-9999., a3out)

def loadSum(Year,Mon,binPr):
  baseDir = "/tank/utsumi/PMM/WNP.261x265"
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




