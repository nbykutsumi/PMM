from numpy import *
from bisect import bisect, bisect_left, bisect_right
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import myfunc.util as util
import calendar
import sys

iYM   = [2014,4]
eYM   = [2015,7]
lYM   = util.ret_lYM(iYM, eYM)
#lYM    = [YM for YM in lYM if YM[1] not in [11,12,1,2,3]]

#clVer = "JMA1"
#clVer = "MyWNP1"
clVer = "MyWNP2"

#rootDir = "/tank/utsumi"
rootDir = "/home/utsumi/mnt/well.share"
if clVer   == "JMA1":
  cl         = CLOUDTYPE.CloudWNP()
  ncltype = 8
  lcltype = range(ncltype)
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.JMA"
  ibaseDirCL = "/tank/utsumi/CLOUDTYPE/WNPAC"

elif clVer[:5] == "MyWNP":
  ver     = int(clVer[5:])
  cl      = CLOUDTYPE.MyCloudWNP(ver=ver)
  ncltype = cl.ncl
  lcltype = cl.licl
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.My%d"%(ver)
  ibaseDirCL = "/home/utsumi/mnt/well.share/CLOUDTYPE/MyWNP%d"%(ver)

Lat   = cl.Lat
Lon   = cl.Lon
ny    = cl.ny
nx    = cl.nx

dclName = cl.dclName
dclShortName = cl.dclShortName

#dclName ={0:"Clear Sky",   1:"Cumulonimbus(Cb)",  2:"High Cloud",3:"Mid Cloud"
#         ,4:"Cumulus(Cu)", 5:"Stratocumulus(Sc)", 6:"Fog/St"    ,7:"Cloudy"}
#
#dclShortName={0:"no", 1:"Cb",  2:"hi",3:"md"
#             ,4:"Cu", 5:"Sc",  6:"St",7:"cw"}

dommask  = "RA"
#dommask  = None

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



da2num = {}
for icl in lcltype:
  da2num[icl]  = zeros([ny,nx],int32)
  for (Year,Mon) in lYM:
    print Year,Mon
    sDir  = ibaseDirCL + "/num"
    sPath  = sDir + "/num.%04d%02d.%s.%dx%d"%(Year,Mon,dclShortName[icl],ny,nx)
    a2tmp  = fromfile(sPath, int32).reshape(ny,nx)
    da2num[icl] = da2num[icl] + a2tmp

da2num[-1] =  array([da2num[icl] for icl in lcltype]).sum(axis=0)

da2frac = { icl: da2num[icl].astype(float32)/da2num[-1].astype(float32) for icl in lcltype }

for icl in lcltype:
  v = ma.masked_where(a2dommask == miss, da2frac[icl]).mean()
  print dclShortName[icl], v

