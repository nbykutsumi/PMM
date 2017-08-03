from numpy import *
import Image
import myfunc.IO.CLOUDTYPE as CLOUDTYPE

#clVer = "JMA1"
#clVer = "MyWNP1"
clVer = "MyWNP2"


ldattype = ["KuPR","GMI","IMERG","IMERG.IR","IMERG.MW","GSMaP","GSMaP.IR","GSMaP.MW"]

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
  dclShortName = cl.dclShortName
  ibaseDir     = rootDir + "/PMM/WNP.261x265/CL.My%d"%(ver)
  ibaseDirCL   = rootDir + "/CLOUDTYPE/MyWNP%d"%(ver)


figDir  = ibaseDir + "/pict"

iy  = 0  # top
ey  = -1  # bottom
ix  = 0   # left
ex  = -1 # right


for dattype in ldattype:
  da2dat = {}
  for icltype in [1,3,4]:
    figPath   = figDir  + "/bar.vsRA.%s.%s.png"%(dattype, dclShortName[icltype])
    iimg      = Image.open(figPath)
    a2array   = asarray(iimg)
    print shape(a2array)
    da2dat[icltype] = a2array[iy:ey, ix:ex]
  
  a2line1  = vstack([da2dat[icltype] for icltype in [1,3,4]])
  a2oarray = a2line1
  oimg     = Image.fromarray(a2oarray)
  
  oPath    = figDir + "/join.vsRA.bar.%s.png"%(dattype)
  oimg.save(oPath)
  print oPath
  
  
