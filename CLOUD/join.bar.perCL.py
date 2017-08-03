from numpy import *
import Image
import myfunc.IO.CLOUDTYPE as CLOUDTYPE

#clVer = "JMA1"
#clVer = "MyWNP1"
#clVer = "MyWNP2"
clVer = "MyWNP.M.3"

rootDir = "/home/utsumi/mnt/well.share"
if clVer   == "JMA1":
  cl         = CLOUDTYPE.CloudWNP()
  ncltype = 8
  lcltype = range(ncltype)
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.JMA"
  ibaseDirCL = "/tank/utsumi/CLOUDTYPE/WNPAC"

elif clVer[:5] == "MyWNP":
  ver     = clVer[5:]
  cl      = CLOUDTYPE.MyCloudWNP(ver=ver)
  ncltype = cl.ncl
  lcltype = cl.licl
  dclShortName = cl.dclShortName
  ibaseDir     = rootDir + "/PMM/WNP.261x265/CL.My%s"%(ver)
  ibaseDirCL   = rootDir + "/CLOUDTYPE/MyWNP%s"%(ver)


figDir  = ibaseDir + "/pict"

lvartype = ["acc","prate","bias"]
#lvartype = ["acc","prate"]
llndsea  = ["any","lnd","sea","cst"]

lcltype = lcltype[1:] + [99]

iy  = 0  # top
ey  = -1  # bottom
ix  = 0   # left
ex  = -1 # right

for vartype in lvartype:
  for lndsea in llndsea:
    da2dat = {}
    for icl in lcltype:
      clName = dclShortName[icl]
      figPath   = figDir  + "/RAdom.bar.%s.perCL.%s.%s.png"%(vartype,lndsea,clName)
      iimg      = Image.open(figPath)
      a2array   = asarray(iimg)
      print shape(a2array)
      da2dat[icl] = a2array[iy:ey, ix:ex]
  
    a2line1  = vstack([da2dat[icl] for icl in lcltype])
    a2oarray = a2line1
    oimg     = Image.fromarray(a2oarray)
  
    oPath    = figDir + "/join.bar.%s.%s.perCL.png"%(vartype,lndsea)
    oimg.save(oPath)
    print oPath
  
  
