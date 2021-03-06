from numpy import *
import Image
import myfunc.IO.CLOUDTYPE as CLOUDTYPE

#clVer = "JMA1"
#clVer = "MyWNP1"
clVer = "MyWNP2"


ldattype = ["RA","KuPR","GMI","IMERG","IMERG.IR","IMERG.MW","GSMaP","GSMaP.IR","GSMaP.MW"]

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

da2dat = {}
for idattype, dattype in enumerate(ldattype):
  figPath   = figDir  + "/RAdom.bar.ratioOfCL.%s.png"%(dattype)
  iimg      = Image.open(figPath)
  a2array   = asarray(iimg)
  print shape(a2array)
  da2dat[idattype] = a2array[iy:ey, ix:ex]

a2line1  = vstack([da2dat[idattype] for idattype in range(len(ldattype))])
a2oarray = a2line1
oimg     = Image.fromarray(a2oarray)

oPath    = figDir + "/join.RAdom.bar.ratioOfCL.png"
oimg.save(oPath)
print oPath


