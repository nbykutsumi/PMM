from numpy import *
import os,sys
import Image
import myfunc.IO.CLOUDTYPE as CLOUDTYPE


#-- JMA-Cloud ---
clVer = "MyWNP.M.3"
rootDir = "/home/utsumi/mnt/well.share"
if clVer[:5] == "MyWNP":
  ver        = clVer[5:]
  cl         = CLOUDTYPE.MyCloudWNP(ver=ver)
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.My%s"%(ver)
  ibaseDirCL = rootDir + "/CLOUDTYPE/%s"%(clVer)
#-- CloudNames --
dclName      = cl.dclName
dclShortName = cl.dclShortName
lclid       = sort(dclName.keys())

iy   = 0   #top
ey   = -1  # bottom
ix   = 0
ex   = -50
da2dat = {}

for i,clid in enumerate(lclid[:-1]):
  figDir  = os.path.join(ibaseDir, "pict")
  figPath = os.path.join(figDir, "prof.CL.%s.png"%(dclShortName[clid]))
  iimg    = Image.open(figPath)
  a2array = asarray(iimg)
  da2dat[i] = a2array[iy:ey,ix:ex]

a2line0 = vstack([da2dat[0],da2dat[1],da2dat[2]])
a2line1 = vstack([da2dat[3],da2dat[4],da2dat[5]])
a2oarray= hstack([a2line0,a2line1])
oimg    = Image.fromarray(a2oarray)

oPath   = os.path.join(figDir, "join.prof.CL.png")
oimg.save(oPath)
print oPath
