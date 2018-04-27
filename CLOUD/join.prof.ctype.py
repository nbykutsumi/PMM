from numpy import *
import os,sys
import Image
import myfunc.IO.CLOUDTYPE as CLOUDTYPE


#-- JMA-Cloud ---
clVer = "MyWNP.M.3"
rootDir = "/home/utsumi/mnt/wellshare"
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
ex   = -1
da2dat = {}

for i,clid in enumerate(lclid[:-1]):
  figDir  = os.path.join(ibaseDir, "pict")
  #figPath = os.path.join(figDir, "prof.CL.%s.png"%(dclShortName[clid]))
  figPath = os.path.join(figDir, "prof.CL.2prob.%s.png"%(dclShortName[clid]))
  iimg    = Image.open(figPath)
  a2array = asarray(iimg)[iy:ey,ix:ex]
  a2top   = (ones(a2array.shape,dtype="uint8")*255)[:20]
  da2dat[i] = vstack([a2top, a2array])
  #da2dat[i] =  a2array
  print da2dat[i].shape, a2array.shape
  print "*"*50
  print a2array

a2line0 = vstack([da2dat[0],da2dat[1],da2dat[2]])
a2line1 = vstack([da2dat[3],da2dat[4],da2dat[5]])
a2oarray= hstack([a2line0,a2line1])
print a2oarray.shape
oimg    = Image.fromarray(a2oarray)

#oPath   = os.path.join(figDir, "join.prof.CL.png")
oPath   = os.path.join(figDir, "join.prof.CL.2prob.png")
oimg.save(oPath)
print oPath
