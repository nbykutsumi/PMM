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
ey   = -65  # bottom
ix   = 0
ex   = -1

llndsea = ["lnd","sea"]
lbinMin = [0, 0.1, 0.5]

for binMin in lbinMin:
  for lndsea in llndsea:
    da2dat = {}
    for i,clid in enumerate([99,1,3,4]):
      figDir  = os.path.join(ibaseDir, "pict")
      figPath = os.path.join(figDir, "TS.rat.min%.1f.%s.%s.png"%(binMin, lndsea, dclShortName[clid]))
      iimg    = Image.open(figPath)
      a2array = asarray(iimg)
      if clid == 4:
        da2dat[clid] = a2array[iy:-1,ix:ex]
      else:
        da2dat[clid] = a2array[iy:ey,ix:ex]
    
    a2line0 = vstack([da2dat[99],da2dat[1],da2dat[3],da2dat[4]])
    a2oarray= a2line0
    oimg    = Image.fromarray(a2oarray)
    
    oPath   = os.path.join(figDir, "join.TS.rat.min%.1f.%s.png"%(binMin, lndsea))
    oimg.save(oPath)
    print oPath
