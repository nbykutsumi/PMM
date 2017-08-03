from numpy import *
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import Image

clVer = "MyWNP.M.3"
lcl_tmp = [99,1,3,4]
llndsea = ["sea","lnd"]


rootDir = "/home/utsumi/mnt/wellshare"
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
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.My%s"%(ver)
  ibaseDirCL = rootDir + "/CLOUDTYPE/MyWNP%s"%(ver)

dclShortName = cl.dclShortName
figDir       = ibaseDir + "/pict"

for lndsea in llndsea:
  da2dat  = {}
  for i,icl in enumerate(lcl_tmp):
    figPath   = figDir + "/lines.mulProd.vsKuPR.%s.%s.png"%(lndsea,dclShortName[icl]) 
    iimg      = Image.open(figPath)
    da2dat[i] = asarray(iimg)

  a2line1 = hstack([da2dat[0], da2dat[1]])
  a2line2 = hstack([da2dat[2], da2dat[3]])
  a2oarray= vstack([a2line1, a2line2])
  oimg    = Image.fromarray(a2oarray)

  oPath   = figDir + "/join.lines.mulProd.vsKuPR.%s.png"%(lndsea)
  oimg.save(oPath)
  print oPath

