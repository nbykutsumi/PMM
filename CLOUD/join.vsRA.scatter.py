from numpy import *
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import Image

#clVer = "JMA1"
#clVer = "MyWNP1"
#clVer = "MyWNP2"
clVer = "MyWNP.M.3"

#ldattype = ["KuPR","GMI","IMERG","IMERG.IR","IMERG.MW","GSMaP","GSMaP.IR","GSMaP.MW"]
ldattype = ["KuPR","IMERG","IMERG.IR","IMERG.MW","GMI","GSMaP","GSMaP.IR","GSMaP.MW"]

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
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.My%s"%(ver)
  ibaseDirCL = rootDir + "/CLOUDTYPE/MyWNP%s"%(ver)

dclShortName = cl.dclShortName
figDir       = ibaseDir + "/pict"

#iy  = 125  # top
iy  = 108  # top
ey  = -1  # bottom
ix  = 1   # left
ex  = -118 # right

for icl in lcltype + [99]:

  da2dat = {}
  for i,dattype in enumerate(ldattype):
    figPath   = figDir  + "/scatter.%s.vsRA.%s.png"%(dattype,dclShortName[icl])
    iimg      = Image.open(figPath)
    a2array   = asarray(iimg)
    print shape(a2array)
    da2dat[i] = a2array[iy:ey, ix:ex]

  a2line1  = hstack([da2dat[0], da2dat[1], da2dat[2], da2dat[3]])
  a2line2  = hstack([da2dat[4], da2dat[5], da2dat[6], da2dat[7]])
  a2oarray = vstack([a2line1, a2line2])
  oimg     = Image.fromarray(a2oarray)

  oPath    = figDir + "/join.pdf.%s.png"%(dclShortName[icl])
  oimg.save(oPath)
  print oPath


