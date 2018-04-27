from numpy import *
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import Image

#clVer = "JMA1"
#clVer = "MyWNP1"
#clVer = "MyWNP2"
clVer = "MyWNP.M.3"

logcont=True

#ldattype = ["KuPR","GMI","IMERG","IMERG.IR","IMERG.MW","GSMaP","GSMaP.IR","GSMaP.MW"]
#ldattype = ["GMI","IMERG","GSMaP.MW","IMERG.MW","GSMaP.IR","IMERG.IR"]
ldattype = ["GMI","GSMaP.MW","GSMaP.IR","IMERG.MW","IMERG.IR"]

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

#iy  = 125  # top
iy  = 1  # top
ey  = -1  # bottom
ix  = 1   # left
ex  = -50 # right

for icl in [99]+[1,3,4]:
  da2dat = {}
  for lndsea in llndsea:
    for i,dattype in enumerate(ldattype):
      if logcont==False:
        figPath   = figDir  + "/scatter.%s.vsKuPR.%s.%s.png"%(dattype,lndsea,dclShortName[icl])
      elif logcont==True:
        figPath   = figDir  + "/scatter.logcont.%s.vsKuPR.%s.%s.png"%(dattype,lndsea,dclShortName[icl])
      else:
        print "check logcont",logcont
        sys.exit()

      iimg      = Image.open(figPath)
      a2array   = asarray(iimg)
      print shape(a2array)
      da2dat[lndsea, i] = a2array[iy:ey, ix:ex]

    da2dat[lndsea, -9999]=da2dat[lndsea,0]*0+255
  
  a2line1  = hstack([da2dat["sea",0],     da2dat["sea",1], da2dat["sea",2]])
  a2line2  = hstack([da2dat["sea",-9999],     da2dat["sea",3], da2dat["sea",4]])
  a2line3  = hstack([da2dat["lnd",0],     da2dat["lnd",1], da2dat["lnd",2]])
  a2line4  = hstack([da2dat["lnd",-9999],     da2dat["lnd",3], da2dat["lnd",4]])

  a2oarray = vstack([a2line1, a2line2, a2line3, a2line4])
  oimg     = Image.fromarray(a2oarray)
 
  if logcont==False: 
    oPath    = figDir + "/join.scatter.vsKuPR.%s.png"%(dclShortName[icl])
  elif logcont==True:
    oPath    = figDir + "/join.scatter.logcont.vsKuPR.%s.png"%(dclShortName[icl])
  oimg.save(oPath)
  print oPath
  
  
