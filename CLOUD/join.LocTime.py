import sys
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
from myfunc import util

icl     = 1

TrackFlag = True

#clVer = "MyWNP3"
clVer = "MyWNP.M.3"

iYM    = [2014,4]
eYM    = [2015,6]
#iYM    = [2015,6]
#eYM    = [2015,6]

lYM    = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [11,12,1,2,3]]

#ldattype = ["RA","KuPR","GMI","GSMaP","GSMaP.IR","GSMaP.MW","IMERG","IMERG.IR","IMERG.MW"]
#ldattype = ["KuPR","GMI","GSMaP","GSMaP.IR","GSMaP.MW","IMERG","IMERG.IR","IMERG.MW"]
dattype = "IMERG.IR"

#rootDir = "/tank/utsumi"
rootDir = "/home/utsumi/mnt/well.share"
if clVer   == "JMA1":
  cl         = CLOUDTYPE.CloudWNP()
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.JMA"
  ibaseDirCL = "/tank/utsumi/CLOUDTYPE/WNPAC"

elif clVer[:5] == "MyWNP":
  ver        = clVer[5:]
  cl         = CLOUDTYPE.MyCloudWNP(ver=ver)
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.My%s"%(ver)
  ibaseDirCL = rootDir + "/CLOUDTYPE/%s"%(clVer)

dclName = cl.dclName
dclShortName = cl.dclShortName
dclid   = cl.dclid

baseDir = ibaseDir
sout = ""
for i,[Year,Mon] in enumerate(lYM):
  if TrackFlag ==True:
    sDir  = baseDir + "/Tr.VsRA.CL.%s/%04d"%(dattype,Year)
  elif TrackFlag ==False:
    sDir  = baseDir + "/VsRA.CL.%s/%04d"%(dattype,Year)
  iPath   = sDir + "/Check.%s.%04d.%02d.%s.csv"%(dattype,Year,Mon,dclShortName[icl])
  f = open(iPath,"r")
  stmp = f.read()
  f.close()
  if i == 0:
    sout = stmp
  else:
    sout = sout +  "\n".join(stmp.split("\n")[1:])

oDir   = sDir[:-4] + "/csv"
oPath  = oDir + "/Check.%s.%s.csv"%(dattype,dclShortName[icl])
util.mk_dir(oDir)
f = open(oPath, "w"); f.write(sout); f.close()
print oPath
