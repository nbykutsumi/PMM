import matplotlib
matplotlib.use("Agg")
from numpy import *
import os, sys
import numpy
import myfunc.util         as util
import myfunc.IO           as IO
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import myfunc.IO.CloudSat  as CloudSat
import myfunc.IO.CloudSat.util as CSutil
import matplotlib.pyplot   as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable


iYM    = [2014,4]
eYM    = [2015,6]
#iYM    = [2014,4]
#eYM    = [2014,4]

lYM    = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [11,12,1,2,3]]


#-- CloudSat ---
prdLev  = "2B"
#prdName = "GEOPROF"
prdName = "CLDCLASS"
prdVer  = "P_R04"
varName = {"GEOPROF" :"Radar_Reflectivity"
          ,"CLDCLASS":"cloud_scenario"
          }[prdName]

cs  = CloudSat.CloudSat(prdLev, prdName, prdVer)
nbin= cs.nbin
liz = arange(nbin)
lz  = arange(-4810, 24950+1, 240) /1000. # [km]

#-- JMA-Cloud ---
clVer = "MyWNP.M.3"
rootDir = "/home/utsumi/mnt/wellshare"
if clVer[:5] == "MyWNP":
  ver        = clVer[5:]
  cl         = CLOUDTYPE.MyCloudWNP(ver=ver)
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.My%s"%(ver)
  ibaseDirCL = rootDir + "/CLOUDTYPE/%s"%(clVer)
  [[lllat,lllon],[urlat,urlon]] = cl.BBox
  dlat       = cl.dLat
  dlon       = cl.dLon
#-- CloudNames --
dclName      = cl.dclName
dclShortName = cl.dclShortName
lclid       = sort(dclName.keys())

dcsName      = CSutil.ret_dclName()
dcsShortName = CSutil.ret_dclShortName()
lcsid        = sort(dcsName.keys())


da1Num = {}
dNum   = {csid:0.0 for csid in lcsid}
nobs   = 0.0
for Year,Mon in lYM:
  print Year,Mon
  srcDir = os.path.join(ibaseDir, "Trc.%s"%(prdName), "%04d.%02d"%(Year,Mon))
  strYM    = "%04d.%02d"%(Year,Mon)
  pathMapcl= os.path.join(srcDir, "mapcl.%s.bn"%(strYM))
  pathProf = os.path.join(srcDir, "ctype.%s.%dlevs.bn"%(strYM,nbin))
  #pathLat  = os.path.join(srcDir, "lat.%s.bn"  %(strYM))
  #pathLon  = os.path.join(srcDir, "lon.%s.bn"  %(strYM))
  #pathTime = os.path.join(srcDir, "tstmp.%s.bn"%(strYM))

  aMapcl   = fromfile(pathMapcl,int16)
  aProf    = fromfile(pathProf, int16).reshape(-1,nbin)
  #aLat     = fromfile(pathLat,  float64)
  #aLon     = fromfile(pathLon,  float64)
  #aDTime   = util.tstmp2dtime(fromfile(pathTime,int32))

  nobs_tmp  = (aProf.shape)[0]
  nobs      = nobs + nobs_tmp
  for csid in lcsid[1:]:
    amsk  = ma.masked_equal(aProf, csid).mask
    if type(amsk)==numpy.bool_:
      anum = zeros(nobs_tmp, float32)
    else:
      anum  = amsk.sum(axis=1)
    da1Num[csid] = ma.masked_greater(anum,0).filled(1.0)

  #-- number of cloud types at each column --
  da1Num[-9999] = zeros(len(anum),float32)   # number of cloud types
  for csid in lcsid[1:]:  # don't count no cloud
    da1Num[-9999] = da1Num[-9999] + da1Num[csid]

  #-- number of each cloud type in entire track---
  for csid in lcsid[1:]:
    anum  = ma.masked_invalid(da1Num[csid] / da1Num[-9999]).filled(0.0)
    dNum[csid] = dNum[csid] + anum.sum()

  #-- number of clear sky --
  dNum[0] = dNum[0] + (ma.masked_equal(da1Num[-9999],0.0).mask).sum()

#-- table ------
outDir = "/home/utsumi/mnt/wellshare/PMM/WNP.261x265/CL.My.M.3/pict"
outPath= outDir + "/table.freq.CloudSat.csv"
#outPath= outDir + "/temp.csv"
sout = "cloud type,frequency\n"
for csid in lcsid:
  freq = float(dNum[csid]) / nobs
  sout = sout + "%s,%s\n"%(dcsShortName[csid], freq)

sout = sout.strip()
f=open(outPath,"w"); f.write(sout); f.close()
print outPath
