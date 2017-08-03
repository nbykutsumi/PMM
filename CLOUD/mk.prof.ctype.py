import os, sys
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import myfunc.IO.CloudSat  as CloudSat
import myfunc.util         as util
import calendar
import PMM.pmm_fsub        as pmm_fsub
from   collections   import deque
from   datetime      import datetime, timedelta
from   numpy         import *


#-- CloudSat ---
prdLev  = "2B"
#prdName = "GEOPROF"
prdName = "CLDCLASS"
prdVer  = "P_R04"
varName = {"GEOPROF" :"Radar_Reflectivity"
          ,"CLDCLASS":"cloud_scenario"
          }[prdName]

BBox    = [[-0.1, 113.875],[52.1, 180.125]]
cs  = CloudSat.CloudSat(prdLev, prdName, prdVer)
nbin= cs.nbin

#-- JMA-Cloud ---
clVer = "MyWNP.M.3"
rootDir = "/home/utsumi/mnt/well.share"
if clVer[:5] == "MyWNP":
  ver        = clVer[5:]
  cl         = CLOUDTYPE.MyCloudWNP(ver=ver)
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.My%s"%(ver)
  ibaseDirCL = rootDir + "/CLOUDTYPE/%s"%(clVer)


iYM    = [2014,4]
eYM    = [2014,4]
lYM    = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [11,12,1,2,3]]
#lYM = [YM for YM in lYM if YM[1] not in [11,12,1,2,3] + [5,7,9]]
#lYM = [YM for YM in lYM if YM[1] not in [11,12,1,2,3] + [4,6,8,10]]

for Year,Mon in lYM:
  iDay   = 1
  eDay   = calendar.monthrange(Year,Mon)[1]
  #eDay   = 2
  iDTime = datetime(Year,Mon,iDay,0,0)
  eDTime = datetime(Year,Mon,eDay,23,0) + timedelta(hours=1)

  print "*"*50
  print iDTime, "-->", eDTime
  print "*"*50
  try:
    csobt  = cs(varName, iDTime, eDTime, BBox=BBox, verbose=False)
  except IOError:
    print "Skip!", DTime,"-",DTime+dDTime
    print "-"*50
    continue

  if len(csobt.data)==0:
    continue 

  a2ctype= cs.Resolve_CloudScenario(csobt.data[0]).astype(int)

  aProf  = a2ctype.flatten()
  aLat   = csobt.lat[0]
  aLon   = csobt.lon[0]   
  aDTime = csobt.dtime[0]

  aTime = [util.dtime2tstmp(DTime)
                    for DTime in aDTime]

  # Save
  oDir = os.path.join(ibaseDir, "Trc.%s"%(prdName), "%04d.%02d"%(Year,Mon))
  util.mk_dir(oDir)

  strYM    = "%04d.%02d"%(Year,Mon)
  pathProf = os.path.join(oDir, "ctype.%s.%dlevs.bn"%(strYM,nbin))
  pathLat  = os.path.join(oDir, "lat.%s.bn"  %(strYM))
  pathLon  = os.path.join(oDir, "lon.%s.bn"  %(strYM))
  pathTime = os.path.join(oDir, "tstmp.%s.bn"%(strYM))

  array(aProf).astype(int16  ).tofile(pathProf)
  array(aLat ).astype(float64).tofile(pathLat )
  array(aLon ).astype(float64).tofile(pathLon )
  array(aTime).astype(int32  ).tofile(pathTime)
  print pathProf

