from   numpy         import *
from   datetime      import datetime, timedelta
from   PMM.pmm_fsub  import *
from   bisect        import bisect_left, bisect_right
import os,sys
import calendar
import myfunc.util         as util
import myfunc.IO           as IO
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import myfunc.IO.CloudSat  as CloudSat

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

miss_out = -9999
#-- JMA-Cloud ---
clVer = "MyWNP.M.3"
rootDir = "/home/utsumi/mnt/well.share"
if clVer[:5] == "MyWNP":
  ver        = clVer[5:]
  cl         = CLOUDTYPE.MyCloudWNP(ver=ver)
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.My%s"%(ver)
  ibaseDirCL = rootDir + "/CLOUDTYPE/%s"%(clVer)
  [[lllat,lllon],[urlat,urlon]] = cl.BBox
  dlat       = cl.dLat
  dlon       = cl.dLon

iYM    = [2015,6]
eYM    = [2015,7]
#iYM    = [2014,5]
#eYM    = [2015,5]

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
  dDTime = timedelta(hours=1)
  lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)

  srcDir = os.path.join(ibaseDir, "Trc.%s"%(prdName), "%04d.%02d"%(Year,Mon))
  strYM    = "%04d.%02d"%(Year,Mon)
  pathLat  = os.path.join(srcDir, "lat.%s.bn"  %(strYM))
  pathLon  = os.path.join(srcDir, "lon.%s.bn"  %(strYM))
  pathTime = os.path.join(srcDir, "tstmp.%s.bn"%(strYM))
  pathProf = os.path.join(srcDir, "ctype.%s.%dlevs.bn"%(strYM,nbin)) # test

  aLat      = fromfile(pathLat, float64)
  aLon      = fromfile(pathLon, float64)
  aDTime    = util.tstmp2dtime(fromfile(pathTime,int32))
  aProf     = fromfile(pathProf, int16).reshape(-1,nbin)   # test

  aCL       = empty(len(aLat)) 
  for DTime in lDTime:
  #for DTime in [datetime(2014,8,31,4,0)]:
    iIdx    = bisect_left (aDTime, DTime)
    eIdx    = bisect_right(aDTime, DTime+timedelta(hours=1))
    alat    = aLat[iIdx:eIdx]
    alon    = aLon[iIdx:eIdx]

    if len(alat) ==0:
      continue
    print DTime

    try:
      a2cl    = cl.loadData(DTime+timedelta(minutes=30)) 
    except IOError:
      continue

    aCL[iIdx:eIdx] = (pmm_fsub.pickup_data(a2cl.T, lllat, lllon, urlat, urlon, dlat, dlon, alon, alat, miss_out).T)[0]


    ## test ------------


    #c1 = ma.masked_not_equal(aCL[iIdx:eIdx], 1.0)

    #ainfo  = zip(c1.compressed(),alon[~c1.mask], alat[~c1.mask])
    #if len(c1.compressed()) >0:
    #   print "-"*50
    #   print DTime
    #   for sinfo in  ainfo:
    #     print sinfo

    #aprof  = aProf[iIdx:eIdx][~c1.mask]

    #ii    = bisect_left (aDTime, datetime(2014,8,31,4,19,39))
    #ei    = bisect_right(aDTime, datetime(2014,8,31,4,22,50))
    #adtime_seg = aDTime[ii:ei]
    #alat_seg   = aLat[ii:ei]
    #alon_seg   = aLon[ii:ei]
    #aprof_seg  = aProf[ii:ei,:]
    #c1_seg= aCL[ii:ei]
    ##------------------


   
  # write to file
  oPath  = os.path.join(srcDir, "mapcl.%s.bn"%(strYM))
  aCL.astype(int16).tofile(oPath)
  print oPath 
