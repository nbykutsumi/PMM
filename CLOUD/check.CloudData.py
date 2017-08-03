import os
from datetime import datetime, timedelta
import myfunc.util as util
import calendar
import myfunc.IO.CLOUDTYPE as CLOUDTYPE

ibaseDir = "/tank/utsumi/CLOUDTYPE/WNPAC"

clVer = "MyWNP2"
clVer = "MyWNP.M.3"
#obaseDir= "/home/utsumi/mnt/well.share/CLOUDTYPE/MyWNP3"
obaseDir= "/home/utsumi/mnt/well.share/CLOUDTYPE/%s"%(clVer)

if clVer[5:8]==".M.":
  MidFlag = True
else:
  MidFlag = False

iYM    = [2014,4]
eYM    = [2015,7]
#eYM    = [2014,4]
lYM    = util.ret_lYM(iYM,eYM)

for YM in lYM:
  Year   = YM[0]
  Mon    = YM[1]
  eDay   = calendar.monthrange(Year,Mon)[1]
  if MidFlag   == True:
    iDTime = datetime(Year,Mon,1,0,30)
    eDTime = datetime(Year,Mon,eDay,23,30)
  elif MidFlag == False:
    iDTime = datetime(Year,Mon,1,0)
    eDTime = datetime(Year,Mon,eDay,23)

  dDTime = timedelta(hours=1)
  lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)
 
  lout   = [] 
  for DTime in lDTime:
    Year  = DTime.year
    Mon   = DTime.month
    Day   = DTime.day
    Hour  = DTime.hour
    Minute= DTime.minute

    if MidFlag ==True:
      DTime0 = DTime - timedelta(minutes=30)
      DTime1 = DTime + timedelta(minutes=30)
      Year0 = DTime0.year
      Mon0  = DTime0.month
      Day0  = DTime0.day
      Hour0 = DTime0.hour
   
      Year1 = DTime1.year
      Mon1  = DTime1.month
      Day1  = DTime1.day
      Hour1 = DTime1.hour

      srcDir0  = ibaseDir + "/%04d%02d/%04d%02d%02d"%(Year0,Mon0,Year0,Mon0,Day0)
      srcDir1  = ibaseDir + "/%04d%02d/%04d%02d%02d"%(Year1,Mon1,Year1,Mon1,Day1)
      srcPath0 = srcDir0  + "/Z__C_RJTD_%04d%02d%02d%02d0000_OBS_SAT_PSclc_RDnwp_grib2.bin"%(Year0,Mon0,Day0,Hour0)
      srcPath1 = srcDir1  + "/Z__C_RJTD_%04d%02d%02d%02d0000_OBS_SAT_PSclc_RDnwp_grib2.bin"%(Year1,Mon1,Day1,Hour1)

      exist0 = os.path.exists(srcPath0)
      exist1 = os.path.exists(srcPath1)
      if exist0*exist1==False:
        lout.append("%04d-%02d-%02d-%02d-%02d"%(Year,Mon,Day,Hour,Minute))
        print Year,Mon,Day,Hour,Minute

    elif MidFlag == False: 
      srcDir  = ibaseDir + "/%04d%02d/%04d%02d%02d"%(Year,Mon,Year,Mon,Day)
      srcPath = srcDir  + "/Z__C_RJTD_%04d%02d%02d%02d0000_OBS_SAT_PSclc_RDnwp_grib2.bin"%(Year,Mon,Day,Hour)
      if not os.path.exists(srcPath):
        lout.append("%04d-%02d-%02d-%02d-%02d"%(Year,Mon,Day,Hour,Minute))
        print Year,Mon,Day,Hour,Minute
 
  # Save monthly file
  oDir = obaseDir + "/%04d%02d"%(Year,Mon)
  oPath= oDir + "/list.nodata.txt"
  util.mk_dir(oDir)
  
  sout = "\n".join(lout).strip()
  f=open(oPath, "w"); f.write(sout); f.close()
  print oPath  
    
  
