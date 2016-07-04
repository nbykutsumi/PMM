import os
from datetime import datetime, timedelta
import myfunc.util as util
import calendar

#baseDir = "/tank/utsumi/CLOUDTYPE/WNPAC/201404/20140408/"
baseDir = "/tank/utsumi/CLOUDTYPE/WNPAC"

iYM    = [2014,4]
eYM    = [2015,12]
lYM    = util.ret_lYM(iYM,eYM)

for YM in lYM:
  Year   = YM[0]
  Mon    = YM[1]
  eDay   = calendar.monthrange(Year,Mon)[1]
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
    srcDir  = baseDir + "/%04d%02d/%04d%02d%02d"%(Year,Mon,Year,Mon,Day)
    srcPath = srcDir  + "/Z__C_RJTD_%04d%02d%02d%02d0000_OBS_SAT_PSclc_RDnwp_grib2.bin"%(Year,Mon,Day,Hour)
    if not os.path.exists(srcPath):
      lout.append("%04d-%02d-%02d-%02d"%(Year,Mon,Day,Hour))
      print Year,Mon,Day,Hour
 
  # Save monthly file
  oDir = baseDir + "/%04d%02d"%(Year,Mon)
  oPath= oDir + "/list.nodata.txt"

  sout = "\n".join(lout).strip()
  f=open(oPath, "w"); f.write(sout); f.close()
  print oPath  
    
  
