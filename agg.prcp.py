from numpy import *
from datetime import datetime, timedelta
import ioPMM
import myfunc.util as util
import calendar

iYM = [2014,12]
eYM = [2014,12]
lYM = util.ret_lYM(iYM, eYM)

iopmm = ioPMM.ioRAdom()
iopmm.init_GSMaP(prj="standard", ver="v6")
ny,nx = iopmm.ny, iopmm.nx
miss  = -999.

a2one = ones([ny,nx],float32)
for Year, Mon in lYM:
  print Year, Mon
  eDay   = calendar.monthrange(Year, Mon)[1]
  iDTime = datetime(Year,Mon,1,0)
  eDTime = datetime(Year,Mon,eDay,23)
  dDTime = timedelta(hours=1)
  lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)

  a2sum  = zeros([ny,nx])
  a2num  = zeros([ny,nx])
  for DTime in lDTime:
    try: #---------------   Caution !!  --------
      a2in  = iopmm.loadRAmmh(DTime, maskflag=True)
    except IOError:
      if DTime == datetime(2014,12,31,23):
        pass   
    a2num = a2num + ma.masked_where(a2in.mask, a2one).filled(0.0)
    a2sum = a2sum + a2in.filled(0.0)

  a2ave = (ma.masked_where(a2num==0.0, a2sum)/a2num).filled(miss)
  #- save ----
  baseDir, sDir, sumPath, numPath, avePath = iopmm.ret_pathMonRAmmh(Year,Mon)
  a2sum.astype(float32).tofile(sumPath)
  a2num.astype(float32).tofile(numPath)
  a2ave.astype(float32).tofile(avePath)
  print avePath 
  
