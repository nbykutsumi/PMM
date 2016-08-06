from numpy import *
from datetime import datetime, timedelta
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import myfunc.util as util
import calendar
import sys

iYM   = [2014,4]
eYM   = [2015,6]
lYM   = util.ret_lYM(iYM, eYM)

cl    = CLOUDTYPE.CloudWNP()
Lat   = cl.Lat
Lon   = cl.Lon
ny    = cl.ny
nx    = cl.nx

lcltype = range(0,7+1)
dclid   = {0:0, 1:1, 2:201, 3:202, 4:4, 5:3, 6:204, 7:200}

dclName ={0:"Clear Sky",   1:"Cumulonimbus(Cb)",  2:"High Cloud",3:"Mid Cloud"
         ,4:"Cumulus(Cu)", 5:"Stratocumulus(Sc)", 6:"Fog/St"    ,7:"Cloudy"}

dclShortName={0:"no", 1:"Cb",  2:"hi",3:"md"
             ,4:"Cu", 5:"Sc",  6:"St",7:"cw"}


a2one = ones([ny,nx],int32)

for YM in lYM:
  Year   = YM[0]
  Mon    = YM[1]
  iDay   = 1
  eDay   = calendar.monthrange(Year,Mon)[1]

  iDTime = datetime(Year,Mon,iDay,0,)
  eDTime = datetime(Year,Mon,eDay,23)

  dDTime = timedelta(hours=1)
  lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)

  # No Data list for Cloud
  f = open(cl.baseDir + "/%04d%02d/list.nodata.txt"%(Year,Mon), "r")
  lNoData = [s.strip() for s in f.readlines()]
  f.close()

  a3num  = zeros([len(lcltype),ny,nx],float32)
  for DTime in lDTime:
    Year = DTime.year
    Mon  = DTime.month
    Day  = DTime.day
    Hour = DTime.hour
    if "%04d-%02d-%02d-%02d"%(Year,Mon,Day,Hour) in lNoData:
      continue
    a2cl = cl.loadData(DTime, DType="clc")
    for icl in lcltype:
      a3num[icl] = a3num[icl] + ma.masked_where(a2cl !=dclid[icl], a2one).filled(0.0)
  #-- Save ----------
  for icl in lcltype:
    sDir   = "/tank/utsumi/CLOUDTYPE/WNPAC/num"
    sPath  = sDir + "/num.%04d%02d.%s.%dx%d"%(Year,Mon,dclShortName[icl],ny,nx)
    util.mk_dir(sDir)

    a3num[icl].astype(int32).tofile(sPath)
    print sPath

    
