from numpy import *
from datetime import datetime, timedelta
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import myfunc.util as util
import calendar
import sys

iYM = [2014,4]
eYM = [2015,7]
lYM = util.ret_lYM(iYM, eYM)

#cltype = "MyWNP1"
#cltype = "MyWNP2"
#cltype = "MyWNP3"
cltype = "MyWNP.M.3"

cl = CLOUDTYPE.CloudWNP()
ny, nx  = cl.ny, cl.nx   # 261, 265
#baseDir = "/home/utsumi/mnt/wellshare/PMM/WNP.261x265/MyCLTYPE"
baseDir = "/home/utsumi/mnt/wellshare/CLOUDTYPE/%s"%(cltype)

a2miss  = ones([ny,nx], int32)*(-9999)

def LoadCloudJust(DTime):
  try:
    a2tot = cl.loadData(DTime, DType="tac")
    a2hig = cl.loadData(DTime, DType="ahc")
    a2cnv = cl.loadData(DTime, DType="cvc")
    a2mlc = a2tot - a2hig - a2cnv
  except IOError:
    Day   = DTime.day
    Hour  = DTime.hour
    if "%04d-%02d-%02d-%02d"%(Year,Mon,Day,Hour) in lNoData:
      print "skip", DTime
      a2tot = None
      a2hig = None
      a2cnv = None
      a2mlc = None
    else:
      print "check cloud data"
      print DTime
      print "stop"
      sys.exit()
  return a2tot, a2hig, a2cnv, a2mlc


def LoadCloudMean(DTime):
  DTime0 = DTime - timedelta(minutes=30)
  DTime1 = DTime + timedelta(minutes=30)
  try:
    a2tot0 = cl.loadData(DTime0, DType="tac")
    a2hig0 = cl.loadData(DTime0, DType="ahc")
    a2cnv0 = cl.loadData(DTime0, DType="cvc")
    a2mlc0 = a2tot0 - a2hig0 - a2cnv0

    a2tot1 = cl.loadData(DTime1, DType="tac")
    a2hig1 = cl.loadData(DTime1, DType="ahc")
    a2cnv1 = cl.loadData(DTime1, DType="cvc")
    a2mlc1 = a2tot1 - a2hig1 - a2cnv1

    a2tot  = (a2tot0 + a2tot1)*0.5
    a2hig  = (a2hig0 + a2hig1)*0.5
    a2cnv  = (a2cnv0 + a2cnv1)*0.5
    a2mlc  = (a2mlc0 + a2mlc1)*0.5
  
  except IOError:
    Year0, Mon0, Day0, Hour0 = DTime0.year, DTime0.month, DTime0.day, DTime0.hour
    Year1, Mon1, Day1, Hour1 = DTime1.year, DTime1.month, DTime1.day, DTime1.hour
    if "%04d-%02d-%02d-%02d"%(Year0,Mon0,Day0,Hour0) in lNoData:
      print "skip", DTime
      a2tot, a2hig, a2cnv, a2mlc = None, None, None, None
    elif "%04d-%02d-%02d-%02d"%(Year1,Mon1,Day1,Hour1) in lNoData:
      print "skip", DTime
      a2tot, a2hig, a2cnv, a2mlc = None, None, None, None
    else:
      print "check cloud data"
      print DTime
      print "stop"
      sys.exit()
  return a2tot, a2hig, a2cnv, a2mlc


def MyClassification1(a2tot, a2hig, a2cnv, a2mlc):
  a2cls = a2miss
  # 1: Deep Convection
  a2cls = ma.masked_where(a2cnv > 20., a2cls).filled(1)
  
  # 2: High Clouds
  a2cls = ma.masked_where((a2cnv <=20.) & (a2hig>50.), a2cls).filled(2)
  
  # 3: Mid and Low Clouds
  a2cls = ma.masked_where(a2mlc>50., a2cls).filled(3)
  
  # 4: Mixed Clouds
  a2cls = ma.masked_where((a2cnv<=20.)&(a2hig<=50.)&(a2mlc<=50.), a2cls).filled(4)
  
  # 0: Clear
  a2cls = ma.masked_where(a2tot <=10., a2cls).filled(0)

  return a2cls 
  
def MyClassification2(a2tot, a2hig, a2cnv, a2mlc):
  a2cls = a2miss
  # 1: Deep Convection I
  a2cls = ma.masked_where(a2cnv > 50., a2cls).filled(1)

  # 2: Deep Convection II
  a2cls = ma.masked_where((a2cnv >20.)&(a2cnv <= 50.), a2cls).filled(2)
 
  # 3: High Clouds
  a2cls = ma.masked_where((a2cnv <=20.) & (a2hig>50.), a2cls).filled(3)
  
  # 4: Mid and Low Clouds
  a2cls = ma.masked_where((a2cnv <=20.) & (a2mlc>50.), a2cls).filled(4)
  
  # 5: Mixed Clouds
  a2cls = ma.masked_where((a2cnv<=20.)&(a2hig<=50.)&(a2mlc<=50.), a2cls).filled(5)
  
  # 0: Clear
  a2cls = ma.masked_where(a2tot <=10., a2cls).filled(0)

  return a2cls 

def MyClassification3(a2tot, a2hig, a2cnv, a2mlc):
  a2cls = a2miss
  # 1: Deep Convection I
  a2cls = ma.masked_where(a2cnv > 50., a2cls).filled(1)

  # 2: Deep Convection II
  a2cls = ma.masked_where((a2cnv >10.)&(a2cnv <= 50.), a2cls).filled(2)
 
  # 3: High Clouds
  a2cls = ma.masked_where((a2cnv <=10.) & (a2hig>50.), a2cls).filled(3)
  
  # 4: Mid and Low Clouds
  a2cls = ma.masked_where((a2cnv <=10.) & (a2mlc>50.), a2cls).filled(4)
  
  # 5: Mixed Clouds
  a2cls = ma.masked_where((a2cnv<=10.)&(a2hig<=50.)&(a2mlc<=50.), a2cls).filled(5)
  
  # 0: Clear
  a2cls = ma.masked_where(a2tot <=10., a2cls).filled(0)

  return a2cls 

#--------------
if   cltype == "MyWNP1":
  MyClassification = MyClassification1
elif cltype == "MyWNP2":
  MyClassification = MyClassification2
elif cltype == "MyWNP3":
  MyClassification = MyClassification3
elif cltype == "MyWNP.M.3":
  MyClassification = MyClassification3

if cltype[5:8]==".M.":
  MeanFlag=True
  LoadCloud = LoadCloudMean
else:
  MeanFlag=False
  LoadCloud = LoadCloudJust


for (Year,Mon) in lYM:
  eDay   = calendar.monthrange(Year,Mon)[1]
  if MeanFlag == False:
    iDTime = datetime(Year,Mon,1,0)
    eDTime = datetime(Year,Mon,eDay,23)
    dDTime = timedelta(hours=1)
    lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)
  else:
    iDTime = datetime(Year,Mon,1,0) + timedelta(minutes=30)
    eDTime = datetime(Year,Mon,eDay,23) + timedelta(minutes=30)
    dDTime = timedelta(hours=1)
    lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)

  # No Data list for Cloud
  f = open(cl.baseDir + "/%04d%02d/list.nodata.txt"%(Year,Mon), "r")
  lNoData = [s.strip() for s in f.readlines()]
  f.close()
  #----------------------- 
 
  for DTime in lDTime:
    # Load data
    a2tot, a2hig, a2cnv, a2mlc = LoadCloud(DTime)
    if a2tot == None:
      continue

    # classification
    a2cls = MyClassification(a2tot, a2hig, a2cnv, a2mlc)
 
    # Write
    Day, Hour, Minute = DTime.day, DTime.hour, DTime.minute
    oDir  = baseDir + "/%04d%02d/%02d"%(Year,Mon,Day)
    sPath = oDir + "/CLTYPE.%04d%02d%02d%02d%02d.%dx%d"%(Year,Mon,Day,Hour,Minute,ny,nx)
    util.mk_dir(oDir)
    a2cls.astype(int32).tofile(sPath)

    #if DTime ==lDTime[0]:
    #  print sPath
    print sPath

