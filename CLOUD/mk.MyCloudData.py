from numpy import *
from datetime import datetime, timedelta
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import myfunc.util as util
import calendar

iYM = [2015,7]
eYM = [2015,7]
lYM = util.ret_lYM(iYM, eYM)

#cltype = "MyWNP1"
cltype = "MyWNP2"

cl = CLOUDTYPE.CloudWNP()
ny, nx  = cl.ny, cl.nx   # 261, 265
#baseDir = "/home/utsumi/mnt/well.share/PMM/WNP.261x265/MyCLTYPE"
baseDir = "/home/utsumi/mnt/well.share/CLOUDTYPE/%s"%(cltype)

a2miss  = ones([ny,nx], int32)*(-9999)

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

#--------------
if   cltype == "MyWNP1":
  MyClassification = MyClassification1
elif cltype == "MyWNP2":
  MyClassification = MyClassification2


for (Year,Mon) in lYM:
 
  eDay   = calendar.monthrange(Year,Mon)[1]
  iDTime = datetime(Year,Mon,1,0)
  eDTime = datetime(Year,Mon,eDay,23)
  dDTime = timedelta(hours=1)
  lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)

  # No Data list for Cloud
  f = open(cl.baseDir + "/%04d%02d/list.nodata.txt"%(Year,Mon), "r")
  lNoData = [s.strip() for s in f.readlines()]
  f.close()
  #----------------------- 
 
  for DTime in lDTime:
    Day   = DTime.day
    Hour  = DTime.hour

    # Load data
    try:
      a2tot = cl.loadData(DTime, DType="tac")
      a2hig = cl.loadData(DTime, DType="ahc")
      a2cnv = cl.loadData(DTime, DType="cvc")
      a2mlc = a2tot - a2hig - a2cnv
  
    except IOError:
      if "%04d-%02d-%02d-%02d"%(Year,Mon,Day,Hour) in lNoData:
        print "skip", DTime
      else:
        print "check cloud data"
        print DTimeCL
        print "stop"
        sys.exit()

    a2cls = MyClassification(a2tot, a2hig, a2cnv, a2mlc)
 
    # Write
    oDir  = baseDir + "/%04d%02d/%02d"%(Year,Mon,Day)
    sPath = oDir + "/CLTYPE.%04d%02d%02d%02d.%dx%d"%(Year,Mon,Day,Hour,ny,nx)
    util.mk_dir(oDir)
    a2cls.astype(int32).tofile(sPath)

    if DTime ==lDTime[0]:
      print sPath


