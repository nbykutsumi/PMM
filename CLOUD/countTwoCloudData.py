from numpy import *
from datetime import datetime, timedelta
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import myfunc.util as util
import calendar
#
iYM  = [2014,4]
eYM  = [2014,10]
lYM  = util.ret_lYM(iYM, eYM)

cl1 = CLOUDTYPE.CloudWNP()
cl2 = CLOUDTYPE.MyCloudWNP()

lcl1 = range(0,7+1)
lcl2 = range(0,4+1)

dclid1   = {0:0, 1:1, 2:201, 3:202, 4:4, 5:3, 6:204, 7:200}
dclid2   = {0:0, 1:1, 2:2, 3:3, 4:4}

dnum = {(icl1, icl2):0 for icl1 in lcl1 for icl2 in lcl2}

for (Year,Mon) in lYM:
  # No Data list for Cloud
  f = open(cl1.baseDir + "/%04d%02d/list.nodata.txt"%(Year,Mon), "r")
  lNoData = [s.strip() for s in f.readlines()]
  f.close()
  #
  eDay   = calendar.monthrange(Year,Mon)[1]
  iDTime = datetime(Year,Mon,1,0)
  eDTime = datetime(Year,Mon,eDay,23)
  dDTime = timedelta(hours=1)
  lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)
  for DTime in lDTime:
    print DTime
    try:
      a2dat1 = cl1.loadData(DTime)
      a2dat2 = cl2.loadData(DTime)
    except IOError:
      Year = DTime.year
      Mon  = DTime.month
      Day  = DTime.day
      Hour = DTime.hour
      if "%04d-%02d-%02d-%02d"%(Year,Mon,Day,Hour) in lNoData:
        print "skip", DTime
      else:
        print "check cloud data"
        print DTime
        print "stop"
        sys.exit()
    #-----------------

    for icl1 in lcl1:
      a  = ma.masked_not_equal(a2dat1, dclid1[icl1])
      for icl2 in lcl2:
        b  = ma.masked_where(a2dat2 != dclid2[icl2], a)
        dnum[icl1, icl2] = dnum[icl1, icl2] + b.count()
  

sout = "" + "," + ",".join(map(str, lcl2)) + "\n"
for icl1 in lcl1:
  sout = sout + "%d"%(icl1) + "," + ",".join(map(str, [dnum[icl1, icl2] for icl2 in lcl2])) + "\n"

rootDir = "/home/utsumi/mnt/well.share"
baseDir = rootDir + "/PMM/WNP.261x265/CL.My1/pict"
oPath = baseDir + "/vsCount.csv"
f = open(oPath, "w"); f.write(sout); f.close() 
print oPath
