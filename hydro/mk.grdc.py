from numpy import *
from datetime import datetime, timedelta
from collections import deque
import myfunc.IO.GRDC as GRDC
import myfunc.util as util
#----------------------------------------
grdc    = GRDC.GRDC()

stnid   = 3625320
#----
srcDir  = '/work/data2/GRDC/dat_day'
srcPath = srcDir + '/%07d.day'%(stnid)
iDTime = datetime(2005,1,4)
eDTime = datetime(2005,1,8)
dDTime = timedelta(days=1)
lDTime = util.ret_lDTime(iDTime,eDTime,dDTime)
ldat   = grdc.mk_grdcDailyArray(srcPath, iDTime, eDTime)
for i in range(len(ldat)):
    print lDTime[i], ldat[i]
