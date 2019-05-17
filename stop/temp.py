from numpy import *
import numpy as np
from datetime import datetime, timedelta
import myfunc.util as util

iDTime = datetime(2017,1,1)
eDTime = datetime(2017,5,31)
dDTime = timedelta(days=1)
lDTime = util.ret_lDTime(iDTime,eDTime,dDTime)

isurf = 3
n = 0
for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]
    #tcDir = '/work/hk01/utsumi/PMM/stop/data/Tc/%04d/%02d/%02d'%(Year,Mon,Day)
    #tcPath = tcDir + '/Tc1.1dy.2dx.%02dsurf.npy'%(isurf)
    #a1tc = np.load(tcPath)
    
    stopDir='/work/hk01/utsumi/PMM/stop/data/stop/%04d/%02d/%02d'%(Year,Mon,Day)
    stopPath = stopDir + '/stop.%02dsurf.npy'%(isurf)
    a1stop = np.load(stopPath)
    n = n + len(a1stop)
    print DTime,len(a1stop)

print n 
