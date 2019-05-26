from numpy import *
import numpy as np
from datetime import datetime, timedelta
import myfunc.util as util

iDTime = datetime(2017,1,1)
eDTime = datetime(2017,1,10)
dDTime = timedelta(days=1)
lDTime = util.ret_lDTime(iDTime,eDTime,dDTime)
isurf = 3
for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]
    tcDir = '/work/hk01/utsumi/PMM/stop/data/Tc/%04d/%02d/%02d'%(Year,Mon,Day)
    tcPath1 = tcDir + '/Tc1.1dy.2dx.%02dsurf.npy'%(isurf)
    tcPath2 = tcDir + '/Tc1.0dy.0dx.%02dsurf.npy'%(isurf)
    a1tc1 = np.load(tcPath1)
    a1tc2 = np.load(tcPath2)
    
    stopDir='/work/hk01/utsumi/PMM/stop/data/stop/%04d/%02d/%02d'%(Year,Mon,Day)
    stopPath1 = stopDir + '/stop.1dy.2dx.%02dsurf.npy'%(isurf)
    stopPath2 = stopDir + '/stop.0dy.0dx.%02dsurf.npy'%(isurf)
    a1stop1 = np.load(stopPath1)
    a1stop2 = np.load(stopPath2)
    
    #print a1tc1.shape, a1tc2.shape, a1stop1.shape, a1stop2.shape
    print a1stop1[101],a1stop2[101]
