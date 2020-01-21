import numpy as np
from numpy import *
import myfunc.util as util
from datetime import datetime, timedelta
import glob, subprocess
import os,sys, socket
import calendar
import h5py
from collections import deque
import myfunc.util as util

myhost = socket.gethostname()
if myhost == 'shui':
    workbaseDir= '/work'
    tankbaseDir= '/tank'
    figDir   = '/home/utsumi/temp/geo'

elif myhost == 'well':
    workbaseDir= '/home/utsumi/mnt/lab_work'
    tankbaseDir= '/home/utsumi/mnt/lab_tank'
    figDir   = '/home/utsumi/temp/stop'

else:
    print 'check hostname',myhost
    sys.exit()

iDTime = datetime(2017,7,1)
eDTime = datetime(2017,7,5)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
ny,nx = 7, 7

aout = None
for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]
    baseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/himawari/obt.ptype/%04d/%02d/%02d'%(Year,Mon,Day)
    ssearch = baseDir + '/*/ptype.npy'

    lsrcPath= glob.glob(ssearch)
    lsrcPath= np.sort(lsrcPath)
    if len(lsrcPath)==0:
        continue
    for srcPath in lsrcPath:
        if aout is None:
            aout = np.load(srcPath)
        else:
            aout = np.concatenate([aout, np.load(srcPath)])

aout = ma.masked_less(aout,0).compressed()
aout = (aout/10000000).astype('int16')

nstra = ma.masked_not_equal(aout,1).count()
nconv = ma.masked_not_equal(aout,2).count()
nothr = ma.masked_not_equal(aout,3).count()
nall  = nstra + nconv + nothr
rstra = nstra/float(nall)
rconv = nconv/float(nall)
rothr = nothr/float(nall)

print nstra, nconv, nothr
print rstra, rconv, rothr


