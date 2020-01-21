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
lch = [['tir',1],['tir',2],['tir',3],['tir',4],['tir',5],['tir',6],['tir',7],['tir',8],['tir',9],['tir',10],['sir',1],['sir',2]]
#lch = [['sir',1]]
miss = -9999.

lout = [['ch','chnum','ave','std','vmin','vmax']]
for [ch,chnum] in lch:
    aout = None
    for DTime in lDTime:
        Year,Mon,Day = DTime.timetuple()[:3]
        baseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/himawari/obt.ptype/%04d/%02d/%02d'%(Year,Mon,Day)
        ssearch = baseDir + '/*/%s.%02d.npy'%(ch, chnum)

        lsrcPath= glob.glob(ssearch)
        lsrcPath= np.sort(lsrcPath)
        if len(lsrcPath)==0:
            continue
        for srcPath in lsrcPath:
            if aout is None:
                aout = np.load(srcPath)
            else:
                aout = np.concatenate([aout, np.load(srcPath)])

    aout = ma.masked_equal(aout,miss)
    ave = aout.mean()
    std = aout.std()
    vmin= aout.min()
    vmax= aout.max()
    print aout.shape
    print ave, std, vmin, vmax
    ltmp = [ch, chnum, ave, std, vmin, vmax]
    lout.append(ltmp)

 
sout = util.list2csv(lout)

csvDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/himawari/para'
csvPath= csvDir + '/norm-para.csv'
util.mk_dir(csvDir)
f=open(csvPath,'w'); f.write(sout); f.close()
print csvPath


 
