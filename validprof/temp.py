import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import glob
import numpy as np
import h5py
import sys
import myfunc.util as util
from datetime import datetime, timedelta
#-- Read DPR-Ku  ----------------------------------------
iDTime = datetime(2014,8,1)
eDTime = datetime(2014,8,10)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))

for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]
    srcDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/validprof/pair/gprof/%04d/%02d/%02d'%(Year,Mon,Day)
    ssearch = srcDir + '/profpmw.*.npy'
    print ssearch
    lsrcTmp = np.sort(glob.glob(ssearch))
    for srcTmp in lsrcTmp:
        oid = int(srcTmp.split('.')[-2])

        surfPath = srcDir + '/surfaceTypeIndex.%06d.npy'%(oid)
        profPath = srcDir + '/profpmw.%06d.npy'%(oid)

        asurf = np.load(surfPath)
        aprof = np.load(profPath)

        amask1 = ma.masked_equal(asurf,2).mask
        amask2 = ma.masked_inside(asurf,8,11).mask
        amask3 = ma.masked_equal(asurf,14).mask
        amask  = amask1 + amask2 + amask3
        for i in range(asurf.shape[0]):
            #if amask[i]==True:
            if asurf[i]==-99:
                print asurf[i], aprof[i]
