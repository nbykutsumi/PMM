import numpy as np
from numpy import *
import myfunc.util as util
from datetime import datetime, timedelta
import glob
import os,sys
import calendar
import h5py
from collections import deque
import myfunc.util as util

varName = 'nltb'
iYM = [2017,1]
eYM = [2017,12]
lYM = util.ret_lYM(iYM,eYM)
#outDir= '/work/hk01/utsumi/PMM/TPCDB/PC_COEF'
outDir= '/work/hk01/utsumi/PMM/stop/data'
verGMI = '05'
subverGMI = 'A'
fullverGMI = '%s%s'%(verGMI,subverGMI)
miss = -9999.
for (Year,Mon) in lYM:
    eDay   = calendar.monthrange(Year,Mon)[1]
    iDTime = datetime(Year,Mon,1)
    eDTime = datetime(Year,Mon,eDay)

    dDTime = timedelta(days=1)
    lDTime = util.ret_lDTime(iDTime,eDTime,dDTime)
    matchBaseDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A'
    #-- Read list -----
    listDir  = '/work/hk01/utsumi/PMM/TPCDB/list'
    listPath = listDir + '/list.1C.V05.%04d%02d.csv'%(Year,Mon)
    f=open(listPath,'r'); lines = f.readlines(); f.close()
    dlorbit = {}
    for line in lines:
        line = map(int, line.split(','))
        oid,Year,Mon,Day,itime,etime = line
        try:
            dlorbit[Day].append(line)
        except KeyError: 
            dlorbit[Day] = [line]
    #-------------------
    for DTime in lDTime:
        Day = DTime.day
        try:
            lorbit = dlorbit[Day]
        except KeyError:
            continue

        #-- Initialize --
        dastop = {}
        for isurf in range(1,15+1):
            dastop[isurf] = deque([])

        #----------------
        #lorbit = lorbit[:2]  # test
        for orbinfo in lorbit:
            oid,Year,Mon,Day,itime,etime = orbinfo
    
            #-- Storm Top Height ----
            stopDir  = matchBaseDir + '/S1.ABp103-117.Ku.V06A.heightStormTop/%04d/%02d/%02d'%(Year,Mon,Day)
            stopPath = stopDir + '/heightStormTop.1.%06d.npy'%(oid)
            a2stop = np.load(stopPath)
            if a2stop.max()<=0: continue
  
            a1stop = a2stop.flatten() 
    
            #-- Surface Type Index --
            surftypeDir = matchBaseDir + '/S1.ABp103-117.GMI.surfaceTypeIndex/%04d/%02d/%02d'%(Year,Mon,Day)
            surftypePath= surftypeDir + '/surfaceTypeIndex.%06d.npy'%(oid)
            a1surftype = np.load(surftypePath).flatten()
    
            #-- Make flag array -----
            a1flagStop = ma.masked_greater(a2stop,0).mask.flatten()
            for isurf in range(1,15+1):
                a1flagSurf = ma.masked_equal(a1surftype,isurf).mask
                a1flag = a1flagStop * a1flagSurf
                if a1flag.any()==False: continue

                a1sc   = a1stop[a1flag]
                dastop[isurf].extend(a1sc)
                    
    
        #******** Save ******************
        for isurf in range(1,15+1):
            aout = dastop[isurf]
            Year,Mon,Day = DTime.timetuple()[:3]
            outDir = '/work/hk01/utsumi/PMM/stop/data/stop/%04d/%02d/%02d'%(Year,Mon,Day)
            outPath= outDir + '/stop.%02dsurf.npy'%(isurf)
            util.mk_dir(outDir)
            np.save(outPath, aout)

            if isurf==0:
                print outPath
            
            

 
