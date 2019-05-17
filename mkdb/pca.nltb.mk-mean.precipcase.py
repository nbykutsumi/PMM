import numpy as np
from numpy import *
import myfunc.util as util
from datetime import datetime, timedelta
import glob
import os,sys
import calendar
import mkdbfunc

#calcmon = True
calcmon = False
calcall = True
#calcall = False
varName = 'nltb'
iYM = [2017,1]
eYM = [2017,12]
lYM = util.ret_lYM(iYM,eYM)
outDir= '/work/hk01/utsumi/PMM/TPCDB/PC_COEF'

for (Year,Mon) in lYM:
    if calcmon != True: continue

    eDay   = calendar.monthrange(Year,Mon)[1]
    iDTime = datetime(Year,Mon,1)
    eDTime = datetime(Year,Mon,eDay)
    #iDTime = datetime(2017,7,12)
    #eDTime = datetime(2017,7,12)

    dDTime = timedelta(days=1)
    lDTime = util.ret_lDTime(iDTime,eDTime,dDTime)
    matchBaseDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A'
   
    #-- initialize ----
    asum = None
    asum2= None
    anum = None


    #atmp = []   # test
    #------------------
    #lDTime = lDTime[:1]  # test  
    for DTime in lDTime:
        print DTime
        Year,Mon,Day = DTime.timetuple()[:3]
        #-- Make KuPR HDF file list ---
        dprDir = matchBaseDir + '/S1.ABp083-137.Ku.V06A.9ave.precipRate/%04d/%02d/%02d'%(Year,Mon,Day)
        #ssearch= dprDir + '/precipRate.016173.npy'
        ssearch= dprDir + '/precipRate.??????.npy'

        ##-- test --
        #dprDir = matchBaseDir + '/S1.ABp103-117.Ku.V06A.precipRate/%04d/%02d/%02d'%(Year,Mon,Day)
        #ssearch= dprDir + '/precipRate.1.??????.npy'
        #print ssearch
        ##----------


        ldprPath= glob.glob(ssearch)
        if len(ldprPath)==0:continue


        #ldprPath = ldprPath[:4] # test
        for dprPath in ldprPath:
            oid = int(dprPath.split('.')[-2]) 

            #-- Read KuPR HDF file -----
            a3dpr = np.load(dprPath)[:,103-83:103-83+15,:]
            #a3dpr = np.load(dprPath)[:]  # test

            #-- Read Tb--------
            tcDir1 = matchBaseDir + '/S1.ABp103-117.GMI.Tc/%04d/%02d/%02d'%(Year,Mon,Day)
            tcDir2 = matchBaseDir + '/S1.ABp103-117.GMI.TcS2/%04d/%02d/%02d'%(Year,Mon,Day)
            
            tcPath1= tcDir1 + '/Tc.%06d.npy'%(oid)
            tcPath2= tcDir2 + '/TcS2.1.%06d.npy'%(oid)
            a3tc1 = np.load(tcPath1)
            a3tc2 = np.load(tcPath2)
            a3tc  = concatenate([a3tc1,a3tc2],axis=2)
            ntc   = a3tc.shape[2]
            #print a3dpr.shape, a3tc1.shape, a3tc2.shape, a3tc.shape
            #-- Extract pixels with precip (in column) --
            if a3dpr.max()<=0: continue

            a1flagtc = ma.masked_inside(a3tc,50,350).mask.all(axis=2).flatten()
            a1flagpr = ma.masked_greater(a3dpr,0).mask.any(axis=2).flatten()
            a1flag = a1flagtc * a1flagpr
            a2tc = a3tc.reshape(-1,ntc)
            a2tc = a2tc[a1flag]

            #-- Make non-linear combination --
            nlen  = a2tc.shape[0]
            adat  = mkdbfunc.mk_nonlin_comb(a2tc.reshape(-1,1,ntc))
            adat  = adat.reshape(nlen,-1)

            print a3tc.shape, adat.shape
            #-- Calc sum, num --
            asumTmp  = adat.sum(axis=0)
            asum2Tmp = np.square(adat).sum(axis=0)
            anumTmp  = ma.count(adat, axis=0)
            try:
                asum  = asum  + asumTmp
                asum2 = asum2 + asum2Tmp
                anum  = anum  + anumTmp
            except TypeError:
                asum  = asumTmp
                asum2 = asum2Tmp
                anum  = anumTmp

            ##-- test --
            #atmp.append(adat)

    ##-- test ---
    #atmp = concatenate(atmp,axis=0)
    #ameantmp = atmp.mean(axis=0)
    #astdtmp  = atmp.std(axis=0) 

    #-- save and calc ---
    sym   = '%04d.%02d'%(Year,Mon)
    amean = asum / anum
    astd  = np.sqrt( (asum2 - 2*amean*asum + anum*np.square(amean))/(anum-1) )
    sumPath  = outDir + '/sum.%s.%s.npy'%(varName, sym)
    sum2Path = outDir + '/sum2.%s.%s.npy'%(varName, sym)
    numPath  = outDir + '/num.%s.%s.npy'%(varName, sym)
    #meanPath = outDir + '/mean.%s.npy'%(varName)
    #stdPath  = outDir + '/std.%s.npy'%(varName)
    np.save(sumPath, asum)
    np.save(sum2Path,asum2)
    np.save(numPath, anum)
    print sumPath
    #--------------------


#--- Make mean and std over the period --
if calcall ==True:
    asum = None
    asum2= None
    anum = None
    for (Year,Mon) in lYM:
        sym   = '%04d.%02d'%(Year,Mon)
        sumPath  = outDir + '/sum.%s.%s.npy'%(varName, sym)
        sum2Path = outDir + '/sum2.%s.%s.npy'%(varName, sym)
        numPath  = outDir + '/num.%s.%s.npy'%(varName, sym)
     
        asumTmp = np.load(sumPath)  
        asum2Tmp= np.load(sum2Path)  
        anumTmp = np.load(numPath)  
        try:
            asum  = asum  + asumTmp
            asum2 = asum2 + asum2Tmp
            anum  = anum  + anumTmp
        except TypeError:
            asum  = asumTmp
            asum2 = asum2Tmp
            anum  = anumTmp

    amean = asum / anum
    astd  = np.sqrt( (asum2 - 2*amean*asum + anum*np.square(amean))/(anum-1) )

    meanPath = outDir + '/mean.%s.npy'%(varName)
    stdPath  = outDir + '/std.%s.npy'%(varName)
    np.save(meanPath, amean)
    np.save(stdPath,  astd)
    print stdPath



