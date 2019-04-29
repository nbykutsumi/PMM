import numpy as np
from numpy import *
import myfunc.util as util
from datetime import datetime, timedelta
import glob
import os,sys
import calendar

iYM = [2017,1]
eYM = [2017,1]
lYM = util.ret_lYM(iYM,eYM)
for (Year,Mon) in lYM:
    eDay   = calendar.monthrange(Year,Mon)[1]
    iDTime = datetime(Year,Mon,1)
    eDTime = datetime(Year,Mon,eDay)
    dDTime = timedelta(days=1)
    lDTime = util.ret_lDTime(iDTime,eDTime,dDTime)
    matchBaseDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A'
   
    #-- initialize ----
    asum = None
    asum2= None
    anum = None
    #------------------
    lDTime = lDTime[:1]  # test  
    for DTime in lDTime:
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

        for dprPath in ldprPath:
            oid = int(dprPath.split('.')[-2]) 

            #-- Read KuPR HDF file -----
            a3dpr = np.load(dprPath)[:,103-83:103-83+15,:]
            #a3dpr = np.load(dprPath)[:]  # test

            #-- Read Tc--------
            tcDir1 = matchBaseDir + '/S1.ABp103-117.GMI.Tc/%04d/%02d/%02d'%(Year,Mon,Day)
            tcDir2 = matchBaseDir + '/S1.ABp103-117.GMI.TcS2/%04d/%02d/%02d'%(Year,Mon,Day)
            
            tcPath1= tcDir1 + '/Tc.%06d.npy'%(oid)
            tcPath2= tcDir2 + '/TcS2.1.%06d.npy'%(oid)
            a3tc1 = np.load(tcPath1)
            a3tc2 = np.load(tcPath2)
            a3tc  = concatenate([a3tc1,a3tc2],axis=2)

            print a3dpr.shape, a3tc1.shape, a3tc2.shape, a3tc.shape
            
            #-- Extract pixels with precip (in column) --
            a2flag = ma.masked_greater(a3dpr,0).mask.any(axis=2)
            a1flag = a2flag.flatten()
            ntc  = a3tc.shape[2]
            a2tc = a3tc.reshape(-1,ntc)
            a2tc = a2tc[a1flag]
            print a2tc.shape
            print ''

            #-- Calc sum, num --
            adatTmp  = ma.masked_equal(a2tc, miss)
            asumTmp  = avar.sum(axis=0)
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

    #-- save and calc ---
    outDir= '/work/hk01/utsumi/PMM/TBPCDB'
    #amean = asum / anum
    #astd  = np.sqrt( (asum2 - 2*amean*asum + anum*np.square(amean))/(anum-1) )
    numPath  = outDir + '/num.%s.npy'%(varName)
    sumPath  = outDir + '/sum.%s.npy'%(varName)
    sum2Path = outDir + '/sum2.%s.npy'%(varName)
    #meanPath = outDir + '/mean.%s.npy'%(varName)
    #stdPath  = outDir + '/std.%s.npy'%(varName)
    np.save(meanPath, amean)
    #--------------------


'''
lvarName= ['S1.ABp103-117.GMI.Tc','S1.ABp103-117.GMI.TcS2']
#lvarName= ['S1.ABp103-117.GMI.Tc']


for varName in lvarName:
    varNameSimple = varName.split('.')[-1]

    if varName in ['S1.ABp103-117.GMI.Tc','S1.ABp103-117.GMI.TcS2']:
        ndim = 3
        miss = -9999.
    else:
        print 'check varName',varName
        sys.exit()
   
    #-- initialize ----
    asum = None
    asum2= None
    anum = None
    #------------------- 
    for DTime in lDTime:
        Year,Mon,Day = DTime.timetuple()[:3]
        varDir = baseDir + '/%s/%04d/%02d/%02d'%(varName,Year,Mon,Day)
        if varName in ['S1.ABp103-117.GMI.Tc']:
            lvarPath = sorted(glob.glob(varDir + '/*.npy'))
        else:
            lvarPath = sorted(glob.glob(varDir + '/%s.1.*.npy'%(varNameSimple)))
    
        for varPath in lvarPath:
            avar = np.load(varPath)
            if ndim ==2:
                ny,nx = avar.shape
                avar  = avar.reshape(ny*nx,1) 
            elif ndim ==3:
                ny,nx,nz = avar.shape
                avar  = avar.reshape(ny*nx, nz)
            avarTmp = ma.masked_equal(avar, miss)
            asumTmp  = avar.sum(axis=0)
            asum2Tmp = np.square(avar).sum(axis=0)
            anumTmp  = ma.count(avar, axis=0)
            try:
                asum  = asum  + asumTmp
                asum2 = asum2 + asum2Tmp
                anum  = anum  + anumTmp
            except TypeError:
                asum  = asumTmp
                asum2 = asum2Tmp
                anum  = anumTmp

    #-- save and calc ---
    amean = asum / anum
    astd  = np.sqrt( (asum2 - 2*amean*asum + anum*np.square(amean))/(anum-1) )
    meanPath = outDir + '/mean.%s.npy'%(varName)
    stdPath  = outDir + '/std.%s.npy'%(varName)
    np.save(meanPath, amean)
'''
