import numpy as np
import os, sys
from numpy import ma
import myfunc.util as util
from numpy import *
import calendar, glob
from datetime import datetime, timedelta
import random

iYM = [2017,1]
eYM = [2017,12]
lYM = util.ret_lYM(iYM,eYM)

srcDir = '/tank/utsumi/validprof/maperror'
lprec = arange(0.25, 20+0.01, 0.5)
lbiaslev = [0,1,2,3]
#dbiasrange = {0:[0,1],1:[1,2],2:[2,5],3:[5,999]}
dbiasrange = {0:[0,1],1:[1,5],2:[5,10],3:[10,999]}
abiasrange = dbiasrange.values()

def taylor_index(a2ref, a2dat, miss):
    print a2ref.shape, a2dat.shape


    a2mask1 = ma.masked_less_equal(a2ref, miss).mask
    a2mask2 = ma.masked_less_equal(a2dat, miss).mask
    a2mask  = a2mask1 + a2mask2
    a2ref = ma.masked_where(a2mask, a2ref)
    a2dat = ma.masked_where(a2mask, a2dat)

    a1num = (~a2mask).sum(axis=1)
    a1stdref = a2ref.std(axis=1)
    a1stddat = a2dat.std(axis=1)
    a1mref   = a2ref.mean(axis=1).reshape(-1,1)
    a1mdat   = a2dat.mean(axis=1).reshape(-1,1)

    a1cov    = ((a2ref - a1mref)*(a2dat - a1mdat)).sum(axis=1)/a1num
    a1corr = a1cov / (a1stdref * a1stddat)
    corrmax= 1.0

    S = 4*(1.0+a1corr)**4 /((a1stdref/a1stddat + a1stddat/a1stdref)**2) / (1.0+corrmax)**4
    S = ma.masked_invalid(S).filled(miss)
    return S



da1sum = {}
da1sum2= {}
da1num = {}
for biaslev in lbiaslev:
    for sgn in [-1,1]:
        da1sum [biaslev,sgn]= np.zeros(len(lprec), float32)
        da1sum2[biaslev,sgn]= np.zeros(len(lprec), float32)
        da1num [biaslev,sgn]= np.zeros(len(lprec), int32)

for Year,Mon in lYM:
    eDay = calendar.monthrange(Year,Mon)[1]
    iDTime = datetime(Year,Mon,1)
    eDTime = datetime(Year,Mon,eDay)
    dDTime = timedelta(days=1)
    lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)
    #lDTime = lDTime[:2]  # test 
    for DTime in lDTime:
        Year,Mon,Day = DTime.timetuple()[:3]

        #if Day !=18: continue 

        srcDir = '/tank/utsumi/validprof/pair.gprof/%04d/%02d/%02d'%(Year,Mon,Day)
        ssearch= srcDir + '/Latitude.*.npy'
        llatPath= sort(glob.glob(ssearch))

        #-- random sampling -----
        if len(llatPath)>=3:
            llatPath = random.sample(llatPath, 3)
        else:
            continue
        #------------------------
        
        for latPath in llatPath:
            print latPath
            oid = int(latPath.split('.')[-2])
            a1lat = np.load(srcDir + '/Latitude.%06d.npy'%(oid))
            a1lon = np.load(srcDir + '/Longitude.%06d.npy'%(oid))

            a1precradTmp= np.load(srcDir+ '/precrad.%06d.npy'%(oid))
            a1precpmwTmp= np.load(srcDir+ '/precpmw.%06d.npy'%(oid))
            a2profradTmp= np.load(srcDir+ '/profrad.%06d.npy'%(oid))
            a2profpmwTmp= np.load(srcDir+ '/profpmw.%06d.npy'%(oid))
            a1qflag   = np.load(srcDir + '/qualityFlag.%06d.npy'%(oid))
            a1surftype= np.load(srcDir + '/surfaceTypeIndex.%06d.npy'%(oid))

            #--- Check quality flag and surface type --
            a1flagQ = ma.masked_equal(a1qflag,0).mask   # Good quality (flag=0)
            a1flagS1 = ma.masked_equal(a1surftype,1).mask   # Ocean
            a1flagS2 = ma.masked_inside(a1surftype,3,7).mask # Vegetation
            a1flagS3 = ma.masked_inside(a1surftype,12,13).mask # Standing water and rivers & Water/Coast boundary
            a1flagS  = a1flagS1 + a1flagS2 + a1flagS3
            a1flag   = a1flagQ * a1flagS 
            
            a1precpmwTmp = a1precpmwTmp[a1flag]
            a1precradTmp = a1precradTmp[a1flag]
            a2profpmwTmp = a2profpmwTmp[a1flag,:]
            a2profradTmp = a2profradTmp[a1flag,:] 

            #--- Mask missing data --
            a1precpmwTmp = ma.masked_less(a1precpmwTmp,0)
            a1precradTmp = ma.masked_less(a1precradTmp,0)

            #a2profpmwTmp = ma.masked_outside(a2profpmwTmp,0,30)
            a2profpmwTmp = ma.masked_less(a2profpmwTmp,0)
            a2profradTmp = ma.masked_less(a2profradTmp,0)


            if a1precpmwTmp.shape[0]==0:
                continue
            #--- Bias -------
            a1prbiasTmp = a1precpmwTmp - a1precradTmp

            #--- Taylor -----
            a1taylorTmp = ma.masked_less(taylor_index(a2profradTmp, a2profpmwTmp, miss=-9999.), 0)

            #--- Classify by bias --
            for biaslev in [0,1,2,3]:
                for sgn in [-1,1]:
                    ibias,ebias = dbiasrange[biaslev]
                    a1biasmask = ma.masked_outside(a1prbiasTmp,sgn*ibias,sgn*ebias).mask
    
                    a1sum = np.zeros(len(lprec), float32)
                    a1sum2= np.zeros(len(lprec), float32)
                    a1num = np.zeros(len(lprec), int32)
                    for iprec, prec in enumerate(lprec):
    
                        a1precmask = ma.masked_outside(a1precradTmp, prec-0.25, prec+0.25).mask
                        a1mask = a1biasmask + a1precmask
                        a1taylorTmp2 = ma.masked_where(a1mask, a1taylorTmp)
    
                        n = a1taylorTmp2.count()
                        if n==0:
                            s = 0
                        else:
                            s = a1taylorTmp2.sum() 
                        a1sum[iprec] = a1sum[iprec] + s
                        a1sum2[iprec]= a1sum2[iprec] + s**2
                        a1num[iprec] = a1num[iprec] + n
    
                    da1sum [biaslev,sgn] = da1sum [biaslev,sgn] + a1sum
                    da1sum2[biaslev,sgn] = da1sum2[biaslev,sgn] + a1sum2
                    da1num [biaslev,sgn] = da1num [biaslev,sgn] + a1num
    
                ##-- test ---
                #iprec = 5
                #cnt   = lprec[iprec]
                #print ''
                #print 'pr=%.1f-0.25 ~ %.1f+0.25'%(cnt,cnt)
                #print 'all count=',ma.masked_outside(a1precradTmp, cnt-0.25, cnt+0.25).count()
                #for biaslev in [0,1,2,3]:
                #    #print 'biaslev=',biaslev,'count=',da1num[biaslev][iprec],'mean=',da1sum[biaslev][iprec]/da1num[biaslev]/da1num[biaslev][iprec]
                #    s = da1sum[biaslev][iprec]
                #    n = da1num[biaslev][iprec]
                #    print 'biaslev=',biaslev,'sum=','%.1f'%s,'num=',n,'mean=','%.2f'%(s/n)
                ##sys.exit()

    #--- Save -----
    for biaslev in lbiaslev:
        for sgn in [-1,1]:
            sym = '%04d.%02d'%(Year,Mon)
            outDir  = '/tank/utsumi/validprof/pr.vs.metrics'
            sumPath = outDir + '/taylor.sum.bias=%dx%d.%s.npy'%(sgn,biaslev, sym)
            sum2Path= outDir + '/taylor.sum2.bias=%dx%d.%s.npy'%(sgn,biaslev, sym)
            numPath = outDir + '/taylor.num.bias=%dx%d.%s.npy'%(sgn,biaslev, sym)
           
            np.save(sumPath,  da1sum [biaslev,sgn])
            np.save(sum2Path, da1sum2[biaslev,sgn])
            np.save(numPath,  da1num [biaslev,sgn])
    
    prbinPath = outDir + '/taylor.by.bias.prbin.%s.npy'%(sym)
    biasrangePath = outDir + '/taylor.by.bias.biasrange.%s.npy'%(sym)
    np.save(prbinPath, lprec)
    np.save(biasrangePath, abiasrange) 
    print prbinPath


