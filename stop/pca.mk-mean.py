from numpy import *
import numpy as np
from datetime import datetime,timedelta
import myfunc.util as util
import random


iDTime = datetime(2017,1,1)
eDTime = datetime(2017,12,31)
dDTime = timedelta(days=1)
lDTime = util.ret_lDTime(iDTime,eDTime,dDTime)

#ldydx = [[dy,dx] for dy in [-1,0,1] for dx in [-3,-2,-1,0,1,2,3]]
ldydx = [[0,0]]

lisurf = range(1,15+1)  # surface type index
#lisurf = [1]  # surface type index
ntc1 = 9
ntc2 = 4

for (dy,dx) in ldydx:
    for isurf in lisurf:
        asx1 = zeros(ntc1).astype(float64)
        asx2 = zeros(ntc2).astype(float64)
        asxx1= zeros(ntc1).astype(float64)
        asxx2= zeros(ntc2).astype(float64)
        anum1= zeros(ntc1).astype(int64)
        anum2= zeros(ntc2).astype(int64)

        for DTime in lDTime:
            print isurf,DTime
            Year,Mon,Day = DTime.timetuple()[:3]
            baseDir = '/work/hk01/utsumi/PMM/stop/data/Tc'
            srcDir  = baseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
            srcPath1 = srcDir + '/Tc1.%ddy.%ddx.%02dsurf.npy'%(dy,dx,isurf)
            srcPath2 = srcDir + '/Tc2.%ddy.%ddx.%02dsurf.npy'%(dy,dx,isurf)
   
            try: 
                a2tmp1 = np.load(srcPath1).reshape(-1,ntc1)
                a2tmp2 = np.load(srcPath2).reshape(-1,ntc2)
            except IOError:
                continue

            a2tmp1 = ma.masked_outside(a2tmp1,50,350)
            a2tmp2 = ma.masked_outside(a2tmp2,50,350)
            asxtmp1 = a2tmp1.sum(axis=0)
            asxtmp2 = a2tmp2.sum(axis=0)

            asxxtmp1= (a2tmp1*a2tmp1).sum(axis=0)
            asxxtmp2= (a2tmp2*a2tmp2).sum(axis=0)

            anumtmp1= a2tmp1.count(axis=0) 
            anumtmp2= a2tmp2.count(axis=0) 

            asx1 = asx1 + asxtmp1
            asx2 = asx2 + asxtmp2
            asxx1= asxx1+ asxxtmp1
            asxx2= asxx2+ asxxtmp2
            anum1=anum1+anumtmp1
            anum2=anum2+anumtmp2

        #-- Mean & Std ----
        amean1 = asx1 / anum1
        amean2 = asx2 / anum2
        astd1  = np.sqrt( (asxx1 - 2*amean1*asx1 + anum1*np.square(amean1))/(anum1-1) )
        astd2  = np.sqrt( (asxx2 - 2*amean2*asx2 + anum2*np.square(amean2))/(anum2-1) )

        amean1 = amean1.filled(-9999.)
        amean2 = amean2.filled(-9999.)
        astd1  = astd1.filled(-9999.)
        astd2  = astd2.filled(-9999.)

        #-- Save ----------
        outDir = '/work/hk01/utsumi/PMM/stop/data/coef'
        meanPath1= outDir + '/mean1.%ddy.%ddx.%02dsurf.npy'%(dy,dx,isurf)
        meanPath2= outDir + '/mean2.%ddy.%ddx.%02dsurf.npy'%(dy,dx,isurf)
        stdPath1 = outDir + '/std1.%ddy.%ddx.%02dsurf.npy'%(dy,dx,isurf)
        stdPath2 = outDir + '/std2.%ddy.%ddx.%02dsurf.npy'%(dy,dx,isurf)

        np.save(meanPath1, amean1)
        np.save(meanPath2, amean2)
        np.save(stdPath1, astd1)
        np.save(stdPath2, astd2)

        print meanPath1
