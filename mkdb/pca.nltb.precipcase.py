import numpy as np
from numpy import *
import myfunc.util as util
from datetime import datetime, timedelta
import glob
import os,sys
import calendar
import mkdbfunc
from sklearn.decomposition import IncrementalPCA

varName = 'nltb'
iYM = [2017,1]
eYM = [2017,12]
#eYM = [2017,1]
lYM = util.ret_lYM(iYM,eYM)
iYear,iMon = iYM
eYear,eMon = eYM
outDir= '/work/hk01/utsumi/PMM/TPCDB/PC_COEF/%04d%02d-%04d%02d'%(iYear,iMon,eYear,eMon)
ncomb = 116  
#-- Read mean and std --
coefDir= '/work/hk01/utsumi/PMM/TPCDB/PC_COEF'
meanPath = coefDir + '/mean.nltb.npy'
stdPath  = coefDir + '/std.nltb.npy'

a1mean = np.load(meanPath)
a1std  = np.load(stdPath)
if ncomb != len(a1mean):
    print 'check ncomb'
    print 'ncomb=',ncomb
    print 'len(mean file)=',len(a1mean)
    sys.exit()

#-----------------------
inc_pca = IncrementalPCA(n_components=ncomb)  # n_components = 195

for (Year,Mon) in lYM:
    eDay   = calendar.monthrange(Year,Mon)[1]
    iDTime = datetime(Year,Mon,1)
    eDTime = datetime(Year,Mon,eDay)
    #iDTime = datetime(2017,7,12)
    #eDTime = datetime(2017,7,12)

    dDTime = timedelta(days=1)
    lDTime = util.ret_lDTime(iDTime,eDTime,dDTime)
    matchBaseDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A'
   
    #------------------
    #lDTime = lDTime[:2]  # test  
    for DTime in lDTime:
        print DTime
        Year,Mon,Day = DTime.timetuple()[:3]

        #-- Make KuPR HDF file list ---
        dprDir = matchBaseDir + '/S1.ABp083-137.Ku.V06A.9ave.precipRate/%04d/%02d/%02d'%(Year,Mon,Day)
        #ssearch= dprDir + '/precipRate.016173.npy'
        ssearch= dprDir + '/precipRate.??????.npy'

        ldprPath= glob.glob(ssearch)
        if len(ldprPath)==0:continue

        #ldprPath = ldprPath[:4] # test
        for dprPath in ldprPath:
            oid = int(dprPath.split('.')[-2]) 

            #-- Read KuPR HDF file -----
            a3dpr = np.load(dprPath)[:,103-83:103-83+15,:]

            #-- Read Tb--------
            tcDir1 = matchBaseDir + '/S1.ABp103-117.GMI.Tc/%04d/%02d/%02d'%(Year,Mon,Day)
            tcDir2 = matchBaseDir + '/S1.ABp103-117.GMI.TcS2/%04d/%02d/%02d'%(Year,Mon,Day)
            
            tcPath1= tcDir1 + '/Tc.%06d.npy'%(oid)
            tcPath2= tcDir2 + '/TcS2.1.%06d.npy'%(oid)
            a3tc1 = np.load(tcPath1)
            a3tc2 = np.load(tcPath2)
            a3tc  = concatenate([a3tc1,a3tc2],axis=2)
            ntc   = a3tc.shape[2]

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

            #-- Normalize --------------------
            adat  = (adat - a1mean)/a1std

            print adat.shape
            #***** Incremental PCA *********** 
            inc_pca.partial_fit(adat) 

#** Make eigen vectors and values --
egvec = inc_pca.components_      # (n-th, nComb)
egval = inc_pca.explained_variance_
varratio= inc_pca.explained_variance_ratio_

#** Save ----------
coefDir = outDir 
util.mk_dir(coefDir)
egvecPath = coefDir + '/egvec.%s.npy'%(varName)
egvalPath = coefDir + '/egval.%s.npy'%(varName)
varratioPath=coefDir+ '/varratio.%s.npy'%(varName)
np.save(egvecPath, egvec)
np.save(egvalPath, egval)
np.save(varratioPath, varratio)
print egvecPath
#print evec
