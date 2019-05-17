import numpy as np
from numpy import *
import myfunc.util as util
from datetime import datetime, timedelta
import glob
import os,sys
import calendar
import mkdbfunc
import random

calcmon = True
#calcmon = False
calcall = True
#calcall = False
sampletype = 'all'
#sampletype = 'precip'

varName = 'nltb'
iYM = [2017,1]
eYM = [2017,12]
#eYM = [2017,1]
lYM = util.ret_lYM(iYM,eYM)
outDir= '/work/hk01/utsumi/PMM/TPCDB/PC_COEF'
ncomb = 116  
#npc_use=4
npc_use=3
#NPCHIST = 25
NPCHIST = 22
#NPCHIST = 10
#-- percentile ranges --
wp  = int(100/(NPCHIST-2))
lp0 = array([0,1])
lp1 = arange(wp,100,wp)
lp2 = array([99,100])
lp  = r_[lp0, lp1, lp2]
print lp
#sys.exit()
#-- Read mean and std --
meanPath = outDir + '/mean.nltb.npy'
stdPath  = outDir + '/std.nltb.npy'

a1mean = np.load(meanPath)
a1std  = np.load(stdPath)
if ncomb != len(a1mean):
    print 'check ncomb'
    print 'ncomb=',ncomb
    print 'len(mean file)=',len(a1mean)
    sys.exit()

#-- Read PC coefficients --
egvecPath = outDir + '/egvec.nltb.npy'
a2egvec   = np.load(egvecPath)  # (n-th, nComb)

#-- Read orbit list  ---
listDir = '/work/hk01/utsumi/PMM/TPCDB/list'
lorbit  = []
for (Year,Mon) in lYM:
    listPath = listDir + '/list.1C.V05.%04d%02d.csv'%(Year,Mon)
    f=open(listPath,'r'); lines = f.readlines(); f.close()
    for line in lines:
        line = map(int, line.split(','))
        #print line
        lorbit.append(line)
#-- Sampling fron orbit list ---
random.seed(0)
samplerate = 0.1
#samplerate = 0.3
nsample = int(len(lorbit)*samplerate)
lorbitTmp = random.sample(lorbit, nsample)
lorbitTmp.sort(key=lambda x:x[0])
#-- Iititalize -----------------
a2pc  = None

#-------------------------------
for orbinfo in lorbitTmp:
    print orbinfo[:3]
    oid,Year,Mon,Day,itime,etime = orbinfo

    matchBaseDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A'
    dprDir = matchBaseDir + '/S1.ABp083-137.Ku.V06A.9ave.precipRate/%04d/%02d/%02d'%(Year,Mon,Day)
    dprPath = dprDir + '/precipRate.%06d.npy'%(oid)
    if not os.path.exists(dprPath): continue
    ##-- Read KuPR HDF file -----
    if sampletype=='precip':
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
    if sampletype=='precip':
        if a3dpr.max()<=0: continue

    a1flagtc = ma.masked_inside(a3tc,50,350).mask.all(axis=2).flatten()

    if sampletype=='precip':
        a1flagpr = ma.masked_greater(a3dpr,0).mask.any(axis=2).flatten()
        a1flag = a1flagtc * a1flagpr
    elif sampletype=='all':
        a1flag = a1flagtc
    else:
        print 'check sampletype',sampletype
        sys.exit()

    a2tc = a3tc.reshape(-1,ntc)
    a2tc = a2tc[a1flag]

    #-- Make non-linear combination --
    nlen  = a2tc.shape[0]
    adat  = mkdbfunc.mk_nonlin_comb(a2tc.reshape(-1,1,ntc))
    adat  = adat.reshape(nlen,-1)

    #-- Normalize --------------------
    adat  = (adat - a1mean)/a1std

    #-- Convert to PC ----------------
    apc   = np.dot(adat, a2egvec[:npc_use,:].T)
    print apc.shape
    try: 
        a2pc = concatenate([a2pc,apc],axis=0)
    except ValueError:
        a2pc = apc


#-- Calc percentiles ----
#ltmp = []  # test
sout = ''
for ipc in range(npc_use):
    vmin = a2pc[:,ipc].min()
    vmax = a2pc[:,ipc].max()

    #lp = np.linspace(0,100,NPCHIST+1)
    lvpercentile = []
    for p in lp:
        vpercentile = np.percentile(a2pc[:,ipc], p)
        lvpercentile.append(vpercentile)


    sline = '%d\t%.5f\t%.5f'%(ipc, vmin, vmax)
    for i in range(len(lvpercentile)-1):
        bnd0 = lvpercentile[i]
        bnd1 = lvpercentile[i+1]
        sline = sline + '\t%.5f\t%.5f'%(bnd0,bnd1)

    sout = sout + sline + '\n'

    #ltmp.append(lvpercentile) # test

#-- Save range file ---
rangePath = outDir + '/PC_MIN_MAX_%d_no_overlap_%scases.txt'%(NPCHIST,sampletype)
f=open(rangePath,'w'); f.write(sout); f.close()
print rangePath

#
##--- test ---------------
#a2bnd = array(ltmp)
#print a2bnd.shape
#a2hist= []
#for i in range(4):
#    a1bnd = a2bnd[i,:]
#    a1hist = np.histogram(a2pc[:,i], bins=a1bnd)[0]
#    a2hist.append(a1hist)
#
#a2hist=array(a2hist)
#


