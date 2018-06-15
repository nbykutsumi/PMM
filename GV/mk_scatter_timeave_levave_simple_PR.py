import matplotlib
matplotlib.use('Agg')
from numpy import *
from datetime import datetime, timedelta
from collections import deque
from gv_fsub import *
import GPMGV
import numpy as np
import myfunc.util as util
import matplotlib.pyplot as plt
import sys, os
from matplotlib import rcParams, cycler


calc = True
#calc = False
iYM = [2005,4]
eYM = [2014,10]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
print lYM

#thdist = 2.5
thdist = 5.0
minNum = 3

prdName = 'L2A25'

gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()

dgName = gv.ret_ddomYM2gName()

ldomain = gv.domains
#ldomain = ['FLORIDA-STJ','FLORIDA-SFL-N','N.Carolina-IPHEx_Duke','N.Carolina-IPHEx_NASA','KWAJALEIN-KWA']
#ldomain = ['FLORIDA-STJ']

offset_bef = 15  # 'bef' should be identical to that in mk_match.py
#offset_aft = 30

#ldt  = [1, 5,10,15,20,25,30,35]
#ldh  = [1, 5,10,15,20,25]

ldt = [25]
ldh = [15]

basepr = 'gv'
#basepr = 'sate'

nt = len(ldt)
nh = len(ldh)
#nh = 5

lprtype = ['heavy']
#lprtype = ['all','light','mod','heavy']
dlthpr = {'all':[-0.1,999],'light':[-0.1,2],'mod':[2,10], 'heavy':[10,50],'extreme':[50,9999]}
ldattype = ['rain','cc','bias','brat','rmse','num','gv']

da2rain = {prtype: empty([nh,nt]) for prtype in lprtype}
da2gv   = {prtype: empty([nh,nt]) for prtype in lprtype}
da2cc   = {prtype: empty([nh,nt]) for prtype in lprtype}
da2rmse = {prtype: empty([nh,nt]) for prtype in lprtype}
da2bias = {prtype: empty([nh,nt]) for prtype in lprtype}
da2brat = {prtype: empty([nh,nt]) for prtype in lprtype}
da2num  = {prtype: empty([nh,nt],int32) for prtype in lprtype}

# for eSurf
de1rain = {prtype: empty([nt]) for prtype in lprtype}
de1gv   = {prtype: empty([nt]) for prtype in lprtype}
de1cc   = {prtype: empty([nt]) for prtype in lprtype}
de1rmse = {prtype: empty([nt]) for prtype in lprtype}
de1bias = {prtype: empty([nt]) for prtype in lprtype}
de1brat = {prtype: empty([nt]) for prtype in lprtype}
de1num  = {prtype: empty([nt],int32) for prtype in lprtype}



a2prof     = deque([]) 
a1esurf    = deque([]) 
a1nsurf    = deque([]) 
a2gv       = deque([]) 
a1nsurfbin = deque([]) 

for domain in ldomain:
    if calc ==False: continue

    for YM in lYM:
        Year, Mon = YM
        if (domain,Year,Mon) not in dgName.keys():
            print 'no obs',domain,Year,Mon
            continue
        # load data
        srcbaseDir = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.L2A25'
        srcDir     = srcbaseDir + '/%.1fkm/%s/%04d%02d'%(thdist,domain, Year,Mon)

        profPath  = srcDir  + '/p_prof.npy'
        eSurfPath = srcDir  + '/p_eSurf.npy'
        nSurfPath = srcDir  + '/p_nSurf.npy'
        gvPath    = srcDir  + '/p_gvprcp.npy'
        ngvPath   = srcDir  + '/p_ngv.npy'
        nSurfBinPath = srcDir  + '/p_nSurfBin.npy'
        

        if not os.path.exists(profPath):
            print 'no file',profPath
            print 'skip', domain, YM
            continue

        aprof  = np.load(profPath)
        aesurf = np.load(eSurfPath)
        ansurf = np.load(nSurfPath)
        agv    = abs(np.load(gvPath))
        angv   = np.load(ngvPath)
        ansurfbin= np.load(nSurfBinPath)

        asate  = aprof

        #-- mean array for masks
        a1idx0    = zeros(aprof.shape[0])
        a1sateAll = gv_fsub.mean_slice_negativemask(aprof.T, a1idx0, ldh[-1])
        a1idx15   = ones(agv.shape[0])*15
        a1gvAll   = gv_fsub.mean_slice_negativemask(agv.T, a1idx15, ldt[-1])
        a1gvAll   = ma.masked_less(agv[:,15:],0).mean(axis=1)
        a1sateMin = asate.min(axis=1)


        #-- mask when both satellite and gv are zero
        amsk1    = ma.masked_equal(a1sateAll,0).mask
        amsk2    = ma.masked_equal(a1gvAll,0).mask
        amskzero = amsk1 * amsk2

        ##-- mask when sate has  missing data in the profile
        #amskmiss = ma.masked_less(a1sateMin,0).mask

        #-- mask when ng  < minNum
        amskN    = ma.masked_less(angv, minNum).mask

        #-- overlay masks
        #amsk     = amskzero + amskmiss amskN
        amsk     = amskzero + amskN

        aidxTmp  = arange(len(amsk)).astype(int32)
        aidxTmp  = ma.masked_where(amsk, aidxTmp).compressed()
    
        asateTmp     = asate[aidxTmp] 
        aesurfTmp    = aesurf[aidxTmp]
        ansurfTmp    = ansurf[aidxTmp]
        agvTmp       = agv[aidxTmp]
        ansurfbinTmp = ansurfbin[aidxTmp]
        #------------------
    
        a2gv.append(agvTmp)
        a2prof.append(asateTmp)
        a1esurf.extend(aesurfTmp)
        a1nsurf.extend(ansurfTmp)
        a1nsurfbin.extend(ansurfbinTmp)

if calc ==True:
    a1esurf= array(a1esurf)
    a1nsurf= array(a1nsurf)*0.01 #mm/h
    a1nsurfbin= array(a1nsurfbin)
    a2gv   = concatenate(a2gv,axis=0)
    a2prof = concatenate(a2prof,axis=0) *0.01
    a2joinprof= concatenate([a1esurf.reshape(-1,1) ,a2prof], axis=1)

#----
dprof = {}
dgv   = {}

ih = -1
for dh in [-99] + ldh: 
    if calc == False: continue

    if dh == -99:
        a1prof = a1esurf
    else:
        ih     = ih + 1
        a1idx0 = zeros(a2prof.shape[0])
        a1prof = gv_fsub.mean_slice_negativemask(a2joinprof.T, a1idx0, dh+1)

    for it,dt in enumerate(ldt):
        a1idx15 = ones(a2gv.shape[0])*15
        a1gv   = gv_fsub.mean_slice_negativemask(a2gv.T, a1idx15, dt)

        # mask events when at least one of prof or gv is missing
        amsk1 = ma.masked_less(a1prof,0).mask
        #amsk2 = ma.masked_less(a1gv,0).mask
        #amskMiss = amsk1 + amsk2
        amskMiss = amsk1

        for prtype in lprtype:

            # mask with nSurf --
            thmin, thmax = dlthpr[prtype]
            if basepr =='sate':
                amskP      = ma.masked_outside(a1nsurf, thmin, thmax).mask
            elif basepr=='gv':
                amskP      = ma.masked_outside(a1gv, thmin, thmax).mask

            # overlay masks
            if (amskMiss is False) or (amskMiss is True):
                amskMiss = array([False]*len(a1prof))
            if (amskP is False) or (amskMiss is True):
                amskP    = array([False]*len(a1prof))

            amsk  = amskMiss + amskP

            a1profTmp = ma.masked_where(amsk, a1prof ).compressed()
            a1gvTmp   = ma.masked_where(amsk, a1gv   ).compressed()

            dprof[prtype,dt,dh] = a1profTmp
            dgv  [prtype,dt,dh] = a1gvTmp
 

for prtype in lprtype:
    #--- Figure -------
    fig  = plt.figure(figsize=(3,3))
    ax   = fig.add_axes([0.1,0.1,0.8,0.8])

    x = dgv  [prtype,25,-99]
    y = dprof[prtype,25,-99]
    ax.scatter(x, y, label='-99')

    x = dgv  [prtype,25,15]
    y = dprof[prtype,25,15]
    ax.scatter(x, y, label='15')

    # legend
    plt.legend()

    # 1:1 line
    ax.plot([0,999],[0,999],'--',color='k')

    plt.xlim([0,50])
    plt.ylim([0,50])

    stitle  = '%s'%(prtype)
    plt.title(stitle)
    figDir  = '/work/a01/utsumi/GPMGV/fig'
    figPath = figDir + '/scatter.simple.%s.%.1fkm.minNum.%d.%s.png'%(prdName, thdist,minNum,prtype)
    plt.savefig(figPath)
    print figPath
    plt.clf()
    
    




