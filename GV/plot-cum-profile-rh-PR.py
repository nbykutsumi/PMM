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

calc = True
#calc = False
iYM = [2005,4]
eYM = [2014,10]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
print lYM

#thdist = 2.5 # km
thdist = 5 # km
minNum = 3
prdName = 'L2A25'
nh = 40
gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()

dgName = gv.ret_ddomYM2gName()

ldomain = gv.domains
#ldomain = ['FLORIDA-STJ','FLORIDA-SFL-N','N.Carolina-IPHEx_Duke','N.Carolina-IPHEx_NASA','KWAJALEIN-KWA']
#ldomain = ['FLORIDA-STJ']

offset_bef = 15  # 'bef' should be identical to what is used in mk_match.py
#offset_aft = 45

#lprtype = ['all','mod','heavy','extreme']
#lprtype = ['light','mod','heavy']
lprtype = ['all']
dlthpr = {'wetevent':[0.1,999],'all':[-0.1,999],'light':[-0.1,2],'mod':[2,10], 'heavy':[10,50],'extreme':[50,9999]}

lrhtype = ['all','dry','hum']
dlthrh  = {'all':[-0.1,9999],'dry':[-0.1,0.7-0.0001],'hum':[0.7,9999]}

miss    = -9999.

#------------------------------------------------
# container
a2prof = deque([])
a1esurf= deque([])
a1rh   = deque([])

for domain in ldomain:
    for YM in lYM:
        Year, Mon = YM
        if (domain,Year,Mon) not in dgName.keys():
            print 'no obs',domain,Year,Mon
            continue
        # load data
        srcbaseDir = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.%s'%(prdName)
        srcDir     = srcbaseDir + '/%.1fkm/%s/%04d%02d'%(thdist, domain, Year,Mon)

        profPath      = srcDir  + '/p_prof.npy'
        eSurfPath     = srcDir + '/p_eSurf.npy'
        gvPath        = srcDir + '/p_gvprcp.npy'
        ngvPath       = srcDir + '/p_ngv.npy'
        rhPath        = srcDir + '/p_rh.npy'



        if not os.path.exists(profPath):
            print 'no file',profPath
            print 'skip', domain, YM
            continue

        aprof    = np.load(profPath)*0.01
        agv      = np.load(gvPath)
        angv     = np.load(ngvPath)
        aesurf   = np.load(eSurfPath)
        arh      = np.load(rhPath)
        aidxtmp  = zeros(len(angv))

        # accumulation arrays for mask
        agroundBin= zeros(len(angv))
        a1sateAll = gv_fsub.mean_slice_negativemask(aprof.T, agroundBin, nh)

        a1idx15   = ones(agv.shape[0])*15
        a1gvAll   = ma.masked_less(agv[:,15:],0).mean(axis=1)


        #-- mask when both satellite and gv are zero
        amsk1    = ma.masked_equal(a1sateAll,0).mask
        amsk2    = ma.masked_equal(a1gvAll,0).mask
        #amsk3    = ma.masked_equal(aesurf,0  ).mask
        amskzero = amsk1 * amsk2


        #-- mask when ngv < minNum
        amskN    = ma.masked_less(angv, minNum).mask

        #-- overlay masks
        amsk     = amskzero + amskN


        aidxTmp  = arange(len(amsk)).astype(int32)
        aidxTmp  = ma.masked_where(amsk, aidxTmp).compressed()

        aprofTmp   = aprof[aidxTmp,:]
        aesurfTmp  = aesurf[aidxTmp]
        arhTmp     = arh[aidxTmp]

        a2prof.append(aprofTmp)
        a1esurf.extend(aesurfTmp) 
        a1rh.extend(arhTmp)


a2prof = concatenate(a2prof, axis=0)
a1esurf= array(a1esurf)
a1rh   = array(a1rh)

for prtype in lprtype:
    dmean = {}
    dstd  = {}
    desurf= {}

    for rhtype in lrhtype:
        thmin,thmax = dlthpr[prtype]
        amskP = ma.masked_outside(a1esurf, thmin, thmax).mask

        thrh0,thrh1 = dlthrh[rhtype]
        amskRH= ma.masked_outside(a1rh, thrh0, thrh1).mask

        amsk  = amskP + amskRH
        a1idx = arange(len(a1esurf))
        a1idx = ma.masked_where(amsk, a1idx).compressed()
        a2profTmp = a2prof[a1idx]
        a1esurfTmp= a1esurf[a1idx]
        a2datTmp  = concatenate([a1esurfTmp.reshape(-1,1), a2profTmp],axis=1)
 
        if len(a1idx)==0:
            a1mean = [nan]*(nh+1)
        else:
            a2cum  = empty(a2datTmp.shape)
            for i in range(nh+1):
                a1tmp      = ma.masked_less(a2datTmp[:,:i+1],0).mean(axis=1).filled(miss)
                a2cum[:,i] = a1tmp

        a1mean = ma.masked_less(a2cum,0).mean(axis=0)
        a1std  = ma.masked_less(a2cum,0).std(axis=0)
        esurf  = a1esurfTmp.mean()
 
        dmean[rhtype] = a1mean
        dstd [rhtype] = a1std
        desurf[rhtype] = esurf

    # figure  -------------
    fig = plt.figure(figsize=(4,3))
    ax  = fig.add_axes([0.1, 0.1, 0.8, 0.7])
    
    a1y = arange(nh+1)*0.25
   
    ax.plot(dmean['all'],a1y,'-',color='k', label='all')
    ax.plot(dstd['all'], a1y,'--',color='k')
    ax.plot(desurf['all'], 0.1, 'v', color='k')
    
    ax.plot(dmean['dry'],a1y,'-',color='r', label='dry')
    ax.plot(dstd['dry'], a1y,'--',color='r')
    ax.plot(desurf['dry'], 0.1, 'v', color='r')
    
    ax.plot(dmean['hum'],a1y,'-',color='b', label='hum')    
    ax.plot(dstd['hum'],a1y,'--',color='b')    
    ax.plot(desurf['hum'], 0.1, 'v', color='b')
   
 
    # legend
    plt.legend()
    
    stitle  = '%s average from surface'%(prdName)
    stitle  = stitle + '\n' + '%04d.%02d-%04d.%02d'%(iYM[0],iYM[1],eYM[0],eYM[1])
    plt.title(stitle)        
    figDir  = '/work/a01/utsumi/GPMGV/fig'
    figPath = figDir + '/plt.aveprof.%s.%.1fkm.minNum.%d.%s.RH.png'%(prdName, thdist,minNum,'prcp')
    plt.savefig(figPath)
    print figPath
    plt.clf()


