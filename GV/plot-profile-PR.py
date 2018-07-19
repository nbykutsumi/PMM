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
#iYM = [2014,4]
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

lraintype = ['all','strat','conv']
dlthrtype = {'all':[-999,999],'strat':[100,170],'conv':[200,297]}


miss    = -9999.

#------------------------------------------------
def ret_a2cumave_negativemask(a2dat):
    miss  = -9999.
    a2cum = empty(a2dat.shape)
    for i in range(a2dat.shape[1]):
        a1tmp = ma.masked_less(a2dat[:,:i+1],0).mean(axis=1).filled(miss)
        a2cum[:,i] = a1tmp
    return a2cum 
#------------------------------------------------
# container
a2prof = deque([])
a1esurf= deque([])
a1raintype= deque([])

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
        rainTypePath  = srcDir + '/p_rainType.npy'



        if not os.path.exists(profPath):
            print 'no file',profPath
            print 'skip', domain, YM
            continue

        aprof    = np.load(profPath)*0.01
        agv      = abs(np.load(gvPath))
        angv     = np.load(ngvPath)
        aesurf   = np.load(eSurfPath)
        araintype= np.load(rainTypePath)
        aidxtmp  = zeros(len(angv))

        # accumulation arrays for mask
        agroundBin= zeros(len(angv))
        a1sateAll = gv_fsub.mean_slice_negativemask(aprof.T, agroundBin, nh)

        a1idx15   = ones(agv.shape[0])*15
        a1gvAll   = ma.masked_less(agv[:,15:15+30],0).mean(axis=1)


        #-- mask when both satellite and gv are zero or miss
        amsk1    = ma.masked_less_equal(a1sateAll,0).mask
        amsk2    = ma.masked_less_equal(a1gvAll,0).mask
        amsk3    = ma.masked_less_equal(aesurf,0).mask
        amskzero = amsk1 * amsk2 * amsk3



        #-- mask when ngv < minNum
        amskN    = ma.masked_less(angv, minNum).mask

        #-- overlay masks
        amsk     = amskzero + amskN

        aidxTmp  = arange(len(amsk)).astype(int32)
        aidxTmp  = ma.masked_where(amsk, aidxTmp).compressed()

        aprofTmp   = aprof[aidxTmp,:]
        aesurfTmp  = aesurf[aidxTmp]
        araintypeTmp=araintype[aidxTmp]

        a2prof.append(aprofTmp)
        a1esurf.extend(aesurfTmp) 
        a1raintype.extend(araintypeTmp)


a2prof = concatenate(a2prof, axis=0)
a1esurf= array(a1esurf)
a1raintype=array(a1raintype)

print a2prof.shape

for prtype in lprtype:
    dmean = {}
    dstd  = {}
    dcount= {}
    desurf= {}

    dcumave= {}
    dcumstd= {}

    for raintype in lraintype:
        thmin,thmax = dlthpr[prtype]
        amskP = ma.masked_outside(a1esurf, thmin, thmax).mask

        thrtype0,thrtype1 = dlthrtype[raintype]
        amskT = ma.masked_outside(a1raintype, thrtype0, thrtype1).mask

        amsk  = amskP + amskT
        a1idx = arange(len(a1esurf))
        a1idx = ma.masked_where(amsk, a1idx).compressed()
        a2profTmp = ma.masked_less(a2prof[a1idx], 0)
        a1esurfTmp= ma.masked_less(a1esurf[a1idx],0)   


        #-- plain profile --
        if len(a1idx)==0:
            a1mean = [nan]*nh
            a1std  = [nan]*nh
        else:

            a1mean  = a2profTmp.mean(axis=0)
            a1std   = a2profTmp.std(axis=0)
            a1count  = a2profTmp.count(axis=0)
            esurf   = a1esurfTmp.mean()
 
        dmean[raintype] = a1mean
        dstd [raintype] = a1std
        dcount[raintype] = a1count
        desurf[raintype] = esurf

        #-- cum-average profile --
        if len(a1idx)==0:
            a1cumave = [nan]*nh
            a1cumstd = [nan]*nh
        else:
            a2joinprof= concatenate([a1esurfTmp.reshape(-1,1), a2profTmp], axis=1)
            a2cum    = ret_a2cumave_negativemask(a2joinprof)
            a1cumave = ma.masked_less(a2cum,0).mean(axis=0)
            a1cumstd = ma.masked_less(a2cum,0).std(axis=0)
    
        dcumave[raintype] = a1cumave
        dcumstd[raintype] = a1cumstd


    # figure  plain profile -------------
    fig = plt.figure(figsize=(4,3))
    ax  = fig.add_axes([0.1, 0.1, 0.8, 0.7])
    
    a1y = (arange(nh) +1)*0.25
   
    imin = 4 
    ax.plot(dmean['all'][imin:],a1y[imin:],'-',color='k', label='all')
    ax.plot(dstd['all'][imin:], a1y[imin:],'--',color='k')
    ax.plot(desurf['all'], 0.1, 'v', color='k')
    

    
    ax.plot(dmean['conv'][imin:],a1y[imin:],'-',color='r', label='conv')
    ax.plot(dstd['conv'][imin:], a1y[imin:],'--',color='r')
    ax.plot(desurf['conv'], 0.1, 'v', color='r')
    
    ax.plot(dmean['strat'][imin:],a1y[imin:],'-',color='b', label='strat')    
    ax.plot(dstd['strat'][imin:],a1y[imin:],'--',color='b')    
    ax.plot(desurf['strat'], 0.1, 'v', color='b')
   
 
    # legend
    plt.legend()
    
    stitle  = '%s %s'%(prdName, 'prcp')
    stitle  = stitle + '\n' + '%04d.%02d-%04d.%02d'%(iYM[0],iYM[1],eYM[0],eYM[1])
    plt.title(stitle)        
    figDir  = '/work/a01/utsumi/GPMGV/fig'
    figPath = figDir + '/plt.prof.%s.%.1fkm.minNum.%d.%s.raintype.png'%(prdName, thdist,minNum,'prcp')
    plt.savefig(figPath)
    print figPath
    plt.clf()



    # figure cumave -------------
    fig = plt.figure(figsize=(4,3))
    ax  = fig.add_axes([0.1, 0.1, 0.8, 0.7])
    
    a1y = arange(nh+1)*0.25
    
    ax.plot(dcumave['all'],a1y,'-',color='k', label='all')
    ax.plot(dcumstd['all'], a1y,'--',color='k')
    ax.plot(desurf['all'], 0.1, 'v', color='k')
    
    ax.plot(dcumave['conv'],a1y,'-',color='r', label='conv')
    ax.plot(dcumstd['conv'], a1y,'--',color='r')
    ax.plot(desurf['conv'], 0.1, 'v', color='r')
    
    ax.plot(dcumave['strat'],a1y,'-',color='b', label='strat')
    ax.plot(dcumstd['strat'],a1y,'--',color='b')
    ax.plot(desurf['strat'], 0.1, 'v', color='b')

    
    # legend
    plt.legend()
    
    stitle  = '%s average from surface'%(prdName)
    stitle  = stitle + '\n' + '%04d.%02d-%04d.%02d'%(iYM[0],iYM[1],eYM[0],eYM[1])
    plt.title(stitle)
    figDir  = '/work/a01/utsumi/GPMGV/fig'
    figPath = figDir + '/plt.aveprof.%s.%.1fkm.minNum.%d.%s.raintype.png'%(prdName, thdist,minNum,'prcp')
    plt.savefig(figPath)
    print figPath
    plt.clf()



    # figure: count  -------------
    fig = plt.figure(figsize=(4,3))
    ax  = fig.add_axes([0.1, 0.1, 0.8, 0.7])
    
    a1y = (arange(nh) +1)*0.25
    
    ax.plot(dcount['all'],a1y,'-',color='k', label='all')
    
    ax.plot(dcount['conv'],a1y,'-',color='r', label='conv')
    
    ax.plot(dcount['strat'],a1y,'-',color='b', label='strat')    
    
    # legend
    plt.legend()
    
    stitle  = '%s %s'%(prdName, 'count')
    stitle  = stitle + '\n' + '%04d.%02d-%04d.%02d'%(iYM[0],iYM[1],eYM[0],eYM[1])
    plt.title(stitle)        
    figDir  = '/work/a01/utsumi/GPMGV/fig'
    figPath = figDir + '/plt.countprof.%s.%.1fkm.minNum.%d.%s.raintype.png'%(prdName, thdist,minNum,'prcp')
    plt.savefig(figPath)
    print figPath
    plt.clf()

    print dcount['all']
