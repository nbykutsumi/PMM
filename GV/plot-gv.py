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
nh = 34
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
# container
a1raintype= deque([])
a2gv   = deque([])

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
        ngvPath       = srcDir + '/p_ngv.npy'
        rainTypePath  = srcDir + '/p_rainType.npy'
        gvPath        = srcDir + '/p_gvprcp.npy'


        if not os.path.exists(profPath):
            print 'no file',profPath
            print 'skip', domain, YM
            continue

        aprof    = np.load(profPath)*0.01
        agv      = np.load(gvPath)
        angv     = np.load(ngvPath)
        aesurf   = np.load(eSurfPath)
        araintype= np.load(rainTypePath)
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

        agvTmp      =agv[aidxTmp,:]
        araintypeTmp=araintype[aidxTmp]

        a2gv.append(agvTmp)
        a1raintype.extend(araintypeTmp)


a2gv       = concatenate(a2gv, axis=0)
a1raintype = array(a1raintype)

dmean = {}
dstd  = {}
dcount= {}
dtimeave= {}

print a2gv.shape

for raintype in lraintype:

    thrtype0,thrtype1 = dlthrtype[raintype]
    amskT = ma.masked_outside(a1raintype, thrtype0, thrtype1).mask

    amsk  = amskT
    a1idx = arange(len(a1raintype))
    a1idx = ma.masked_where(amsk, a1idx).compressed()
    a2gvTmp = ma.masked_less(a2gv[a1idx,:], 0)

    if len(a1idx)==0:
        a1mean = [nan]*31
        a1std  = [nan]*31
    else:
        a1mean  = a2gvTmp.mean(axis=0)[15:15+31]
        a1std   = a2gvTmp.std(axis=0)[15:15+31]
        a1count  = a2gvTmp.count(axis=0)[15:15+31]
        a2timeave= empty([a2gvTmp.shape[0],31])

        for i in range(31):
            a1tmp = ma.masked_less(a2gvTmp[:,15:15+i+1],0).mean(axis=1).filled(miss)

            a2timeave[:,i] = a1tmp
        a1timeave = ma.masked_less(a2timeave, 0).mean(axis=0)

    dmean[raintype] = a1mean
    dstd [raintype] = a1std
    dcount[raintype] = a1count
    dtimeave[raintype] = a1timeave

# figure  -------------
for raintype in lraintype:
    fig = plt.figure(figsize=(4,3))
    ax  = fig.add_axes([0.2, 0.1, 0.75, 0.7])
    
    
    ax.plot(arange(31),dmean[raintype],'-',color='k', label='rain(mean)')
    ax.plot(arange(31),dstd[raintype], '--',color='k',label='rain(std)')
    ax.plot(arange(31),dtimeave[raintype], '-', color='r',label='rain(time ave)')
    
    # legend
    plt.legend()
   
    plt.ylim([0,12])
 
    stitle  = 'gauge rain %s'%(raintype)
    stitle  = stitle + '\n' + '%04d.%02d-%04d.%02d'%(iYM[0],iYM[1],eYM[0],eYM[1])
    plt.title(stitle)        
    figDir  = '/work/a01/utsumi/GPMGV/fig'
    figPath = figDir + '/plt.gv.%s.%s.%.1fkm.minNum.%d.%s.png'%(prdName, raintype, thdist,minNum,'prcp')
    plt.savefig(figPath)
    print figPath
    plt.clf()
    
# figure: count  -------------
fig = plt.figure(figsize=(4,3))
ax  = fig.add_axes([0.2, 0.1, 0.75, 0.7])


ax.plot(range(31),dcount['all'],'-',color='k', label='all')

ax.plot(range(31),dcount['conv'],'-',color='r', label='conv')

ax.plot(range(31),dcount['strat'],'-',color='b', label='strat')    

# legend
plt.legend()

plt.ylim(ymin=0)
stitle  = '%s %s'%(prdName, 'count')
stitle  = stitle + '\n' + '%04d.%02d-%04d.%02d'%(iYM[0],iYM[1],eYM[0],eYM[1])
plt.title(stitle)        
figDir  = '/work/a01/utsumi/GPMGV/fig'
figPath = figDir + '/plt.countgv.%s.%.1fkm.minNum.%d.%s.png'%(prdName, thdist,minNum,'prcp')
plt.savefig(figPath)
print figPath
plt.clf()

print dcount['all']
