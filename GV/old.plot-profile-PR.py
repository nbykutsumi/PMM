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
lprtype = ['wetevent']
dlthpr = {'wetevent':[0.1,999],'all':[-0.1,999],'light':[-0.1,2],'mod':[2,10], 'heavy':[10,50],'extreme':[50,9999]}

miss    = -9999.

#------------------------------------------------
# container
a2prof = deque([])
a1esurf= deque([])

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

        if not os.path.exists(profPath):
            print 'no file',profPath
            print 'skip', domain, YM
            continue

        aprof  = np.load(profPath)*0.01
        angv   = np.load(ngvPath)
        aesurf = np.load(eSurfPath)
        aidxtmp= zeros(len(angv))

        # accumulation arrays for mask
        agroundBin= zeros(len(angv))
        asateAcc = gv_fsub.mean_slice_negativemask(aprof.T, agroundBin, nh)

        #-- mask when ngv < minNum
        amskN    = ma.masked_less(angv, minNum).mask

        #-- mask when satellite is zero
        amskzero = ma.masked_equal(asateAcc,0).mask

        #-- mask when at least one of the data is missing or low quality
        amskmiss = ma.masked_less(asateAcc,0).mask

        #-- mask when sate has too large value
        amskL    = ma.masked_greater(asateAcc, 100).mask

        #-- overlay masks
        amsk     = amskN + amskzero + amskmiss + amskL

        aidxTmp  = arange(len(amsk)).astype(int32)
        aidxTmp  = ma.masked_where(amsk, aidxTmp).compressed()

        aprofTmp   = aprof[aidxTmp,:]
        aprofTmp   = ma.masked_less(aprofTmp,0)
        aesurfTmp  = aesurf[aidxTmp]

        a2prof.append(aprofTmp)
        a1esurf.extend(aesurfTmp) 

a2prof = concatenate(a2prof, axis=0)
a1esurf= array(a1esurf)

dmean = {}
dstd  = {}
for prtype in lprtype:
    thmin,thmax = dlthpr[prtype]
    a1idx = arange(len(a1esurf))
    amskP = ma.masked_outside(a1esurf, thmin, thmax).mask
    a1idx = ma.masked_where(amskP, a1idx).compressed()
    a2profTmp = ma.masked_less(a2prof[a1idx], 0)

    if len(a1idx)==0:
        a1mean = [nan]*nh
        a1std  = [nan]*nh
    else:
        a1mean  = a2profTmp.mean(axis=0)
        a1std   = a2profTmp.std(axis=0)

    dmean[prtype] = a1mean
    dstd [prtype] = a1std

# figure  -------------
fig = plt.figure(figsize=(4,3))
ax  = fig.add_axes([0.1, 0.1, 0.8, 0.7])

a1y = (arange(nh) +1)*0.25

ax.plot(dmean['light'],a1y,'-',color='b', label='light')
ax.plot(dstd['light'], a1y,'--',color='b', label='light')

ax.plot(dmean['mod'],a1y,'-',color='g', label='mod')
ax.plot(dstd['mod'], a1y,'--',color='g', label='mod')

ax.plot(dmean['heavy'],a1y,'-',color='r', label='heavy')    
ax.plot(dstd['heavy'],a1y,'--',color='r', label='heavy')    

# legend
plt.legend()

stitle  = '%s %s'%(prdName, 'prcp')
stitle  = stitle + '\n' + '%04d.%02d-%04d.%02d'%(iYM[0],iYM[1],eYM[0],eYM[1])
plt.title(stitle)        
figDir  = '/work/a01/utsumi/GPMGV/fig'
figPath = figDir + '/plt.prof.%s.%.1fkm.minNum.%d.%s.png'%(prdName, thdist,minNum,'prcp')
plt.savefig(figPath)
print figPath
plt.clf()

