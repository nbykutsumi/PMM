import matplotlib
matplotlib.use('Agg')
from numpy import *
import numpy as np
import myfunc.IO.GPM as GPM
from gv_fsub import *
import myfunc.util as util
import GPMGV
from collections import deque
import sys, os
import matplotlib.pyplot as plt


iYM = [2005,4]
eYM = [2014,10]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
lYM = lYM[::-1]
print lYM
thdist = 5
minNum = 3
gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()
ldomain = gv.domains
basepr = 'gv'
dt = 0
prdName = 'L2A25'
#prdName = '2A-CLIM'

#lprtype = ['heavy','extreme','mod']
lprtype = ['all','light','mod','heavy','extreme']
dlthpr = {'all':[-0.1,999],'light':[-0.1,2],'mod':[2,10], 'heavy':[10,50],'extreme':[50,9999]}
ldattype = ['rain','cc','bias','brat','rmse','num','gv']


a2gvprcp = deque([])
a1eSurf  = deque([])
for domain in ldomain:
    for YM in lYM:
        Year,Mon = YM
        srcDir = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.%s/%.1fkm/%s/%04d%02d'%(prdName, thdist, domain, Year,Mon)
    
        if not os.path.exists(srcDir):
            print 'no directory: ',srcDir
            continue
        print domain, YM
        ngvTmp    = np.load(srcDir + '/p_ngv.npy')
        a2gvprcpTmp  = abs(np.load(srcDir + '/p_gvprcp.npy'))
        eSurfTmp   = np.load(srcDir + '/p_eSurf.npy')
        if prdName =='2A-CLIM':
            qFlag    = np.load(srcDir + '/p_qFlag.npy')


        # mask where ngv < minNum
        mskN  = ma.masked_less(ngvTmp,minNum).mask

        # mask where both are zero
        gvprcpTmp = abs(a2gvprcpTmp[:,15:15:5]).mean(axis=1)
        mskA  = ma.masked_equal(gvprcpTmp,0).mask
        mskB  = ma.masked_equal(eSurfTmp,0).mask
        msk0  = mskA * mskB

        # mask where satellite is low quality
        mskQ  = ma.masked_less(eSurfTmp,0).mask 

        msk   = msk0 + mskN + mskQ

        a1idxTmp  = range(len(eSurfTmp))
        a1idxTmp  = ma.masked_where(msk, a1idxTmp).compressed()

        gvprcpTmp = a2gvprcpTmp[a1idxTmp,:]
        eSurfTmp  = eSurfTmp[a1idxTmp]

        a2gvprcp.extend(gvprcpTmp)
        a1eSurf.extend(eSurfTmp) 

a1eSurf  = np.array(a1eSurf)
a2gvprcp = np.array(a2gvprcp)
#a1gvprcp = a2gvprcp[:,15:15+5].mean(axis=1)
a1gvprcp = a2gvprcp[:,15]


# calc bias
for prtype in lprtype:
    vmin, vmax = dlthpr[prtype]
    a1idxMsk = ma.masked_outside(a1gvprcp, vmin, vmax).mask
    a1idxTmp = range(len(a1eSurf))
    a1idxTmp = ma.masked_where(a1idxMsk, a1idxTmp).compressed()

    a1eSurfTmp = a1eSurf[a1idxTmp]
    a1gvprcpTmp= a1gvprcp[a1idxTmp]

    esurf = a1eSurfTmp.mean(axis=0)
    gv    = a1gvprcpTmp.mean(axis=0)
    if (esurf is nan) or (gv is nan):
        rat = nan
    else:
        rat   = esurf / gv

    print prtype, esurf, gv, rat, 1./rat

'''
figDir = '/work/a01/utsumi/GPMGV/fig'
figPath= figDir + '/scatter.%s_gv.%.1fkm.ndom.%d.png'%(prdName,thdist, minNum)
#plt.loglog(gvprcp,eSurf,'o', color='k')
plt.plot(gvprcp,eSurf,'o', color='k')
stitle = '%s dist=%.1fkm minNum=%d'%(prdName, thdist, minNum)
stitle = stitle + '\n' + '%04d.%02d-%04d.%02d'%(iYM[0],iYM[1],eYM[0],eYM[1])
plt.title(stitle)
vmax  = 100
vmin  = 0.1
plt.ylim([vmin,vmax])
plt.xlim([vmin,vmax])

#-- 1:1 line -
plt.plot([vmin,vmax],[vmin,vmax],'-',color='k')

plt.savefig(figPath)
print figPath
plt.clf()
print len(eSurf)
'''


