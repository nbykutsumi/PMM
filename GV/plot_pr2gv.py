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


iYM = [2014,4]
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

dt = 0
prdName = 'L2A25'
#prdName = '2A-CLIM'

gvprcpAll = deque([])
eSurfAll  = deque([])
for domain in ldomain:
    gvprcp = deque([])
    eSurf  = deque([])

    for YM in lYM:
        Year,Mon = YM
        srcDir = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.%s/%.1fkm/%s/%04d%02d'%(prdName, thdist, domain, Year,Mon)
    
        if not os.path.exists(srcDir):
            print 'no directory: ',srcDir
            continue
        print domain, YM
        ngvTmp    = np.load(srcDir + '/p_ngv.npy')
        #gvprcpTmp  = np.load(srcDir + '/p_gvprcp.npy')[:,dt+15]
        a2gvprcpTmp  = np.load(srcDir + '/p_gvprcp.npy')[:,dt+15-2:dt+15+2]
        gvprcpTmp = abs(a2gvprcpTmp).mean(axis=1)
        eSurfTmp   = np.load(srcDir + '/p_eSurf.npy')
        if prdName =='2A-CLIM':
            qFlag    = np.load(srcDir + '/p_qFlag.npy')


        # mask where ngv < minNum
        mskN  = ma.masked_less(ngvTmp,minNum).mask
        print 'mskN',len(mskN) 
        # mask where both are zero
        mskA  = ma.masked_equal(gvprcpTmp,0).mask
        mskB  = ma.masked_equal(eSurfTmp,0).mask
        msk0  = mskA * mskB

        # mask where gv is low quality
        mskG  = ma.masked_less(gvprcpTmp,0).mask

        # mask where satellite is low quality
        if prdName == '2A-CLIM':
            mskQ  = ma.masked_greater(qFlag, 0).mask
            msk   = msk0 + mskG + mskQ + mskN
        else:
            msk   = msk0 + mskG

        gvprcpTmp = ma.masked_where(msk, gvprcpTmp).compressed()
        eSurfTmp  = ma.masked_where(msk, eSurfTmp).compressed()

        gvprcp.extend(gvprcpTmp)
        eSurf.extend(eSurfTmp) 
        gvprcpAll.extend(gvprcpTmp)
        eSurfAll.extend(eSurfTmp) 



# figure for all data
gvprcp = np.array(gvprcpAll)
eSurf  = np.array(eSurfAll)
 

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






