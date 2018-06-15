import matplotlib
matplotlib.use('Agg')
from numpy import *
import numpy as np
import myfunc.util as util
import GPMGV
import sys, os
import matplotlib.pyplot as plt

iYM = [2014,8]
eYM = [2014,10]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
lYM = lYM[::-1]
print lYM
gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()
ldomain = gv.domains

thdist = 2.5
minNum = 2
dt     = 0
nlev   = 20

prdName = 'L2A25'
for domain in ldomain:
    for YM in lYM:
        Year,Mon = YM
        srcDir = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.%s/%.1fkm/%s/%04d%02d'%(prdName, thdist, domain, Year,Mon)

        if not os.path.exists(srcDir):
            #print 'no directory: ',srcDir
            continue
        print domain, YM
        ngvTmp    = np.load(srcDir + '/p_ngv.npy')
        gvprcpTmp = np.load(srcDir + '/p_gvprcp.npy')[:,dt+15:dt+15+nlev]
        profTmp   = np.load(srcDir + '/p_prof.npy')[:,:nlev]
        a1nsurfbinTmp= np.load(srcDir + '/p_nSurfBin.npy')
        eSurfTmp  = np.load(srcDir + '/p_eSurf.npy')

        if prdName =='2A-CLIM':
            qFlag    = np.load(srcDir + '/p_qFlag.npy')

        # mask where ngv <minNum
        mskN  = ma.masked_less(ngvTmp,minNum).mask

        # mask where both are zero or negative
        a1idxtmp = arange(profTmp.shape[0])
        mskA  = ma.masked_less_equal(gvprcpTmp[:,0],0).mask
        mskB  = ma.masked_less_equal(profTmp[a1idxtmp,a1nsurfbinTmp],0).mask
        msk0  = mskA * mskB

        # mask where gv is low quality
        mskG  = ma.masked_less(gvprcpTmp[:,0],0).mask

        # mask where satellite is low quality
        if prdName == '2A-CLIM':
            mskQ  = ma.masked_greater(qFlag, 0).mask

        # overlay masks
        if prdName == '2A-CLIM':
            msk   = msk0 + mskG + mskQ + mskN
        else:
            msk   = msk0 + mskG + mskN

        a1idx = arange(len(ngvTmp))
        a1idx = ma.masked_where(msk, a1idx).compressed()
        gvprcpTmp = gvprcpTmp[a1idx,:]
        profTmp   = profTmp[a1idx,:]
        eSurfTmp  = eSurfTmp[a1idx]

        print gvprcpTmp.shape
        print profTmp.shape

        for i in range(profTmp.shape[0]):
            a1gv   = abs(gvprcpTmp[i,:])
            a1prof = ma.masked_less(profTmp[i,:],0) * 0.01

            fig = plt.figure(figsize=(3,3))
            ax1 = fig.add_axes([0.2,0.1,0.7,0.3])
            ax2 = fig.add_axes([0.2,0.5,0.7,0.3])
            ax1.plot(a1gv,'-')
            ax2.plot(a1prof,'-')

            # eSurf
            ax2.plot(0,eSurfTmp[i], 'o', color='r')

            ax1.set_xlim([0,nlev])
            ax2.set_xlim([0,nlev])



            plt.title('%s %04d-%02d %d'%(domain,Year,Mon,i))
            figDir  = '/work/a01/utsumi/GPMGV/fig'
            figPath = figDir + '/plot.comp.gv-prof.%s.%04d%02d.%d.png'%(domain,Year,Mon,i)
            plt.savefig(figPath)
            print figPath


