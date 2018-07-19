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
import scipy.stats

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
#corrFlag= 'CR'
corrFlag= 'NO'

prdName = 'L2A25'

gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()

dgName = gv.ret_ddomYM2gName()

ldomain = gv.domains
#ldomain = ['FLORIDA-STJ','FLORIDA-SFL-N','N.Carolina-IPHEx_Duke','N.Carolina-IPHEx_NASA','KWAJALEIN-KWA']
#ldomain = ['FLORIDA-STJ']

offset_bef = 15  # 'bef' should be identical to that in mk_match.py
#offset_aft = 30

ldt  = [1, 15,30]
ldh  = [10,20,30]
basepr = 'gv'
#basepr = 'sate'

nt = len(ldt)
nh = len(ldh)
#nh = 5


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


        ##-- mask when both satellite and gv are zero
        #amsk1    = ma.masked_equal(a1sateAll,0).mask
        #amsk2    = ma.masked_equal(a1gvAll,0).mask
        #amskzero = amsk1 * amsk2

        ##-- mask when sate has  missing data in the profile
        #amskmiss = ma.masked_less(a1sateMin,0).mask

        #-- mask when ng  < minNum
        amskN    = ma.masked_less(angv, minNum).mask

        #-- overlay masks
        #amsk     = amskzero + amskN
        amsk     = amskN

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

#-- bias correction factor --
a1gv    = a2gv[:,15]
amsk    = range(len(a1gv))
amsk    = ma.masked_where(a1gv<0, amsk)
amsk    = ma.masked_where(a1esurf<0, amsk)
amsk    = amsk.mask
bfactor = ma.masked_where(amsk, a1gv).mean() / ma.masked_where(amsk, a1esurf).mean()

#--- averaging time loop ---
for it,dt in enumerate(ldt):
    a1idx15 = ones(a2gv.shape[0])*15
    a1gv   = gv_fsub.mean_slice_negativemask(a2gv.T, a1idx15, dt)
    
    #--- vertical layer loop ---
    ih = -1
   
    dprof   = {}
    dprofCR = {}
    for dh in [-99] + ldh: 
        if calc == False: continue
    
        if dh == -99:
            a1prof = a1esurf
        else:
            ih     = ih + 1
            a1idx0 = zeros(a2prof.shape[0])
            a1prof = gv_fsub.mean_slice_negativemask(a2joinprof.T, a1idx0, dh+1)
    
    
        # mask events when sate is missing
        a1profTmp = ma.masked_less(a1prof,0).compressed()
   
        dprof[dh]   = a1profTmp 
        dprofCR[dh] = a1profTmp * bfactor

        print dt, dh, a1gv.mean(), dprof[dh].mean(), dprofCR[dh].mean(), bfactor
         
    #-- Figure -------
    fig = plt.figure(figsize=(6,4))
    ax  = fig.add_axes([0.15, 0.1, 0.8, 0.78])


    #-- x ----
    xBnd = np.linspace(0,20,40)
    x    = 0.5*(xBnd[:-1]+xBnd[1:])
    # gv 
    #kde = scipy.stats.gaussian_kde(a1gv, bw_method='scott')
    #ax.plot(x, kde(x), '-', color='k', label='gauge')

    y, abin = np.histogram( a1gv, xBnd, density=True)
    ax.plot(x, y, '-', color='k', label='gauge')


    # sate
    cmap = plt.get_cmap('coolwarm')

    for ih,dh in enumerate([-99] + ldh):
        y,   abin = np.histogram( dprof[dh], xBnd, density=True)
        yCR, abin = np.histogram( dprofCR[dh],xBnd, density=True)

        mycm  = cmap(float(ih)/(len(ldh)+1))
        if dh ==-99:
            slabel = 'eSurf'
            ax.plot(x, y,   '-', color=mycm, label=slabel)
        else:
            slabel ='%.1fkm'%(dh*0.25)
            ax.plot(x, y,   '-', color=mycm, label=slabel)
            ax.plot(x, yCR, '--',color=mycm, label=slabel+' CR')



    plt.ylim([0, 0.02])

    # legend
    plt.legend()



    stitle  = 'basepr=%s dt=%dh'%(basepr, dt)
    stitle  = stitle + '\n' + '%04d.%02d-%04d.%02d'%(iYM[0],iYM[1],eYM[0],eYM[1])
    plt.title(stitle)
    figDir  = '/work/a01/utsumi/GPMGV/fig'
    figPath = figDir + '/pdf.nt-nlev.dt.%dh.%s.%.1fkm.minNum.%d.base.%s.png'%(dt, prdName, thdist,minNum, basepr)
    plt.savefig(figPath)
    print figPath
    plt.clf()






