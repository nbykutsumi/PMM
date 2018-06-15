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
import scipy.stats

calc = True
#calc = False
iYM = [2014,4]
eYM = [2014,10]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
lYM = lYM[::-1]
print lYM
#thdist  = 7.5 # km
thdist  = 15 # km
minNum  = 5
#basepr = 'sate'
basepr = 'gv'

prdName = '2A-CLIM'
gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()

dgName = gv.ret_ddomYM2gName()

ldomain = gv.domains
#ldomain = ['FLORIDA-STJ','FLORIDA-SFL-N','N.Carolina-IPHEx_Duke','N.Carolina-IPHEx_NASA','KWAJALEIN-KWA']
#ldomain = ['FLORIDA-STJ']

offset_bef = 15  # 'bef' should be identical to that in mk_match.py
offset_aft = 30
ldt= range(-offset_bef, offset_aft+1)

nt = offset_aft + offset_bef +1
nh = 20  # 0.5 - 10km, at most.
#nh = 5

lprtype = ['all','light','mod','heavy']
dlthpr = {'all':[-0.1,999],'light':[-0.1,2],'mod':[2,10], 'heavy':[10,50],'extreme':[50,9999]}

ldattype = ['rain','cc','bias','brat','rmse','num','gv']


 
coef_def = 4.17 * 1000/(60*60)  #=1.158:  g/m3 * m/s *1000/(60*60) = mm/h

#-- cluserProfile --

clustbaseDir = "/home/utsumi/mnt/wellshare/GPMGV/%s"%(prdName)
clusterProfPath = clustbaseDir + "/clusterProf.npy"
a4clusterProf = np.load(clusterProfPath)

#------------------------------------------------
def ret_aprof(a4clusterProf, a1tIndex, a2profNum, a2profScale, a1groundBin):
    lspecies = [0,2,3,4]
    nh_in    = 28
    nh_out   = nh
    a3out = empty([len(lspecies),len(a1tIndex),nh_in]).astype(float32)

    for i,species in enumerate(lspecies):
        a1profNum  = a2profNum[:,species]
        a1profScale= a2profScale[:,species]
        a2prof = a1profScale.reshape(-1,1) * a4clusterProf[a1profNum-1,:,a1tIndex-1,species]
        a3out[i] = a2prof

    a2prof = a3out.sum(axis=0)
    a2out  = gv_fsub.extract_slice_clusterprof(a2prof.T, agroundBin, nh_out).T 

    return a2out

#------------------------------------------------

a2prof  = deque([])
a2gv    = deque([])
a1esurf = deque([])
a1coef  = deque([])

for domain in ldomain:
    if calc != True:
        continue

    for YM in lYM:
        Year, Mon = YM
        if (domain,Year,Mon) not in dgName.keys():
            print 'no obs',domain,Year,Mon
            continue
        print domain, YM
        # load data
        srcbaseDir = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.%s'%(prdName)
        srcDir     = srcbaseDir + '/%.1fkm/%s/%04d%02d'%(thdist, domain, Year,Mon)
    
        gvPath        = srcDir + '/p_gvprcp.npy'
        qFlagPath     = srcDir + '/p_qFlag.npy'
        profNumPath   = srcDir + '/p_profNum.npy'
        profScalePath = srcDir + '/p_profScale.npy'
        tIndexPath    = srcDir + '/p_tIndex.npy'
    
        eSurfPath     = srcDir + '/p_eSurf.npy'
        mlPrcpPath    = srcDir + '/p_mlPrecip.npy'
        groundBinPath = srcDir + '/p_groundBin.npy'
        ngvPath       = srcDir + '/p_ngv.npy' 


        if not os.path.exists(profNumPath):
            print 'no file',profNumPath
            print 'skip', domain, YM
            continue

        aqFlag     = np.load(qFlagPath) 
        aprofNum   = np.load(profNumPath)
        aprofScale = np.load(profScalePath)
        atIndex    = np.load(tIndexPath)
        agroundBin = np.load(groundBinPath)
        angv       = np.load(ngvPath)

        aprof     = ret_aprof(a4clusterProf, atIndex, aprofNum, aprofScale, agroundBin) 

        aesurf = np.load(eSurfPath)

        agv    = abs(np.load(gvPath))

        acoef  = ma.masked_invalid(aesurf/aprof[:,0]).filled(coef_def)

        #---- mask where ngv < minNum
        amskN  = ma.masked_less(angv, minNum).mask
 
        #---- mask where qFlag >0 -
        amskQ  = ma.masked_greater(aqFlag, 0).mask

        #-- mask when both satellite and gv are zero
        a1sateAll= ma.masked_less(aprof, 0).mean(axis=1)
        a1gvAll  = ma.masked_less(agv,   0).mean(axis=1)

        amsk1    = ma.masked_equal(a1sateAll, 0).mask
        amsk2    = ma.masked_equal(a1gvAll,   0).mask
        amskzero = amsk1 * amsk2       

        #-- overlay --
        print amskN.shape, amskQ.shape, amskzero.shape
        amsk     = amskN + amskQ + amskzero
 

        aidxTmp  = arange(len(amsk)).astype(int32)
        aidxTmp  = ma.masked_where(amsk, aidxTmp).compressed()
    
        aprofTmp     = aprof[aidxTmp] 
        agvTmp       = agv[aidxTmp]
        aesurfTmp    = aesurf[aidxTmp]
        acoefTmp     = acoef[aidxTmp]
        #------------------
    
        a2prof.append(aprofTmp)
        a2gv.append(agvTmp)
        a1esurf.extend(aesurfTmp)
        a1coef.extend(acoefTmp)

a2gv   = concatenate(a2gv,axis=0)
a2prof = concatenate(a2prof,axis=0)
a1esurf= array(a1esurf)
a1coef = array(a1coef)


for prtype in lprtype:
    # mask with precipitation intensity--
    thmin, thmax = dlthpr[prtype]
    if   basepr =='sate':
        amskP   = ma.masked_outside(a1esurf, thmin, thmax).mask
    elif basepr =='gv':
        a1gvNow = a2gv[:,15]
        amskP   = ma.masked_outside(a1gvNow, thmin, thmax).mask

    amsk = amskP
    a1idxTmp  = ma.masked_where(amsk, range(len(a1coef))).compressed()
    a1coefTmp = a1coef[a1idxTmp]
    a1vel     = a1coefTmp / (1000./(60*60))
    a2profTmp = a2prof[a1idxTmp,:] 

    #-- mask zero velocity -
    a1vel = ma.masked_equal(a1vel, 0).compressed()

    #-- figure -------
    fig = plt.figure(figsize=(4,8))
    ax  = fig.add_axes([0.15, 0.1, 0.8, 0.8])

    hist, binBnd = np.histogram(a1vel,bins=arange(0,150,1.0),density=False)
    x   = 0.5*(binBnd[:-1] + binBnd[1:])
    ax.plot(x, hist, '-')

    # title
    stitle = 'v[m/s] %.1fkm minNum=%d basepr=%s'%(thdist, minNum, basepr)
    stitle = stitle + '\n' + '%04d.%02d-%04d.%02d'%(iYM[0],iYM[1],eYM[0],eYM[1])

    figDir  = '/work/a01/utsumi/GPMGV/fig'
    figPath = figDir + '/histo.velocity.%s.%.1fkm.minNum.%d.base.%s.%s.png'%(prdName, thdist, minNum, basepr, prtype)
    plt.savefig(figPath)
    print figPath
    plt.clf()
    

