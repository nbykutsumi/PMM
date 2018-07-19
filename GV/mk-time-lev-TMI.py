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
#nh = 20  # 0.5 - 10km, at most.
nh = 13  # 0.5 - 10km, at most.
#nh = 5

#lprtype = ['all']
lprtype = ['all','light','mod','heavy']
dlthpr = {'all':[-0.1,999],'light':[-0.1,2],'mod':[2,10], 'heavy':[10,50],'extreme':[50,9999]}

ldattype = ['rain','cc','bias','brat','rmse','num','gv']


 
#coef_def = 4.17 * 1000/(60*60)  # g/m3 * m/s *1000/(60*60) = mm/h
coef_def = 0.0  # g/m3 * m/s *1000/(60*60) = mm/h

#-- cluserProfile --

clustbaseDir = "/home/utsumi/mnt/wellshare/GPMGV/%s"%(prdName)
clusterProfPath = clustbaseDir + "/clusterProf.npy"
a4clusterProf = np.load(clusterProfPath)

#------------------------------------------------
def ret_aprof(a4clusterProf, a1tIndex, a2profNum, a2profScale, a1groundBin, lspecies=[0,2,3,4]):
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

da2rain = {prtype: empty([nh,nt],float32) for prtype in lprtype}
da2gv   = {prtype: empty([nh,nt],float32) for prtype in lprtype}
da2cc   = {prtype: empty([nh,nt],float32) for prtype in lprtype}
da2rmse = {prtype: empty([nh,nt],float32) for prtype in lprtype}
da2bias = {prtype: empty([nh,nt],float32) for prtype in lprtype}
da2brat = {prtype: empty([nh,nt],float32) for prtype in lprtype}
da2num  = {prtype: empty([nh,nt],int32) for prtype in lprtype}

a2prof     = deque([])
a2gv       = deque([])
a1esurf    = deque([])
a1coef     = deque([])

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

if calc==True: 
    a2gv   = concatenate(a2gv,axis=0)
    a2prof = concatenate(a2prof,axis=0)
    a1esurf= array(a1esurf)
    a1coef = array(a1coef)


for prtype in lprtype:
    if calc != True: continue

    # mask with eSurf --
    thmin, thmax = dlthpr[prtype]
    if   basepr =='sate':
        amskP   = ma.masked_outside(a1esurf, thmin, thmax).mask
    elif basepr =='gv':
        a1gvNow = a2gv[:,15]
        amskP   = ma.masked_outside(a1gvNow, thmin, thmax).mask

    for ih in range(nh): 
        #-------------
        a1prof = a2prof[:,ih]   # g/m3

        #-- mask when sate is missing
        amskMiss  = ma.masked_less(a1prof, 0).mask

        #-- overlay masks --
        amsk = amskP + amskMiss

        a1idxTmp  = ma.masked_where(amsk, range(len(a1prof))).compressed()
        a1profTmp = a1prof[a1idxTmp] 
        a1esurfTmp= a1esurf[a1idxTmp]
        a1coefTmp = a1coef[a1idxTmp]
    
        a1sateTmp = a1profTmp * a1coefTmp

        #print '---',ih,'--------'
        #print 'prof',a1profTmp.max(), a1prof.max(), a2prof.max()
        #print 'coef',a1coefTmp.max()
        #print 'sate',a1sateTmp.min(), a1sateTmp.max()
        for it, dt in enumerate(ldt):
            a1gv   = a2gv[:,it] 
            a1gvTmp= a1gv[a1idxTmp]

            #print ih,it,'gvTmp',a1gvTmp.max(),a2gv.max()
            #-- calc indices --
            cc     = np.corrcoef(a1sateTmp, a1gvTmp)[0,1]
            rmse   = np.sqrt( ((a1sateTmp- a1gvTmp)**2).mean() )
           
            bias   = (a1sateTmp - a1gvTmp).mean()
            brat   = (a1sateTmp - a1gvTmp).mean() / a1gvTmp.mean()*100
            num    = len(a1profTmp)
            rain   = a1sateTmp.mean()
            gv     = a1gvTmp.mean()
    
            da2cc  [prtype][ih, it] = cc
            da2rmse[prtype][ih, it] = rmse
            da2bias[prtype][ih, it] = bias
            da2brat[prtype][ih, it] = brat
            da2num [prtype][ih, it] = num       
            da2rain[prtype][ih, it] = rain
            da2gv  [prtype][ih, it] = gv



#-- save data ----
for prtype in lprtype:
    if calc != True:
        continue
    
    outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-%s/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])
    util.mk_dir(outDir)
    
    for dattype in ldattype:  
        datPath= outDir + '/dt-lev.base.%s.%s.%s.%s.npy'%(basepr, prdName, dattype, prtype)
    
        if dattype=='cc':
            a2dat = da2cc[prtype]
        elif dattype=='rmse':
            a2dat = da2rmse[prtype]
        elif dattype=='bias':
            a2dat = da2bias[prtype]
        elif dattype=='brat':
            a2dat = da2brat[prtype]
        elif dattype=='num':
            a2dat = da2num[prtype]
        elif dattype=='rain':
            a2dat = da2rain[prtype]
        elif dattype=='gv':
            a2dat = da2gv[prtype]
    
        np.save(datPath, a2dat)
        print datPath

#-- figure -------
for prtype in lprtype:
    outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-%s/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])
    
    for dattype in ldattype:
    
        datPath=outDir + '/dt-lev.base.%s.%s.%s.%s.npy'%(basepr, prdName, dattype, prtype)
        a2dat  = np.load(datPath)
    
        fig = plt.figure(figsize=(4,8))
        ax  = fig.add_axes([0.15, 0.1, 0.8, 0.8])
    
        a1t  = range(-offset_bef, offset_aft+1)
    
        cmap = plt.cm.coolwarm
        rcParams['axes.prop_cycle'] = cycler(color=cmap(np.linspace(0,1,nh)))
        cmap = plt.get_cmap('coolwarm')
        for ih in range(nh):
            a1y = a2dat[ih,:]
            a1y = ma.masked_invalid(a1y)
    
            if ih==0:
                ax.plot(a1t, a1y, '--', c='k', zorder=10, label='%d'%(ih))
            else:
                ax.plot(a1t, a1y, '-', c=cmap(float(ih)/nh), label='%d'%(ih))
        # legend
        ax.legend()
    
        # ylim
        if dattype in ['cc']:
            ymin = 0.0
            ymax = 0.5
            plt.ylim([ymin,ymax])
    
        elif dattype in ['rain','gv','rmse','num']:
            ymin = 0
            ymax = ma.masked_invalid(a2dat).max()*1.1
            plt.ylim([ymin,ymax])
    
        elif dattype in ['bias','brat']:
            ymax = abs(ma.masked_invalid(a2dat)).max()*1.1
            ymin = -ymax
            plt.ylim([ymin,ymax])
    
    
        #- zero line ---
        if dattype in ['cc','bias','brat']:
            plt.plot([-10,50],[0,0],'--',color='k', linewidth=0.5)
    
        plt.xlim([-6,31])
    
    
        # title
        stitle  = '%s %.1fkm minNum=%d %s'%(prtype, thdist, minNum, dattype)
        stitle  = stitle + '\n' + 'base:%s %04d.%02d-%04d.%02d'%(basepr, iYM[0],iYM[1],eYM[0],eYM[1])
        plt.title(stitle)        


        figDir  = '/work/a01/utsumi/GPMGV/fig'
        figPath = figDir + '/plot.dt-lev.onebox.%s.%.1fkm.minNum.%d.base.%s.%s.%s.png'%(prdName, thdist, minNum, basepr, dattype, prtype)
        plt.savefig(figPath)
        print figPath
        plt.clf()
    




