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
iYM = [2005,5]
eYM = [2014,10]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
lYM = lYM[::-1]
print lYM

thdist  = 7.5
minNum  = 5
prdName = '2A-CLIM'
gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()

dgName = gv.ret_ddomYM2gName()

ldomain = gv.domains
#ldomain = ['FLORIDA-STJ','FLORIDA-SFL-N','N.Carolina-IPHEx_Duke','N.Carolina-IPHEx_NASA','KWAJALEIN-KWA']
#ldomain = ['FLORIDA-STJ']

offset_bef = 15  # 'bef' should be identical to that in mk_match.py
offset_aft = 30


nt = offset_aft + offset_bef +1
nh = 18  # 0.5 - 10km, at most.
#nh = 5

#lprtype = ['heavy']
lprtype = ['all','mod','heavy','extreme']
dlthpr = {'all':[-0.1,9999],'mod':[-0.1,10], 'heavy':[10,50],'extreme':[50,9999]}
ldattype = ['rain','cc','bias','brat','rmse','num','gv']

coef_def = 4.17 * 1000/(60*60)  # g/m3 * m/s *1000/(60*60) = mm/h

#-- cluserProfile --

obaseDir = "/home/utsumi/mnt/wellshare/GPMGV/%s"%(prdName)
clusterProfPath = obaseDir + "/clusterProf.npy"
a4clusterProf = np.load(clusterProfPath)

#------------------------------------------------

da1rain = {prtype: empty([nt]) for prtype in lprtype}
da1gv   = {prtype: empty([nt]) for prtype in lprtype}
da1cc   = {prtype: empty([nt]) for prtype in lprtype}
da1rmse = {prtype: empty([nt]) for prtype in lprtype}
da1bias = {prtype: empty([nt]) for prtype in lprtype}
da1brat = {prtype: empty([nt]) for prtype in lprtype}
da1num  = {prtype: empty([nt],int32) for prtype in lprtype}

for idt, dt in enumerate(range(-offset_bef, offset_aft+1)):
    if calc != True:
        continue

    daesurf    = deque([])
    dagv       = deque([])
    dacoef     = deque([])

    for domain in ldomain:
        for YM in lYM:
            Year, Mon = YM
            if (domain,Year,Mon) not in dgName.keys():
                print 'no obs',domain,Year,Mon
                continue
            print domain, 'dt=',dt,YM
            # load data
            srcbaseDir = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.%s'%(prdName)
            srcDir     = srcbaseDir + '/%.1fkm/%s/%04d%02d'%(thdist, domain, Year,Mon)
    
            gvPath      = srcDir  + '/p_gvprcp.npy'
            qFlagPath   = srcDir  + '/p_qFlag.npy'
    
            eSurfPath   = srcDir + '/p_eSurf.npy'
            ngvPath     = srcDir + '/p_ngv.npy'

            if not os.path.exists(qFlagPath):
                print 'no file',qFlagPath
                print 'skip', domain, YM
                continue

            aqFlag = np.load(qFlagPath) 
            aesurf = np.load(eSurfPath)
            agv    = np.load(gvPath)[:,idt]
            angv   = np.load(ngvPath)

            #-- mask where ngv < minNum
            amskN    = ma.masked_less(angv,minNum).mask

            #-- mask where qFlag >0 ---
            amskQ    = ma.masked_greater(aqFlag,0).mask

            #-- mask when both satellite and gv are zero
            amsk1    = ma.masked_equal(aesurf,0).mask
            amsk2    = ma.masked_equal(agv,0).mask
            amskzero = amsk1 * amsk2

            #-- mask when at least one of the data is missing or low quality
            amsk3    = ma.masked_less(aesurf,0).mask
            amsk4    = ma.masked_less(agv,0).mask
            amskmiss = amsk3 + amsk4

            #-- overlay masks
            amsk     = amskN + amskQ + amskzero + amskmiss

            aidxTmp  = arange(len(amsk)).astype(int32)
            aidxTmp  = ma.masked_where(amsk, aidxTmp).compressed()
        
            aesurfTmp    = aesurf[aidxTmp]
            agvTmp       = agv[aidxTmp]
            #------------------
            daesurf.extend(aesurfTmp)
            dagv.extend(agvTmp)
 
    #-------------
    aesurf    = array(daesurf)  # mm/h

    agv   = array(dagv)

    #- eSurf vs gv --
    for prtype in lprtype:
        # mask with eSurf --
        thmin, thmax = dlthpr[prtype]
        amsk      = ma.masked_outside(aesurf, thmin, thmax).mask
        #amsk      = ma.masked_outside(agv, thmin, thmax).mask
        agvTmp    = ma.masked_where(amsk, agv   ).compressed()
        aesurfTmp = ma.masked_where(amsk, aesurf).compressed()
    
         
        cc     = np.corrcoef(aesurfTmp, agvTmp)[0,1]
        rmse   = np.sqrt( ((aesurfTmp- agvTmp)**2).mean() )
       
        bias   = (aesurfTmp - agvTmp).mean()
        brat   = (aesurfTmp - agvTmp).mean() / agvTmp.mean()*100
        num    = len(aesurfTmp)
        rain   = aesurfTmp.mean()
        gv     = agvTmp.mean()

        da1cc  [prtype][idt] = cc
        da1rmse[prtype][idt] = rmse
        da1bias[prtype][idt] = bias
        da1brat[prtype][idt] = brat
        da1num [prtype][idt] = num       
        da1rain[prtype][idt] = rain
        da1gv  [prtype][idt] = gv

        #----- test --------------
        if ((dt ==0)&(prtype=='all')):
            figDir = '/work/a01/utsumi/GPMGV/fig'
            figPath= figDir + '/temp.%s.png'%(prtype)
            #plt.plot(agvTmp, aesurfTmp, 'o')
            plt.loglog(agvTmp, aesurfTmp, 'o')
            
            plt.ylim([0.01,100])
            plt.xlim([0.01,100])
            plt.savefig(figPath)
            print figPath
            plt.clf()


#--- save data ----
for prtype in lprtype:
    if calc != True:
        continue

    outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-eSurf-%s/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])
    util.mk_dir(outDir)

    for dattype in ldattype:  
        datPath= outDir + '/dt-eSurf.%s.%s.%s.npy'%(prdName, dattype, prtype)

        if dattype=='cc':
            a1dat = da1cc[prtype]
        elif dattype=='rmse':
            a1dat = da1rmse[prtype]
        elif dattype=='bias':
            a1dat = da1bias[prtype]
        elif dattype=='brat':
            a1dat = da1brat[prtype]
        elif dattype=='num':
            a1dat = da1num[prtype]
        elif dattype=='rain':
            a1dat = da1rain[prtype]
        elif dattype=='gv':
            a1dat = da1gv[prtype]


        np.save(datPath, a1dat)
        print datPath


#--- figure -------

for prtype in lprtype:
    outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-eSurf-%s/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])

    for dattype in ldattype:

        datPath=outDir + '/dt-eSurf.%s.%s.%s.npy'%(prdName, dattype, prtype)
        a1dat  = np.load(datPath)
    
        fig = plt.figure(figsize=(4,3))
        a1t  = range(-offset_bef, offset_aft+1)


        print '-'*50
        ax = fig.add_axes([0.15, 0.1, 0.8, 0.8])
        a1y = a1dat
        a1y = ma.masked_invalid(a1y)

        ax.plot(a1t, a1y, '-', zorder=10)


        # ylim
        if dattype in ['cc']:
            ymin = -1.2
            ymax = 1.2
            plt.ylim([ymin,ymax])

        elif dattype in ['rain','gv','rmse','num']:
            ymin = 0
            ymax = ma.masked_invalid(a1dat).max()*1.1
            plt.ylim([ymin,ymax])

        elif dattype in ['bias','brat']:
            ymax = abs(ma.masked_invalid(a1dat)).max()*1.1
            ymin = -ymax
            plt.ylim([ymin,ymax])


        #- zero line ---
        if dattype in ['cc','bias','brat']:
            plt.plot([-10,50],[0,0],'--',color='k', linewidth=0.5)

        plt.xlim([-6,31])


        stitle  = '%s %s'%(prtype, dattype)
        plt.title(stitle)        
        figDir  = '/work/a01/utsumi/GPMGV/fig'
        figPath = figDir + '/plot.dt-eSurf.%s.%s.%s.png'%(prdName, dattype, prtype)
        plt.savefig(figPath)
        print figPath
        plt.clf()
    



