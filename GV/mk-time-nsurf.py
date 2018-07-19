import matplotlib
matplotlib.use('Agg')
from numpy import *
from datetime import datetime, timedelta
from collections import deque
import GPMGV
import numpy as np
import myfunc.util as util
import matplotlib.pyplot as plt
import sys, os

calc = True
#calc = False
iYM = [2010,4]
eYM = [2014,10]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]

thdist = 5.0 # km
minNum = 3

gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()

dgName = gv.ret_ddomYM2gName()

ldomain = gv.domains
#ldomain = ['FLORIDA-STJ','FLORIDA-SFL-N','N.Carolina-IPHEx_Duke','N.Carolina-IPHEx_NASA','KWAJALEIN-KWA']
#ldomain = ['FLORIDA-STJ']

offset_bef = 15  # 'bef' should be identical to what is used in mk_match.py
offset_aft = 30


lrdt     = range(-offset_bef, offset_aft+1)
nt = len(lrdt)
nh = 20

#lprtype = ['heavy','extreme','mod']
lprtype = ['all','mod','heavy']
dlthpr = {'all':[-0.1,9999],'mod':[-0.1,10], 'heavy':[10,50],'extreme':[50,9999]}
ldattype = ['rain','cc','bias','brat','rmse','num']
#ldattype = ['rmse']


dansurf    = {rdt:deque([]) for rdt in lrdt} 
dagv       = {rdt:deque([]) for rdt in lrdt}
dansurfbin = {rdt:deque([]) for rdt in lrdt}

for idt, dt in enumerate(range(-offset_bef, offset_aft+1)):
    if calc != True:
        continue

    for domain in ldomain:

        
        for YM in lYM:
            Year, Mon = YM
            if (domain,Year,Mon) not in dgName.keys():
                print 'no obs',domain,Year,Mon
                continue
            print domain, 'dt=',dt,YM
            # load data
            srcbaseDir = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.L2A25'

            srcDir     = srcbaseDir + '/%.1fkm/%s/%04d%02d'%(thdist, domain, Year,Mon)
    
            profPath  = srcDir  + '/p_prof.npy'
            eSurfPath = srcDir  + '/p_eSurf.npy'
            nSurfPath = srcDir  + '/p_nSurf.npy'
            gvPath    = srcDir  + '/p_gvprcp.npy'
            ngvPath   = srcDir  + '/p_ngv.npy'
            nSurfBinPath = srcDir  + '/p_nSurfBin.npy'
    
            if not os.path.exists(profPath):
                print 'no file',profPath
                print 'for',domain,YM
                print 'skip'
                continue


            aprof  = np.load(profPath)
            #aesurf = np.load(eSurfPath)
            ansurf = np.load(nSurfPath)
            agv    = np.load(gvPath)[:,idt]
            angv   = np.load(ngvPath)
            ansurfbin= np.load(nSurfBinPath)

            asate  = ansurf


            #-- mask when ngv < minNum
            amskN    = ma.masked_less(angv, minNum).mask

            #-- mask when both satellite and gv are zero
            amsk1    = ma.masked_equal(asate,0).mask
            amsk2    = ma.masked_equal(agv,0).mask
            amskzero = amsk1 * amsk2

            #-- mask when at least one of the data is missing or low quality
            amsk3    = ma.masked_less(asate,0).mask
            amsk4    = ma.masked_less(agv,0).mask
            amskmiss = amsk3 + amsk4

            #-- overlay masks
            amsk     = amskN + amskzero + amskmiss

            aidxTmp  = arange(len(amsk)).astype(int32)

            aidxTmp  = ma.masked_where(amsk, aidxTmp).compressed()
        
            #aesurfTmp    = aesurf[aidxTmp]
            ansurfTmp    = asate[aidxTmp]
            agvTmp       = agv[aidxTmp]
            ansurfbinTmp = ansurfbin[aidxTmp]
            artime       = ansurfbinTmp + 1
            ardt         = dt - artime 


            #------------------
            for irdt, rdt in enumerate(lrdt):

                aidxTmp2  = arange(len(ansurfTmp))
                aidxTmp2  = ma.masked_where(ardt != rdt, aidxTmp2).compressed()

                dansurf[rdt].extend( ansurfTmp[aidxTmp2] )
                dagv[rdt] .extend( agvTmp [aidxTmp2] )
  


da1rain = {prtype: empty([len(lrdt)]) for prtype in lprtype}
da1cc   = {prtype: empty([len(lrdt)]) for prtype in lprtype}
da1rmse = {prtype: empty([len(lrdt)]) for prtype in lprtype}
da1bias = {prtype: empty([len(lrdt)]) for prtype in lprtype}
da1brat = {prtype: empty([len(lrdt)]) for prtype in lprtype}
da1num  = {prtype: empty([len(lrdt)],int32) for prtype in lprtype}

for irdt, rdt in enumerate(lrdt): 
    if calc != True: continue

    ansurf= array(dansurf[rdt]) *0.01  # mm/h
    agv   = array(dagv[rdt])
    #- nsurf vs gv --
    for prtype in lprtype:
        # mask with nSurf --
        thmin, thmax = dlthpr[prtype]
        amsk      = ma.masked_outside(ansurf, thmin, thmax).mask
        #amsk      = ma.masked_outside(agv, thmin, thmax).mask
        ansurfTmp  = ma.masked_where(amsk, ansurf ).compressed()
        agvTmp    = ma.masked_where(amsk, agv   ).compressed()

         
        cc     = np.corrcoef(ansurfTmp, agvTmp)[0,1]
        rmse   = np.sqrt( ((ansurfTmp- agvTmp)**2).mean() )
       
        bias   = (ansurfTmp - agvTmp).mean()
        brat   = (ansurfTmp - agvTmp).mean() / agvTmp.mean()*100
        num    = len(ansurfTmp)
        rain   = ansurfTmp.mean()

        da1rain[prtype][irdt] = rain
        da1cc  [prtype][irdt] = cc
        da1rmse[prtype][irdt] = rmse
        da1bias[prtype][irdt] = bias
        da1brat[prtype][irdt] = brat
        da1num [prtype][irdt] = num       


#--- save data ----
for prtype in lprtype:
    if calc != True:
        continue

    outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-PR/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])
    util.mk_dir(outDir)

    for dattype in ['cc','rmse','bias','brat','num','rain']:  
        datPath= outDir + '/dt-nsurf.%s.%s.npy'%(dattype, prtype)

        if dattype=='cc':
            a1dat = da1cc[prtype]
        elif dattype=='rmse':
            a1dat = da1rmse[prtype]
            print 'rmse max',da1rmse[prtype].max()
        elif dattype=='bias':
            a1dat = da1bias[prtype]
        elif dattype=='brat':
            a1dat = da1brat[prtype]
        elif dattype=='num':
            a1dat = da1num[prtype]
        elif dattype=='rain':
            a1dat = da1rain[prtype]

        np.save(datPath, a1dat)
        print datPath

#--- figure -------
dylim = {}

for prtype in lprtype:
    outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-PR/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])

    for dattype in ldattype:

        datPath=outDir + '/dt-nsurf.%s.%s.npy'%(dattype, prtype)
        a1dat  = np.load(datPath)
        fig = plt.figure(figsize=(4,3))
        a1t  = lrdt

        ax = fig.add_axes([0.15, 0.1, 0.8, 0.8])
        a1y = ma.masked_invalid(a1dat)

        if len(a1y.compressed())==0:
            a1y = zeros(len(a1t))
        
        ax.plot(a1t, a1y, '-', zorder=10)


        # ylim
        if dattype in ['cc']:
            ymin = -1.2
            ymax = 1.2
            plt.ylim([ymin,ymax])

        elif dattype in ['rain','rmse','num']:
            ymin = 0
            ymax = ma.masked_invalid(a1dat).max()*1.1
            plt.ylim([ymin,ymax])

        elif dattype in ['bias','brat']:
            ymax = abs(ma.masked_invalid(a1dat)).max()*1.1
            ymin = -ymax
            plt.ylim([ymin,ymax])


        #- zero line ---
        if dattype in ['cc','bias','brat']:
            plt.plot([-50,50],[0,0],'--',color='k', linewidth=0.5)

        plt.xlim([min(lrdt)-1, max(lrdt)+1])


        stitle  = '%s %s'%(prtype, dattype)
        plt.title(stitle)        
        figDir  = '/work/a01/utsumi/GPMGV/fig'
        figPath = figDir + '/plot.dt-nsurf.%s.%s.png'%(dattype, prtype)
        plt.savefig(figPath)
        print figPath
        plt.clf()
    



