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
iYM = [2010,4]
eYM = [2014,10]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
print lYM

thdist = 5.0 # km
minNum = 3
prdName= 'L2A25'
gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()

dgName = gv.ret_ddomYM2gName()

ldomain = gv.domains
#ldomain = ['FLORIDA-STJ','FLORIDA-SFL-N','N.Carolina-IPHEx_Duke','N.Carolina-IPHEx_NASA','KWAJALEIN-KWA']
#ldomain = ['FLORIDA-STJ']

offset_bef = 15  # 'bef' should be identical to what is used in mk_match.py
#offset_aft = 45

ldt  = [1, 5,10,15,20,25,30,35]
ldh  = [1, 5,10,15,20,25]

nt  = len(ldt)
nh  = len(ldh)

lprtype = ['all','light','mod','heavy']
dlthpr = {'all':[-0.1,999],'light':[-0.1,2],'mod':[2,10], 'heavy':[10,50],'extreme':[50,9999]}
ldattype = ['rain','cc','bias','brat','rmse','num','gv']



miss    = -9999.

da2rain = {prtype: empty([nh,nt]) for prtype in lprtype}
da2gv   = {prtype: empty([nh,nt]) for prtype in lprtype}
da2cc   = {prtype: empty([nh,nt]) for prtype in lprtype}
da2rmse = {prtype: empty([nh,nt]) for prtype in lprtype}
da2bias = {prtype: empty([nh,nt]) for prtype in lprtype}
da2brat = {prtype: empty([nh,nt]) for prtype in lprtype}
da2num  = {prtype: empty([nh,nt],int32) for prtype in lprtype}

for idt, dt in enumerate(ldt):
    if calc != True:
        continue

    daprof     = {ih:deque([]) for ih in range(nh)}
    daesurf    = {ih:deque([]) for ih in range(nh)}
    dansurf    = {ih:deque([]) for ih in range(nh)}
    dagv       = {ih:deque([]) for ih in range(nh)}
    dansurfbin = {ih:deque([]) for ih in range(nh)}

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
                print 'skip', domain, YM
                continue
 
            aprof  = np.load(profPath)
            aesurf = np.load(eSurfPath)
            ansurf = np.load(nSurfPath)
            ansurfbin= np.load(nSurfBinPath)
            angv   = np.load(ngvPath)

            a2gv   = abs(np.load(gvPath))
            #agv   = gv_fsub.mean_slice_negativemask(a2gv.T, ansurfbin, dt)
            agv    = a2gv[:,15:15+dt].mean(axis=1)

            #-- accumulation for mask -
            #asateAcc  = gv_fsub.mean_slice_negativemask(aprof.T, ansurfbin, ldh[-1])
            asateAcc  = ma.masked_less(aprof[:,:ldh[-1]], 0).mean(axis=1)
            #agvAcc    = gv_fsub.mean_slice_negativemask(a2gv.T, ansurfbin, ldt[-1])
            agvAcc    = a2gv[:,15:15+ldt[-1]].mean(axis=1)

            #-- mask when ngv < minNum
            amskN    = ma.masked_less(angv, minNum).mask

            #-- mask when both satellite and gv are zero
            amsk1    = ma.masked_equal(asateAcc,0).mask
            amsk2    = ma.masked_equal(agvAcc,0).mask
            amskzero = amsk1 * amsk2

            #-- mask when at least one of the data is missing or low quality
            amsk3    = ma.masked_less(asateAcc,0).mask
            amsk4    = ma.masked_less(agvAcc,0).mask
            amskmiss = amsk3 + amsk4

            #-- overlay masks
            amsk     = amskN + amskzero + amskmiss

 
            for ih,dh in enumerate(ldh):
                #asate    = gv_fsub.mean_slice_negativemask(aprof.T, ansurfbin, dh)
                asate    = ma.masked_less(aprof[:,:dh], 0).mean(axis=1)

                aidxTmp  = arange(len(amsk)).astype(int32)
                aidxTmp  = ma.masked_where(amsk, aidxTmp).compressed()
        
                asateTmp     = asate[aidxTmp] 
                aesurfTmp    = aesurf[aidxTmp]
                ansurfTmp    = ansurf[aidxTmp]
                agvTmp       = agv[aidxTmp]
                ansurfbinTmp = ansurfbin[aidxTmp]
                #------------------
        
                daprof[ih].extend(asateTmp)
                daesurf[ih].extend(aesurfTmp)
                dansurf[ih].extend(ansurfTmp)
                dagv[ih].extend(agvTmp)
                dansurfbin[ih].extend(ansurfbinTmp)
   
    #----
    for ih in range(nh): 
        aprof = array(daprof[ih])  *0.01  # mm/h
        aesurf= array(daesurf[ih])
        ansurf= array(dansurf[ih]) *0.01  # mm/h
        agv   = array(dagv[ih])
        ansurfbin= array(dansurfbin[ih])


        print ih, '-'*50
        print ansurfbin.shape


        #- profile vs gv --
        for prtype in lprtype:
            # mask with nSurf --
            thmin, thmax = dlthpr[prtype]
            #amsk      = ma.masked_outside(ansurf, thmin, thmax).mask
            amsk      = ma.masked_outside(agv, thmin, thmax).mask
            aprofTmp  = ma.masked_where(amsk, aprof ).compressed()
            agvTmp    = ma.masked_where(amsk, agv   ).compressed()
            ansurfTmp = ma.masked_where(amsk, ansurf).compressed()
       
             
            cc     = np.corrcoef(aprofTmp, agvTmp)[0,1]
            rmse   = np.sqrt( ((aprofTmp- agvTmp)**2).mean() )
           
            bias   = (aprofTmp - agvTmp).mean()
            brat   = (aprofTmp - agvTmp).mean() / agvTmp.mean()*100
            num    = len(aprofTmp)
            rain   = aprofTmp.mean()
            gv     = agvTmp.mean()

            da2cc  [prtype][ih, idt] = cc
            da2rmse[prtype][ih, idt] = rmse
            da2bias[prtype][ih, idt] = bias
            da2brat[prtype][ih, idt] = brat
            da2num [prtype][ih, idt] = num       
            da2rain[prtype][ih, idt] = rain
            da2gv  [prtype][ih, idt] = gv


#--- save data ----
for prtype in lprtype:
    if calc != True:
        continue

    outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-%s/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])
    util.mk_dir(outDir)

    for dattype in ldattype:  
        datPath= outDir + '/nt-nlev.%s.%s.npy'%(dattype, prtype)

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


#--- figure -------
for prtype in lprtype:
    outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-%s/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])

    for dattype in ldattype:

        datPath=outDir + '/nt-nlev.%s.%s.npy'%(dattype, prtype)
        a2dat  = np.load(datPath)
    
        fig = plt.figure(figsize=(8,4))
        ax  = fig.add_axes([0.15, 0.1, 0.8, 0.8])
        a1t  = ldt

        wbar = 0.7/nh
        print '-'*50
        for ih, dh in enumerate(ldh):
            a1y = a2dat[ih,:]
            a1y = ma.masked_invalid(a1y)
            a1x = arange(nt) -0.35 +ih*wbar
            ax.bar(a1x, a1y, width=wbar, tick_label=ldt, align='center', label='%d bins'%(dh))
            #if ih !=0:
            #    ax.tick_params(labelbottom='off') 

        # legend
        plt.legend()
        # ylim
        if dattype in ['cc']:
            ymin = -1.2
            ymax = 1.2
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

        plt.xlim([a1x[0]-1,a1x[-1]+2])


        stitle  = '%s %s'%(dattype, prtype)
        plt.title(stitle)        
        figDir  = '/work/a01/utsumi/GPMGV/fig'
        figPath = figDir + '/bar.nt-nlev.%s.%.1fkm.minNum.%d.%s.%s.png'%(prdName, thdist,minNum,dattype, prtype)
        plt.savefig(figPath)
        print figPath
        plt.clf()


'''
#--- figure -------
for prtype in lprtype:
    outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-PR/ndom.%02d.%04d.%02d-%04d.%02d'%(len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])

    for dattype in ldattype:

        datPath=outDir + '/nt-nlev.%s.%s.npy'%(dattype, prtype)
        a2dat  = np.load(datPath)
    
        fig = plt.figure(figsize=(4,8))
        a1t  = ldt


        print '-'*50
        for ih in range(nh):
            ax = fig.add_axes([0.15, 0.1 + 0.8/nh*ih, 0.8, 0.8/nh])
            a1y = a2dat[ih,:]
            a1y = ma.masked_invalid(a1y)

            ax.plot(a1t, a1y, '-', zorder=10)


            # ylim
            if dattype in ['cc']:
                ymin = -1.2
                ymax = 1.2
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

            plt.xlim([ldt[0]-1,ldt[-1]+1])



            if ih !=0:
                ax.tick_params(labelbottom='off') 


        stitle  = '%s %s'%(prtype, dattype)
        plt.title(stitle)        
        figDir  = '/work/a01/utsumi/GPMGV/fig'
        figPath = figDir + '/plot.nt-nlev.%s.%s.png'%(dattype, prtype)
        plt.savefig(figPath)
        print figPath
        plt.clf()
    
'''
