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
iYM = [2011,6]
eYM = [2014,10]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
lYM = lYM[::-1]
print lYM

#thdist  = 7.5 # km
thdist  = 15 # km
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
nh = 20  # 0.5 - 10km, at most.
#nh = 5

#lprtype = ['heavy']
lprtype = ['all','mod','heavy','extreme']
dlthpr = {'all':[-0.1,10], 'mod':[-0.1,10], 'heavy':[10,50],'extreme':[50,9999]}
ldattype = ['rain','cc','bias','brat','rmse','num','gv']

coef_def = 4.17 * 1000/(60*60)  # g/m3 * m/s *1000/(60*60) = mm/h

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

da2rain = {prtype: empty([nh,nt]) for prtype in lprtype}
da2gv   = {prtype: empty([nh,nt]) for prtype in lprtype}
da2cc   = {prtype: empty([nh,nt]) for prtype in lprtype}
da2rmse = {prtype: empty([nh,nt]) for prtype in lprtype}
da2bias = {prtype: empty([nh,nt]) for prtype in lprtype}
da2brat = {prtype: empty([nh,nt]) for prtype in lprtype}
da2num  = {prtype: empty([nh,nt],int32) for prtype in lprtype}

for idt, dt in enumerate(range(-offset_bef, offset_aft+1)):
    if calc != True:
        continue

    daprof     = {ih:deque([]) for ih in range(nh)}
    daesurf    = {ih:deque([]) for ih in range(nh)}
    dagv       = {ih:deque([]) for ih in range(nh)}
    dacoef     = {ih:deque([]) for ih in range(nh)}

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
    
            gvPath        = srcDir  + '/p_gvprcp.npy'
            qFlagPath     = srcDir  + '/p_qFlag.npy'
            profNumPath   = srcDir  + '/p_profNum.npy'
            profScalePath = srcDir  + '/p_profScale.npy'
            tIndexPath    = srcDir  + '/p_tIndex.npy'
    
            eSurfPath    = srcDir + '/p_eSurf.npy'
            mlPrcpPath   = srcDir + '/p_mlPrecip.npy'
            groundBinPath= srcDir + '/p_groundBin.npy'
            ngvPath      = srcDir + '/p_ngv.npy' 


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

            agv    = np.load(gvPath)[:,idt]

            acoef  = ma.masked_invalid(aesurf/aprof[:,0]).filled(coef_def)

            #---- mask where ngv < minNum
            amskN  = ma.masked_less(angv, minNum).mask
 
            #---- mask where qFlag >0 -
            amskQ  = ma.masked_greater(aqFlag, 0).mask
 
            #-------------------------- 
            for ih in range(nh):
                asate  = aprof[:,ih]
                #-- mask when both satellite and gv are zero
                amsk1    = ma.masked_equal(asate,0).mask
                amsk2    = ma.masked_equal(agv,0).mask
                amskzero = amsk1 * amsk2

                #-- mask when at least one of the data is missing or low quality
                amsk3    = ma.masked_less(asate,0).mask
                amsk4    = ma.masked_less(agv,0).mask
                amskmiss = amsk3 + amsk4

                #-- overlay masks
                amsk     = amskN + amskQ + amskzero + amskmiss

                aidxTmp  = arange(len(amsk)).astype(int32)
                aidxTmp  = ma.masked_where(amsk, aidxTmp).compressed()
        
                asateTmp     = asate[aidxTmp] 
                aesurfTmp    = aesurf[aidxTmp]
                agvTmp       = agv[aidxTmp]
                acoefTmp     = acoef[aidxTmp]
                #------------------
        
                daprof[ih].extend(asateTmp)
                daesurf[ih].extend(aesurfTmp)
                dagv[ih].extend(agvTmp)
                dacoef[ih].extend(acoefTmp)

 
    for ih in range(nh): 
        #-------------
        # coeffitiont for aprof kg/m3 --> mm/h
        #-------------
        aprof_tmp = array(daprof[ih])   # g/m3
        aesurf    = array(daesurf[ih])  # mm/h
        a1coef    = array(dacoef[ih])  

        aprof = aprof_tmp * a1coef
        agv   = array(dagv[ih])

        #- profile vs gv --
        for prtype in lprtype:
            # mask with eSurf --
            thmin, thmax = dlthpr[prtype]
            amsk      = ma.masked_outside(aesurf, thmin, thmax).mask
            #amsk      = ma.masked_outside(agv, thmin, thmax).mask
            aprofTmp  = ma.masked_where(amsk, aprof ).compressed()
            agvTmp    = ma.masked_where(amsk, agv   ).compressed()
            aesurfTmp = ma.masked_where(amsk, aesurf).compressed()
       
             
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
        datPath= outDir + '/dt-lev.%s.%s.%s.npy'%(prdName, dattype, prtype)

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
#-- Multi rows ----
'''
for prtype in lprtype:
    outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-%s/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])

    for dattype in ldattype:

        datPath=outDir + '/dt-lev.%s.%s.%s.npy'%(prdName, dattype, prtype)
        a2dat  = np.load(datPath)
    
        fig = plt.figure(figsize=(4,8))
        a1t  = range(-offset_bef, offset_aft+1)


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

            plt.xlim([-6,31])



            if ih !=0:
                ax.tick_params(labelbottom='off') 


        stitle  = '%s %.1fkm minNum=%d %s'%(prtype, thdist, minNum, dattype)
        plt.title(stitle)        
        figDir  = '/work/a01/utsumi/GPMGV/fig'
        figPath = figDir + '/plot.dt-lev.%s.%.1fkm.minNum.%d.%s.%s.png'%(prdName, thdist, minNum, dattype, prtype)
        plt.savefig(figPath)
        print figPath
        plt.clf()
    
'''

for prtype in lprtype:
    outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-%s/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])

    for dattype in ldattype:

        datPath=outDir + '/dt-lev.%s.%s.%s.npy'%(prdName, dattype, prtype)
        a2dat  = np.load(datPath)
    
        fig = plt.figure(figsize=(4,8))
        ax  = fig.add_axes([0.15, 0.1, 0.8, 0.8])

        a1t  = range(-offset_bef, offset_aft+1)

        cmap = plt.cm.coolwarm
        rcParams['axes.prop_cycle'] = cycler(color=cmap(np.linspace(0,1,nh)))

        print '-'*50
        lines = []
        for ih in range(nh):
            a1y = a2dat[ih,:]
            a1y = ma.masked_invalid(a1y)

            ax.plot(a1t, a1y, '-', zorder=10, label='%d'%(ih))
        # legend
        ax.legend()



        # ylim
        if dattype in ['cc']:
            ymin = -0.6
            ymax = 0.6
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
        plt.title(stitle)        
        figDir  = '/work/a01/utsumi/GPMGV/fig'
        figPath = figDir + '/plot.dt-lev.onebox.%s.%.1fkm.minNum.%d.%s.%s.png'%(prdName, thdist, minNum, dattype, prtype)
        plt.savefig(figPath)
        print figPath
        plt.clf()
    




