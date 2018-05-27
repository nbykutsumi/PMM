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
iYM = [2013,4]
eYM = [2014,9]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
print lYM


gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()

dgName = gv.ret_ddomYM2gName()

ldomain = gv.domains
#ldomain = ['FLORIDA-STJ','FLORIDA-SFL-N','N.Carolina-IPHEx_Duke','N.Carolina-IPHEx_NASA','KWAJALEIN-KWA']
#ldomain = ['FLORIDA-STJ']

offset_bef = 15  # 'bef' should be identical to that in mk_match.py
offset_aft = 30


nt = offset_aft + offset_bef +1
nh = 20
#nh = 5

#lprtype = ['heavy','extreme','mod']
lprtype = ['mod','heavy','extreme']
dlthpr = {'mod':[-0.1,10], 'heavy':[10,50],'extreme':[50,9999]}
ldattype = ['rain','cc','bias','brat','rmse','num','gv']

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
            srcDir     = srcbaseDir + '/%s/%04d%02d'%(domain, Year,Mon)
    
            profPath  = srcDir  + '/prof.npy'
            eSurfPath = srcDir  + '/eSurf.npy'
            nSurfPath = srcDir  + '/nSurf.npy'
            gvPath    = srcDir  + '/gvprcp.npy'
            nSurfBinPath = srcDir  + '/nSurfBin.npy'
   
            if not os.path.exists(profPath):
                print 'no file',profPath
                print 'skip', domain, YM
                continue
 
            aprof  = np.load(profPath)
            aesurf = np.load(eSurfPath)
            ansurf = np.load(nSurfPath)
            agv    = np.load(gvPath)[:,idt]
            ansurfbin= np.load(nSurfBinPath)
  
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
                amsk     = amskzero + amskmiss

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
            amsk      = ma.masked_outside(ansurf, thmin, thmax).mask
            #amsk      = ma.masked_outside(agv, thmin, thmax).mask
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

    outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-PR/ndom.%02d.%04d.%02d-%04d.%02d'%(len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])
    util.mk_dir(outDir)

    for dattype in ldattype:  
        datPath= outDir + '/dt-lev.%s.%s.npy'%(dattype, prtype)

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
    outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-PR/ndom.%02d.%04d.%02d-%04d.%02d'%(len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])

    for dattype in ldattype:

        datPath=outDir + '/dt-lev.%s.%s.npy'%(dattype, prtype)
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


        stitle  = '%s %s'%(prtype, dattype)
        plt.title(stitle)        
        figDir  = '/work/a01/utsumi/GPMGV/fig'
        figPath = figDir + '/plot.dt-lev.%s.%s.png'%(dattype, prtype)
        plt.savefig(figPath)
        print figPath
        plt.clf()
    



