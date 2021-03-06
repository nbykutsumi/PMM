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
from matplotlib import rcParams, cycler

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

prdName = 'L2A25'
basepr = 'gv'
#basepr = 'sate'
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
lprtype = ['all','light','mod','heavy']
dlthpr = {'all':[-0.1,999],'light':[-0.1,2],'mod':[2,10], 'heavy':[10,50],'extreme':[50,9999]}
ldattype = ['rain','cc','bias','brat','rmse','num','gv']
#ldattype = ['rain']

da2rain = {prtype: empty([nh,nt]) for prtype in lprtype}
da2gv   = {prtype: empty([nh,nt]) for prtype in lprtype}
da2cc   = {prtype: empty([nh,nt]) for prtype in lprtype}
da2rmse = {prtype: empty([nh,nt]) for prtype in lprtype}
da2bias = {prtype: empty([nh,nt]) for prtype in lprtype}
da2brat = {prtype: empty([nh,nt]) for prtype in lprtype}
da2num  = {prtype: empty([nh,nt],int32) for prtype in lprtype}

# for eSurf
de1rain = {prtype: empty([nt]) for prtype in lprtype}
de1gv   = {prtype: empty([nt]) for prtype in lprtype}
de1cc   = {prtype: empty([nt]) for prtype in lprtype}
de1rmse = {prtype: empty([nt]) for prtype in lprtype}
de1bias = {prtype: empty([nt]) for prtype in lprtype}
de1brat = {prtype: empty([nt]) for prtype in lprtype}
de1num  = {prtype: empty([nt],int32) for prtype in lprtype}



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

        asate  = aprof[:,:nh]

        #-- mean array for masks
        a1sateAll = ma.masked_less(asate, 0).mean(axis=1)
        a1gvAll   = ma.masked_less(agv,0).mean(axis=1)
        a1sateMin = asate.min(axis=1)

        #-- mask when both satellite and gv are zero
        amsk1    = ma.masked_equal(a1sateAll,0).mask
        amsk2    = ma.masked_equal(a1gvAll,0).mask
        amskzero = amsk1 * amsk2

        ##-- mask when sate has  missing data in the profile
        #amskmiss = ma.masked_less(a1sateMin,0).mask

        #-- mask when ng  < minNum
        amskN    = ma.masked_less(angv, minNum).mask

        #-- overlay masks
        #amsk     = amskzero + amskmiss amskN
        amsk     = amskzero + amskN

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
    a2gv   = concatenate(a2gv,axis=0)
    a2prof = concatenate(a2prof,axis=0)
    a1esurf= array(a1esurf)
    a1nsurf= array(a1nsurf)*0.01 #mm/h
    a1nsurfbin= array(a1nsurfbin)
        
#----

for prtype in lprtype:

    # mask with eSurf --
    thmin, thmax = dlthpr[prtype]

    if basepr=='sate':
        amskP = ma.masked_outside(a1esurf, thmin, thmax).mask

    elif basepr=='gv':
        a1gvNow= a2gv[:,15]
        amskP = ma.masked_outside(a1gvNow, thmin, thmax).mask


    for ih in [-99] + range(nh): 
        if calc == False: continue
    
        if ih == -99:
            a1prof = a1esurf
        else:
            a1prof = a2prof[:,ih]  *0.01  # mm/h
    
        # mask when sate is missing
        amsk1 = ma.masked_less(a1prof,0).mask
        amskMiss = amsk1

        # overlay masks
        if (amskMiss is False) or (amskMiss is True):
            amskMiss = array([False]*len(a1prof))
        if (amskP is False) or (amskMiss is True):
            amskP    = array([False]*len(a1prof))

        amsk  = amskMiss + amskP
  
        for it in range(nt):
            a1gv = a2gv[:,it]

            a1profTmp = ma.masked_where(amsk, a1prof ).compressed()
            a1gvTmp   = ma.masked_where(amsk, a1gv   ).compressed()
       
             
            cc     = np.corrcoef(a1profTmp, a1gvTmp)[0,1]
            rmse   = np.sqrt( ((a1profTmp- a1gvTmp)**2).mean() )
           
            bias   = (a1profTmp - a1gvTmp).mean()
            brat   = (a1profTmp - a1gvTmp).mean() / a1gvTmp.mean()*100
            num    = len(a1profTmp)
            rain   = a1profTmp.mean()
            gv     = a1gvTmp.mean()
   
            if ih==-99:
                de1cc  [prtype][it] = cc
                de1rmse[prtype][it] = rmse
                de1bias[prtype][it] = bias
                de1brat[prtype][it] = brat
                de1num [prtype][it] = num       
                de1rain[prtype][it] = rain
                de1gv  [prtype][it] = gv

            else: 
                da2cc  [prtype][ih, it] = cc
                da2rmse[prtype][ih, it] = rmse
                da2bias[prtype][ih, it] = bias
                da2brat[prtype][ih, it] = brat
                da2num [prtype][ih, it] = num       
                da2rain[prtype][ih, it] = rain
                da2gv  [prtype][ih, it] = gv
    

#- save data ----
for prtype in lprtype:
    if calc == False: continue

    outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-%s/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])
    util.mk_dir(outDir)
    
    for dattype in ldattype:  
        datPath= outDir + '/dt-lev.base.%s.%s.%s.npy'%(basepr, dattype, prtype)
        esPath = outDir + '/dt-eSurf.base.%s.%s.%s.npy'%(basepr, dattype, prtype)
    
        if dattype=='cc':
            a2dat = da2cc[prtype]
            e1dat = de1cc[prtype]
        elif dattype=='rmse':
            a2dat = da2rmse[prtype]
            e1dat = de1rmse[prtype]
        elif dattype=='bias':
            a2dat = da2bias[prtype]
            e1dat = de1bias[prtype]
        elif dattype=='brat':
            a2dat = da2brat[prtype]
            e1dat = de1brat[prtype]
        elif dattype=='num':
            a2dat = da2num[prtype]
            e1dat = de1num[prtype]
        elif dattype=='rain':
            a2dat = da2rain[prtype]
            e1dat = de1rain[prtype]
        elif dattype=='gv':
            a2dat = da2gv[prtype]
            e1dat = de1gv[prtype]
    
    
        np.save(datPath, a2dat)
        np.save(esPath, e1dat)
        print datPath
    

#-- Figure in one panel ------------
for prtype in lprtype:
    outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-%s/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])

    for dattype in ldattype:

        datPath=outDir + '/dt-lev.base.%s.%s.%s.npy'%(basepr,dattype, prtype)
        esPath =outDir + '/dt-eSurf.base.%s.%s.%s.npy'%(basepr, dattype, prtype)

        a2dat  = np.load(datPath)
        e1dat  = np.load(esPath)

        fig = plt.figure(figsize=(4,8))
        ax  = fig.add_axes([0.15, 0.1, 0.8, 0.8])

        a1t  = range(-offset_bef, offset_aft+1)

        cmap = plt.get_cmap('coolwarm')

        lh = range(nh)[2:]
        for itmp,ih in enumerate(lh):
            a1y = a2dat[ih,:]
            a1y = ma.masked_invalid(a1y)

            ax.plot(a1t, a1y, '-', c=cmap(float(itmp)/len(lh)), zorder=10, label='%d'%(ih))

        # eSurf
        a1y = ma.masked_invalid(e1dat)
        ax.plot(a1t, a1y, '--', zorder=10, label='eSurf', color='k')
        # legend
        ax.legend()


        # ylim
        if dattype in ['cc']:
            ymin = -0.1
            ymax = 0.8
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



'''
#--- figure : multi-rows -------

for prtype in lprtype:
    outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-%s/ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])

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
        figPath = figDir + '/plot.dt-lev.simple.%s.%s.%s.png'%(prdName, dattype, prtype)
        plt.savefig(figPath)
        print figPath
        plt.clf()
    
'''

