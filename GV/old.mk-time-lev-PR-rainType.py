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
nh = 30

#lprtype = ['heavy','extreme','mod']
#lprtype = ['all','light','mod','heavy']
lprtype = ['all']
dlthpr = {'all':[-0.1,999],'light':[-0.1,2],'mod':[2,10], 'heavy':[10,50],'extreme':[50,9999]}
ldattype = ['rain','cc','bias','brat','rmse','num','gv']
#ldattype = ['rain']

lraintype = ['alltype','strat','conv']
dlthrtype = {'alltype':[0,999],'strat':[100,170],'conv':[200,297]}

#-----------------------------------------------------
#def running_mean(a1dat):
#    n = len(a1dat)
#    vfirst= a1dat[0]*2./3  + a1dat[1]*1./3
#    vlast = a1dat[-1]*2./3 + a1dat[-2]*1./3
#    a1out = empty(n, float32)
#    for i in range(1,n-1):
#        a1out[i] = a1dat[i-1:i+1].mean()
#
#    a1out[0] = vfirst
#    a1out[-1] = vlast
#    return a1out

def running_mean(a1dat):
    ''' 1,2,3,2,1 karnel '''

    n = len(a1dat)
    a1out = empty(n, float32)
    for i in range(n):
        a1msk = ma.masked_outside(range(n),i-2,i+2).mask
        a1wt  = 3-abs( i - arange(n))
        a1wt  = a1wt.astype(float32)
        a1datTmp = ma.masked_where(a1msk, a1dat)
        a1wtTmp  = ma.masked_where(a1msk, a1wt)
        a1wtTmp  = a1wtTmp / a1wtTmp.sum()
        a1out[i] = (a1datTmp * a1wtTmp).sum()

    return a1out




#-----------------------------------------------------
da2rain = {(prtype,raintype): empty([nh,nt]) for prtype in lprtype for raintype in lraintype}
da2gv   = {(prtype,raintype): empty([nh,nt]) for prtype in lprtype for raintype in lraintype}
da2cc   = {(prtype,raintype): empty([nh,nt]) for prtype in lprtype for raintype in lraintype}
da2rmse = {(prtype,raintype): empty([nh,nt]) for prtype in lprtype for raintype in lraintype}
da2bias = {(prtype,raintype): empty([nh,nt]) for prtype in lprtype for raintype in lraintype}
da2brat = {(prtype,raintype): empty([nh,nt]) for prtype in lprtype for raintype in lraintype}
da2num  = {(prtype,raintype): empty([nh,nt],int32) for prtype in lprtype for raintype in lraintype}

# for eSurf
de1rain = {(prtype,raintype): empty([nt]) for prtype in lprtype for raintype in lraintype}
de1gv   = {(prtype,raintype): empty([nt]) for prtype in lprtype for raintype in lraintype}
de1cc   = {(prtype,raintype): empty([nt]) for prtype in lprtype for raintype in lraintype}
de1rmse = {(prtype,raintype): empty([nt]) for prtype in lprtype for raintype in lraintype}
de1bias = {(prtype,raintype): empty([nt]) for prtype in lprtype for raintype in lraintype}
de1brat = {(prtype,raintype): empty([nt]) for prtype in lprtype for raintype in lraintype}
de1num  = {(prtype,raintype): empty([nt],int32) for prtype in lprtype for raintype in lraintype}



a2prof     = deque([]) 
a1esurf    = deque([]) 
a1nsurf    = deque([]) 
a2gv       = deque([]) 
a1nsurfbin = deque([]) 
a1raintype = deque([])

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
        
        rainTypePath = srcDir  + '/p_rainType.npy'

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
        araintype= np.load(rainTypePath)

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
        araintypeTmp = araintype[aidxTmp]
        #------------------
    
        a2gv.append(agvTmp)
        a2prof.append(asateTmp)
        a1esurf.extend(aesurfTmp)
        a1nsurf.extend(ansurfTmp)
        a1nsurfbin.extend(ansurfbinTmp)

        a1raintype.extend(araintypeTmp)

if calc ==True:
    a2gv   = concatenate(a2gv,axis=0)
    a2prof = concatenate(a2prof,axis=0)
    a1esurf= array(a1esurf)
    a1nsurf= array(a1nsurf)*0.01 #mm/h
    a1nsurfbin= array(a1nsurfbin)
    a1raintype= array(a1raintype) 
#----

for raintype in lraintype:
    if calc != True: continue

    raintypemin, raintypemax = dlthrtype[raintype]
    amskRainType = ma.masked_outside(a1raintype, raintypemin, raintypemax).mask


    for prtype in lprtype:
    
        # mask with eSurf --
        thmin, thmax = dlthpr[prtype]
    
        if basepr=='sate':
            amskP = ma.masked_outside(a1esurf, thmin, thmax).mask
    
        elif basepr=='gv':
            a1gvNow= a2gv[:,15]
            amskP = ma.masked_outside(a1gvNow, thmin, thmax).mask
    
    
        #for ih in [-99] + range(nh): 
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
    
            amsk  = amskMiss + amskP + amskRainType
      
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
                    de1cc  [prtype,raintype][it] = cc
                    de1rmse[prtype,raintype][it] = rmse
                    de1bias[prtype,raintype][it] = bias
                    de1brat[prtype,raintype][it] = brat
                    de1num [prtype,raintype][it] = num       
                    de1rain[prtype,raintype][it] = rain
                    de1gv  [prtype,raintype][it] = gv
    
                else: 
                    da2cc  [prtype,raintype][ih, it] = cc
                    da2rmse[prtype,raintype][ih, it] = rmse
                    da2bias[prtype,raintype][ih, it] = bias
                    da2brat[prtype,raintype][ih, it] = brat
                    da2num [prtype,raintype][ih, it] = num       
                    da2rain[prtype,raintype][ih, it] = rain
                    da2gv  [prtype,raintype][ih, it] = gv
        
    
#- save data ----
for raintype in lraintype:
    for prtype in lprtype:
        if calc == False: continue
    
        outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-%s/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])
        util.mk_dir(outDir)
        
        for dattype in ldattype:  
            datPath= outDir + '/dt-lev.base.%s.%s.%s.%s.npy'%(basepr, dattype, raintype, prtype)
            esPath = outDir + '/dt-eSurf.base.%s.%s.%s.%s.npy'%(basepr, dattype, raintype, prtype)
        
            if dattype=='cc':
                a2dat = da2cc[prtype,raintype]
                e1dat = de1cc[prtype,raintype]
            elif dattype=='rmse':
                a2dat = da2rmse[prtype,raintype]
                e1dat = de1rmse[prtype,raintype]
            elif dattype=='bias':
                a2dat = da2bias[prtype,raintype]
                e1dat = de1bias[prtype,raintype]
            elif dattype=='brat':
                a2dat = da2brat[prtype,raintype]
                e1dat = de1brat[prtype,raintype]
            elif dattype=='num':
                a2dat = da2num[prtype,raintype]
                e1dat = de1num[prtype,raintype]
            elif dattype=='rain':
                a2dat = da2rain[prtype,raintype]
                e1dat = de1rain[prtype,raintype]
            elif dattype=='gv':
                a2dat = da2gv[prtype,raintype]
                e1dat = de1gv[prtype,raintype]
        
        
            np.save(datPath, a2dat)
            np.save(esPath, e1dat)
            print datPath

'''    
#-- Figure imshow ------------
for raintype in lraintype:
    for prtype in lprtype:
        outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-%s/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])
        for dattype in ['cc']:
            datPath=outDir + '/dt-lev.base.%s.%s.%s.%s.npy'%(basepr,dattype, raintype, prtype)
            esPath =outDir + '/dt-eSurf.base.%s.%s.%s.%s.npy'%(basepr, dattype, raintype, prtype)
    
            a2dat  = np.load(datPath)
            e1dat  = np.load(esPath)
    
            fig = plt.figure(figsize=(4,4))
            ax  = fig.add_axes([0.15, 0.1, 0.8, 0.8])
    
            a1t  = range(-offset_bef, offset_aft+1)
    
            cmap = plt.get_cmap('coolwarm')

            a2fig= concatenate([e1dat.reshape(1,-1), a2dat], axis=0)

            im = ax.imshow(a2fig, origin='lower', vmin=0,vmax=0.6) 
            plt.colorbar(im)
 
            # title
            stitle  = '%s %s %.1fkm minNum=%d %s'%(raintype, prtype, thdist, minNum, dattype)
            stitle  = stitle + '\n' + 'base:%s %04d.%02d-%04d.%02d'%(basepr, iYM[0],iYM[1],eYM[0],eYM[1])
            plt.title(stitle)
            figDir  = '/work/a01/utsumi/GPMGV/fig'
            figPath = figDir + '/grid.dt-lev.rainType.runmean.%s.%.1fkm.minNum.%d.base.%s.%s.%s.%s.png'%(prdName, thdist, minNum, basepr, dattype, raintype, prtype)

            plt.savefig(figPath)
            print figPath
            plt.clf()

'''    


#-- Figure in one panel ------------
for raintype in lraintype:
    for prtype in lprtype:
        outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-%s/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])
    
        for dattype in ldattype:
    
            datPath=outDir + '/dt-lev.base.%s.%s.%s.%s.npy'%(basepr,dattype, raintype, prtype)
            esPath =outDir + '/dt-eSurf.base.%s.%s.%s.%s.npy'%(basepr, dattype, raintype, prtype)
    
            a2dat  = np.load(datPath)
            e1dat  = np.load(esPath)
    
            #fig = plt.figure(figsize=(4,8))
            fig = plt.figure(figsize=(4,5))
            ax  = fig.add_axes([0.15, 0.1, 0.8, 0.8])
    
            a1t  = range(-offset_bef, offset_aft+1)
    
            cmap = plt.get_cmap('coolwarm')
    
            # eSurf
            a1y = ma.masked_invalid(e1dat)
            a1y = running_mean(a1y)
            ax.plot(a1t, a1y, '--', zorder=10, label='eSurf', color='k', linewidth=2)
            #-- find peak or bottom for eSurf ---
            imin = a1y.argmin()
            imax = a1y.argmax()

            if dattype=='cc':
                #ax.plot(a1t[imax],0.01,'D',c='k', markersize=12)
                ax.plot(a1t[imax],a1y[imax],'D',c='k', markersize=12)
 

    
            #lh = range(nh)[2:]
            #lh = [2,4,8,12,16]
            lh = [4,8,12,16,20,24]
    
            for itmp,ih in enumerate(lh):
                #a1yA = a2dat[ih,:]
                #a1yB = a2dat[ih+1,:]
                #a1y  = (a1yA + a1yB)*0.5
                a1y = a2dat[ih,:]
                a1y = ma.masked_invalid(a1y)
                a1y = running_mean(a1y)
   
                mycm= cmap(float(itmp)/len(lh)) 
                ax.plot(a1t, a1y, '-', c=mycm, label='%.1f'%(ih*0.25), linewidth=2)


                #-- find peak or bottom ---
                imin = a1y.argmin()
                imax = a1y.argmax()

                if dattype=='cc':
                    #ax.plot(a1t[imax],0.01,'v',c=mycm, markersize=10)
                    ax.plot(a1t[imax],a1y[imax],'v',c=mycm, markersize=10)
    
            # legend
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles[::-1], labels[::-1])
    
    
            # ylim
            if dattype in ['cc']:
                ymin = 0.2
                #ymax = 1.0
                ymax = 0.7
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
            stitle  = '%s %s %.1fkm minNum=%d %s'%(raintype, prtype, thdist, minNum, dattype)
            stitle  = stitle + '\n' + 'base:%s %04d.%02d-%04d.%02d'%(basepr, iYM[0],iYM[1],eYM[0],eYM[1])
            plt.title(stitle)
            figDir  = '/work/a01/utsumi/GPMGV/fig'
            figPath = figDir + '/plot.dt-lev.rainType.runmean.%s.%.1fkm.minNum.%d.base.%s.%s.%s.%s.png'%(prdName, thdist, minNum, basepr, dattype, raintype, prtype)
            plt.savefig(figPath)
            print figPath
            plt.clf()


