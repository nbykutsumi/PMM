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
#iYM = [2014,4]
eYM = [2014,10]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
print lYM
figtype = 'bar'
#figtype = 'cont'
#corrFlag= 'CR'
corrFlag= 'NO'



#thdist = 2.5
thdist = 5.0
minNum = 3
prdName = 'L2A25'

gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()

dgName = gv.ret_ddomYM2gName()

ldomain = gv.domains
#ldomain = ['FLORIDA-STJ','FLORIDA-SFL-N','N.Carolina-IPHEx_Duke','N.Carolina-IPHEx_NASA','KWAJALEIN-KWA']
#ldomain = ['FLORIDA-STJ']

offset_bef = 15  # 'bef' should be identical to that in mk_match.py
#offset_aft = 30

if   figtype=='bar':
    ldt  = [1, 5,10,15,20,25,30]
    ldh  = range(2,34+1)[::2]
elif figtype=='cont':
    ldt  = range(1,30+1)
    ldh  = range(1,34+1)


nt = len(ldt)
nh = len(ldh)
#nh = 5

lrhtype = ['all','dry','hum']
dlthrh  = {'all':[-0.1,9999],'dry':[-0.1,0.7-0.0001],'hum':[0.7,9999]}

ldattype = ['rain','cc','bias','brat','rmse','gv','num']

da2rain = {rhtype: empty([nh,nt]) for rhtype in lrhtype}
da2gv   = {rhtype: empty([nh,nt]) for rhtype in lrhtype}
da2cc   = {rhtype: empty([nh,nt]) for rhtype in lrhtype}
da2rmse = {rhtype: empty([nh,nt]) for rhtype in lrhtype}
da2bias = {rhtype: empty([nh,nt]) for rhtype in lrhtype}
da2brat = {rhtype: empty([nh,nt]) for rhtype in lrhtype}
da2num  = {rhtype: empty([nh,nt],int32) for rhtype in lrhtype}

# for eSurf
de1rain = {rhtype: empty([nt]) for rhtype in lrhtype}
de1gv   = {rhtype: empty([nt]) for rhtype in lrhtype}
de1cc   = {rhtype: empty([nt]) for rhtype in lrhtype}
de1rmse = {rhtype: empty([nt]) for rhtype in lrhtype}
de1bias = {rhtype: empty([nt]) for rhtype in lrhtype}
de1brat = {rhtype: empty([nt]) for rhtype in lrhtype}
de1num  = {rhtype: empty([nt],int32) for rhtype in lrhtype}



a2prof     = deque([]) 
a1esurf    = deque([]) 
a1nsurf    = deque([]) 
a2gv       = deque([]) 
a1rh       = deque([])

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
        rhPath    = srcDir  + '/p_rh.npy'
        

        if not os.path.exists(profPath):
            print 'no file',profPath
            print 'skip', domain, YM
            continue

        aprof  = np.load(profPath)
        aesurf = np.load(eSurfPath)
        ansurf = np.load(nSurfPath)
        agv    = abs(np.load(gvPath))
        angv   = np.load(ngvPath)
        arh    = np.load(rhPath)

        asate  = aprof

        #-- mean array for masks
        a1idx0    = zeros(aprof.shape[0])
        a1sateAll = gv_fsub.mean_slice_negativemask(aprof.T, a1idx0, ldh[-1])
        a1idx15   = ones(agv.shape[0])*15
        a1gvAll   = ma.masked_less(agv[:,15:],0).mean(axis=1)
        a1sateMin = asate.min(axis=1)


        #-- mask when both satellite and gv are zero
        amsk1    = ma.masked_equal(a1sateAll,0).mask
        amsk2    = ma.masked_equal(a1gvAll,0).mask
        amskzero = amsk1 * amsk2

        #-- mask when ng  < minNum
        amskN    = ma.masked_less(angv, minNum).mask

        #-- overlay masks
        amsk     = amskzero + amskN 

        aidxTmp  = arange(len(amsk)).astype(int32)
        aidxTmp  = ma.masked_where(amsk, aidxTmp).compressed()
    
        asateTmp     = asate[aidxTmp] 
        aesurfTmp    = aesurf[aidxTmp]
        ansurfTmp    = ansurf[aidxTmp]
        agvTmp       = agv[aidxTmp]
        arhTmp       = arh[aidxTmp]
        #------------------
    
        a2gv.append(agvTmp)
        a2prof.append(asateTmp)
        a1esurf.extend(aesurfTmp)
        a1nsurf.extend(ansurfTmp)
        a1rh.extend(arhTmp)

if calc ==True:
    a1esurf= array(a1esurf)
    a1nsurf= array(a1nsurf)*0.01 #mm/h
    a1rh   = array(a1rh)
    a2gv   = concatenate(a2gv,axis=0)
    a2prof = concatenate(a2prof,axis=0) *0.01
    a2joinprof= concatenate([a1esurf.reshape(-1,1) ,a2prof], axis=1)

#----
for rhtype in lrhtype:
    if calc==False: continue

    # mask with RH --
    thmin, thmax = dlthrh[rhtype]
    amskRH = ma.masked_outside(a1rh, thmin, thmax).mask


    #--- vertical layer loop ---
    ih = -1
    for dh in [-99] + ldh: 
        if calc == False: continue
    
        if dh == -99:
            a1prof = a1esurf
        else:
            ih     = ih + 1
            a1idx0 = zeros(a2prof.shape[0])
            #a1prof = gv_fsub.mean_slice_negativemask(a2joinprof.T, a1idx0, dh+1)
            a1prof = ma.masked_less(a2joinprof[:,:dh+1], 0).mean(axis=1)
    

        # mask events when sate is missing
        if ma.isMaskedArray(a1prof):
            amskMiss = a1prof.mask
        else:
            amskMiss = False


        # overlay masks
        if (amskMiss is False) or (amskMiss is True):
            amskMiss = array([False]*len(a1prof))
        if (amskRH is False) or (amskMiss is True):
            amskRH    = array([False]*len(a1prof))

        amsk  = amskMiss + amskRH

        #-- bias correction factor --
        if corrFlag=='CR':
            print dh,'in CR'
            if dh==-99:
                a1gvTmp = a2gv[:,15]
                a1gvTmp  = ma.masked_where(amsk, a1gvTmp).compressed()
                a1profTmp= ma.masked_where(amsk, a1prof ).compressed()
                bfactor = a1gvTmp.mean() / a1profTmp.mean()
        #--- averating time loop ---
        for it,dt in enumerate(ldt):
            #a1idx15 = ones(a2gv.shape[0])*15
            #a1gv   = gv_fsub.mean_slice_negativemask(a2gv.T, a1idx15, dt)
            a1gv   = ma.masked_less(a2gv[:,15:15+dt],0).mean(axis=1)
    
            a1profTmp = ma.masked_where(amsk, a1prof ).compressed()
            a1gvTmp   = ma.masked_where(amsk, a1gv   ).compressed()

            #-- check miss --
            amskMiss1 = ma.masked_less(a1profTmp,0).mask
            amskMiss2 = ma.masked_less(a1gvTmp, 0).mask
            amskMiss  = amskMiss1 + amskMiss2
            a1profTmp = ma.masked_where(amskMiss, a1profTmp)
            a1gvTmp   = ma.masked_where(amskMiss, a1gvTmp)

            #-- bias correction --
            if corrFlag=='CR':
                print rhtype, dh,it, bfactor
                a1profTmp = a1profTmp * bfactor
            #--------------------- 
             
            cc     = np.corrcoef(a1profTmp, a1gvTmp)[0,1]
            rmse   = np.sqrt( ((a1profTmp- a1gvTmp)**2).mean() )
           
            bias   = (a1profTmp - a1gvTmp).mean()
            brat   = (a1profTmp - a1gvTmp).mean() / a1gvTmp.mean()*100
            num    = len(a1profTmp)
            rain   = a1profTmp.mean()
            gv     = a1gvTmp.mean()
  
            if dh==-99:
                de1cc  [rhtype][it] = cc
                de1rmse[rhtype][it] = rmse
                de1bias[rhtype][it] = bias
                de1brat[rhtype][it] = brat
                de1num [rhtype][it] = num       
                de1rain[rhtype][it] = rain
                de1gv  [rhtype][it] = gv

            else: 
                da2cc  [rhtype][ih, it] = cc
                da2rmse[rhtype][ih, it] = rmse
                da2bias[rhtype][ih, it] = bias
                da2brat[rhtype][ih, it] = brat
                da2num [rhtype][ih, it] = num       
                da2rain[rhtype][ih, it] = rain
                da2gv  [rhtype][ih, it] = gv
    
#- save data ----
for rhtype in lrhtype:
    if calc == False: continue

    outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-%s/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])
    util.mk_dir(outDir)
    
    for dattype in ldattype:  
        datPath= outDir + '/nt-nlev.%s.%s.RH.%s.npy'%(dattype, corrFlag, rhtype)
        esPath = outDir + '/nt-eSurf.%s.%s.RH.%s.npy'%(dattype, corrFlag, rhtype)
    
        if dattype=='cc':
            a2dat = da2cc[rhtype]
            e1dat = de1cc[rhtype]
        elif dattype=='rmse':
            a2dat = da2rmse[rhtype]
            e1dat = de1rmse[rhtype]
        elif dattype=='bias':
            a2dat = da2bias[rhtype]
            e1dat = de1bias[rhtype]
        elif dattype=='brat':
            a2dat = da2brat[rhtype]
            e1dat = de1brat[rhtype]
        elif dattype=='num':
            a2dat = da2num[rhtype]
            e1dat = de1num[rhtype]
        elif dattype=='rain':
            a2dat = da2rain[rhtype]
            e1dat = de1rain[rhtype]
        elif dattype=='gv':
            a2dat = da2gv[rhtype]
            e1dat = de1gv[rhtype]
    
    
        np.save(datPath, a2dat)
        np.save(esPath, e1dat)
        print datPath
        print datPath
 


#-- Figure: contour -----
for rhtype in lrhtype:
    if figtype !='cont': continue

    outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-%s/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])

    for dattype in ldattype:
        if dattype=='num':
            print 'skip num for cont figure'
            continue

        datPath=outDir + '/nt-nlev.%s.%s.RH.%s.npy'%(dattype, corrFlag, rhtype)
        esPath =outDir + '/nt-eSurf.%s.%s.RH.%s.npy'%(dattype, corrFlag, rhtype)


        a1es   = np.load(esPath)
        a2dat  = np.load(datPath)

        fig = plt.figure(figsize=(6,4))
        ax  = fig.add_axes([0.15, 0.1, 0.8, 0.78])
        a1t  = ldt

        a2fig= concatenate([a1es.reshape(1,-1), a2dat], axis=0)

        a1x  = ldt
        a1y  = array([0]+ldh)*0.25
        X,Y  = meshgrid(a1x, a1y)

        #im = ax.imshow(a2fig, origin='lower')
        #plt.colorbar(im)
        im = ax.contour(X,Y,a2fig, colors='k')
        plt.clabel(im, inline=1, fontisze=10)

        stitle  = '%s %s'%(dattype, rhtype)
        stitle  = stitle + '\n' + '%04d.%02d-%04d.%02d BiasCorr=%s'%(iYM[0],iYM[1],eYM[0],eYM[1], corrFlag)

        plt.title(stitle)
        figDir  = '/work/a01/utsumi/GPMGV/fig'
        figPath = figDir + '/cont.nt-nlev.%s.%.1fkm.minNum.%d.s.%s.RH.%s.png'%(prdName, thdist,minNum, dattype, corrFlag, rhtype)
        plt.savefig(figPath)
        print figPath
        plt.clf()



 

#--- Figure -------
for rhtype in lrhtype:
    if figtype !='bar': continue

    outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-%s/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])

    for dattype in ldattype:

        datPath=outDir + '/nt-nlev.%s.%s.RH.%s.npy'%(dattype, corrFlag, rhtype)
        esPath =outDir + '/nt-eSurf.%s.%s.RH.%s.npy'%(dattype, corrFlag, rhtype)
        a1es   = np.load(esPath)
        a2dat  = np.load(datPath)

        fig = plt.figure(figsize=(8,4))
        ax  = fig.add_axes([0.10, 0.22, 0.85, 0.66])
        a1t  = ldt

        wbar = 0.7/(nh+1)
        print '-'*50
        # plot eSurf ---
        a1y  = a1es
        a1x  = arange(nt) -0.7*0.5
        ax.bar(a1x, a1y, width=wbar, tick_label=ldt, align='center', label='eSurf')

        # plot average --
        for ih, dh in enumerate(ldh):
            a1y = a2dat[ih,:]
            a1y = ma.masked_invalid(a1y)
            a1x = arange(nt) -0.7*0.5 +(ih+1)*wbar
            ax.bar(a1x, a1y, width=wbar, tick_label=ldt, align='center', label='%.1fkm'%(dh*0.25))

        # x-ticks
        plt.xticks(range(nt), ldt)

        # ylim
        if dattype in ['cc']:
            if rhtype in ['all','hum','dry']:
                ymin = 0.38
                ymax = 0.68
            else:
                ymin = 0.38
                ymax = 0.6
            plt.ylim([ymin,ymax])

        elif dattype in ['rain','gv','num']:
            ymin = 0
            ymax = ma.masked_invalid(a2dat).max()*1.1
            plt.ylim([ymin,ymax])

        elif dattype in ['rmse']:
            ymax = 3.8
            ymin = 2
            plt.ylim([ymin,ymax])



        elif dattype in ['bias']:
            ymax = 0.1
            ymin = -0.3
            plt.ylim([ymin,ymax])

        elif dattype in ['brat']:
            ymax = 10
            ymin = -38 
            plt.ylim([ymin,ymax])


        #- zero line ---
        if dattype in ['cc','bias','brat']:
            plt.plot([-10,50],[0,0],'--',color='k', linewidth=0.5)

        #plt.xlim([a1x[0]-1,a1x[-1]+3.5])
        plt.xlim([a1x[0]-1,a1x[-1]+0.5])

        #- legend
        plt.legend(bbox_to_anchor=(-0.05, -0.25, 1.1,  0.1), loc='upper left',
                    ncol=9, mode='expand', handletextpad=0.)

        stitle  = '%s %s %s'%(dattype, corrFlag, rhtype)
        stitle  = stitle + '\n' + '%04d.%02d-%04d.%02d'%(iYM[0],iYM[1],eYM[0],eYM[1])
        plt.title(stitle)
        figDir  = '/work/a01/utsumi/GPMGV/fig'
        figPath = figDir + '/bar.nt-nlev.%s.%.1fkm.minNum.%d.%s.%s.RH.%s.png'%(prdName, thdist,minNum, dattype, corrFlag, rhtype)
        plt.savefig(figPath)
        print figPath
        plt.clf()


#-- table ----
for rhtype in lrhtype:
    if figtype !='bar': continue

    outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-%s/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])

    for dattype in ldattype:

        datPath=outDir + '/nt-nlev.%s.%s.RH.%s.npy'%(dattype, corrFlag, rhtype)
        esPath =outDir + '/nt-eSurf.%s.%s.RH.%s.npy'%(dattype, corrFlag, rhtype)
        a1es   = np.load(esPath)
        a2dat  = np.load(datPath)

        slabel = ','+','.join(map(str,ldt))
        sout   = slabel + '\n'
        sout   = sout + 'eSurf,' + ','.join(map(str, a1es)) + '\n'
        for ih in range(nh):
            sout = sout + '%.1fkm,'%(ldh[ih]*0.25) + ','.join(map(str,a2dat[ih])) + '\n'


        figDir  = '/work/a01/utsumi/GPMGV/fig'
        csvPath = figDir + '/table.nt-nlev.%s.%.1fkm.minNum.%d.%s.%s.RH.%s.csv'%(prdName, thdist,minNum, dattype, corrFlag, rhtype)

        f = open(csvPath,'w'); f.write(sout); f.close()
        print csvPath
