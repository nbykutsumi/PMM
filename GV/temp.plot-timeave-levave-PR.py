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
#iYM = [2005,9]
iYM = [2005,4]
eYM = [2014,10]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
print lYM

#thdist = 2.5
thdist = 5.0
minNum = 3

cls = 'RH'
#cls = 'RainType'

prdName = 'L2A25'
gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()

dgName = gv.ret_ddomYM2gName()

ldomain = gv.domains
#ldomain = ['FLORIDA-STJ','FLORIDA-SFL-N','N.Carolina-IPHEx_Duke','N.Carolina-IPHEx_NASA','KWAJALEIN-KWA']
#ldomain = ['FLORIDA-STJ']

offset_bef = 15  # 'bef' should be identical to that in mk_match.py
offset_aft = 29


nt = 30
#nh = 20
nh = 30

#lprtype = ['heavy','extreme','mod']
#lprtype = ['all','light','mod','heavy']
lprtype = ['all']
dlthpr = {'all':[-0.1,999],'light':[-0.1,2],'mod':[2,10], 'heavy':[10,50],'extreme':[50,9999]}
ldattype = ['rain','cc','bias','brat','rmse','num','gv']
#ldattype = ['rain']

if cls=='RH':
    lclstype = ['all','dry','hum']
elif cls=='RainType':
    lclstype = ['alltype','strat','conv']
else:
    print 'check cls',cls
    sys.exit()
#----------------------------------
def load_csv(srcPath):
    f=open(srcPath,'r'); lines=f.readlines(); f.close()
    ny = len(lines)
    nx = len(lines[0].strip().split(','))
    a2dat = empty([ny,nx],float32)
    for i,line in enumerate(lines):
        a1dat = map(float, line.strip().split(','))
        a2dat[i] = a1dat

    return a2dat

#----------------------------------
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
for clstype in lclstype:

    a2rain =  empty([nh+1,nt]) 
    a2gv   =  empty([nh+1,nt]) 
    a2cc   =  empty([nh+1,nt]) 
    a2rmse =  empty([nh+1,nt]) 
    a2bias =  empty([nh+1,nt]) 
    a2brat =  empty([nh+1,nt]) 
    a2num  =  empty([nh+1,nt],int32) 


    csvDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-%s/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])

    profavePath  = csvDir + '/%s.profave.%s.csv'%(cls, clstype)
    gvavePath    = csvDir + '/%s.gvave.%s.csv'%(cls, clstype)

    a2profave  = load_csv(profavePath)
    a2gvave    = load_csv(gvavePath)

    for ih in range(nh+1): 
        a1prof = a2profave[:,ih] # mm/h
    
        for it in range(nt):
            a1gv = a2gvave[:,it]
    
            a1profTmp = ma.masked_less(a1prof, 0 )
            a1gvTmp   = ma.masked_less(a1gv, 0 )

            a1profTmp = ma.masked_invalid(a1profTmp)      
            a1gvTmp   = ma.masked_invalid(a1gv)
 
            cc     = np.ma.corrcoef(a1profTmp, a1gvTmp, allow_masked=True)[0,1]
            rmse   = np.sqrt( ((a1profTmp- a1gvTmp)**2).mean() )
           
            bias   = (a1profTmp - a1gvTmp).mean()
            brat   = (a1profTmp - a1gvTmp).mean() / a1gvTmp.mean()
            num    = len(a1profTmp)
            rain   = a1profTmp.mean()
            gv     = a1gvTmp.mean()
    
            a2cc  [ih, it] = cc
            a2rmse[ih, it] = rmse
            a2bias[ih, it] = bias
            a2brat[ih, it] = brat
            a2num [ih, it] = num       
            a2rain[ih, it] = rain
            a2gv  [ih, it] = gv
        
   
    #-- Figure in one panel ------------
    for dattype in ldattype:
        if dattype == 'cc'  : a2dat = a2cc
        if dattype == 'rmse': a2dat = a2rmse
        if dattype == 'bias': a2dat = a2bias
        if dattype == 'num' : a2dat = a2num
        if dattype == 'rain': a2dat = a2rain
        if dattype == 'gv'  : a2dat = a2gv 

        fig = plt.figure(figsize=(4,5))
        ax  = fig.add_axes([0.15, 0.1, 0.8, 0.8])
    
        a1t  = range(nt)
    
        cmap = plt.get_cmap('coolwarm')
    
        # eSurf
        a1y = ma.masked_invalid(a2dat[0,:])
        a1y = running_mean(a1y)
        ax.plot(a1t, a1y, '--', zorder=10, label='eSurf', color='k', linewidth=2)
        #-- find peak or bottom for eSurf ---
        imin = a1y.argmin()
        imax = a1y.argmax()

        if dattype=='cc':
            #ax.plot(a1t[imax],0.01,'D',c='k', markersize=12)
            ax.plot(a1t[imax],a1y[imax],'D',c='k', markersize=12)
 

    
        #lh = range(nh)[2:]
        lh = [2,4,8,12,16,20,24,28,30]
        #lh = [4,8,12,16,20,24]
    
        for itmp,ih in enumerate(lh):
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
            if clstype=='dry':
                ymin=0.2
                ymax=0.7
            else:
                ymin = 0.2
                ymax = 0.7
            plt.ylim([ymin,ymax])
    
        elif dattype in ['rain','gv','rmse','num']:
            ymin = 0
            ymax = ma.masked_invalid(a2dat).max()*1.1
            plt.ylim([ymin,ymax])
    
        elif dattype in ['bias']:
            ymax = abs(ma.masked_invalid(a2dat)).max()*1.1
            ymin = -ymax
            plt.ylim([ymin,ymax])
        elif dattype in ['brat']:
            ymax = 1.8
            ymin = -1.8
            plt.ylim([ymin,ymax])
    
    
        #- zero line ---
        if dattype in ['cc','bias','brat']:
            plt.plot([-10,50],[0,0],'--',color='k', linewidth=0.5)
    
        plt.xlim([-6,31])
    
    
        # title
        stitle  = 'TAve-LevAve %s %s %.1fkm minNum=%d'%(dattype, clstype, thdist, minNum)
        stitle  = stitle + '\n' + '%04d.%02d-%04d.%02d'%(iYM[0],iYM[1],eYM[0],eYM[1])
        plt.title(stitle)
        figDir  = '/work/a01/utsumi/GPMGV/fig'
        figPath = figDir + '/plot.nt-nlev.%s.runmean.%s.%.1fkm.minNum.%d.%s.%s.png'%(cls, prdName, thdist, minNum, dattype, clstype)
        plt.savefig(figPath)
        print figPath
        plt.close()


