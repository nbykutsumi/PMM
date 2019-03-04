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
corrFlag= 'CR'
#corrFlag= 'NO'

cls = 'RH'
#cls = 'RainType'

#thdist = 2.5
thdist = 5.0
minNum = 3
prdName = 'L2A25'

nozero = 'nozero'
#nozero = 'withzero'

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
    ldh  = [0] + range(2,34+1,2)
elif figtype=='cont':
    ldt  = range(1,30+1)
    ldh  = range(1,34+1)


nt = len(ldt)
nh = len(ldh)
#nh = 5


if cls=='RH':
    #lclstype = ['all','dry','hum']
    lclstype = ['all']
elif cls=='RainType':
    lclstype = ['alltype','strat','conv']
else:
    print 'check cls',cls
    sys.exit()


#ldattype = ['rain','cc','bias','brat','rmse','gv','num']
ldattype = ['brat','rmse']
#ldattype = ['cc','brat','rmse']

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
for clstype in lclstype:

    csvDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-%s/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])

    joinprofPath = csvDir + '/%s.%s.joinprof.%s.csv'%(cls, nozero, clstype)
    gvPath       = csvDir + '/%s.%s.gv.%s.csv'%(cls, nozero, clstype)
    profavePath  = csvDir + '/%s.%s.profave.%s.csv'%(cls, nozero, clstype)
    gvavePath    = csvDir + '/%s.%s.gvave.%s.csv'%(cls, nozero, clstype)
   
    a2joinprof = load_csv(joinprofPath)
    a2gv       = load_csv(gvPath)[:,15:]
    a2profave  = load_csv(profavePath)
    a2gvave    = load_csv(gvavePath)


    #-- bias correction factor --
    if corrFlag=='CR':
        amskTmp1 = ma.masked_less(a2gv[:,0], 0).mask
        amskTmp2 = ma.masked_less(a2joinprof[:,0],0 ).mask
        amskTmp  = amskTmp1 + amskTmp2

        a1gvTmp  = ma.masked_where(amskTmp, a2gv[:,0])
        a1profTmp= ma.masked_where(amskTmp, a2joinprof[:,0])
        bfactor = a1gvTmp.mean() / a1profTmp.mean()

        a2joinprof = a2joinprof * bfactor
        a2profave  = a2profave * bfactor
    #----------------------------

    a2rain = empty([nh,nt]) 
    a2gv   = empty([nh,nt]) 
    a2cc   = empty([nh,nt]) 
    a2rmse = empty([nh,nt]) 
    a2bias = empty([nh,nt]) 
    a2brat = empty([nh,nt]) 
    a2num  = empty([nh,nt],int32)

    for ih,dh in enumerate(ldh): 
        for it,dt in enumerate(ldt):
            a1profTmp = a2profave[:,dh]
            a1gvTmp   = a2gvave[:,dt-1]

            a1profTmp = ma.masked_less(a1profTmp,0)
            a1gvTmp   = ma.masked_less(a1gvTmp, 0)

            cc     = np.ma.corrcoef(a1profTmp, a1gvTmp, allow_masked=True)[0,1]  # allow_masked_True is important
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
    

    #--- Figure -------
    for dattype in ldattype:
        if dattype =='rain': a2dat = a2rain
        if dattype =='cc':   a2dat = a2cc
        if dattype =='bias': a2dat = a2bias
        if dattype =='brat': a2dat = a2brat
        if dattype =='rmse': a2dat = a2rmse
        if dattype =='gv':   a2dat = a2gv
        if dattype =='num':  a2dat = a2num
        

        #-- One panel for each dt ---
        fig = plt.figure(figsize=(8,2.3))

        for idt,dt in enumerate(ldt):
            wax= 0.8/len(ldt)
            ax = fig.add_axes([0.1+wax*idt, 0.19, wax, 0.71])
            #-- plot --
            a1y  = range(1,len(ldh))
            a1x  = a2dat[1:,idt]
            #ax.plot(a1x,a1y, marker='o',color='k',markerfacecolor='k', markersize=2,linewidth=1)
            ax.plot(a1x,a1y, marker='o',color='k',markerfacecolor='k', markersize=1.5,linewidth=0)
            #ax.plot(a1x,a1y, '-', color='k', linewidth=1)

            #-- plot best --
            if dattype in ['cc']:
                ibest  = np.argmax(a1x)
                xbest  = a1x[ibest]
                ybest  = a1y[ibest]
                ax.plot(xbest, ybest, marker='o',color='k',markerfacecolor='k', markersize=6)

            elif dattype in ['brat']:
                # find top-2 absolute minima
                aidx   = range(len(a1x))
                a2tmp  = zip( aidx, abs(a1x) )
                a2sort = array(sorted( a2tmp, key=lambda x:x[1]))
                aidx_sort=a2sort[:,0].astype(int32)
                ibest0 = aidx_sort[0]
                ibest1 = aidx_sort[1]
                                    
                xbest0 = a1x[ibest0]
                xbest1 = a1x[ibest1]
                ybest0 = a1y[ibest0]
                ybest1 = a1y[ibest1]

                ax.plot(xbest0, ybest0, marker='o',color='k',markerfacecolor='k', markersize=6)
                ax.plot(xbest1, ybest1, marker='o',color='k',markerfacecolor='k', markersize=6)

            elif dattype in ['rmse']:
                ibest  = np.argmin(a1x)
                xbest  = a1x[ibest]
                ybest  = a1y[ibest]
                ax.plot(xbest, ybest, marker='o',color='k',markerfacecolor='k', markersize=6)


            # plot eSurf ---
            y  = 0
            x  = a2dat[0,idt]
            ax.plot(x,y, 'o', color='k',markerfacecolor='none',markeredgewidth=2)
   
            # xlim
            if dattype in ['cc']:
                xmax = 0.60
                xmin = 0.48
                plt.xticks([0.5,0.55])
                plt.xlim([xmin,xmax])
    
            elif dattype in ['rain','gv','num']:
                xmax = ma.masked_invalid(a2dat).max()*1.1
                xmin = 0
                plt.xlim([xmin,xmax])
    
            elif dattype in ['rmse']:
                xmax = 4.8
                xmin = 3.8
                plt.xticks([4.0,4.5])
                plt.xlim([xmin,xmax])
    
            elif dattype in ['bias']:
                xmax = 0.3
                xmin = -1.0
                plt.xticks([-0.5,0])
                plt.xlim([xmin,xmax])
    
            elif dattype in ['brat']:
                if clstype in ['conv','dry']:
                    xmax = 0.45
                    xmin = -0.45
                    plt.xticks([-0.3,0,0.3])
                else:
                    xmax = 0.15
                    xmin = -0.45
                    plt.xticks([-0.3, 0])

                plt.xlim([xmin,xmax])
 
            
            # y-ticks
            ax.set_yticks([0]+range(2,20,2))
            ax.yaxis.set_ticks_position('both')
            if idt==0:
                ax.set_yticklabels(['eSurf']+list(arange(2,20,2)*0.5))
            elif idt==len(ldt)-1:
                ax.set_yticklabels([])

            else:
                ax.yaxis.set_ticks_position('both')
                ax.set_yticklabels([])
 
            #- zero line ---
            if dattype in ['cc','bias','brat']:
                plt.plot([0,0],[-1,len(ldh)+4],'--',color='k', linewidth=0.8)
   
            #-- ylim -------
            plt.ylim([-1,len(ldh)+2]) 

            #-- text --
            stext = '%d-min ave.'%(dt)
            ax.text(0.08,0.92, stext, transform=ax.transAxes)
   
        #-- Title --
        if dattype=='cc':
            dattype_long='Correlation Coefficient'
        elif dattype=='brat':
            dattype_long='Bias Ratio'
        elif dattype=='rmse':
            dattype_long='RMSE'
        else:
            dattype_long=dattype

        if corrFlag=='CR':
            corr_long= 'with bias correction'
        elif corrFlag=='NO':
            corr_long= 'no bias correction'

        stitle = '%s (%s)'%(dattype_long, corr_long)
        fig.suptitle(stitle)

        #-- common axis labels --
        fig.text(0.01, 0.55, 'Averaging height [km]', va='center', rotation='vertical', fontsize=12)
        fig.text(0.5, 0.02, '%s'%(dattype_long), ha='center', fontsize=12)

        #-- Save file ---
        figDir  = '/work/a01/utsumi/GPMGV/fig'
        figPath = figDir + '/line.nt-nlev.%s.%.1fkm.minNum.%d.%s.%s.%s.%s.png'%(prdName, thdist,minNum, dattype, corrFlag, cls, clstype)
        plt.savefig(figPath)
        print figPath
        plt.clf()
        
    #-- table ----
    for dattype in ldattype:
        if dattype =='rain': a2dat = a2rain
        if dattype =='cc':   a2dat = a2cc
        if dattype =='bias': a2dat = a2bias
        if dattype =='brat': a2dat = a2brat
        if dattype =='rmse': a2dat = a2rmse
        if dattype =='gv':   a2dat = a2gv
        if dattype =='num':  a2dat = a2num
 

        slabel = ','+','.join(map(str,ldt))
        sout   = slabel + '\n'
        sout   = sout + 'eSurf,' + ','.join(map(str, a2dat[0,:])) + '\n'
        for ih in range(nh):
            if ih == 0: continue
            sout = sout + '%.1fkm,'%(ldh[ih]*0.25) + ','.join(map(str,a2dat[ih])) + '\n'


        figDir  = '/work/a01/utsumi/GPMGV/fig'
        csvPath = figDir + '/table.nt-nlev.%s.%.1fkm.minNum.%d.%s.%s.%s.%s.csv'%(prdName, thdist,minNum, dattype, corrFlag, cls, clstype)

        f = open(csvPath,'w'); f.write(sout); f.close()
        print csvPath
