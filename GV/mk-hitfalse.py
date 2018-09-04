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

#corrFlag= 'CR'
corrFlag= 'NO'

nozero = 'withzero'

dt  = 30

thpr = 0.1   # mm/h
#thpr = 0.254   # mm/h
#thpr = 20.0   # mm/h

cls = 'RH'
#cls = 'RainType'

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

#ldh  = [0] + range(2,34+1,2)
ldh  = [0] + [2,4,6,8,10,12,14,16,18]


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


ldattype = ['rain','cc','bias','brat','rmse','gv','num']

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


    #----------------------------
    if thpr !=0:
        a2profave = ma.masked_less(a2joinprof, thpr).filled(0.0)
        a2gvave   = ma.masked_less(a2gvave, thpr).filled(0.0) 
        


    #----------------------------

    a1miss = empty([nh])
    a1hit  = empty([nh])
    a1false= empty([nh])

    for ih,dh in enumerate(ldh): 
        a1profTmp = a2profave[:,dh]
        a1gvTmp   = a2gvave[:,dt-1]

        amskTmp1 = ma.masked_less(a1profTmp,0).mask
        amskTmp2 = ma.masked_less(a1gvTmp, 0).mask
        amskTmp  = amskTmp1 + amskTmp2

        a1profTmp = ma.masked_where(amskTmp, a1profTmp).compressed()
        a1gvTmp   = ma.masked_where(amskTmp, a1gvTmp).compressed()

        #-- hit --
        a1tmp  = arange(len(a1profTmp))
        a1flag1= ma.masked_where(a1profTmp>0, a1tmp).mask
        a1flag2= ma.masked_where(a1gvTmp>0, a1tmp).mask
        a1flag = a1flag1*a1flag2
        hit    = a1flag.sum()

        #-- false --
        a1tmp  = arange(len(a1profTmp))
        a1flag1= ma.masked_where(a1profTmp>0, a1tmp).mask
        a1flag2= ma.masked_where(a1gvTmp==0, a1tmp).mask
        a1flag = a1flag1*a1flag2
        false  = a1flag.sum()

        #-- miss --
        a1tmp  = arange(len(a1profTmp))
        a1flag1= ma.masked_where(a1profTmp==0, a1tmp).mask
        a1flag2= ma.masked_where(a1gvTmp>0, a1tmp).mask
        a1flag = a1flag1*a1flag2
        miss   = a1flag.sum()

        #-- negc --
        a1tmp  = arange(len(a1profTmp))
        a1flag1= ma.masked_where(a1profTmp==0, a1tmp).mask
        a1flag2= ma.masked_where(a1gvTmp==0, a1tmp).mask
        a1flag = a1flag1*a1flag2
        negc   = a1flag.sum()     

        #-- total --
        total  = float(len(a1tmp))
        #total  = float(hit + false + miss)

        rhit   = hit  / total * 100
        rfalse = false/ total * 100
        rmiss  = miss / total * 100
        rnegc  = negc / total * 100

        if dh==0:
            hit0   = hit
            false0 = false
            miss0  = miss
            negc0  = negc
            total0 = total

            rhit0  = rhit
            rfalse0= rfalse
            rmiss0 = rmiss
            rnegc0 = rnegc


        print ''
        print 'thpr=',thpr,'dh=',dh*0.25
        print 'hit miss false neg0'
        #print hit,false,miss,negc
        print hit+false+miss+negc
        #print '%d %d %d %d'%(hit,miss, false, negc)
        print '%.4f %.4f %.4f %.4f'%(rhit,rmiss, rfalse, rnegc)
        print '%.4f %.4f %.4f %.4f'%(rhit-rhit0,rmiss-rmiss0,rfalse-rfalse0,rnegc-rnegc0)


    sys.exit()
    #--- Figure -------
    for dattype in ldattype:
        if dattype =='rain': a2dat = a2rain
        if dattype =='cc':   a2dat = a2cc
        if dattype =='bias': a2dat = a2bias
        if dattype =='brat': a2dat = a2brat
        if dattype =='rmse': a2dat = a2rmse
        if dattype =='gv':   a2dat = a2gv
        if dattype =='num':  a2dat = a2num
        
        fig = plt.figure(figsize=(8,4))
        ax  = fig.add_axes([0.10, 0.22, 0.85, 0.66])
        a1t  = ldt

        wbar = 0.7/(nh+1)
        # plot eSurf ---
        a1y  = a2dat[0,:]
        a1x  = arange(nt) -0.7*0.5
        ax.bar(a1x, a1y, width=wbar, tick_label=ldt, align='center', label='eSurf')

        # plot average --
        for ih, dh in enumerate(ldh):
            if dh==0: continue

            a1y = a2dat[ih,:]
            a1y = ma.masked_invalid(a1y)
            a1x = arange(nt) -0.7*0.5 +(ih+1)*wbar
            ax.bar(a1x, a1y, width=wbar, tick_label=ldt, align='center', label='%.1fkm'%(dh*0.25))

        # x-ticks
        plt.xticks(range(nt), ldt)

        # ylim
        if dattype in ['cc']:
            ymax = 0.67
            ymin = 0.33
            plt.ylim([ymin,ymax])

        elif dattype in ['rain','gv','num']:
            ymax = ma.masked_invalid(a2dat).max()*1.1
            ymin = 0
            plt.ylim([ymin,ymax])

        elif dattype in ['rmse']:
            ymax = 5.9
            ymin = 3.2
            plt.ylim([ymin,ymax])

        elif dattype in ['bias']:
            ymax = 0.3
            ymin = -1.0
            plt.ylim([ymin,ymax])

        elif dattype in ['brat']:
            ymax = 0.15
            ymin = -0.45
            plt.ylim([ymin,ymax])


        #- zero line ---
        if dattype in ['cc','bias','brat']:
            plt.plot([-10,50],[0,0],'--',color='k', linewidth=0.5)

        #plt.xlim([a1x[0]-1,a1x[-1]+3.5])
        plt.xlim([a1x[0]-1,a1x[-1]+0.5])

        #- legend
        plt.legend(bbox_to_anchor=(-0.05, -0.25, 1.1,  0.1), loc='upper left',
                    ncol=9, mode='expand', handletextpad=0.)

        stitle  = '%s %s %s'%(dattype, corrFlag, clstype)
        stitle  = stitle + '\n' + '%04d.%02d-%04d.%02d'%(iYM[0],iYM[1],eYM[0],eYM[1])
        plt.title(stitle)
        figDir  = '/work/a01/utsumi/GPMGV/fig'
        figPath = figDir + '/bar.nt-nlev.%s.%.1fkm.minNum.%d.%s.%s.%s.%s.png'%(prdName, thdist,minNum, dattype, corrFlag, cls, clstype)
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
