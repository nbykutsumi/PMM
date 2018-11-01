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
#iYM = [2005,4]
iYM = [2014,9]
eYM = [2014,10]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
print lYM
figtype = 'bar'
#figtype = 'cont'
corrFlag= 'CR'
#corrFlag= 'NO'

#cls = 'RH'
#cls = 'RainType'
cls = 'ALL'

#thdist = 2.5
thdist = 5.0
minNum = 3
prdName = 'L2A25'

nozero = 'nozero'
#nozero = 'withzero'
miss   = -9999.
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


#if cls=='RH':
#    lclstype = ['all','dry','hum']
#elif cls=='RainType':
#    lclstype = ['alltype','strat','conv']
#else:
#    print 'check cls',cls
#    sys.exit()


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
csvDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-%s/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])

clstype      = 'ALL'
joinprofPath = csvDir + '/ALL.%s.joinprof.csv'%(nozero)
gvPath       = csvDir + '/ALL.%s.gv.csv'%(nozero)
profavePath  = csvDir + '/ALL.%s.profave.csv'%(nozero)
gvavePath    = csvDir + '/ALL.%s.gvave.csv'%(nozero)
rhPath       = csvDir + '/ALL.%s.RH.csv'%(nozero)
rtypePath    = csvDir + '/ALL.%s.rainType.csv'%(nozero)
stormHPath   = csvDir + '/ALL.%s.stormH.csv'%(nozero)

a2joinprof = load_csv(joinprofPath)
a2gv       = load_csv(gvPath)[:,15:]
a2profave  = load_csv(profavePath)
a2gvave    = load_csv(gvavePath)
a1rh       = load_csv(rhPath)
a1rtype    = load_csv(rtypePath)
a1stormH   = load_csv(stormHPath)

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
#-- profile adjustment based on climatology --
lhup     = [2,3,4,5,6,7,8]
llstormH = [[2,3],[3,4],[4,6],[6,8],[8,15]]
ratDir   = '/home/utsumi/mnt/wellshare/GPMGV/DOM.L2A25/cbin.11/ratio'
llRH     = [[0.0,0.6],[0.6,0.7],[0.7,0.8],[0.8,0.9],[0.9,1.0]]

a1height = arange(a2joinprof.shape[1])*0.25
dratio   = {}
for raintype in ['all','strat','conv']:
    for istormH,lstormH in enumerate(llstormH):
        stormH0,stormH1 = lstormH
        ratPath = ratDir + '/ratio.%s.H%d-%d.npy'%(raintype,stormH0,stormH1)
        a2tmp   = np.load(ratPath)
        for iRH,lRH in enumerate(llRH):
            a1ratioTmp =  ones(a2joinprof.shape[1])*miss
            a1idx      =  arange(a2joinprof.shape[1])
            for ihup,hup in enumerate(lhup):
                if ihup == len(lhup):
                    hup0 = hup - 0.5
                    hup1 = 99
                else:
                    hup0 = hup - 0.5
                    hup1 = hup + 0.5

                a1msk   = ma.masked_less(a1height,hup0)
                a1msk   = ma.masked_greater_equal(a1msk, hup1)
                a1idxTmp= ma.masked_where(a1msk.mask, a1idx).compressed()
                a1ratioTmp[a1idxTmp] = a2tmp[ihup,iRH]


            #a1ratioTmp  = ma.masked_equal(a1ratioTmp,miss).filled(1.0)
            #-- fill lower and higher ends --
            a1idxTmp = ma.masked_where(a1ratioTmp==miss,a1idx).compressed()
            idxmin   = a1idxTmp.min()
            idxmax   = a1idxTmp.max()
            a1ratioTmp[:idxmin] = 1.0
            a1ratioTmp[idxmax+1:] = a1ratioTmp[idxmax]


            dratio[raintype,istormH,iRH] = a1ratioTmp
            if (raintype=='conv')&(istormH==2)&(iRH==0):
                print ''
                print '***********************'
                print dratio[raintype,istormH,iRH] 
                print ''
                print a2tmp[:,0]
                print llstormH
                print lhup
                sys.exit()
sys.exit()
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
        if clstype in ['conv','dry']:
            ymax = 0.45
            ymin = -0.45
        else:
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
