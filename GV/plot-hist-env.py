import matplotlib
matplotlib.use('Agg')
from numpy import *
import numpy as np
import myfunc.util as util
from collections import deque
import sys, os
import matplotlib.pyplot as plt
from   datetime import datetime, timedelta
import scipy.stats as stats

#iYM = [2001,4]
#eYM = [2004,8]

iYM = [1998,4]
eYM = [2004,10]



lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
#lYM = lYM[::-1]
print lYM

ldomain = ['FLORIDA-KSC', 'KWAJALEIN-KWA', 'TEXAS-HAR', 'FLORIDA-SFL-N', 'N.Carolina-IPHEx_Duke', 'FLORIDA-STJ', 'N.Carolina-IPHEx_NASA']
#ldomain = ['FLORIDA-KSC', 'KWAJALEIN-KWA', 'TEXAS-HAR', 'FLORIDA-SFL-N']
#ldomain = ['FLORIDA-KSC']

lrtype    = ['all','strat','conv']
dlthrtype = {'all':[-999,999],'strat':[100,170],'conv':[200,297]}


nh  = 25
hcutoff = 1 # km
#thdist = 5 # km
thdist = 10 # km

thcount = 50 # for drawing lower levels

wfig, hfig = 4,3

baseDir = '/home/utsumi/mnt/wellshare/GPMGV/GLOC.L2A25/%.1fkm'%(thdist)
for domain in ldomain:

    a1rh    = array([])
    a1div850= array([])
    a1div10m= array([])
    a2prof  = array([])
    a1rtype = array([])

    for YM in lYM:
        Year,Mon = YM
        srcDir = baseDir + '/%s/%04d%02d'%(domain, Year,Mon)
        rhPath = srcDir  + '/rh.npy'
        div850Path= srcDir  + '/div850.npy'
        div10mPath= srcDir  + '/div10m.npy'
        profPath =srcDir + '/prof.npy'
        rtypePath=srcDir + '/rainType.npy'

        a1rhTmp     = np.load(rhPath)
        a1div850Tmp = np.load(div850Path)
        a1div10mTmp = np.load(div10mPath)
        a2profTmp   = np.load(profPath)
        a1rtypeTmp  = np.load(rtypePath)

        a1max       = a2profTmp.max(axis=1)
        a1idx       = arange(len(a1rhTmp)).astype(int32)
        a1idx       = ma.masked_where(a1max<=0, a1idx).compressed()
        
        a1rhTmp     = a1rhTmp [a1idx]
        a1div850Tmp = a1div850Tmp[a1idx]
        a1div10mTmp = a1div10mTmp[a1idx]
        a2profTmp   = a2profTmp[a1idx,:]
        a1rtypeTmp  = a1rtypeTmp[a1idx]


        a1rh        = concatenate([a1rh,    a1rhTmp    ])
        a1div850    = concatenate([a1div850,a1div850Tmp])
        a1div10m    = concatenate([a1div10m,a1div10mTmp])
        a1rtype     = concatenate([a1rtype, a1rtypeTmp ])

        try:
            a2prof  = concatenate([a2prof,  a2profTmp ], axis=0)
        except ValueError:
            a2prof  = a2profTmp


    a1rh     = ma.masked_less(a1rh, 0)
    a1div850 = ma.masked_equal(a1div850,-9999.)*1e5
    a1div10m = ma.masked_equal(a1div10m,-9999.)*1e5
    a2prof   = a2prof[:,:nh]
    #*********************************************************
    #--- by div 850hPa -------
    #'''
    ddens = {}
    dhist = {}
    dx    = {}
    for rtype in lrtype:
        thrtype0,thrtype1 = dlthrtype[rtype]
        amskT = ma.masked_outside(a1rtype, thrtype0, thrtype1).mask
        a1idx = arange(len(a1div850)).astype(int32)
        a1idx = ma.masked_where(amskT, a1idx)
        a1idx = a1idx.compressed()
        a1dat = a1div850[a1idx]
        a1dat = ma.masked_equal(a1dat,-9999).compressed()

        ddens[rtype] = stats.gaussian_kde(a1dat)
        dhist[rtype], dx[rtype] = np.histogram(a1dat, density=True)
    #-- draw -----------
    lcolor = ['k','b','r']
    #-- draw precip ----
    fig    = plt.figure(figsize=(wfig,hfig))
    ax = fig.add_subplot(111)
    for i, rtype in enumerate(lrtype):
        slabel = lrtype[i]
        a1x    = (array(dx[rtype][:-1]) + array(dx[rtype][1:]))*0.5
        density = ddens[rtype]
        a1hist  = dhist[rtype]
        #ax.plot(a1x, density(a1x), '-',linewidth=2, color=lcolor[i], label=slabel)
        ax.plot(a1x, a1hist, '-',linewidth=2, color=lcolor[i], label=slabel)

    plt.legend(loc='upper right')
    ssuptitle = '%s (div850hPa e-5)'%(domain)
    plt.suptitle(ssuptitle)

    figDir  = '/work/a01/utsumi/GPMGV/fig/env-prof'
    figPath = figDir + '/hist.PR.div850.%s.%.1fkm.png'%(domain,thdist)
    util.mk_dir(figDir)
    plt.savefig(figPath)
    plt.clf()
    print figPath


