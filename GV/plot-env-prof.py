import matplotlib
matplotlib.use('Agg')
from numpy import *
import numpy as np
import myfunc.util as util
from collections import deque
import sys, os
import matplotlib.pyplot as plt
from   datetime import datetime, timedelta

#iYM = [2002,4]
#eYM = [2003,5]

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

cbins   = 11 # number of angle bins
thcount = 50 # for drawing lower levels

wfig, hfig = 6,4

baseDir = '/home/utsumi/mnt/wellshare/GPMGV/DOM.L2A25/cbin.%d'%(cbins)
for domain in ldomain:

    a1esurf = array([])
    a1rh    = array([])
    a1div850= array([])
    a1div10m= array([])
    a1stormH= array([])
    a2prof  = array([])
    a1rtype = array([])

    for YM in lYM:
        Year,Mon = YM
        srcDir = baseDir + '/%s/%04d%02d'%(domain, Year,Mon)
        rhPath = srcDir  + '/rh.npy'
        div850Path= srcDir  + '/div850.npy'
        div10mPath= srcDir  + '/div10m.npy'
        stormHPath= srcDir  + '/stormH.npy'
        esurfPath=srcDir + '/eSurf.npy'
        profPath =srcDir + '/prof.npy'
        rtypePath=srcDir + '/rainType.npy'

        a1esurfTmp  = np.load(esurfPath)
        a1rhTmp     = np.load(rhPath)
        a1div850Tmp = np.load(div850Path)
        a1div10mTmp = np.load(div10mPath)
        a1stormHTmp = np.load(stormHPath)
        a2profTmp   = np.load(profPath)
        a1rtypeTmp  = np.load(rtypePath)

        a1max       = a2profTmp.max(axis=1)
        a1idx       = arange(len(a1esurfTmp)).astype(int32)
        a1idx       = ma.masked_where(a1max<=0, a1idx).compressed()
        
        a1esurfTmp  = a1esurfTmp[a1idx]
        a1rhTmp     = a1rhTmp [a1idx]
        a1div850Tmp = a1div850Tmp[a1idx]
        a1div10mTmp = a1div10mTmp[a1idx]
        a1stormHTmp = a1stormHTmp[a1idx]
        a2profTmp   = a2profTmp[a1idx,:]
        a1rtypeTmp  = a1rtypeTmp[a1idx]


        a1esurf     = concatenate([a1esurf, a1esurfTmp ])
        a1rh        = concatenate([a1rh,    a1rhTmp    ])
        a1div850    = concatenate([a1div850,a1div850Tmp])
        a1div10m    = concatenate([a1div10m,a1div10mTmp])
        a1stormH    = concatenate([a1stormH,a1stormHTmp])
        a1rtype     = concatenate([a1rtype, a1rtypeTmp ])

        try:
            a2prof  = concatenate([a2prof,  a2profTmp ], axis=0)
        except ValueError:
            a2prof  = a2profTmp


    a1esurf  = ma.masked_less(a1esurf,0)
    a1rh     = ma.masked_less(a1rh, 0)
    a1div850 = ma.masked_equal(a1div850,-9999.)*1e5
    a1div10m = ma.masked_equal(a1div10m,-9999.)*1e5
    a1stormH = ma.masked_equal(a1stormH,-9999.)*0.001  # m->km
    a2prof   = a2prof[:,:nh]
    #*********************************************************
    #--- by div 850hPa -------
    #'''

    ldivrange = [[-999,-1],[-2,-1],[-1,1],[1,2],[2,999]]
    da1ave    = {}
    da1std    = {}
    da1num    = {}

    for rtype in lrtype:

        thrtype0,thrtype1 = dlthrtype[rtype]
        amskT = ma.masked_outside(a1rtype, thrtype0, thrtype1).mask

        for i,divrange in enumerate(ldivrange):
            divmin, divmax = divrange
            a1idx = arange(len(a1div850)).astype(int32)
            a1idx = ma.masked_where(a1div850 <divmin, a1idx)
            a1idx = ma.masked_where(divmax <= a1div850, a1idx)
            a1idx = ma.masked_where(amskT, a1idx)
            a1idx = a1idx.compressed()
            a2profTmp = a2prof[a1idx]
            a2profTmp = ma.masked_less(a2profTmp,0)
    
            da1ave[rtype,i] = a2profTmp.mean(axis=0)*0.01
            da1std[rtype,i] = a2profTmp.std(axis=0)*0.01
            da1num[rtype,i] = a2profTmp.count(axis=0)


    #-- draw -----------
    lcolor = ['r','r','k','b','b']
    llinewidth=[2,0.5,1,0.5,2]

    #-- draw precip ----
    fig    = plt.figure(figsize=(wfig,hfig))
    for irtype, rtype in enumerate(lrtype):
        axloc = 100 + 10*len(lrtype) + irtype + 1
        ax = fig.add_subplot( axloc )

        stitle = rtype
        plt.title(stitle)

        for i,divrange in enumerate(ldivrange):
            divmin, divmax = divrange
            slabel = '%d ~ %d'%(divmin, divmax)
            a1y = (arange(nh)+1) * 0.25
            a1x = da1ave[rtype, i]
            a1x = ma.masked_where(da1num[rtype,i]<thcount, a1x)
            a1x = ma.masked_where(a1y<hcutoff, a1x)

            ax.plot(a1x, a1y, '-',linewidth=llinewidth[i], color=lcolor[i], label=slabel)

    plt.legend(loc='upper right')
    ssuptitle = 'rain profile %s (div850hPa)'%(domain)
    plt.suptitle(ssuptitle)

    figDir  = '/work/a01/utsumi/GPMGV/fig/env-prof'
    figPath = figDir + '/prof.PR.div850.prcp.%s.cbin.%d.png'%(domain,cbins)
    util.mk_dir(figDir)
    plt.savefig(figPath)
    plt.clf()
    print figPath

    #-- draw count ----
    fig    = plt.figure(figsize=(wfig,hfig))
    for irtype, rtype in enumerate(lrtype):
        axloc = 100 + 10*len(lrtype) + irtype + 1
        ax = fig.add_subplot( axloc )

        stitle = rtype
        plt.title(stitle)

        for i,divrange in enumerate(ldivrange):
            divmin, divmax = divrange
            slabel = '%d ~ %d'%(divmin, divmax)
            a1y = (arange(nh)+1) * 0.25
            a1x = da1num[rtype, i]
            a1x = ma.masked_where(a1y<hcutoff, a1x)
            ax.plot(a1x, a1y, '-',linewidth=llinewidth[i], color=lcolor[i], label=slabel)

    plt.legend(loc='upper right')
    ssuptitle = 'count %s (div850hPa)'%(domain)
    plt.suptitle(ssuptitle)
    figDir  = '/work/a01/utsumi/GPMGV/fig/env-prof'
    figPath = figDir + '/prof.PR.div850.count.%s.cbin.%d.png'%(domain,cbins)
    util.mk_dir(figDir)
    plt.savefig(figPath)
    plt.clf()
    
    print figPath
    #'''


    #*********************************************************
    #--- by div 10m -------
    #'''

    ldivrange = [[-999,-1],[-2,-1],[-1,1],[1,2],[2,999]]
    da1ave    = {}
    da1std    = {}
    da1num    = {}

    for rtype in lrtype:

        thrtype0,thrtype1 = dlthrtype[rtype]
        amskT = ma.masked_outside(a1rtype, thrtype0, thrtype1).mask

        for i,divrange in enumerate(ldivrange):
            divmin, divmax = divrange
            a1idx = arange(len(a1div10m)).astype(int32)
            a1idx = ma.masked_where(a1div10m <divmin, a1idx)
            a1idx = ma.masked_where(divmax <= a1div10m, a1idx)
            a1idx = ma.masked_where(amskT, a1idx)
            a1idx = a1idx.compressed()
            a2profTmp = a2prof[a1idx]
            a2profTmp = ma.masked_less(a2profTmp,0)
    
            da1ave[rtype,i] = a2profTmp.mean(axis=0)*0.01
            da1std[rtype,i] = a2profTmp.std(axis=0)*0.01
            da1num[rtype,i] = a2profTmp.count(axis=0)


    #-- draw -----------
    lcolor = ['r','r','k','b','b']
    llinewidth=[2,0.5,1,0.5,2]

    #-- draw precip ----
    fig    = plt.figure(figsize=(wfig,hfig))
    for irtype, rtype in enumerate(lrtype):
        axloc = 100 + 10*len(lrtype) + irtype + 1
        ax = fig.add_subplot( axloc )

        stitle = rtype
        plt.title(stitle)

        for i,divrange in enumerate(ldivrange):
            divmin, divmax = divrange
            slabel = '%d ~ %d'%(divmin, divmax)
            a1y = (arange(nh)+1) * 0.25
            a1x = da1ave[rtype, i]
            a1x = ma.masked_where(da1num[rtype,i]<thcount, a1x)
            a1x = ma.masked_where(a1y<hcutoff, a1x)

            ax.plot(a1x, a1y, '-',linewidth=llinewidth[i], color=lcolor[i], label=slabel)

    plt.legend(loc='upper right')
    ssuptitle = 'rain profile %s (div10m)'%(domain)
    plt.suptitle(ssuptitle)

    figDir  = '/work/a01/utsumi/GPMGV/fig/env-prof'
    figPath = figDir + '/prof.PR.div10m.prcp.%s.cbin.%d.png'%(domain,cbins)
    util.mk_dir(figDir)
    plt.savefig(figPath)
    plt.clf()
    print figPath

    #-- draw count ----
    fig    = plt.figure(figsize=(wfig,hfig))
    for irtype, rtype in enumerate(lrtype):
        axloc = 100 + 10*len(lrtype) + irtype + 1
        ax = fig.add_subplot( axloc )

        stitle = rtype
        plt.title(stitle)

        for i,divrange in enumerate(ldivrange):
            divmin, divmax = divrange
            slabel = '%d ~ %d'%(divmin, divmax)
            a1y = (arange(nh)+1) * 0.25
            a1x = da1num[rtype, i]
            a1x = ma.masked_where(a1y<hcutoff, a1x)
            ax.plot(a1x, a1y, '-',linewidth=llinewidth[i], color=lcolor[i], label=slabel)

    plt.legend(loc='upper right')
    ssuptitle = 'count %s (div10m)'%(domain)
    plt.suptitle(ssuptitle)
    figDir  = '/work/a01/utsumi/GPMGV/fig/env-prof'
    figPath = figDir + '/prof.PR.div10m.count.%s.cbin.%d.png'%(domain,cbins)
    util.mk_dir(figDir)
    plt.savefig(figPath)
    plt.clf()
    
    print figPath
    #'''


    #*********************************************************
    #--- by RH -------
    #'''

    #lrhrange  = [[0,0.6],[0.6,0.7],[0.7,0.8],[0.8,0.9],[0.9,1.0]]
    lrhrange  = [[0,0.6],[0.6,0.8],[0.8,1.0]]
    da1ave    = {}
    da1std    = {}
    da1num    = {}

    for rtype in lrtype:
        thrtype0,thrtype1 = dlthrtype[rtype]
        amskT = ma.masked_outside(a1rtype, thrtype0, thrtype1).mask

        for i,rhrange in enumerate(lrhrange):
            rhmin,rhmax = rhrange
            a1idx = arange(len(a1rh)).astype(int32)
            a1idx = ma.masked_where(a1rh <rhmin, a1idx)
            a1idx = ma.masked_where(rhmax <= a1rh, a1idx)
            a1idx = ma.masked_where(amskT, a1idx)
            a1idx = a1idx.compressed()
            a2profTmp = a2prof[a1idx]
            a2profTmp = ma.masked_less(a2profTmp,0)
    
            da1ave[rtype,i] = a2profTmp.mean(axis=0)*0.01
            da1std[rtype,i] = a2profTmp.std(axis=0)*0.01
            da1num[rtype,i] = a2profTmp.count(axis=0)


    #-- draw -----------
    lcolor = ['b','k','r']
    llinewidth=[2,1,2]


    #-- draw precip ----
    fig    = plt.figure(figsize=(wfig,hfig))
    for irtype, rtype in enumerate(lrtype):
        axloc = 100 + 10*len(lrtype) + irtype + 1
        ax = fig.add_subplot( axloc )

        stitle = rtype
        plt.title(stitle)

        for i,rhrange in enumerate(lrhrange):
            rhmin,rhmax = rhrange
            slabel = '%d ~ %d%%'%(rhmin*100, rhmax*100)
            a1y = (arange(nh)+1) * 0.25
            a1x = da1ave[rtype, i]
            a1x = ma.masked_where(da1num[rtype,i]<thcount, a1x)
            a1x = ma.masked_where(a1y<hcutoff, a1x)

            ax.plot(a1x, a1y, '-',linewidth=llinewidth[i], color=lcolor[i], label=slabel)

    plt.legend(loc='upper right')
    ssuptitle = 'rain profile %s (RH)'%(domain)
    plt.suptitle(ssuptitle)

    figDir  = '/work/a01/utsumi/GPMGV/fig/env-prof'
    figPath = figDir + '/prof.PR.RH.prcp.%s.cbin.%d.png'%(domain,cbins)
    util.mk_dir(figDir)
    plt.savefig(figPath)
    plt.clf()
    print figPath

    #-- draw count ----
    fig    = plt.figure(figsize=(wfig,hfig))
    for irtype, rtype in enumerate(lrtype):
        axloc = 100 + 10*len(lrtype) + irtype + 1
        ax = fig.add_subplot( axloc )

        stitle = rtype
        plt.title(stitle)

        for i,rhrange in enumerate(lrhrange):
            rhmin,rhmax = rhrange
            slabel = '%d ~ %d%%'%(rhmin*100, rhmax*100)
            a1y = (arange(nh)+1) * 0.25
            a1x = da1num[rtype, i]
            #a1x = ma.masked_where(da1num[rtype,i]<50, a1x)
            a1x = ma.masked_where(a1y<hcutoff, a1x)
            ax.plot(a1x, a1y, '-',linewidth=llinewidth[i], color=lcolor[i], label=slabel)

    plt.legend(loc='upper right')
    ssuptitle = 'count %s (RH)'%(domain)
    plt.suptitle(ssuptitle)

    figDir  = '/work/a01/utsumi/GPMGV/fig/env-prof'
    figPath = figDir + '/prof.PR.RH.count.%s.cbin.%d.png'%(domain,cbins)
    util.mk_dir(figDir)
    plt.savefig(figPath)
    plt.clf()
    print figPath

    #'''


    #*********************************************************
    #--- by stormH -------
    #'''
    lhrange  = [[0,3],[3,5],[5,15]]
    da1ave    = {}
    da1std    = {}
    da1num    = {}

    for rtype in lrtype:
        thrtype0,thrtype1 = dlthrtype[rtype]
        amskT = ma.masked_outside(a1rtype, thrtype0, thrtype1).mask

        for i,hrange in enumerate(lhrange):
            hmin,hmax = hrange
            a1idx = arange(len(a1stormH)).astype(int32)
            a1idx = ma.masked_where(a1stormH <hmin, a1idx)
            a1idx = ma.masked_where(hmax <= a1stormH, a1idx)

            a1idx = ma.masked_where(amskT, a1idx)
            a1idx = a1idx.compressed()
            a2profTmp = a2prof[a1idx]
            a2profTmp = ma.masked_less(a2profTmp,0)
    
            da1ave[rtype,i] = a2profTmp.mean(axis=0)*0.01
            da1std[rtype,i] = a2profTmp.std(axis=0)*0.01
            da1num[rtype,i] = a2profTmp.count(axis=0)


    #-- draw -----------
    lcolor = ['b','k','r']
    llinewidth=[2,1,2]


    #-- draw precip ----
    fig    = plt.figure(figsize=(wfig,hfig))
    for irtype, rtype in enumerate(lrtype):
        axloc = 100 + 10*len(lrtype) + irtype + 1
        ax = fig.add_subplot( axloc )

        stitle = rtype
        plt.title(stitle)

        for i,hrange in enumerate(lhrange):
            hmin,hmax = hrange
            slabel = '%d ~ %d'%(hmin, hmax)
            a1y = (arange(nh)+1) * 0.25
            a1x = da1ave[rtype, i]
            a1x = ma.masked_where(da1num[rtype,i]<thcount, a1x)
            a1x = ma.masked_where(a1y<hcutoff, a1x)

            ax.plot(a1x, a1y, '-',linewidth=llinewidth[i], color=lcolor[i], label=slabel)

    plt.legend(loc='upper right')
    ssuptitle = 'rain profile %s (stormH)'%(domain)
    plt.suptitle(ssuptitle)

    figDir  = '/work/a01/utsumi/GPMGV/fig/env-prof'
    figPath = figDir + '/prof.PR.stormH.prcp.%s.cbin.%d.png'%(domain,cbins)
    util.mk_dir(figDir)
    plt.savefig(figPath)
    plt.clf()
    print figPath

    #-- draw count ----
    fig    = plt.figure(figsize=(wfig,hfig))
    for irtype, rtype in enumerate(lrtype):
        axloc = 100 + 10*len(lrtype) + irtype + 1
        ax = fig.add_subplot( axloc )

        stitle = rtype
        plt.title(stitle)

        for i,hrange in enumerate(lhrange):
            hmin,hmax = hrange
            slabel = '%d ~ %d%%'%(hmin, hmax)
            a1y = (arange(nh)+1) * 0.25
            a1x = da1num[rtype, i]
            #a1x = ma.masked_where(da1num[rtype,i]<50, a1x)
            a1x = ma.masked_where(a1y<hcutoff, a1x)
            ax.plot(a1x, a1y, '-',linewidth=llinewidth[i], color=lcolor[i], label=slabel)

    plt.legend(loc='upper right')
    ssuptitle = 'count %s (stormH)'%(domain)
    plt.suptitle(ssuptitle)

    figDir  = '/work/a01/utsumi/GPMGV/fig/env-prof'
    figPath = figDir + '/prof.PR.stormH.count.%s.cbin.%d.png'%(domain,cbins)
    util.mk_dir(figDir)
    plt.savefig(figPath)
    plt.clf()
    print figPath
    #'''


