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

#ldomain = ['FLORIDA-SFL-N']
ldomain = ['FLORIDA-KSC', 'KWAJALEIN-KWA', 'TEXAS-HAR', 'FLORIDA-SFL-N', 'N.Carolina-IPHEx_Duke', 'FLORIDA-STJ', 'N.Carolina-IPHEx_NASA']
#ldomain = ['KWAJALEIN-KWA']

lrtype    = ['all','strat','conv']
dlthrtype = {'all':[-999,999],'strat':[100,170],'conv':[200,297]}

cbins   = 11 # number of angle bins
baseDir = '/home/utsumi/mnt/wellshare/GPMGV/DOM.L2A25/cbin.%d'%(cbins)

wfig,hfig= 7,9
miss     = -9999.
thnum    = 50
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


    #***************************************************
    # RH and stormH
    #---------------------------------------------------
    llRH     = [[0.0,0.6],[0.6,0.7],[0.7,0.8],[0.8,0.9],[0.9,1.0]]
    llstormH = [[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,20]]
    lhup     = [2,3,4]
    dmean = {}
    dstd  = {}
    drpos = {}
    dnum  = {}
    for rtype in lrtype:
        thrtype0,thrtype1 = dlthrtype[rtype]
        a1mskType = ma.masked_outside(a1rtype, thrtype0, thrtype1).mask
        for iRH,lRH in enumerate(llRH):
            RHmin,RHmax = lRH
            a1mskRH = ma.masked_outside(a1rh, RHmin, RHmax).mask

            for istormH,lstormH in enumerate(llstormH):
                stormHmin,stormHmax = lstormH
                a1mskstormH = ma.masked_outside(a1stormH, stormHmin, stormHmax).mask

                a1msk = a1mskType + a1mskRH + a1mskstormH
                a1idx = arange(a2prof.shape[0])

                a1idx = ma.masked_where(a1msk, a1idx).compressed()

                a2profTmp = a2prof[a1idx]
                a1stormHTmp=a1stormH[a1idx]

                a1rainlw = a2profTmp[:,4]  # rain rate @1km
                a1rainlw = ma.masked_less(a1rainlw,0)

                for hup in lhup:
                    iup  = int(hup/0.25)
                    a1rainup = a2profTmp[:,iup]
                    a1rainup = ma.masked_less(a1rainup,0)
                    a1grad   = (a1rainlw - a1rainup)*0.01 / (hup-1)

                    grad_mean= a1grad.mean()
                    grad_std = a1grad.std()
                    grad_rpos= ma.masked_less(a1grad,0).count() / float( ma.masked_equal(a1grad,0).count())
                    grad_num = ma.masked_equal(a1grad,miss).count()

                    dmean[rtype,iRH,istormH,hup] = grad_mean
                    dstd [rtype,iRH,istormH,hup] = grad_std
                    drpos[rtype,iRH,istormH,hup] = grad_rpos
                    dnum [rtype,iRH,istormH,hup] = grad_num


                    #if (rtype=='all')&(hup==2)&(istormH==6)&(iRH==1):
                    #    fig = plt.figure()
                    #    ax  = fig.add_subplot(111)
                    #    a1tmp= ma.masked_equal(a1grad,miss).compressed()
                    #    ax.hist(a1tmp,bins=arange(-20,20))
                    #    figDir  = '/work/a01/utsumi/GPMGV/fig/env-prof'
                    #    figPath = figDir + '/hist.png'
                    #    plt.savefig(figPath)
                    #    plt.clf()
                    #    print '-------- hist --------------'
                    #    print 'mean',a1tmp.mean()
                    #    print 'len=',a1tmp.shape
                    #    print 'pos=',ma.masked_less_equal(a1tmp,0).count()
                    #    print 'neg=',ma.masked_greater_equal(a1tmp,0).count()
                    #    print 'zero',ma.masked_not_equal(a1tmp,0).count()


    #-- Figure ------------------
    for dattype in ['mean','std','positive','num']:
        if   dattype =='mean':
            ddat = dmean
            vmin,vmax = -0.8, 0.8
            mycm= 'RdBu_r'
        elif dattype == 'std':
            ddat = dstd
            vmin,vmax = 0, 1
            mycm= 'jet'
        elif dattype == 'positive':
            ddat = drpos
            vmin,vmax = 0,1.0
            mycm= 'RdBu_r'
        elif dattype == 'num':
            ddat = dnum
            vmin,vmax = 0,1000
            mycm= 'jet'

 
        fig =  plt.figure(figsize=(wfig,hfig))
        iloc= 0

        for ihup,hup in enumerate(lhup[::-1]):
            for irtype,rtype in enumerate(lrtype):
                iloc  = iloc + 1
                axloc = 100*len(lhup) + 10*len(lrtype) + iloc
                ax = fig.add_subplot(axloc)
        
                a2dat = ones([len(llstormH),len(llRH)])*miss
                a2num = ones(a2dat.shape)*miss
                a2stormH = ones(a2dat.shape)*miss

                for iRH,lRH in enumerate(llRH):
                    for istormH,lstormH in enumerate(llstormH):
                        a2dat[istormH,iRH] = ddat[rtype,iRH,istormH,hup]
                        a2num[istormH,iRH] = dnum[rtype,iRH,istormH,hup]
                        a2stormH[istormH,iRH] = lstormH[1]

                
                a2msk = ma.masked_where(a2stormH>hup, ones(a2dat.shape))
                #a2msk = ma.masked_where(a2num>=thnum,  a2msk)
                a2dat = ma.masked_equal(a2dat,miss)
                a2dat = ma.masked_where(a2msk.mask==False, a2dat)

                a1x = arange(len(llRH)+1)
                a1y = arange(len(llstormH)+1)
                X,Y = meshgrid(a1x,a1y)
                #im  = ax.pcolormesh(X,X,a2msk, cmap='Set1_r')
                #im  = ax.pcolormesh(X,Y,a2dat,vmin=vmin,vmax=vmax,cmap=mycm)
                ax.imshow(a2msk, cmap='Set1_r',origin='lower')
                im  = ax.imshow(a2dat, vmin=vmin, vmax=vmax, cmap=mycm, origin='lower')


                plt.colorbar(im,orientation='vertical',fraction=0.056,pad=0.04)
                stitle = 'Hup=%d %s'%(hup,rtype)
                plt.title(stitle)

                #-- x ticks ------
                lRHbnd = list(zip(*llRH)[0]) + [llRH[-1][1]]
                lxBnd     = arange(len(llRH)+1) -0.5
                #plt.xticks([-0.5,0.5,1.5,2.5,3.5,4.5], lRHbnd)
                plt.xticks(lxBnd, lRHbnd)

                #-- y ticks ------
                lstormHbnd = list(zip(*llstormH)[0]) + [llstormH[-1][1]]
                lyBnd      = arange(len(llstormH)+1) -0.5
                plt.yticks(lyBnd, lstormHbnd)


        plt.tight_layout()
        plt.suptitle('%s %s'%(domain,dattype))
        figDir  = '/work/a01/utsumi/GPMGV/fig/env-prof'
        util.mk_dir(figDir)
        figPath = figDir + '/slope.%s.%s.png'%(dattype,domain)

        plt.savefig(figPath)
        print figPath
        plt.clf()
