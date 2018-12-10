import matplotlib
matplotlib.use('Agg')
from numpy import *
import numpy as np
import myfunc.util as util
from collections import deque
import sys, os
import matplotlib.pyplot as plt
from   datetime import datetime, timedelta


iYM = [1998,4]
eYM = [2004,10]
#eYM = [2004,10]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
#lYM = lYM[::-1]
print lYM

#figflag  = False # True / False
figflag  = True # True / False
ldomain = ['FLORIDA-KSC', 'KWAJALEIN-KWA', 'TEXAS-HAR', 'FLORIDA-SFL-N', 'N.Carolina-IPHEx_Duke', 'FLORIDA-STJ', 'N.Carolina-IPHEx_NASA']
#ldomain = ['FLORIDA-KSC']

lrtype    = ['alltype','strat','conv']
dlthrtype = {'alltype':[-999,999],'strat':[100,170],'conv':[200,297]}

cbins   = 11 # number of angle bins

#baseDir = '/home/utsumi/mnt/wellshare/GPMGV/DOM.L2A25/cbin.%d'%(cbins)
baseDir = '/media/disk2/share/GPMGV/DOM.L2A25/cbin.%d'%(cbins)

#wfig,hfig= 7,12
wfig,hfig= 7,10
miss     = -9999.
thnum    = 50

a1esurf = array([])
a1rh    = array([])
a1div850= array([])
a1div10m= array([])
a1stormH= array([])
a1freezH= array([])
a1elev  = array([])
a2prof  = array([])
a1rtype = array([])


for domain in ldomain:
    for YM in lYM:
        Year,Mon = YM
        srcDir = baseDir + '/%s/%04d%02d'%(domain, Year,Mon)
        rhPath = srcDir  + '/rh.npy'
        div850Path= srcDir  + '/div850.npy'
        div10mPath= srcDir  + '/div10m.npy'
        stormHPath= srcDir  + '/stormH.npy'
        freezHPath= srcDir  + '/freezH.npy'

        elevPath  = srcDir  + '/elev.npy'

        esurfPath=srcDir + '/eSurf.npy'
        profPath =srcDir + '/prof.npy'
        rtypePath=srcDir + '/rainType.npy'

        a1esurfTmp  = np.load(esurfPath)
        a1rhTmp     = np.load(rhPath)
        a1div850Tmp = np.load(div850Path)
        a1div10mTmp = np.load(div10mPath)
        a1stormHTmp = np.load(stormHPath)
        #a1freezHtmp = np.load(freezHPath)
        a1elevTmp   = np.load(elevPath)

        a2profTmp   = np.load(profPath)
        a1rtypeTmp  = np.load(rtypePath)


        a1max       = a2profTmp.max(axis=1)
        a1idx       = arange(len(a1esurfTmp)).astype(int32)
        a1idx       = ma.masked_where(a1elevTmp<0, a1idx)    # screen ocean
        a1idx       = ma.masked_where(a1max<=0, a1idx).compressed()

        a1esurfTmp  = a1esurfTmp[a1idx]
        a1rhTmp     = a1rhTmp [a1idx]
        a1div850Tmp = a1div850Tmp[a1idx]
        a1div10mTmp = a1div10mTmp[a1idx]
        a1stormHTmp = a1stormHTmp[a1idx]
        #a1freezHTmp = a1freezHTmp[a1idx]
        a1elevTmp   = a1elevTmp  [a1idx]

        a2profTmp   = a2profTmp[a1idx,:]
        a1rtypeTmp  = a1rtypeTmp[a1idx]


        #-- convert stormH and freezH: relative to the ground surface --
        a1stormHTmp = a1stormHTmp - a1elevTmp
        #a1freezHTmp = a1freezHTmp - a1elevTmp

        #------------------
        a1esurf     = concatenate([a1esurf, a1esurfTmp ])
        a1rh        = concatenate([a1rh,    a1rhTmp    ])
        a1div850    = concatenate([a1div850,a1div850Tmp])
        a1div10m    = concatenate([a1div10m,a1div10mTmp])
        a1stormH    = concatenate([a1stormH,a1stormHTmp])
        #a1freezH    = concatenate([a1freezH,a1freezHTmp])
        a1elev      = concatenate([a1elev,  a1elevTmp  ])
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
#a1freezH = ma.masked_equal(a1freezH,-9999.)*0.001  # m->km


#***************************************************
# RH and stormH
#---------------------------------------------------
llRH     = [[0.0,0.6],[0.6,0.7],[0.7,0.8],[0.8,0.9],[0.9,1.0]]
#llstormH = [[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,20]]
llstormH = [[2,3],[3,4],[4,6],[6,8],[8,15]]
lhup     = arange(1,10+0.1,0.25)
llstormHFig = [[2,3],[3,4],[4,6],[6,8],[8,15]]
#lhupFig  = [2,3,4,5,6,7,8]
lhupFig  = [2,3,4,5]

dmean = {}
dstd  = {}
drpos = {}
dnum  = {}
drat  = {}
dnormprof = {}
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

            a1rainlw = a2profTmp[:,3]  # rain rate @1km
            a1rainlw = ma.masked_less(a1rainlw,0)

            for hup in lhup:
                iup  = int(hup/0.25) -1
                a1rainup = a2profTmp[:,iup]
                a1rainup = ma.masked_less(a1rainup,0)
                a1grad   = (a1rainlw - a1rainup)*0.01 / (hup-1)
                grad_mean= a1grad.mean()
                grad_std = a1grad.std()
                grad_rpos= ma.masked_less(a1grad,0).count() / float( ma.masked_equal(a1grad,0).count())
                grad_num = ma.masked_equal(a1grad,miss).count()
                grad_rat = a1rainlw.mean() / a1rainup.mean()
                grad_normprof = a1rainup.mean() / a1rainlw.mean()

                dmean[rtype,iRH,istormH,hup] = grad_mean
                dstd [rtype,iRH,istormH,hup] = grad_std
                drpos[rtype,iRH,istormH,hup] = grad_rpos
                dnum [rtype,iRH,istormH,hup] = grad_num
                drat [rtype,iRH,istormH,hup] = grad_rat
                dnormprof[rtype,iRH,istormH,hup] = grad_normprof


#-- Figure ------------------
for dattype in ['normprof','rat','mean','std','positive','num']:
    if figflag==False: continue

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
    elif dattype == 'rat':
        ddat = drat
        vmin,vmax = 0,2
        mycm= 'RdBu_r'
    elif dattype == 'normprof':
        ddat = dnormprof
        vmin,vmax = 0,2
        mycm= 'RdBu_r'


    fig =  plt.figure(figsize=(wfig,hfig))
    iloc= 0

    da2rat = {}

    for istormH,lstormH in enumerate(llstormHFig):
        for irtype,rtype in enumerate(lrtype):
            nrow  = len(llstormHFig)
            ncol  = len(lrtype)
            iloc  = iloc + 1
            ax = fig.add_subplot(nrow,ncol,iloc)
    
            a2dat = ones([len(lhupFig),len(llRH)])*miss
            a2num = ones(a2dat.shape)*miss
            a2hup = ones(a2dat.shape)*miss

            for ihup,hup in enumerate(lhupFig):
                for iRH,lRH in enumerate(llRH):
                    a2dat[ihup,iRH] = ddat[rtype,iRH,istormH,hup]
                    a2num[ihup,iRH] = dnum[rtype,iRH,istormH,hup]
                    a2hup[ihup,iRH] = hup           
 
            #a2msk = ma.masked_where(a2hup<=lstormH[0], ones(a2dat.shape))
            a2msk = ma.masked_where(a2hup<lstormH[1], ones(a2dat.shape))
            a2dat = ma.masked_equal(a2dat,miss)
            a2dat = ma.masked_where(a2msk.mask==False, a2dat)

            if dattype=='rat':
                da2rat[istormH] = a2dat.filled(miss)

            a1x = arange(len(llRH)+1)
            a1y = arange(len(llstormHFig)+1)
            X,Y = meshgrid(a1x,a1y)
            #im  = ax.pcolormesh(X,X,a2msk, cmap='Set1_r')
            #im  = ax.pcolormesh(X,Y,a2dat,vmin=vmin,vmax=vmax,cmap=mycm)
            ax.imshow(a2msk, cmap='Set1_r',origin='lower',interpolation='nearest')
            im  = ax.imshow(a2dat, vmin=vmin, vmax=vmax, cmap=mycm, origin='lower',interpolation='nearest')


            plt.colorbar(im,orientation='vertical',fraction=0.056,pad=0.1)
            stitle = 'H=%d-%d %s'%(lstormH[0],lstormH[1],rtype)
            plt.title(stitle)

            #-- x ticks ------
            lRHbnd = list(zip(*llRH)[0]) + [llRH[-1][1]]
            lxBnd     = arange(len(llRH)+1) -0.5
            #plt.xticks([-0.5,0.5,1.5,2.5,3.5,4.5], lRHbnd)
            plt.xticks(lxBnd, lRHbnd)

            #-- y ticks ------
            lylabel = lhupFig
            lyBnd   = arange(len(lhupFig))
            plt.yticks(lyBnd, lylabel)

    plt.suptitle('all domain %s'%(dattype))
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    #figDir  = '/work/a01/utsumi/GPMGV/fig/env-prof'
    figDir  = baseDir + '/ratio'
    util.mk_dir(figDir)
    figPath = figDir + '/slope.%s.%s.png'%(dattype,'ALL')

    plt.savefig(figPath)
    print figPath
    plt.clf()

#--- save ratio --------------------
dattype=='rat'
ddat   = drat
for istormH,lstormH in enumerate(llstormH):
    stormH0,stormH1 = lstormH
    for irtype,rtype in enumerate(lrtype):
        
        a2dat = ones([len(lhup)+int(lhup[0]/0.25),len(llRH)])*miss # height=0, 0.25, 0.5, 0.75, ..., 10  # 41 in total
        a2hup = ones(a2dat.shape)*miss

        for hup in lhup:
            iy = int(hup/0.25)
            for iRH,lRH in enumerate(llRH):
                a2dat[iy,iRH] = ddat[rtype,iRH,istormH,hup]
                a2hup[iy,iRH] = hup           
 
        a2dat = ma.masked_where(a2hup>=lstormH[1], a2dat).filled(miss)

        #-- cut-off at stormH or 8km-km (take lower one) ------
        if stormH1 <= 8:
            iytop = int(stormH1/0.25)-1
        else:
            iytop = int(hup/0.25)

        for iRH, lRH in enumerate(llRH):
            a2dat[iytop+1:,iRH] = a2dat[iytop,iRH] 

        a2dat  = ma.masked_equal(a2dat,miss).filled(1.0)
        outDir = baseDir + '/ratio'
        outPath= outDir + '/ratio.%s.H%d-%d.npy'%(rtype,stormH0,stormH1)
        np.save(outPath, a2dat)
        print outPath



