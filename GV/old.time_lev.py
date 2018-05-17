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

iYM = [2014,1]
eYM = [2014,2]
lYM = util.ret_lYM(iYM, eYM)

gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()

dgName = gv.ret_ddomYM2gName()

ldomain = gv.domains
#ldomain = ['FLORIDA-STJ','FLORIDA-SFL-N','N.Carolina-IPHEx_Duke','N.Carolina-IPHEx_NASA','KWAJALEIN-KWA']
ldomain = ['FLORIDA-STJ']

offset_bef = -5
offset_aft = 30


nt = offset_aft - offset_bef +1
nh = 20

lprtype = ['heavy','extreme','mod']
dlthpr = {'mod':[-0.1,10], 'heavy':[10,50],'extreme':[50,9999]}

for domain in ldomain:
    for YM in lYM:
        Year, Mon = YM
        if (domain,Year,Mon) not in dgName.keys():
            print 'no obs',domain,Year,Mon
            continue

        # load data
        srcbaseDir = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.L2A25'
        srcDir     = srcbaseDir + '/%s/%04d%02d'%(domain, Year,Mon)
    
        profPath  = srcDir  + '/prof.npy'
        eSurfPath = srcDir  + '/eSurf.npy'
        nSurfPath = srcDir  + '/nSurf.npy'
        gvAllPath = srcDir  + '/gvAll.npy'
        gvGdPath  = srcDir  + '/gvGd.npy'
        nSurfBinPath = srcDir  + '/nSurfBin.npy'
    
        a2prof   = np.load(profPath)
        aesurf   = np.load(eSurfPath)
        ansurf   = np.load(nSurfPath)
        #a2gvall  = np.load(gvAllPath)
        a2gvgd   = np.load(gvGdPath)
        ansurfbin= np.load(nSurfBinPath)

        # container empty array

        da2sx  = {prtype:empty([nh,nt]).astype(float32) for prtype in lprtype}
        da2sy  = {prtype:empty([nh,nt]).astype(float32) for prtype in lprtype}
        da2syy = {prtype:empty([nh,nt]).astype(float32) for prtype in lprtype}
        da2sxx = {prtype:empty([nh,nt]).astype(float32) for prtype in lprtype}
        da2sxy = {prtype:empty([nh,nt]).astype(float32) for prtype in lprtype}
        da2n   = {prtype:empty([nh,nt]).astype(int32  ) for prtype in lprtype}

 
        for ih in range(nh):
            asate  = a2prof[:,ih]

            for idt, dt in enumerate(range(offset_bef, offset_aft+1)):

                agvgd = a2gvgd[:,idt]

                #-- mask when both satellite and gv are zero
                amsk1    = ma.masked_equal(asate,0).mask
                amsk2    = ma.masked_equal(agvgd,0).mask
                amskzero = amsk1 * amsk2

                #-- mask when at least one of the data is missing
                amsk3    = ma.masked_less(asate,0).mask
                amsk4    = ma.masked_less(agvgd,0).mask
                amskmiss = amsk3 + amsk4

                #-- overlay masks
                amsk     = amskzero + amskmiss

                aidxTmp  = arange(len(amsk)).astype(int32)
                aidxTmp  = ma.masked_where(amsk, aidxTmp).compressed()
        
                asateTmp     = asate[aidxTmp]*0.01 # mm/h 
                aesurfTmp    = aesurf[aidxTmp]
                ansurfTmp    = ansurf[aidxTmp]*0.01 # mm/h
                agvgdTmp     = agvgd[aidxTmp]
                ansurfbinTmp = ansurfbin[aidxTmp]
                #------------------
                for prtype in lprtype:
                    # mask with nSurf --
                    thmin, thmax = dlthpr[prtype]
                    amsk      = ma.masked_outside(ansurfTmp, thmin, thmax).mask
                    y  = ma.masked_where(amsk, asateTmp).compressed()
                    x  = ma.masked_where(amsk, agvgdTmp).compressed()

                    sx = x.sum()
                    sy = y.sum()
                    sxx= (x*x).sum()
                    syy= (y*y).sum()
                    sxy= (x*y).sum()
                    n  = len(x) 
                    
                    da2sx [prtype][ih,idt] = sx
                    da2sy [prtype][ih,idt] = sy 
                    da2sxx[prtype][ih,idt] = sxx
                    da2syy[prtype][ih,idt] = syy 
                    da2sxy[prtype][ih,idt] = sxy 
                    da2n  [prtype][ih,idt] = n



        #--- save data ----
        for prtype in lprtype:
            outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-PR/%s'%(domain)
            util.mk_dir(outDir)

            satePath    = outDir + '/ssate.%04d%02d.%s.npy'%(Year,Mon,prtype)
            sate2Path   = outDir + '/ssate2.%04d%02d.%s.npy'%(Year,Mon,prtype)
            gvPath      = outDir + '/sgv.%04d%02d.%s.npy'%(Year,Mon,prtype)
            gv2Path     = outDir + '/sgv2.%04d%02d.%s.npy'%(Year,Mon,prtype)
            sategvPath  = outDir + '/sategv.%04d%02d.%s.npy'%(Year,Mon,prtype)
            numPath     = outDir + '/num.%04d%02d.%s.npy'%(Year,Mon,prtype)

            np.save(gvPath    ,da2sx[prtype])
            np.save(gv2Path   ,da2sxx[prtype])
            np.save(satePath  ,da2sy[prtype])
            np.save(sate2Path ,da2syy[prtype])
            np.save(sategvPath,da2sxy[prtype])
            np.save(numPath   ,da2n[prtype])

            print gvPath



'''
#--- figure -------
for prtype in lprtype:
    outDir = '/home/utsumi/mnt/wellshare/GPMGV/out'
    ccPath= outDir + '/dt-lev.cc.%s.npy'%(prtype)
    a2cc= np.load(ccPath)

    fig = plt.figure(figsize=(4,4))
    nh  = 20
    at  = range(offset_bef, offset_aft+1)
    for ih in range(nh):
        ax = fig.add_axes([0.1, 0.1 + 0.8/nh*ih, 0.8, 0.8/nh])
        ay = a2cc[ih,:]
        if np.isnan(ay.sum()):
            ay = zeros(len(at))
    
        ax.plot(at, ay, '-')
        plt.ylim([0,1.2])
        if ih !=0:
            ax.tick_params(labelbottom='off') 
    
    figDir  = '/work/a01/utsumi/GPMGV/fig'
    figPath = figDir + '/plot.dt-lev.%s.png'%(prtype)
    plt.savefig(figPath)
    print figPath
    plt.clf()
'''


