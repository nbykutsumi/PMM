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

nh  = 20
nt  = offset_aft - offset_bef +1


lprtype = ['heavy','extreme','mod']

def ret_cc(sx, sxx, sy, syy, sxy, n):
    mx = ma.masked_invalid(sx / n)
    my = ma.masked_invalid(sy / n)

    SXY= sxy - my*sx - mx*sy + n*mx*my

    SX = sxx - mx*sx + n*mx*mx
    SY = syy - my*sy + n+my+my

    return SXY/np.square(SX*SY)


for prtype in lprtype:
    asx = zeros([nh,nt])
    asy = zeros([nh,nt])
    asxx= zeros([nh,nt])
    asyy= zeros([nh,nt])
    asxy= zeros([nh,nt])
    anum= zeros([nh,nt])

    for domain in ldomain:
        for YM in lYM:
            Year, Mon = YM
            if (domain,Year,Mon) not in dgName.keys():
                print 'no obs',domain,Year,Mon
                continue

            datDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-PR/%s'%(domain)
            util.mk_dir(datDir)

            satePath    = datDir + '/ssate.%04d%02d.%s.npy'%(Year,Mon,prtype)
            sate2Path   = datDir + '/ssate2.%04d%02d.%s.npy'%(Year,Mon,prtype)
            gvPath      = datDir + '/sgv.%04d%02d.%s.npy'%(Year,Mon,prtype)
            gv2Path     = datDir + '/sgv2.%04d%02d.%s.npy'%(Year,Mon,prtype)
            sategvPath  = datDir + '/sategv.%04d%02d.%s.npy'%(Year,Mon,prtype)
            numPath     = datDir + '/num.%04d%02d.%s.npy'%(Year,Mon,prtype)

            asxTmp  = np.load(gvPath    )
            asxxTmp = np.load(gv2Path   )
            asyTmp  = np.load(satePath  )
            asyyTmp = np.load(sate2Path )
            asxyTmp = np.load(sategvPath)
            anumTmp = np.load(numPath   )

            asx  = asx  + asxTmp           
            asxx = asxx + asxxTmp
            asy  = asy  + asyTmp 
            asyy = asyy + asyyTmp
            asxy = asxy + asxyTmp
            anum = anum + anumTmp

    a2cc = ret_cc(asx, asxx, asy, asyy, asxy, anum)
    #--- figure -------
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
    figPath = figDir + '/tmp.plot.dt-lev.%s.png'%(prtype)
    plt.savefig(figPath)
    print figPath
    plt.clf()

