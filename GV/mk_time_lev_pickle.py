import matplotlib
matplotlib.use('Agg')

import pickle
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

lh = range(nh)
ldt= range(offset_bef, offset_aft+1)

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

        for idt, dt in enumerate(ldt):
            agvgd = a2gvgd[:,idt]

            dsate  = {}
            desurf = {}
            dnsurf = {}
            dgv    = {}
            dnsurfbin= {}

            for ih in range(nh):
                asate  = a2prof[:,ih]
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
                #-----------------
                dsate[ih]    = asateTmp
                desurf[ih]   = aesurfTmp
                dnsurf[ih]   = ansurfTmp
                dgv[ih]      = agvgdTmp
                dnsurfbin[ih]= ansurfbinTmp
                

            #--- save data ----
            for prtype in lprtype:
                outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-PR/%s/%04d%02d'%(domain,Year,Mon)
                util.mk_dir(outDir)

                satePath    = outDir + '/sate.dt%02d.pickle'%(dt)
                esurfPath   = outDir + '/esurf.dt%02d.pickle'%(dt)
                nsurfPath   = outDir + '/nsurf.dt%02d.pickle'%(dt)
                gvPath      = outDir + '/gv.dt%02d.pickle'%(dt)
                nsurfbinPath= outDir + '/nsurfBin.dt%02d.pickle'%(dt)

                f = open(satePath,    'wb'); pickle.dump(dsate,  f); f.close() 
                f = open(esurfPath,   'wb'); pickle.dump(desurf, f); f.close() 
                f = open(nsurfPath,   'wb'); pickle.dump(dnsurf, f); f.close() 
                f = open(gvPath,      'wb'); pickle.dump(dgv,    f); f.close() 
                f = open(nsurfbinPath,'wb'); pickle.dump(dnsurfbin, f); f.close() 

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


