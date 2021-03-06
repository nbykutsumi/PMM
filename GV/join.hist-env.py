from numpy import *
import os, sys
import Image
import myfunc.util as util

ldomain = ['FLORIDA-KSC', 'KWAJALEIN-KWA', 'TEXAS-HAR', 'FLORIDA-SFL-N', 'N.Carolina-IPHEx_Duke', 'FLORIDA-STJ', 'N.Carolina-IPHEx_NASA']
iy = 0   # top
ey = -1
ix = 10
ex = -30
thdist = 10 # km

#for env in ['RH','div850','div10m']:
for env in ['div850']:
    da2dat = {}
    for i,domain in enumerate(ldomain):
        figDir  = '/work/a01/utsumi/GPMGV/fig/env-prof'
        figPath = figDir + '/hist.PR.%s.%s.%.1fkm.png'%(env, domain,thdist)
        iimg    = Image.open(figPath)
        a2array = asarray(iimg)[iy:ey, ix:ex]

        da2dat[i] = a2array

    da2dat[-999] = da2dat[0]*0 + 255

    a2line0 = concatenate([da2dat[0],da2dat[1],da2dat[2],da2dat[3]],axis=1)
    a2line1 = concatenate([da2dat[4],da2dat[5],da2dat[6],da2dat[-999]],axis=1)
    a2oarray= concatenate([a2line0, a2line1], axis=0)

    oimg    = Image.fromarray(a2oarray)
    outDir  = figDir
    oPath   = outDir + '/join.hist.PR.%s.%.1fkm.png'%(env, thdist)

    oimg.save(oPath)
    print oPath


