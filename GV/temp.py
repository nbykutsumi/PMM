import matplotlib
matplotlib.use('Agg')
from numpy import *
import numpy as np
import myfunc.util as util
from collections import deque
import sys, os
import matplotlib.pyplot as plt
from   datetime import datetime, timedelta

iYM = [2001,4]
eYM = [2002,8]

lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
#lYM = lYM[::-1]
print lYM


#ldomain = ['FLORIDA-KSC', 'KWAJALEIN-KWA', 'TEXAS-HAR', 'FLORIDA-SFL-N', 'N.Carolina-IPHEx_Duke', 'FLORIDA-STJ', 'N.Carolina-IPHEx_NASA']
ldomain = ['FLORIDA-KSC', 'KWAJALEIN-KWA', 'TEXAS-HAR', 'FLORIDA-SFL-N']
#ldomain = ['FLORIDA-KSC']

a1esurf = array([])
a1rh    = array([])
a1div   = array([])


baseDir = '/home/utsumi/mnt/wellshare/GPMGV/GLOC.L2A25/5.0km'
for domain in ldomain:

    for YM in lYM:
        Year,Mon = YM
        srcDir = baseDir + '/%s/%04d%02d'%(domain, Year,Mon)
        rhPath = srcDir  + '/rh.npy'
        divPath= srcDir  + '/div850.npy'
        esurfPath=srcDir + '/eSurf.npy'
    
        a1esurfTmp  = np.load(esurfPath)
        a1rhTmp     = np.load(rhPath)
        a1divTmp    = np.load(divPath)
   
        a1esurf     = concatenate([a1esurf, a1esurfTmp])
        a1rh        = concatenate([a1rh,    a1rhTmp   ])
        a1div       = concatenate([a1div,   a1divTmp  ])

print a1esurf.min()
print a1rh.min()
print a1div.min()
#--
outDir = '/work/a01/utsumi/temp'

a1esurf = ma.masked_less(a1esurf,0)
a1rh    = ma.masked_less(a1rh, 0)
a1div   = a1div*1e5
#a1div   = ma.masked_less(a1div,0)

figPath= outDir + '/rh-esurf.png'    
plt.plot(a1rh, a1esurf,'o')
plt.savefig(figPath)
plt.clf()

figPath= outDir + '/div-esurf.png'    
plt.plot(a1div, a1esurf,'o')
plt.savefig(figPath)
plt.clf()

print figPath
print len(a1rh)
 
    
    
