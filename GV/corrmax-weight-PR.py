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
from gv_fsub import *

calc = True
#calc = False
iYM = [2005,4]
#iYM = [2014,4]
eYM = [2014,10]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
print lYM

#corrFlag= 'CR'
corrFlag= 'NO'

#nozero = 'withzero'
nozero = 'nozero'

dt  = 30

cls = 'RH'
#cls = 'RainType'

#thdist = 2.5
thdist = 5.0
minNum = 3
prdName = 'L2A25'
miss   = -9999.

gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()

dgName = gv.ret_ddomYM2gName()

ldomain = gv.domains
#ldomain = ['FLORIDA-STJ','FLORIDA-SFL-N','N.Carolina-IPHEx_Duke','N.Carolina-IPHEx_NASA','KWAJALEIN-KWA']
#ldomain = ['FLORIDA-STJ']

offset_bef = 15  # 'bef' should be identical to that in mk_match.py
#offset_aft = 30

if cls=='RH':
    #lclstype = ['all','dry','hum']
    #lclstype = ['dry','hum']
    lclstype = ['all']
elif cls=='RainType':
    lclstype = ['alltype','strat','conv']
else:
    print 'check cls',cls
    sys.exit()


ldattype = ['rain','cc','bias','brat','rmse','gv','num']


#----------------------------------
def findMaxCorr(a2m, a1y):
    # shape of a2m = (nobs,nh)
    # shape of a1y = (nobs)

    #----------------------
    # A = inv(Cm) (Cm,x) inv(Cx) (Cm,x).T  # npc x npc matrix
    #----------------------
    nobs = (a2m.shape)[0]
    npc  = (a2m.shape)[1]
    if nobs <=1:
        print "nobs=1, exit"
        return None,None,None,None

    Cmx = array([ sum( (a2m[:,i] - mean(a2m[:,i]))*(a1y-mean(a1y)))/(nobs-1) for i in range(npc)]).reshape(npc,-1)   # (npc x 1)

    Cm    = np.cov(a2m, rowvar=0, bias=0)     # (npc x npc)
    Cx  = np.var(a1y, ddof=1)   # scalar


    #-- check if Cm is invertable --
    if not linalg.cond(Cm) < 1/sys.float_info.epsilon:
        print "Cm is not invertable: skip"
        return None,None,None,None
    #-------------------------------

    A   = np.dot(linalg.inv(Cm), Cmx)/Cx
    A   = np.dot(A, Cmx.T)
    leval, levect= linalg.eig(A)

    corrmax = 0
    for i in range(npc):
        a1ktmp  = levect[:,i]
        if np.iscomplex(a1ktmp).sum()>0: continue

        a1kmtmp = np.dot(a1ktmp.real, a2m.T)
        corr= np.corrcoef( a1kmtmp, a1y)[0,1]
        if abs(corrmax) < abs(corr):
            corrmax = corr
            a1k     = a1ktmp
            a1km    = a1kmtmp
            ikmax   = argmax(abs(a1k))
    #print "%.2f"%(corrmax),nobs
    return a1k, a1km, ikmax, corrmax

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
a1esurf= deque([])
a1prof = deque([])
a1gv   = deque([])
a1gv0  = deque([])
a2joinprofAll= deque([])
a2profaveAll = deque([])

for clstype in lclstype:
    csvDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-%s/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])

    joinprofPath = csvDir + '/%s.%s.joinprof.%s.csv'%(cls, nozero, clstype)
    gvPath       = csvDir + '/%s.%s.gv.%s.csv'%(cls, nozero, clstype)
    profavePath  = csvDir + '/%s.%s.profave.%s.csv'%(cls, nozero, clstype)
    gvavePath    = csvDir + '/%s.%s.gvave.%s.csv'%(cls, nozero, clstype)

    a2joinprof = load_csv(joinprofPath)
    a2profave  = load_csv(profavePath)
    a2gvave    = load_csv(gvavePath)
    a1gvave    = a2gvave[:,dt-1]

    #a2joinprof = ma.masked_less(a2joinprof, 0).filled(miss)
    #a2joinprof = gv_fsub.fill_esurf_interp(a2joinprof.T, miss).T
    a2joinprof = gv_fsub.fill_esurf(a2joinprof.T, miss).T
    #a2joinprof = ma.masked_less(a2joinprof, 0).filled(0)
    #----------------------------

    a1k, a1km, ikmax, corrmax = findMaxCorr(a2joinprof, a1gvave)
    a1k = a1k.real
    print corrmax
    print a1k

    a1profTmp = ma.masked_less(a2profave[:,18],0)
    a1gvave   = ma.masked_less(a1gvave, 0)
    corrTmp = np.ma.corrcoef(a1profTmp, a1gvave, allow_masked=True)[0,1]
    print corrmax, corrTmp
    print a2joinprof.min(), a1gvave.min()

    for x in a1k:
        print x
