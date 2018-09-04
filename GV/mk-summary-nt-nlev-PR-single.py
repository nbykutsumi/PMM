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
import math
import scipy.stats

calc = True
#calc = False
iYM = [2005,4]
#iYM = [2014,4]
eYM = [2014,10]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
print lYM

corrFlag= 'CR'
#corrFlag= 'NO'

#nozero = 'withzero'
nozero = 'nozero'

#dh  = 18
dt  = 30
#lbnd = [[0.1,1],[1,5],[5,10],[10,30],[30,50],[50,70],[70,999],[0,999]]
lbnd = [[0.1,1],[1,5],[5,10],[10,30],[30,999],[0,999]]
#lbnd = [[30,999]]

cls = 'RH'
#cls = 'RainType'

#thdist = 2.5
thdist = 5.0
minNum = 3
prdName = 'L2A25'

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
    lclstype = ['dry','hum']
elif cls=='RainType':
    lclstype = ['alltype','strat','conv']
else:
    print 'check cls',cls
    sys.exit()

#----------------------------------


def test_cc(n1,n2,cc1,cc2):
    '''
    http://www.sigmath.es.osaka-u.ac.jp/~kano/old/lecture/faq/q1.html
    '''
    z1 = 0.5*math.log((1+cc1)/(1-cc1))
    z2 = 0.5*math.log((1+cc2)/(1-cc2))
    z  = np.sqrt((n1-3)*(n2-3)/(n1+n2-6))*(z1-z2)
    if abs(z)>1.96:
        #return 1
        return z
    else:
        #return 0
        return z

def test_cc_paired(n,ccesurf,ccprof, ccesurfprof):
    '''
    http://www.sigmath.es.osaka-u.ac.jp/~kano/old/lecture/faq/q1.html
    '''
    z1 = 0.5*math.log((1+ccesurf)/(1-ccesurf))
    z2 = 0.5*math.log((1+ccprof)/(1-ccprof))
    bunshi =np.sqrt(n-3)*(z1-z2)    
    
    r12 = ccesurf
    r13 = ccesurfprof
    r14 = ccesurf
    r23 = ccprof
    r24 = 1.0
    r34 = ccprof  

    a1  = r13*r24 + r14*r23
    a2  = -r34*(r13*r23 + r14*r24)
    a3  = -r12*(r13*r14 + r23*r24)
    a4  = r12*r34*(r13**2 + r14**2 + r23**2 + r24**2)/2.0
    d   = (1-r12**2)*(1-r34**2)

    bunbo= np.sqrt(2-2*(a1+a2+a3+a4)/d)
    z    = bunshi / bunbo

    if abs(z)>1.96:
        #return 1
        return z
    else:
        #return 0
        return z

#def test_bias_paired(x1, x2):
#    n  = len(x1)
#    d  = x1 - x2
#    dm = d.mean()
#    sd = np.sqrt( ((d-dm)**2).sum()   /(n-1))
#    t  = dm / (sd/np.sqrt(n))
#    return t

def test_bias_paired(x1, x2):
    t, p = scipy.stats.ttest_rel(x1, x2)
    return t


#def test_rmse(a1esurf, a1prof, a1gv):
#    a1msk1 = ma.masked_less(a1esurf,0).mask
#    a1msk2 = ma.masked_less(a1esurf,0).mask
#    a1msk  = a1msk1 + a1msk2
#    a1esurf = ma.masked_where(a1msk, a1esurf).compressed()
#    a1prof  = ma.masked_where(a1msk, a1prof ).compressed()
#    a1gv    = ma.masked_where(a1msk, a1gv   ).compressed()
#    
#    ntimes = 100
#    ndat   = len(a1esurf)
#    for i in range(ntimes):
#        a1idx = np.random.choice(ndat,ndat,replace=True)
#        a1esurfTmp = a1esurf[a1idx]
#        a1profTmp  = a1prof[a1idx]
#        a1gvTmp    = a1gv[a1idx]
        

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

for clstype in lclstype:

    csvDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-%s/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])

    joinprofPath = csvDir + '/%s.%s.joinprof.%s.csv'%(cls, nozero, clstype)
    gvPath       = csvDir + '/%s.%s.gv.%s.csv'%(cls, nozero, clstype)
    profavePath  = csvDir + '/%s.%s.profave.%s.csv'%(cls, nozero, clstype)
    gvavePath    = csvDir + '/%s.%s.gvave.%s.csv'%(cls, nozero, clstype)
   
    a2joinprof = load_csv(joinprofPath)
    a2gv       = load_csv(gvPath)[:,15:]
    a2profave  = load_csv(profavePath)
    a2gvave    = load_csv(gvavePath)

    #----------------------------
    #if   clstype == 'dry': dh = int(6.0/0.25)
    #elif clstype == 'hum': dh = int(3.5/0.25)

    #if   clstype == 'dry': dh = int(6.0/0.25)
    #elif clstype == 'hum': dh = int(5.5/0.25)

    if   clstype == 'all': dh = int(4.5/0.25)
    elif clstype == 'dry': dh = int(4.5/0.25)
    elif clstype == 'hum': dh = int(4.5/0.25)


    a1esurfTmp = a2profave[:,0]
    a1profTmp  = a2profave[:,dh]
    a1gvTmp    = a2gvave[:,dt-1]

    a1gvTmp0   = a2gvave[:,0]

    a1esurf.extend(a1esurfTmp)
    a1prof.extend(a1profTmp)
    a1gv.extend(a1gvTmp)
    a1gv0.extend(a1gvTmp0)

a1esurf = array(a1esurf)
a1prof  = array(a1prof)
a1gv    = array(a1gv)
a1gv0   = array(a1gv0)


#-- correction ---
if corrFlag == 'CR':
    a1mskTmp1 = ma.masked_less(a1gv0,  0).mask
    a1mskTmp2 = ma.masked_less(a1esurf,0).mask
    a1mskTmp  = a1mskTmp1 + a1mskTmp2 

    gvTmp    = ma.masked_where(a1mskTmp, a1gv0  ).mean()
    esurfTmp = ma.masked_where(a1mskTmp, a1esurf).mean()
    bfactor  = gvTmp / esurfTmp

    a1esurf = (ma.masked_less(a1esurf, 0) * bfactor).data
    a1prof  = (ma.masked_less(a1prof , 0) * bfactor).data
    a1gv    = (ma.masked_less(a1gv   , 0) * bfactor).data
    a1gv0   = (ma.masked_less(a1gv0  , 0) * bfactor).data





for bnd in lbnd:
    prmin, prmax = bnd
    #a1mskTmp = ma.masked_outside(a1gv, prmin, prmax)
    a1mskTmp = ma.masked_outside(a1esurf, prmin, prmax)
    a1mskTmp = ma.masked_where(a1esurf<0, a1mskTmp)
    a1mskTmp = ma.masked_where(a1prof<0,  a1mskTmp)
    a1mskTmp = a1mskTmp.mask

    a1esurfTmp = ma.masked_where(a1mskTmp, a1esurf)
    a1profTmp  = ma.masked_where(a1mskTmp, a1prof)
    a1gvTmp    = ma.masked_where(a1mskTmp, a1gv)
    a1gvTmp0   = ma.masked_where(a1mskTmp, a1gv0)

    esurfTmp = a1esurfTmp.mean()
    profTmp  = a1profTmp.mean()
    gvTmp    = a1gvTmp.mean()
    numTmp   = a1esurfTmp.count()


    biasesurf = esurfTmp - gvTmp
    biasprof  = profTmp  - gvTmp
    rbiasesurf= biasesurf / gvTmp
    rbiasprof = biasprof  / gvTmp

    #-- metrics --
    nesurf    = a1esurfTmp.count()
    nprof     = a1profTmp.count()

    ccesurf   = np.ma.corrcoef(a1esurfTmp, a1gvTmp, allow_masked=True)[0,1]
    ccprof    = np.ma.corrcoef(a1profTmp,  a1gvTmp, allow_masked=True)[0,1]

    rmseesurf = np.sqrt( ((a1esurfTmp - a1gvTmp)**2).mean() )
    rmseprof  = np.sqrt( ((a1profTmp  - a1gvTmp)**2).mean() )

    #-- statistical test --
    ccesurfprof= np.ma.corrcoef(a1esurfTmp,  a1profTmp, allow_masked=True)[0,1]
    cctest  = test_cc(nesurf, nprof, ccesurf, ccprof)
    cctest2 = test_cc_paired(nesurf, ccesurf, ccprof, ccesurfprof)

    a1biasesurf = (a1esurfTmp - a1gvTmp).compressed()
    a1biasprof  = (a1profTmp  - a1gvTmp).compressed()

    biastest= test_bias_paired( a1biasesurf, a1biasprof)
    #t, p = test_bias_paired2( a1biasesurf, a1biasprof)
    
     


    sout ='%.1f %.1f %05d %.4f %.4f %.4f %.4f %.4f %.4f'%(prmin, prmax, numTmp, ccesurf, ccprof, rbiasesurf, rbiasprof, rmseesurf, rmseprof)

    #sout = sout + ' %.3f %.3f'%(cctest2, biastest)
    print sout
sys.exit()

