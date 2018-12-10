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
#iYM = [2005,4]
iYM = [2005,4]
eYM = [2014,10]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
print lYM

#corrFlag= 'CR'
corrFlag= 'NO'

#nozero = 'withzero'
nozero = 'nozero'
nboot = 1000

#dh  = 18
dt  = 30
lbnd = [[0.1,1],[1,5],[5,10],[10,20],[20,999],[0,999]]
#lbnd = [[0,999]]

cls = 'RainType'

#-- load profile adjustment factors based on climatology --
llstormH = [[2,3],[3,4],[4,6],[6,8],[8,15]]
ratDir   = '/home/utsumi/mnt/wellshare/GPMGV/DOM.L2A25/cbin.11/ratio'
llRH     = [[0.0,0.6],[0.6,0.7],[0.7,0.8],[0.8,0.9],[0.9,1.0]]


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

dlthrtype = {'alltype':[-999,999],'strat':[100,170],'conv':[200,297]}

miss  = -9999

if cls=='RH':
    #lclstype = ['all','dry','hum']
    #lclstype = ['dry','hum']
    #lclstype = ['all']
    #lclstype = ['dry']
    lclstype = ['hum']
elif cls=='RainType':
    #lclstype = ['alltype','strat','conv']
    lclstype = ['alltype']
    #lclstype = ['strat']
    #lclstype = ['conv']
else:
    print 'check cls',cls
    sys.exit()

#----------------------------------
def ret_a2cumave_negativemask(a2dat):
    miss  = -9999.
    a2cum = empty(a2dat.shape)
    for i in range(a2dat.shape[1]):
        a1tmp = ma.masked_less(a2dat[:,:i+1],0).mean(axis=1).filled(miss)
        a2cum[:,i] = a1tmp
    return a2cum
#------------------------------------------------

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

a1profAdj  = deque([])

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



    #** load profile adjustment factors based on climatology **
    llstormH = [[2,3],[3,4],[4,6],[6,8],[8,15]]
    ratDir   = '/home/utsumi/mnt/wellshare/GPMGV/DOM.L2A25/cbin.11/ratio'
    llRH     = [[0.0,0.6],[0.6,0.7],[0.7,0.8],[0.8,0.9],[0.9,1.0]]
    
    a1height = arange(a2joinprof.shape[1])*0.25
    dratio   = {}

    #for raintype in ['alltype','strat','conv']:
    for raintype in ['alltype']:
        for istormH,lstormH in enumerate(llstormH):
            stormH0,stormH1 = lstormH
            ratPath = ratDir + '/ratio.%s.H%d-%d.npy'%(raintype,stormH0,stormH1)
            a2tmp   = np.load(ratPath)
            for iRH,lRH in enumerate(llRH):
                dratio[raintype,istormH,iRH] = a2tmp[:,iRH]
    

    #******* adjustment ********
    csvDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-%s/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])
    
    joinprofAllPath = csvDir + '/ALL.%s.joinprof.csv'%(nozero)
    gvAllPath       = csvDir + '/ALL.%s.gv.csv'%(nozero)
    gvaveAllPath    = csvDir + '/ALL.%s.gvave.csv'%(nozero)
    rhAllPath       = csvDir + '/ALL.%s.RH.csv'%(nozero)
    rtypeAllPath    = csvDir + '/ALL.%s.rainType.csv'%(nozero)
    stormHAllPath   = csvDir + '/ALL.%s.stormH.csv'%(nozero)
    
    a2joinprofAll = load_csv(joinprofAllPath)
    a2gvAll       = load_csv(gvAllPath)[:,15:]
    a2gvaveAll    = load_csv(gvaveAllPath)
    a1rhAll       = load_csv(rhAllPath).flatten()
    a1rtypeAll    = load_csv(rtypeAllPath).flatten()
    a1stormHAll   = load_csv(stormHAllPath).flatten() * 0.001

    thrtype0,thrtype1 = dlthrtype[clstype]
    a1mskT  = ma.masked_outside(a1rtypeAll, thrtype0,thrtype1).mask

    #-- profile adjustment ------
    a1idx = arange(len(a1stormHAll))
    for istormH,lstormH in enumerate(llstormH):
        if istormH == len(llstormH):
            stormH0,stormH1 = lstormH[0],999
        else:
            stormH0,stormH1 = lstormH

        a1mskSH = ma.masked_less(a1stormHAll,stormH0)
        a1mskSH = ma.masked_greater_equal(a1mskSH, stormH1).mask

        for iRH,lRH in enumerate(llRH):
            RH0,RH1 = lRH
            a1mskRH = ma.masked_less(a1rhAll, RH0)
            a1mskRH = ma.masked_greater_equal(a1mskRH, RH1).mask

            a1msk   = a1mskT + a1mskSH + a1mskRH
            a1idxTmp = ma.masked_where(a1msk, a1idx).compressed()

            a1ratio  = dratio[clstype,istormH,iRH]

            a2joinprofAll[a1idxTmp,:] = ma.masked_equal(a2joinprofAll[a1idxTmp,:],miss)*a1ratio

    a1idxT = ma.masked_where(a1mskT, a1idx).compressed()
    a2joinprofAdj = a2joinprofAll[a1idxT]
    #a2gvaveAdj    = a2gvaveAll[a1idxT]

    a2profaveAdj = ret_a2cumave_negativemask(a2joinprofAdj)

    #print a2joinprof.shape, a2joinprofAdj.shape
    #sys.exit()


    #----------------------------
    if   clstype == 'all': dh = int(4.5/0.25)
    elif clstype == 'dry': dh = int(4.5/0.25)
    elif clstype == 'hum': dh = int(4.5/0.25)
    elif clstype == 'alltype': dh = int(4.5/0.25)
    elif clstype == 'strat'  : dh = int(4.5/0.25)
    elif clstype == 'conv'   : dh = int(4.5/0.25)


    a1esurfTmp = a2profave[:,0]
    a1profTmp  = a2profave[:,dh]
    a1gvTmp    = a2gvave[:,dt-1]

    a1gvTmp0   = a2gvave[:,0]

    a1profAdjTmp  = a2profaveAdj[:,dh]

    #-- mask missing ----
    a1mskTmp1 = ma.masked_less(a1esurfTmp, 0).mask
    a1mskTmp2 = ma.masked_less(a1profTmp,  0).mask
    a1mskTmp3 = ma.masked_less(a1gvTmp,    0).mask
    a1mskTmp  = a1mskTmp1 + a1mskTmp2 + a1mskTmp3

    a1esurfTmp = ma.masked_where(a1mskTmp, a1esurfTmp).compressed()
    a1profTmp  = ma.masked_where(a1mskTmp, a1profTmp).compressed()
    a1gvTmp    = ma.masked_where(a1mskTmp, a1gvTmp).compressed()
    a1gvTmp0   = ma.masked_where(a1mskTmp, a1gvTmp0).compressed()


    a1profAdjTmp  = ma.masked_where(a1mskTmp, a1profAdjTmp).compressed()
    #------------------

    a1esurf.extend(a1esurfTmp)
    a1prof.extend(a1profTmp)
    a1gv.extend(a1gvTmp)
    a1gv0.extend(a1gvTmp0)

    a1profAdj.extend(a1profAdjTmp)


a1esurf = array(a1esurf)
a1prof  = array(a1prof)
a1gv    = array(a1gv)
a1gv0   = array(a1gv0)

a1profAdj  = array(a1profAdj)

#-- correction ---
if corrFlag == 'CR':
    gvTmp    = a1gv0.mean()
    esurfTmp = a1esurf.mean()
    bfactor  = gvTmp / esurfTmp

    a1esurf = a1esurf * bfactor
    a1prof  = a1prof  * bfactor
    a1gv    = a1gv    * bfactor
    a1gv0   = a1gv0   * bfactor

    a1profAdj= a1profAdj * bfactor


#-- bnd loop -----
for bnd in lbnd:
    prmin, prmax = bnd
    #a1mskTmp = ma.masked_outside(a1gv, prmin, prmax)
    a1mskTmp = ma.masked_outside(a1esurf, prmin, prmax)
    a1mskTmp = a1mskTmp.mask

    a1esurfTmp = ma.masked_where(a1mskTmp, a1esurf).compressed()
    a1profTmp  = ma.masked_where(a1mskTmp, a1prof).compressed()
    a1gvTmp    = ma.masked_where(a1mskTmp, a1gv).compressed()
    a1gvTmp0   = ma.masked_where(a1mskTmp, a1gv0).compressed()

    a1profAdjTmp=ma.masked_where(a1mskTmp, a1profAdj).compressed()


    esurfTmp = a1esurfTmp.mean()
    profTmp  = a1profTmp.mean()
    gvTmp    = a1gvTmp.mean()
    #numTmp   = a1esurfTmp.count()
    numTmp   = len(a1esurfTmp)

    profAdjTmp=a1profAdjTmp.mean()
    #----------------
    biasesurf = esurfTmp - gvTmp
    biasprof  = profTmp  - gvTmp
    rbiasesurf= biasesurf / gvTmp
    rbiasprof = biasprof  / gvTmp

    biasprofAdj =profAdjTmp - gvTmp
    rbiasprofAdj=biasprofAdj / gvTmp
    #-- metrics --
    nesurf    = len(a1esurfTmp)
    nprof     = len(a1profTmp)

    ccesurf   = np.ma.corrcoef(a1esurfTmp, a1gvTmp, allow_masked=True)[0,1]
    ccprof    = np.ma.corrcoef(a1profTmp,  a1gvTmp, allow_masked=True)[0,1]
    ccprofAdj = np.ma.corrcoef(a1profAdjTmp,a1gvTmp, allow_masked=True)[0,1]

    rmseesurf = np.sqrt( ((a1esurfTmp - a1gvTmp)**2).mean() )
    rmseprof  = np.sqrt( ((a1profTmp  - a1gvTmp)**2).mean() )
    rmseprofAdj=np.sqrt( ((a1profAdjTmp-a1gvTmp)**2)
.mean())

    changermse= (rmseprof - rmseesurf)/rmseesurf
    changerbias=(abs(rbiasprof)-abs(rbiasesurf))/abs(rbiasesurf)
    changecc   =(ccprof - ccesurf)/abs(ccesurf)

    changermseAdj= (rmseprofAdj - rmseesurf)/rmseesurf
    changerbiasAdj=(abs(rbiasprofAdj)-abs(rbiasesurf))/abs(rbiasesurf)
    changeccAdj   =(ccprofAdj - ccesurf)/abs(ccesurf)

    ##-- statistical test --
    #ccesurfprof= np.ma.corrcoef(a1esurfTmp,  a1profTmp, allow_masked=True)[0,1]
    #cctest  = test_cc(nesurf, nprof, ccesurf, ccprof)
    #cctest2 = test_cc_paired(nesurf, ccesurf, ccprof, ccesurfprof)

    #a1biasesurf = (a1esurfTmp - a1gvTmp).compressed()
    #a1biasprof  = (a1profTmp  - a1gvTmp).compressed()

    #biastest= test_bias_paired( a1biasesurf, a1biasprof)
    #t, p = test_bias_paired2( a1biasesurf, a1biasprof)

    sout ='%.1f %.1f %05d %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f'%(prmin, prmax, numTmp, ccesurf, ccprof, ccprofAdj, rbiasesurf, rbiasprof, rbiasprofAdj, rmseesurf, rmseprof, rmseprofAdj)


    #-- bootstrap ---------------
    lnum       = []
    lrbiasprof = []
    lccprof    = []
    lrmseprof  = []
    lrbiasesurf= []
    lccesurf   = []
    lrmseesurf = []

    lrbiasprofAdj = []
    lccprofAdj    = []
    lrmseprofAdj  = []


    lchangecc  = []
    lchangerbias=[]
    lchangermse =[] 

    lchangeccAdj   = []
    lchangerbiasAdj=[]
    lchangermseAdj =[] 

    numAll = len(a1esurf)
    for iboot in range(nboot):
        a1idx     = np.random.choice(numAll, numAll, replace=True)
        a1esurfTmp2 = a1esurf[a1idx]  
        a1profTmp2  = a1prof [a1idx]
        a1gvTmp2    = a1gv   [a1idx]

        a1profAdjTmp2 = a1profAdj[a1idx]

        #- bin --
        prmin, prmax = bnd
        a1mskTmp = ma.masked_outside(a1esurfTmp2, prmin, prmax)
        a1mskTmp = a1mskTmp.mask
    
        a1esurfTmp2 = ma.masked_where(a1mskTmp, a1esurfTmp2).compressed()
        a1profTmp2  = ma.masked_where(a1mskTmp, a1profTmp2).compressed()
        a1gvTmp2    = ma.masked_where(a1mskTmp, a1gvTmp2).compressed()

        a1profAdjTmp2  = ma.masked_where(a1mskTmp, a1profAdjTmp2).compressed()
        #----------------
        num       = len(a1esurfTmp2)
        esurfTmp2 = a1esurfTmp2.mean()
        profTmp2  = a1profTmp2.mean()
        gvTmp2    = a1gvTmp2.mean()

        profAdjTmp2 = a1profAdjTmp2.mean()

        biasesurf2 = esurfTmp2 - gvTmp2
        biasprof2  = profTmp2  - gvTmp2
        biasprofAdj2 = profAdjTmp2 - gvTmp2

        rbiasesurf2= biasesurf2 / gvTmp2
        rbiasprof2 = biasprof2  / gvTmp2
        rbiasprofAdj2 = biasprofAdj2 / gvTmp2 

        #-- metrics for bootstrap samples--
        ccesurf2   = np.ma.corrcoef(a1esurfTmp2, a1gvTmp2, allow_masked=True)[0,1]
        ccprof2    = np.ma.corrcoef(a1profTmp2,  a1gvTmp2, allow_masked=True)[0,1]
        ccprofAdj2 = np.ma.corrcoef(a1profAdjTmp2, a1gvTmp2, allow_masked=True)[0,1]
    
        rmseesurf2 = np.sqrt( ((a1esurfTmp2 - a1gvTmp2)**2).mean() )
        rmseprof2  = np.sqrt( ((a1profTmp2  - a1gvTmp2)**2).mean() )
        rmseprofAdj2 = np.sqrt( ((a1profAdjTmp2 - a1gvTmp2)**2).mean() )
    
        changecc2   = (ccprof2   -ccesurf2   )/ccesurf2
        changeccAdj2= (ccprofAdj2 - ccesurf2 )/ccesurf2

        changerbias2= (rbiasprof2-rbiasesurf2)/rbiasesurf2
        changerbiasAdj2 = (rbiasprofAdj2-rbiasesurf2)/rbiasesurf2

        changermse2 = (rmseprof2 -rmseesurf2 )/rmseesurf2
        changermseAdj2 = (rmseprofAdj2 - rmseesurf2)/rmseesurf2
        #-- 
        lnum        .append(num)
        lrbiasprof  .append(rbiasprof2)
        lccprof     .append(ccprof2)
        lrmseprof   .append(rmseprof2)
        lrbiasesurf .append(rbiasesurf2)
        lccesurf    .append(ccesurf2)
        lrmseesurf  .append(rmseesurf2)

        lrbiasprofAdj.append(rbiasprofAdj2)
        lccprofAdj  .append(ccprofAdj2)
        lrmseprofAdj.append(rmseprofAdj2)

        lchangecc   .append(changecc2)
        lchangerbias.append(changerbias2)
        lchangermse .append(changermse2)

        lchangeccAdj   .append(changeccAdj2)
        lchangerbiasAdj.append(changerbiasAdj2)
        lchangermseAdj .append(changermseAdj2)

    #--- metrics for bootstrap erro range --
    lrbiasprof  = ma.masked_invalid(lrbiasprof)
    lrbiasesurf = ma.masked_invalid(lrbiasesurf)
    lrbiasprofAdj  = ma.masked_invalid(lrbiasprofAdj)

    lccprof     = ma.masked_invalid(lccprof)
    lccesurf    = ma.masked_invalid(lccesurf)
    lccprofAdj  = ma.masked_invalid(lccprofAdj)

    lrmseprof   = ma.masked_invalid(lrmseprof)
    lrmseesurf  = ma.masked_invalid(lrmseesurf)
    lrmseprofAdj= ma.masked_invalid(lrmseprofAdj)

    lchangecc   = ma.masked_invalid(lchangecc)
    lchangerbias= ma.masked_invalid(lchangerbias)
    lchangermse = ma.masked_invalid(lchangermse)

    lchangeccAdj   = ma.masked_invalid(lchangeccAdj)
    lchangerbiasAdj= ma.masked_invalid(lchangerbiasAdj)
    lchangermseAdj = ma.masked_invalid(lchangermseAdj)

    stdrbiasprof  = np.std(lrbiasprof)
    stdrbiasesurf = np.std(lrbiasesurf)
    stdrbiasprofAdj = np.std(lrbiasprofAdj)

    stdccprof     = np.std(lccprof)
    stdccesurf    = np.std(lccesurf)
    stdccprofAdj  = np.std(lccprofAdj)

    stdrmseprof   = np.std(lrmseprof)
    stdrmseesurf  = np.std(lrmseesurf)
    stdrmseprofAdj= np.std(lrmseprofAdj)

    stdchangecc   = np.std(lchangecc)
    stdchangerbias= np.std(lchangerbias)
    stdchangermse = np.std(lchangermse)

    stdchangeccAdj   = np.std(lchangeccAdj)
    stdchangerbiasAdj= np.std(lchangerbiasAdj)
    stdchangermseAdj = np.std(lchangermseAdj)


    sout ='%.1f %.1f %05d %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f'%(prmin, prmax, numTmp, ccesurf, ccprof, ccprofAdj, changecc, changeccAdj, stdccesurf, stdccprof, stdccprofAdj, stdchangecc, stdchangeccAdj, rbiasesurf, rbiasprof, rbiasprofAdj, changerbias, changerbiasAdj, stdrbiasesurf, stdrbiasprof, stdrbiasprofAdj, stdchangerbias, stdchangerbiasAdj, rmseesurf, rmseprof, rmseprofAdj, changermse, changermseAdj, stdrmseesurf, stdrmseprof, stdrmseprofAdj, stdchangermse, stdchangermseAdj)
    print sout

print 'prmin, prmax, numTmp, ccesurf, ccprof, ccprofAdj, changecc, changeccAdj, stdccesurf, stdccprof, stdccprofAdj, stdchangecc, stdchangeccAdj, rbiasesurf, rbiasprof, rbiasprofAdj, changerbias, changerbiasAdj, stdrbiasesurf, stdrbiasprof, stdrbiasprofAdj, stdchangerbias, stdchangerbiasAdj, rmseesurf, rmseprof, rmseprofAdj, changermse, changermseAdj, stdrmseesurf, stdrmseprof, stdrmseprofAdj, stdchangermse, stdchangermseAdj'
