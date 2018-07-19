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


calc = True
#calc = False
iYM = [2005,4]
#iYM = [2014,10]
eYM = [2014,10]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
print lYM
figtype = 'bar'
#figtype = 'cont'
#corrFlag= 'CR'
corrFlag= 'NO'



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

if   figtype=='bar':
    ldt  = [1, 5,10,15,20,25,30]
    ldh  = range(2,34+1)[::2]
elif figtype=='cont':
    ldt  = range(1,30+1)
    ldh  = range(1,34+1)


nt = len(ldt)
nh = len(ldh)
#nh = 5

lrhtype = ['all','dry','hum']
#lrhtype = ['hum']
dlthrh  = {'all':[-0.1,9999],'dry':[-0.1,0.7-0.0001],'hum':[0.7,9999]}

ldattype = ['rain','cc','bias','brat','rmse','gv','num']

a2prof     = deque([]) 
a1esurf    = deque([]) 
a1nsurf    = deque([]) 
a2gv       = deque([]) 
a1rh       = deque([])

for domain in ldomain:
    if calc ==False: continue

    for YM in lYM:
        Year, Mon = YM
        if (domain,Year,Mon) not in dgName.keys():
            print 'no obs',domain,Year,Mon
            continue
        # load data
        srcbaseDir = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.L2A25'
        srcDir     = srcbaseDir + '/%.1fkm/%s/%04d%02d'%(thdist,domain, Year,Mon)

        profPath  = srcDir  + '/p_prof.npy'
        eSurfPath = srcDir  + '/p_eSurf.npy'
        nSurfPath = srcDir  + '/p_nSurf.npy'
        gvPath    = srcDir  + '/p_gvprcp.npy'
        ngvPath   = srcDir  + '/p_ngv.npy'
        nSurfBinPath = srcDir  + '/p_nSurfBin.npy'
        rhPath    = srcDir  + '/p_rh.npy'
        

        if not os.path.exists(profPath):
            print 'no file',profPath
            print 'skip', domain, YM
            continue

        aprof  = np.load(profPath)
        aesurf = np.load(eSurfPath)
        ansurf = np.load(nSurfPath)
        agv    = abs(np.load(gvPath))
        angv   = np.load(ngvPath)
        arh    = np.load(rhPath)

        asate  = aprof

        #-- mean array for masks
        agroundBin= zeros(len(angv))
        a1sateAll = gv_fsub.mean_slice_negativemask(aprof.T, agroundBin, ldh[-1])

        a1idx15   = ones(agv.shape[0])*15
        a1gvAll   = ma.masked_less(agv[:,15:15+30],0).mean(axis=1)

        ##-- mask when both satellite and gv are zero or miss
        #amsk1    = ma.masked_equal(a1sateAll,0).mask
        #amsk2    = ma.masked_equal(a1gvAll,0).mask
        #amskzero = amsk1 * amsk2

        #-- mask when both satellite and gv are zero or miss
        amsk1    = ma.masked_less_equal(a1sateAll,0).mask
        amsk2    = ma.masked_less_equal(a1gvAll,0).mask
        amsk3    = ma.masked_less_equal(aesurf,0).mask
        amskzero = amsk1 * amsk2 * amsk3


        #-- mask when ng  < minNum
        amskN    = ma.masked_less(angv, minNum).mask

        #-- overlay masks
        amsk     = amskzero + amskN 

        aidxTmp  = arange(len(amsk)).astype(int32)
        aidxTmp  = ma.masked_where(amsk, aidxTmp).compressed()

        asateTmp     = asate[aidxTmp] 
        aesurfTmp    = aesurf[aidxTmp]
        ansurfTmp    = ansurf[aidxTmp]
        agvTmp       = agv[aidxTmp]
        arhTmp       = arh[aidxTmp]
        #------------------
    
        a2gv.append(agvTmp)
        a2prof.append(asateTmp)
        a1esurf.extend(aesurfTmp)
        a1nsurf.extend(ansurfTmp)
        a1rh.extend(arhTmp)

if calc ==True:
    a1esurf= array(a1esurf)
    a1nsurf= array(a1nsurf)*0.01 #mm/h
    a1rh   = array(a1rh)
    a2gv   = concatenate(a2gv,axis=0)
    a2prof = concatenate(a2prof,axis=0) *0.01
    a2joinprof= concatenate([a1esurf.reshape(-1,1) ,a2prof], axis=1)


#----
for rhtype in lrhtype:
    if calc==False: continue

    # mask with RH --
    thmin, thmax = dlthrh[rhtype]
    amskRH = ma.masked_outside(a1rh, thmin, thmax).mask

    aidxTmp = arange(a2joinprof.shape[0])
    aidxTmp = ma.masked_where(amskRH, aidxTmp).compressed()
    a2joinprofTmp = a2joinprof[aidxTmp,:]
    a2gvTmp = a2gv[aidxTmp,:]
    a1esurfTmp = a1esurf[aidxTmp]


    a2profave = empty([len(aidxTmp), len(ldh)+1])
    a2gvave   = empty([len(aidxTmp), len(ldt)])

    #--- vertical layer loop ---
    ih = -1
    for dh in [-99] + ldh: 
        if calc == False: continue
        ih     = ih + 1
    
        if dh == -99:
            a1profTmp = a1esurfTmp
        else:
            #a1prof = gv_fsub.mean_slice_negativemask(a2joinprof.T, a1idx0, dh+1)
            a1profTmp = ma.masked_less(a2joinprofTmp[:,:dh+1], 0).mean(axis=1)
    

        a2profave[:,ih] = a1profTmp

    #--- averating time loop ---
    for it,dt in enumerate(ldt):
        #a1idx15 = ones(a2gv.shape[0])*15
        #a1gv   = gv_fsub.mean_slice_negativemask(a2gv.T, a1idx15, dt)
        a1gv   = ma.masked_less(a2gv[:,15:15+dt],0).mean(axis=1)

    
        a1gvTmp   = ma.masked_where(amskRH, a1gv   ).compressed()

        a2gvave[:,it] =a1gvTmp

        '''
        #-- check miss --
        amskMiss1 = ma.masked_less(a1profTmp,0).mask
        amskMiss2 = ma.masked_less(a1gvTmp, 0).mask
        amskMiss  = amskMiss1 + amskMiss2
        a1profTmp = ma.masked_where(amskMiss, a1profTmp)
        a1gvTmp   = ma.masked_where(amskMiss, a1gvTmp)

        #-- bias correction --
        if corrFlag=='CR':
            print rhtype, dh,it, bfactor
            a1profTmp = a1profTmp * bfactor
        #--------------------- 
         
        cc     = np.corrcoef(a1profTmp, a1gvTmp)[0,1]
        rmse   = np.sqrt( ((a1profTmp- a1gvTmp)**2).mean() )
       
        bias   = (a1profTmp - a1gvTmp).mean()
        brat   = (a1profTmp - a1gvTmp).mean() / a1gvTmp.mean()*100
        num    = len(a1profTmp)
        rain   = a1profTmp.mean()
        gv     = a1gvTmp.mean()
        '''
 

    
    figDir  = '/work/a01/utsumi/GPMGV/fig'
    joinprofPath = figDir + '/temp.table.joinprof.%s.csv'%(rhtype)
    gvoutPath    = figDir + '/temp.table.gv.%s.csv'%(rhtype)
    profavePath  = figDir + '/temp.table.profave.%s.csv'%(rhtype)
    gvavePath    = figDir + '/temp.table.gvave.%s.csv'%(rhtype)

    a2joinprofTmp = ma.masked_less(a2joinprofTmp,0).filled(-9999.)
    a2profave     = ma.masked_less(a2profave,0).filled(-9999.)

    sjoinprof  = util.array2csv(a2joinprofTmp)
    sgv        = util.array2csv(a2gvTmp[:,15:15+30])
    sprofave   = util.array2csv(a2profave)
    sgvave     = util.array2csv(a2gvave)
    f = open(joinprofPath, 'w'); f.write(sjoinprof); f.close()
    f = open(gvoutPath, 'w'); f.write(sgv); f.close()
    f = open(profavePath, 'w'); f.write(sprofave); f.close()
    f = open(gvavePath, 'w'); f.write(sgvave); f.close()
    print joinprofPath



