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
#iYM = [2010,4]
eYM = [2014,10]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
print lYM

#thdist = 2.5
thdist = 5.0
minNum = 3
prdName = 'L2A25'

nozero  = 'nozero'
#nozero  = 'withzero'

gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()

dgName = gv.ret_ddomYM2gName()

ldomain = gv.domains
#ldomain = ['FLORIDA-STJ','FLORIDA-SFL-N','N.Carolina-IPHEx_Duke','N.Carolina-IPHEx_NASA','KWAJALEIN-KWA']
#ldomain = ['FLORIDA-STJ']

offset_bef = 15  # 'bef' should be identical to that in mk_match.py
#offset_aft = 30

nh   = 40
nt   = 30
#ldt  = [1, 5,10,15,20,25,30]
#ldh  = range(2,34+1)[::2]

#nt = len(ldt)
#nh = len(ldh)

lrhtype = ['all','dry','hum']
#lrhtype = ['hum']
dlthrh  = {'all':[-0.1,9999],'dry':[-0.1,0.7-0.0001],'hum':[0.7,9999]}

miss = -9999.
#------------------------------------------------
def ret_a2cumave_negativemask(a2dat):
    miss  = -9999.
    a2cum = empty(a2dat.shape)
    for i in range(a2dat.shape[1]):
        a1tmp = ma.masked_less(a2dat[:,:i+1],0).mean(axis=1).filled(miss)
        a2cum[:,i] = a1tmp
    return a2cum
#------------------------------------------------



a2prof     = deque([]) 
a1esurf    = deque([]) 
a1nsurf    = deque([]) 
a2gv       = deque([]) 
a1rh       = deque([])
a1rtype    = deque([])
a1stormH   = deque([])

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
        rhPath    = srcDir  + '/p_rh.npy'
        rainTypePath = srcDir  + '/p_rainType.npy'
        stormHPath = srcDir  + '/p_stormH.npy'

        
        if not os.path.exists(profPath):
            print 'no file',profPath
            print 'skip', domain, YM
            continue

        aprof  = np.load(profPath)[:,:nh]
        aesurf = np.load(eSurfPath)
        ansurf = np.load(nSurfPath)
        agv    = abs(np.load(gvPath))[:,:15+nt]
        angv   = np.load(ngvPath)
        arh    = np.load(rhPath)
        artype = np.load(rainTypePath)
        astormH = np.load(stormHPath)


        #-- mean array for masks
        #a1sateAll = ma.masked_less(aprof[:,:nh],0).mean(axis=1).filled(miss)
        a1sateAll = ma.masked_less(aprof,0).mean(axis=1).filled(miss)

        #a1gvAll   = ma.masked_less(agv[:,15:15+nt],0).mean(axis=1)
        a1gvAll   = ma.masked_less(agv,0).mean(axis=1)

        #-- mask when both satellite and gv are zero or miss
        if nozero=='nozero':
            amsk1    = ma.masked_less_equal(a1sateAll,0).mask
            amsk2    = ma.masked_less_equal(a1gvAll,0).mask
            amsk3    = ma.masked_less_equal(aesurf,0).mask
            amskzero = amsk1 * amsk2 * amsk3

        elif nozero=='withzero':
            amsk1    = ma.masked_less(a1sateAll,0).mask
            amsk2    = ma.masked_less(a1gvAll,0).mask
            amsk3    = ma.masked_less(aesurf,0).mask
            amskzero = amsk1 * amsk2 * amsk3

        else:
            print 'check nozero',nozero
            sys.exit()


        #-- mask when ng  < minNum
        amskN    = ma.masked_less(angv, minNum).mask

        #-- overlay masks
        amsk     = amskzero + amskN 

        aidxTmp  = arange(len(amsk)).astype(int32)
        aidxTmp  = ma.masked_where(amsk, aidxTmp).compressed()

        aprofTmp     = aprof[aidxTmp] 
        aesurfTmp    = aesurf[aidxTmp]
        ansurfTmp    = ansurf[aidxTmp]
        agvTmp       = agv[aidxTmp]
        arhTmp       = arh[aidxTmp]
        artypeTmp    = artype[aidxTmp]
        astormHTmp   = astormH[aidxTmp]
        #------------------
    
        a2gv.append(agvTmp)
        a2prof.append(aprofTmp)
        a1esurf.extend(aesurfTmp)
        a1nsurf.extend(ansurfTmp)
        a1rh.extend(arhTmp)
        a1rtype.extend(artypeTmp)
        a1stormH.extend(astormHTmp)

a1esurf= array(a1esurf)
a1nsurf= array(a1nsurf)*0.01 #mm/h
a1rh   = array(a1rh)
a1rtype= array(a1rtype)
a1stormH=array(a1stormH)
a2gv   = concatenate(a2gv,axis=0)
a2prof = concatenate(a2prof,axis=0) *0.01
a2joinprof= concatenate([a1esurf.reshape(-1,1) ,a2prof], axis=1)


#----
a2profave = ret_a2cumave_negativemask(a2joinprof)

a2gvave   = ret_a2cumave_negativemask(a2gv[:,15:])

# save file
outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-%s/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])
util.mk_dir(outDir)


joinprofPath = outDir + '/ALL.%s.joinprof.csv'%(nozero)
gvoutPath    = outDir + '/ALL.%s.gv.csv'%(nozero)
profavePath  = outDir + '/ALL.%s.profave.csv'%(nozero)
gvavePath    = outDir + '/ALL.%s.gvave.csv'%(nozero)

rhPath       = outDir + '/ALL.%s.RH.csv'%(nozero)
rtypePath    = outDir + '/ALL.%s.rainType.csv'%(nozero)
stormHPath   = outDir + '/ALL.%s.stormH.csv'%(nozero)


a2joinprof   = ma.masked_less(a2joinprof,0).filled(-9999.)
a2profave    = ma.masked_less(a2profave,0).filled(-9999.)

sjoinprof  = util.array2csv(a2joinprof)
sgv        = util.array2csv(a2gv)
sprofave   = util.array2csv(a2profave)
sgvave     = util.array2csv(a2gvave)
sRH        = util.array2csv(a1rh)
srainType  = util.array2csv(a1rtype)
sstormH    = util.array2csv(a1stormH)



f = open(joinprofPath, 'w'); f.write(sjoinprof); f.close()
f = open(gvoutPath, 'w'); f.write(sgv); f.close()
f = open(profavePath, 'w'); f.write(sprofave); f.close()
f = open(gvavePath, 'w'); f.write(sgvave); f.close()
f = open(rhPath, 'w'); f.write(sRH); f.close()
f = open(rtypePath, 'w'); f.write(srainType); f.close()
f = open(stormHPath, 'w'); f.write(sstormH); f.close()


print joinprofPath



