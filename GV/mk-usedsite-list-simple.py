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

lraintype = ['alltype','strat','conv']
dlthrtype = {'alltype':[-1,999],'strat':[100,170],'conv':[200,297]}

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

dnumGauge   = {}
dnumObsSate = {}
dgauge= {}
for domain in ldomain:
    if calc ==False: continue

    lgNameTmp = deque([])
    numObsSate = 0

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

        gnamePath = srcDir  + '/p_gName.pickle'
        
        if not os.path.exists(profPath):
            print 'no file',profPath
            print 'skip', domain, YM
            continue

        aprof  = np.load(profPath)[:,:nh]
        aesurf = np.load(eSurfPath)
        ansurf = np.load(nSurfPath)
        agv    = abs(np.load(gvPath))[:,:15+nt]
        angv   = np.load(ngvPath)

        agName = np.load(gnamePath)

        #-- mean array for masks
        #a1sateAll = ma.masked_less(aprof[:,:nh],0).mean(axis=1).filled(miss)
        a1sateAll = ma.masked_less(aprof,0).mean(axis=1).filled(miss)

        #a1gvAll   = ma.masked_less(agv[:,15:15+nt],0).mean(axis=1)
        a1gvAll   = ma.masked_less(agv,0).mean(axis=1)

        #-- mask when both satellite and gv are zero or miss
        if nozero == 'nozero':
            amsk1    = ma.masked_less_equal(a1sateAll,0).mask
            amsk2    = ma.masked_less_equal(a1gvAll,0).mask
            amsk3    = ma.masked_less_equal(aesurf,0).mask
            amskzero = amsk1 * amsk2 * amsk3

        elif nozero == 'withzero':
            amsk1    = ma.masked_less(a1sateAll,0).mask
            amsk2    = ma.masked_less(a1gvAll,0).mask
            amsk3    = ma.masked_less(aesurf,0).mask
            amskzero = amsk1 * amsk2 * amsk3
        else:
            print 'check nozero', nozero
            sys.exit()

        #-- mask when ng  < minNum
        amskN    = ma.masked_less(angv, minNum).mask

        #-- overlay masks
        amsk     = amskzero + amskN 

        aidxTmp  = arange(len(amsk)).astype(int32)
        aidxTmp  = ma.masked_where(amsk, aidxTmp).compressed()

        agNameTmp= deque([])
        for idx in aidxTmp:
            agNameTmp.extend(agName[idx]) 

        lgNameTmp.extend(set(agNameTmp))
        numObsSate = numObsSate  + len(aidxTmp)

    lgNameTmp = set(lgNameTmp)
    dgauge[domain]= lgNameTmp
    dnumGauge[domain]  = len(lgNameTmp)
    dnumObsSate[domain]= numObsSate


#---- load site list ------
listDir = '/work/a01/utsumi/data/GPMGV/sitelist'
listPath= listDir + '/sitelist_reclassified.csv'
f=open(listPath,'r'); lines=f.readlines(); f.close()

for domain in sort(dgauge.keys()):
    if dnumGauge[domain]==0: continue

    llat = []
    llon = []
    for gName in dgauge[domain]:
        lat, lon = gv.dlatlon[domain, gName]
        llat.append(lat)
        llon.append(lon)

    latmax = max(llat)
    latmin = min(llat)
    lonmax = max(llon)
    lonmin = min(llon) 

    print domain, '%.2f %.2f %.2f %.2f'%(latmin,latmax,lonmin,lonmax), dnumGauge[domain], dnumObsSate[domain]


