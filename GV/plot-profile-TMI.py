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

calc = True
#calc = False
iYM = [2010,4]
eYM = [2014,10]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
print lYM

thdist = 7.5 # km
#thdist = 15 # km
minNum = 3
prdName = '2A-CLIM'
nh = 36
gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()

dgName = gv.ret_ddomYM2gName()

ldomain = gv.domains
#ldomain = ['FLORIDA-STJ','FLORIDA-SFL-N','N.Carolina-IPHEx_Duke','N.Carolina-IPHEx_NASA','KWAJALEIN-KWA']
#ldomain = ['FLORIDA-STJ']

offset_bef = 15  # 'bef' should be identical to what is used in mk_match.py
#offset_aft = 45

#lprtype = ['all','mod','heavy','extreme']
lprtype = ['mod','heavy']
#lprtype = ['mod']
dlthpr = {'all':[-0.1,9999],'mod':[-0.1,10], 'heavy':[10,50],'extreme':[50,9999]}

miss    = -9999.

coef_def = 4.17 * 1000/(60*60)  # g/m3 * m/s *1000/(60*60) = mm/h
#------------------------------------------------
def ret_aprof(a4clusterProf, a1tIndex, a2profNum, a2profScale, a1groundBin):
    lspecies = [0,2,3,4]
    #lspecies = [0]
    '''
    0 Rain Water Content
    1 Cloud Water Content
    2 Ice Water Content
    3 Snow Water Content
    4 Grauple/Hail Content
    Hydrometeor Unit:[g/m3] see filespec.TRMM.V7.2A12
    '''

    nh_in    = 36  # interpolated profile levels
    nh_out   = nh_in
    a3out = empty([len(lspecies),len(a1tIndex),nh_in]).astype(float32)

    for i,species in enumerate(lspecies):
        a1profNum  = a2profNum[:,species]
        a1profScale= a2profScale[:,species]
        a2prof = a1profScale.reshape(-1,1) * a4clusterProf[a1profNum-1,:,a1tIndex-1,species]
        a3out[i] = a2prof

    a2prof = a3out.sum(axis=0)
    a2out  = gv_fsub.extract_slice_clusterprof(a2prof.T, agroundBin, nh_out).T

    return a2out


#------------------------------------------------
def prof_highlevel_interp(a4clusterProf):
    '''
    default level [km]: [0.5, 1.0, 1.5, ...9.5, 10, 11, 12, ..., 17, 18] # 28-level
    output  level [km]: [0.5, 1.0, 1.5, ...9.5, 10, 10.5, 11, 11.5, ..., 17.5, 18] # 36-level
    simply duplicate higher levels
    '''
    a4clusterProfOut = empty([80,36,12,5],dtype='float32')
    for ispec in range(80):
        for iout in range(0,19+1):
            iin = iout
            a4clusterProfOut[ispec,iout,:,:] = a4clusterProf[ispec,iin,:,:]

        for iout in range(20,35+1):
            iin = 20 + int((iout-20)/2)
            a4clusterProfOut[ispec,iout,:,:] = a4clusterProf[ispec,iin,:,:]
    return a4clusterProfOut 


#-- cluserProfile --

clustbaseDir = "/home/utsumi/mnt/wellshare/GPMGV/%s"%(prdName)
clusterProfPath = clustbaseDir + "/clusterProf.npy"
a4clusterProfOrg = np.load(clusterProfPath)
a4clusterProf   = prof_highlevel_interp(a4clusterProfOrg)
#------------------------------------------------

#da1sx     = {prtype: zeros(36) for prtype in lprtype}
#da1sxx    = {prtype: zeros(36) for prtype in lprtype}
#da1sx_mmh  = {prtype: zeros(36) for prtype in lprtype}
#da1sxx_mmh = {prtype: zeros(36) for prtype in lprtype}

# container
a2prof = deque([])
a1coef = deque([])
a1esurf= deque([])

for domain in ldomain:
    for YM in lYM:
        Year, Mon = YM
        if (domain,Year,Mon) not in dgName.keys():
            print 'no obs',domain,Year,Mon
            continue
        # load data
        srcbaseDir = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.%s'%(prdName)
        srcDir     = srcbaseDir + '/%.1fkm/%s/%04d%02d'%(thdist, domain, Year,Mon)

        gvPath        = srcDir  + '/p_gvprcp.npy'
        qFlagPath     = srcDir  + '/p_qFlag.npy'
        profNumPath   = srcDir  + '/p_profNum.npy'
        profScalePath = srcDir  + '/p_profScale.npy'
        tIndexPath    = srcDir  + '/p_tIndex.npy'

        eSurfPath     = srcDir + '/p_eSurf.npy'
        groundBinPath = srcDir + '/p_groundBin.npy'
        ngvPath       = srcDir + '/p_ngv.npy'



        if not os.path.exists(profNumPath):
            print 'no file',profNumPath
            print 'skip', domain, YM
            continue

        aqFlag     = np.load(qFlagPath)
        aprofNum   = np.load(profNumPath)
        aprofScale = np.load(profScalePath)
        atIndex    = np.load(tIndexPath)
        agroundBin = np.load(groundBinPath)
        angv       = np.load(ngvPath)

        aprof     = ret_aprof(a4clusterProf, atIndex, aprofNum, aprofScale, agroundBin)

        aesurf = np.load(eSurfPath)
        aidxtmp= zeros(len(angv))
        acoef  = ma.masked_invalid(aesurf/aprof[:,0]).filled(coef_def)

        # accumulation arrays for mask
        asateAcc = gv_fsub.mean_slice_negativemask(aprof.T, agroundBin, nh)

        #-- mask when ngv < minNum
        amskN    = ma.masked_less(angv, minNum).mask

        #-- mask when satellite is zero
        amskzero = ma.masked_equal(asateAcc,0).mask

        #-- mask when at least one of the data is missing or low quality
        amskmiss = ma.masked_less(asateAcc,0).mask

        #-- mask when sate has too large value
        amskL    = ma.masked_greater(asateAcc*acoef, 100).mask

        #-- overlay masks
        amsk     = amskN + amskzero + amskmiss + amskL

        aidxTmp  = arange(len(amsk)).astype(int32)
        aidxTmp  = ma.masked_where(amsk, aidxTmp).compressed()

        aprofTmp   = aprof[aidxTmp,:]
        aprofTmp   = ma.masked_less(aprofTmp,0)
        aesurfTmp  = aesurf[aidxTmp]
        acoefTmp   = acoef[aidxTmp]

        a2prof.append(aprofTmp)
        a1coef.extend(acoefTmp)
        a1esurf.extend(aesurfTmp) 

a2prof = concatenate(a2prof, axis=0)
a1coef = array(a1coef)
a1esurf= array(a1esurf)

for dattype in ['hydro','prcp']:
    dmean = {}
    dstd  = {}
    for prtype in lprtype:
        thmin,thmax = dlthpr[prtype]
        a1idx = arange(len(a1coef))
        amskP = ma.masked_outside(a2prof[:,0]*a1coef, thmin, thmax).mask
        a1idx = ma.masked_where(amskP, a1idx).compressed()
        a2profTmp = ma.masked_less(a2prof[a1idx], 0)
        a1coefTmp = a1coef[a1idx]

        if dattype=='hydro':
            a1mean  = a2profTmp.mean(axis=0)
            a1std   = a2profTmp.std(axis=0)

        elif dattype=='prcp':    
            a1mean  = (a2profTmp*a1coefTmp.reshape(-1,1)).mean(axis=0)
            a1std   = (a2profTmp*a1coefTmp.reshape(-1,1)).std(axis=0)

        else:
            print 'check dattype',dattype
            sys.exit()    
       
        dmean[prtype] = a1mean
        dstd [prtype] = a1std

    # figure  -------------
    fig = plt.figure(figsize=(4,3))
    ax  = fig.add_axes([0.1, 0.1, 0.8, 0.7])
  
    a1y = (arange(nh) +1)*0.5

    ax.plot(dmean['mod'],a1y,'-',color='r', label='mod')
    ax.plot(dstd['mod'], a1y,'--',color='r', label='mod')

    ax.plot(dmean['heavy'],a1y,'-',color='b', label='heavy')    
    ax.plot(dstd['heavy'],a1y,'--',color='b', label='heavy')    
    
    # legend
    plt.legend()
    
    stitle  = '%s %s'%(prdName, dattype)
    stitle  = stitle + '\n' + '%04d.%02d-%04d.%02d'%(iYM[0],iYM[1],eYM[0],eYM[1])
    plt.title(stitle)        
    figDir  = '/work/a01/utsumi/GPMGV/fig'
    figPath = figDir + '/plt.prof.%s.%.1fkm.minNum.%d.%s.png'%(prdName, thdist,minNum,dattype)
    plt.savefig(figPath)
    print figPath
    plt.clf()
    

