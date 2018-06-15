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
gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()

dgName = gv.ret_ddomYM2gName()

ldomain = gv.domains
#ldomain = ['FLORIDA-STJ','FLORIDA-SFL-N','N.Carolina-IPHEx_Duke','N.Carolina-IPHEx_NASA','KWAJALEIN-KWA']
#ldomain = ['FLORIDA-STJ']

offset_bef = 15  # 'bef' should be identical to what is used in mk_match.py
#offset_aft = 45

ldt  = [1, 5,10,15,20,25,30,35]
ldh  = [1,3,5,7,9,11,13,15,17,19,21,23,25]

nt  = len(ldt)
nh  = len(ldh)

#lprtype = ['all','mod','heavy','extreme']
lprtype = ['all','mod','heavy']
#lprtype = ['mod']
dlthpr = {'all':[-0.1,9999],'mod':[-0.1,10], 'heavy':[10,50],'extreme':[50,9999]}
ldattype = ['rain','cc','bias','brat','rmse','num','gv']
#ldattype = ['rmse']

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

da2rain = {prtype: empty([nh,nt]) for prtype in lprtype}
da2gv   = {prtype: empty([nh,nt]) for prtype in lprtype}
da2cc   = {prtype: empty([nh,nt]) for prtype in lprtype}
da2rmse = {prtype: empty([nh,nt]) for prtype in lprtype}
da2bias = {prtype: empty([nh,nt]) for prtype in lprtype}
da2brat = {prtype: empty([nh,nt]) for prtype in lprtype}
da2num  = {prtype: empty([nh,nt],int32) for prtype in lprtype}

for idt, dt in enumerate(ldt):
    if calc != True:
        continue

    daprof     = {ih:deque([]) for ih in range(nh)}
    daesurf    = {ih:deque([]) for ih in range(nh)}
    dagv       = {ih:deque([]) for ih in range(nh)}

    for domain in ldomain:
        for YM in lYM:
            Year, Mon = YM
            if (domain,Year,Mon) not in dgName.keys():
                print 'no obs',domain,Year,Mon
                continue
            print domain, 'dt=',dt,YM
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
            a2gv   = abs(np.load(gvPath))
            aidxtmp= zeros(a2gv.shape[0])
            agv    = gv_fsub.mean_slice_negativemask(a2gv.T, aidxtmp, dt)
            acoef  = ma.masked_invalid(aesurf/aprof[:,0]).filled(coef_def)

            # accumulation arrays for mask
            asateAcc = gv_fsub.mean_slice_negativemask(aprof.T, agroundBin, ldh[-1])
            agvAcc   = gv_fsub.mean_slice_negativemask(a2gv.T, aidxtmp, ldt[-1])

            #-- mask when ngv < minNum
            amskN    = ma.masked_less(angv, minNum).mask

            #-- mask when both satellite and gv are zero
            amsk1    = ma.masked_equal(asateAcc,0).mask
            amsk2    = ma.masked_equal(agvAcc,0).mask
            amskzero = amsk1 * amsk2

            #-- mask when at least one of the data is missing or low quality
            amsk3    = ma.masked_less(asateAcc,0).mask
            amsk4    = ma.masked_less(agvAcc,0).mask
            amskmiss = amsk3 + amsk4

            #-- mask when sate has too large value
            amsk5    = ma.masked_greater(asateAcc*acoef, 100).mask
            amsk6    = ma.masked_greater(agvAcc, 100).mask
            amskL    = amsk5 + amsk6

            #-- overlay masks
            amsk     = amskN + amskzero + amskmiss + amskL



            for ih,dh in enumerate(ldh):
                asate    = gv_fsub.mean_slice_negativemask(aprof.T, agroundBin, dh)

                aidxTmp  = arange(len(amsk)).astype(int32)
                aidxTmp  = ma.masked_where(amsk, aidxTmp).compressed()

                agroundBinTmp= agroundBin[aidxTmp]        
                asateTmp     = asate[aidxTmp] 
                aesurfTmp    = aesurf[aidxTmp]
                agvTmp       = agv[aidxTmp]
                acoefTmp     = acoef[aidxTmp]
                asate_mmh    = asateTmp * acoefTmp # g/m3 --> mm/h
                #------------------
                daprof[ih].extend(asate_mmh)  
                daesurf[ih].extend(aesurfTmp)
                dagv[ih].extend(agvTmp)

                #if (dt==1)&(dh==25):
                #    aprofTmp = aprof[aidxTmp]

                #    aprofNumTmp   = aprofNum[aidxTmp]
                #    aprofScaleTmp = aprofScale[aidxTmp]
                #    atIndexTmp    = atIndex[aidxTmp]

                #    for i in range(aprofTmp.shape[0]):
                #        coef = acoefTmp[i]

                #        if asate_mmh[i]>100:
                #            print ''
                #            print '-'*50
                #            print domain, i, asate_mmh[i]
                #            print 'profOrf'
                #            print aprofTmp[i]
                #            print 'prof_mmh'
                #            print aprofTmp[i]*coef
                #            print aprofNumTmp.shape, aprofScaleTmp.shape,atIndexTmp.shape
                #            for isp in [0,2,3,4]:
                #                profNum   = aprofNumTmp[i,isp]
                #                profScale = aprofScaleTmp[i,isp]
                #                tIndex    = atIndexTmp[i]

                #                print ''
                #                print 'isp=',isp,'profNum=profNum','tIndex=',tIndex, 'scale=',profScale
                #                print a4clusterProf[profNum-1,:,tIndex-1,isp]

                #            print  
                #            print ''
                #            sys.exit()

    #----
    for ih in range(nh): 
        aprof = array(daprof[ih])  # mm/h
        aesurf= array(daesurf[ih])
        agv   = array(dagv[ih])


        #- profile vs gv --
        for prtype in lprtype:
            # mask negative values
            amsk1 = ma.masked_less(aprof, 0).mask
            amsk2 = ma.masked_less(agv, 0).mask
            amskNega = amsk1 + amsk2
            # mask with nSurf --
            thmin, thmax = dlthpr[prtype]
            #amskP      = ma.masked_outside(ansurf, thmin, thmax).mask
            amskP      = ma.masked_outside(agv, thmin, thmax).mask

            amsk = amskNega + amskP

            aprofTmp  = ma.masked_where(amsk, aprof ).compressed()
            agvTmp    = ma.masked_where(amsk, agv   ).compressed()
    
            # mask negative values 
            
            cc     = np.corrcoef(aprofTmp, agvTmp)[0,1]
            rmse   = np.sqrt( ((aprofTmp- agvTmp)**2).mean() )
           
            bias   = (aprofTmp - agvTmp).mean()
            brat   = (aprofTmp - agvTmp).mean() / agvTmp.mean()*100
            num    = len(aprofTmp)
            rain   = aprofTmp.mean()
            gv     = agvTmp.mean()

            da2cc  [prtype][ih, idt] = cc
            da2rmse[prtype][ih, idt] = rmse
            da2bias[prtype][ih, idt] = bias
            da2brat[prtype][ih, idt] = brat
            da2num [prtype][ih, idt] = num       
            da2rain[prtype][ih, idt] = rain
            da2gv  [prtype][ih, idt] = gv


#--- save data ----
for prtype in lprtype:
    if calc != True:
        continue

    outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-%s/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])
    util.mk_dir(outDir)

    for dattype in ldattype:  
        datPath= outDir + '/nt-nlev.%s.%s.npy'%(dattype, prtype)

        if dattype=='cc':
            a2dat = da2cc[prtype]
        elif dattype=='rmse':
            a2dat = da2rmse[prtype]
        elif dattype=='bias':
            a2dat = da2bias[prtype]
        elif dattype=='brat':
            a2dat = da2brat[prtype]
        elif dattype=='num':
            a2dat = da2num[prtype]
        elif dattype=='rain':
            a2dat = da2rain[prtype]
        elif dattype=='gv':
            a2dat = da2gv[prtype]


        np.save(datPath, a2dat)
        print datPath


#--- figure -------
for prtype in lprtype:
    outDir = '/home/utsumi/mnt/wellshare/GPMGV/dt-lev-%s/dist.%.1fkm.ndom.%02d.%04d.%02d-%04d.%02d'%(prdName, thdist, len(ldomain), iYM[0], iYM[1], eYM[0], eYM[1])

    for dattype in ldattype:

        datPath=outDir + '/nt-nlev.%s.%s.npy'%(dattype, prtype)
        a2dat  = np.load(datPath)
    
        fig = plt.figure(figsize=(8,4))
        ax  = fig.add_axes([0.15, 0.1, 0.8, 0.8])
        a1t  = ldt

        wbar = 0.7/nh
        for ih, dh in enumerate(ldh):
            a1y = a2dat[ih,:]
            a1y = ma.masked_invalid(a1y)
            a1x = arange(nt) -0.35 +ih*wbar
            ax.bar(a1x, a1y, width=wbar, tick_label=ldt, align='center', label='%d bins'%(dh))
            #if ih !=0:
            #    ax.tick_params(labelbottom='off') 

        # legend
        plt.legend()
        # ylim
        if dattype in ['cc']:
            ymin = -1.2
            ymax = 1.2
            plt.ylim([ymin,ymax])

        elif dattype in ['rain','gv','rmse','num']:
            ymin = 0
            ymax = ma.masked_invalid(a2dat).max()*1.1
            plt.ylim([ymin,ymax])

        elif dattype in ['bias','brat']:
            ymax = abs(ma.masked_invalid(a2dat)).max()*1.1
            ymin = -ymax
            plt.ylim([ymin,ymax])


        #- zero line ---
        if dattype in ['cc','bias','brat']:
            plt.plot([-10,50],[0,0],'--',color='k', linewidth=0.5)

        plt.xlim([a1x[0]-1,a1x[-1]+2])


        stitle  = '%s %s %s'%(prdName, dattype, prtype)
        plt.title(stitle)        
        figDir  = '/work/a01/utsumi/GPMGV/fig'
        figPath = figDir + '/bar.nt-nlev.%s.%.1fkm.minNum.%d.%s.%s.png'%(prdName, thdist,minNum,dattype, prtype)
        plt.savefig(figPath)
        print figPath
        plt.clf()



