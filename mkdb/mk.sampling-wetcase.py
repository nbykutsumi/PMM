from numpy import *
import myfunc.util as util
from datetime import datetime, timedelta
import random
import numpy as np
import os, sys
import socket

hostname = socket.gethostname()
if hostname == 'shui':
    tankbaseDir = '/tank'
    workbaseDir = '/work'
    figDir = '/home/utsumi/temp/stop'
elif hostname == 'well':
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    workbaseDir = '/home/utsumi/mnt/lab_work'
    figDir = '/home/utsumi/temp/stop'

else:
    print 'check hostname',hostname
    sys.exit()

#lvarName = ['tqv']
#lvarName = ['t2m']
#lvarName = ['epc']
#lvarName = ['gtopo']
#lvarName = ['Ku_NS_heightStormTop']
#lvarName = ['epc','Ku_NS_heightStormTop']
#lvarName = ['surfaceTypeIndex','tqv','t2m','gtopo']
lvarName = ['flagInvalidTc']


NREC  = 5000
ibaseDir = tankbaseDir + '/utsumi/PMM/EPCDB/GMI.V05A.S1.ABp103-117'
obaseDir = tankbaseDir + '/utsumi/PMM/EPCDB/wetcase.samp.%d.GMI.V05A.S1.ABp103-117'%(NREC)
lepcid = range(0,29*29*29)
#lepcid = range(8630,29*29*29)
tcmin, tcmax = 50, 350

lYear = [2017]
lMon  = range(1,12+1)
lYM   = [[Year,Mon] for Year in lYear for Mon in lMon]
iMon,eMon = lMon[0],lMon[-1]

#-- Read nrec files --
drec = {}
for (Year,Mon) in lYM:
    srcPath = tankbaseDir + '/utsumi/PMM/EPCDB/list/nrec-wetcase.%04d%02d.csv'%(Year,Mon)
    f=open(srcPath,'r'); lines=f.readlines(); f.close()
    lnrecTmp = []
    for line in lines:
        line  = line.strip().split(',')
        nrec  = int(line[1])
        lnrecTmp.append(nrec)

    drec[Year,Mon] = array(lnrecTmp)

drec[0] = zeros(len(drec[lYear[0],lMon[0]]),int32)
for (Year,Mon) in lYM:
    drec[0] = drec[0] + drec[Year,Mon]


#------ Sampling of entries ---------
for varName in lvarName:
    if varName == 'flagInvalidTc':
        varNameTmp = 'Tc'
    else:
        varnameTmp
    #for epcid in [4150]:
    for epcid in lepcid:
        #-- Sample from monthly data ---
        astack = []
        for (Year,Mon) in lYM:
            nrec = drec[Year,Mon][epcid]

            #-- Calc nuse for each month --
            if drec[0][epcid] >= NREC:
                nuse = int(NREC * nrec / drec[0][epcid])
            else:
                nuse = nrec

            print drec[0][epcid],nrec, nuse
            #-- Read data ---
            srcDir = ibaseDir + '/%s/%04d%02d'%(varNameTmp, Year,Mon)
            srcPath= srcDir + '/%s.%05d.npy'%(varNameTmp,epcid)

            prcDir = ibaseDir + '/DPRGMI_NS_surfPrecipTotRate/%04d%02d'%(Year,Mon)
            prcPath= prcDir + '/DPRGMI_NS_surfPrecipTotRate.%05d.npy'%(epcid)

            #print Year,Mon,srcPath
            if os.path.exists(srcPath):
                aTmp   = np.load(srcPath)
                aPrc   = np.load(prcPath)
            else:
                #aTmp   = array([]).astype(dtype=dtype)
                continue

            #-- flagInvalidTc --
            if varName == 'flagInvalidTc':
                nrec = aTmp.shape[0]
                aTmp = ma.masked_outside(aTmp, tcmin, tcmax).mask

                if type(aTmp) is np.bool_:
                    aTmp = np.array([aTmp]*nrec)
                else:
                    aTmp = aTmp.any(axis=1)

                aTmp = aTmp.astype('int16')

                #if aTmp.sum()>0:
                #    if epcid >0:
                #        sys.exit()
            #-- Wet cases ---
            amaskP1 = ma.masked_invalid(aPrc).mask
            amaskP2 = ma.masked_less_equal(aPrc, 0).mask
            amask   = amaskP1 + amaskP2
            aidx = np.arange(aPrc.shape[0]).astype('int32')
            aidx = ma.masked_where(amask, aidx).compressed()
            aidx = aidx[:nuse]

            #print Year,Mon,aTmp.shape
            aTmp = aTmp[:nuse]
            astack.append(aTmp)


        if len(astack)==0: continue
        astack = concatenate(astack, axis=0)
        ##-- Random sampling from stacked data --
        #random.seed(epcid,0)
        #ntot = len(astack)
        #aidx = random.sample(range(ntot), k=ntot)
        #astack = astack[aidx].astype(dtype=dtype)

        #-- Save data ---
        outDir  = obaseDir + '/%s'%(varName)
        util.mk_dir(outDir)
        outPath = outDir + '/%05d.npy'%(epcid)
        np.save(outPath, astack)
        print outPath
        #print astack.shape

