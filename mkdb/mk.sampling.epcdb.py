from numpy import *
import myfunc.util as util
from datetime import datetime, timedelta
import random
import numpy as np
import os, sys

lYear = [2017]
lMon  = range(1,12+1)
lYM   = [[Year,Mon] for Year in lYear for Mon in lMon]
iMon,eMon = lMon[0],lMon[-1]

#lvarName = ['epc','Ka_MS_zFactorMeasured','Tc','gNum','t2m']
#lvarName = ['DPRGMI_NS_precipTotWaterCont','DPRGMI_MS_surfPrecipTotRate','DPRGMI_NS_surfPrecipTotRate','Ku_NS_precipRateNearSurface']
#lvarName = ['Ka_MS_precipRateNearSurface']
lvarName = ['gtopo']
NREC  = 20000
#NREC  = 5000
ibaseDir = '/work/hk01/utsumi/PMM/EPCDB/GMI.V05A.S1.ABp103-117'
obaseDir = '/work/hk01/utsumi/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.%02d-%02d'%(NREC,iMon,eMon)
lepcid = range(0,29*29*29)


'''
int8 : -128 ~ +127
int16: -32768 ~ +32767
int32: -2147483648 ~ +2147483647
'''



#-- Read nrec files --
drec = {}
for (Year,Mon) in lYM:
    srcPath = '/work/hk01/utsumi/PMM/EPCDB/list/nrec.%04d%02d.csv'%(Year,Mon)
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

#------ Ramdom sampling of entries ---------
for varName in lvarName:
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
    
            #-- Read epcdb --
            srcDir = ibaseDir + '/%s/%04d%02d'%(varName, Year,Mon)
            srcPath= srcDir + '/%s.%05d.npy'%(varName,epcid) 

            #print Year,Mon,srcPath
            if os.path.exists(srcPath):
                aTmp   = np.load(srcPath)
            else:
                #aTmp   = array([]).astype(dtype=dtype)
                continue 

            #print Year,Mon,aTmp.shape
            #-- Sample entries --
            #random.seed(epcid,Year,Mon)
            #aidx = random.sample(range(nrec),k=nuse)
            #aTmp = aTmp[aidx] 
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
