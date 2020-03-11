from numpy import *
import numpy as np
import os, sys
import socket
import myfunc.util as util

hostname = socket.gethostname()
if hostname == 'shui':
    tankbaseDir = '/tank'
    stopbaseDir= '/tank/utsumi/PMM/stop'
    figDir = '/home/utsumi/temp/stop'
elif hostname == 'well':
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    stopbaseDir= '/home/utsumi/mnt/lab_tank/utsumi/PMM/stop'
    figDir = '/home/utsumi/temp/stop'

else:
    print 'check hostname',hostname
    sys.exit()


#lvarName = ['tqv']
#lvarName = ['t2m']
#lvarName = ['epc']
#lvarName = ['gtopo']
lvarName = ['flagInvalidTc']
#lvarName = ['Ku_NS_heightStormTop']
#lvarName = ['epc','Ku_NS_heightStormTop']
#lvarName = ['surfaceTypeIndex','tqv','t2m']


baseDir = tankbaseDir + '/utsumi/PMM/EPCDB/wetcase.samp.5000.GMI.V05A.S1.ABp103-117'
lidx_db = np.arange(29*29*29)
#lidx_db = np.arange(3000)
#************************
# Read nrec file
#------------------------
nrecDir = tankbaseDir + '/utsumi/PMM/EPCDB/wetcase.samp.5000.GMI.V05A.S1.ABp103-117/list'
nrecPath= nrecDir + '/nrec.csv'

f=open(nrecPath,'r'); lines=f.readlines(); f.close()
a1nrec = []
for line in lines:
    a1nrec.append( int(line.strip().split(',')[1]) )
a1nrec = np.array(a1nrec)[:len(lidx_db)]
nrecAll = int(a1nrec.sum())
#************************
for varName in lvarName:

    avar = None
    for idx_db in lidx_db:
        if a1nrec[idx_db]==0: continue
        print idx_db, varName

        srcPath = baseDir + '/%s/%05d.npy'%(varName,idx_db)
        avarTmp = np.load(srcPath)

        if avar is None: 
            avar = avarTmp
        else:
            avar = np.concatenate([avar, avarTmp], axis=0)


    aidx = np.arange(nrecAll)
    np.random.seed(0)
    np.random.shuffle(aidx)

    nbatch  = 10
    nrecTmp = int(nrecAll/nbatch)+1
    for i in range(nbatch):
        aidxTmp = aidx[i*nrecTmp:(i+1)*nrecTmp]

        avarTmp = avar[aidxTmp]

        outDir = baseDir + '/batch.%s'%(varName)
        util.mk_dir(outDir)
        outPath= outDir  + '/%02d.npy'%(i)
   
        np.save(outPath, avarTmp)
        print outPath
        print nrecAll, avarTmp.shape
#************************


