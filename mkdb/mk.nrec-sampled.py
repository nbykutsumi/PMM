from numpy import *
import myfunc.util as util
import numpy as np
import os,sys

iYM = [2017,1]
eYM = [2017,12]
lYM = util.ret_lYM(iYM,eYM)
varName = 'Ku_NS_heightStormTop' 
srcDir = '/tank/utsumi/PMM/EPCDB/wetcase.samp.5000.GMI.V05A.S1.ABp103-117/%s'%(varName)

listDir = '/tank/utsumi/PMM/EPCDB/wetcase.samp.5000.GMI.V05A.S1.ABp103-117/list'

#lepcid  = range(0,25*25*25) # 25*25*25=15625
lepcid  = range(0,29*29*29) # 29*29*29 = 24389
#lepcid = range(0,9000)

lout   = []
for epcid in lepcid:
    print epcid
    srcPath = srcDir + '/%05d.npy'%(epcid)
    if os.path.exists(srcPath):
        avar = np.load(srcPath)
        nrec = avar.shape[0]
    else:
        nrec = 0

    lout.append([epcid,nrec])
#-- Output data --
sout = util.list2csv(lout)

#-- Save to file --
util.mk_dir(listDir)
outPath = listDir + '/nrec.csv'
f=open(outPath,'w'); f.write(sout); f.close()
print outPath
