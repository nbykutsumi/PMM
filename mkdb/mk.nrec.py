from numpy import *
import myfunc.util as util
import numpy as np
import os,sys

iYM = [2017,2]
eYM = [2017,6]
lYM = util.ret_lYM(iYM,eYM)
varName = 'DPRGMI_NS_surfPrecipTotRate'  # Any 1-dimensional variable
baseDir = '/work/hk01/utsumi/PMM/EPCDB/GMI.V05A.S1.ABp103-117/%s'%(varName)

listDir = '/work/hk01/utsumi/PMM/EPCDB/list'

#lepcid  = range(0,25*25*25) # 25*25*25=15625
lepcid  = range(0,29*29*29) # 29*29*29 = 24389

for (Year,Mon) in lYM:
    lout   = []
    for epcid in lepcid:
        srcDir  = baseDir + '/%04d%02d'%(Year,Mon)
        srcPath = srcDir + '/%s.%05d.npy'%(varName, epcid)
        print Year,Mon,epcid
        if os.path.exists(srcPath):
            avar    = np.load(srcPath)
            nrec    = len(avar)
        else:
            nrec    = 0

        lout.append([epcid,nrec])
    #-- Output data --
    sout = util.list2csv(lout)
    
    #-- Save to file --
    outPath = listDir + '/nrec.%04d%02d.csv'%(Year,Mon)
    f=open(outPath,'w'); f.write(sout); f.close()
    print outPath
