from numpy import *
import myfunc.util as util
import numpy as np
import os,sys

iYM = [2017,1]
eYM = [2017,12]
lYM = util.ret_lYM(iYM,eYM)
#varName = 'DPRGMI_NS_surfPrecipTotRate'  # Any 1-dimensional variable
#baseDir = '/work/hk01/utsumi/PMM/EPCDB/GMI.V05A.S1.ABp103-117/%s'%(varName)
baseDir = '/work/hk01/utsumi/PMM/EPCDB/GMI.V05A.S1.ABp103-117/dprx'

listDir = '/work/hk01/utsumi/PMM/EPCDB/list'

#lepcid  = range(0,25*25*25) # 25*25*25=15625
lepcid  = range(0,29*29*29) # 29*29*29 = 24389

for (Year,Mon) in lYM:
    lout   = []
    #for epcid in lepcid:
    for epcid in lepcid[1:]:  # test
        srcDir  = baseDir + '/%04d%02d'%(Year,Mon)
        srcPath = srcDir + '/%s.%05d.npy'%('dprx', epcid)
        #print Year,Mon,epcid
        if os.path.exists(srcPath):
            adprx    = np.load(srcPath)
            nrec    = len(adprx)
        else:
            continue

        #amask = ma.masked_equal(adprx, 24)  # Nadir(Xpy=24)
        amask = ma.masked_inside(adprx, 23,25)  # Nadir(Xpy=24) and 23, 26
        amask = ma.masked_outside(amask, 21,27)  # -, 21, 22, 23, -, 25, 26, 27, -
        if amask.count()==0: continue
        airec = np.arange(adprx.shape[0]).astype(int32)
        airec = ma.masked_where(amask.mask, airec).compressed()

        print epcid, amask.shape[0],amask.count()
        if amask.count()>100:sys.exit()
        lout.append([epcid,nrec])
    ##-- Output data --
    #sout = util.list2csv(lout)
    #
    ##-- Save to file --
    #outPath = listDir + '/nrec.%04d%02d.csv'%(Year,Mon)
    #f=open(outPath,'w'); f.write(sout); f.close()
    #print outPath
