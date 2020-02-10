import numpy as np
from numpy import *
import myfunc.util as util
import sys, os
import socket
import glob

myhost = socket.gethostname()
if myhost =='shui':
    tankbaseDir = '/tank'
elif myhost == 'well':
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
else:
    print 'check myhost'
    sys.exit()

nz = 60
iwatDir    = tankbaseDir + '/utsumi/PMM/EPCDB/samp.10000.GMI.V05A.S1.ABp103-117.01-12/DPRGMI_NS_precipTotWaterCont'
surfbinDir = tankbaseDir + '/utsumi/PMM/EPCDB/samp.10000.GMI.V05A.S1.ABp103-117.01-12/DPRGMI_NS_surfaceRangeBin'
clutbinDir = tankbaseDir + '/utsumi/PMM/EPCDB/samp.10000.GMI.V05A.S1.ABp103-117.01-12/DPRGMI_NS_lowestClutterFreeBin'

outDir    = tankbaseDir + '/utsumi/PMM/EPCDB/samp.10000.GMI.V05A.S1.ABp103-117.01-12/DPRGMI_NS_precipTotWaterContRelSurf'
util.mk_dir(outDir)


ssearch = iwatDir + '/?????.npy'
liwatPath = np.sort(glob.glob(ssearch))
for iwatPath in liwatPath:
    #print iwatPath
    idx_db = int(iwatPath.split('/')[-1].split('.')[-2])

    surfbinPath = surfbinDir + '/%05d.npy'%(idx_db)
    clutbinPath = clutbinDir + '/%05d.npy'%(idx_db)

    a2iwat = np.load(iwatPath)   # top to bottom, 250m resol, 60-ranges
    a1surfbin = np.load(surfbinPath)  # surface bin, out of 88 (250m) from top to bottom (1=top, 88=bottom)
    a1clutbin = np.load(clutbinPath)  # clutter bin, out of 88 (250m) from top to bottom (1=top, 88=bottom)

    #-- Replace -9999 to 88 --
    a1surfbin = ma.masked_less(a1surfbin,0).filled(88)
    a1clutbin = ma.masked_less(a1clutbin,0).filled(88)


    #-- convert bins ---
    a1surfbin = a1surfbin -29  # 88-->None, ... 28-->None, 29-->0, 30-->1, ... 88-->59 (python indexing)
    a1clutbin = a1clutbin -29  # 88-->None, ... 28-->None, 29-->0, 30-->1, ... 88-->59 (python indexing)



    nl  = a1surfbin.shape[0]
    a1y = np.arange(nl).astype(int32)

    #-- Fill clutter range --
    a1valid = a2iwat[a1y, a1clutbin].copy()
    setclutbin = list(set(a1clutbin))

    for clutbin in setclutbin:
       
        a1yTmp = ma.masked_where(a1clutbin !=clutbin,  a1y).compressed() 
        a2iwat[a1yTmp, clutbin:] = a1valid[a1yTmp].reshape(-1,1)

    #-- Make double-size array ---

    a2iwat = np.concatenate([np.zeros([nl, nz], float32), a2iwat], axis=1)
    a1surfbin = a1surfbin + nz



    #-- Relative to surface ---------
    a2owat = np.zeros([nl,nz],float32)

    for i in range(nz):
        a1k = a1surfbin - nz+1 + i  # i=0-59
        a2owat[:,i] = a2iwat[a1y,a1k]

    #if (nl >1)and(idx_db !=0):
    #    print a2iwat
    #    print a1surfbin
    #    print a1clutbin
    #    print a2iwat.shape
    #    print a2owat
    #    print a2owat.shape
    #    sys.exit()

    #-- Save -----------
    outPath = outDir + '/%05d.npy'%(idx_db)
    np.save(outPath, a2owat.astype('float32'))
    print outPath


