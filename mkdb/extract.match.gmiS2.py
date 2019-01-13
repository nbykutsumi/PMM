import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import myfunc.util as util
import myfunc.IO.GPM.l2_dpr as l2_dpr
import myfunc.IO.GPM.l1_gmi as l1_gmi
import glob
from datetime import datetime, timedelta
import numpy as np
import sys
from f_match_fov import *


gmi  = l1_gmi.L1_GMI()
dpr  = l2_dpr.L2_DPR()
mwscan= 'S1'
radar = 'Ku'

iDTime = datetime(2016,12,31)
eDTime = datetime(2017,1,5)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
lDTime = [DTime for DTime in lDTime if not (datetime(2017,9,26)<=DTime)&(datetime(2017,9,29)]

ix0 = 83   # in python indexing. GMI angle bins= 0, 1, 2, ..., 220 : in total=221
ex0 = 137  # in python indexing. GMI angle bins= 0, 1, 2, ..., 220 : in total=221
cx  = 110  # GMI center angle bin (py-idx)
cw  = 15    # extract this width around center
w   = int(cw/2)


verGMI = '05'
subverGMI = 'A'
fullverGMI = '%s%s'%(verGMI,subverGMI)

baseDirGMI = '/work/hk01/PMM/NASA/GPM.GMI/1C/V%s'%(verGMI)
idxbaseDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V%s/%s.ABp%03d-%03d.GMI.S2idx'%(fullverGMI, mwscan, ix0, ex0)

outrootDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V%s'%(fullverGMI)

for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]

    baseDirGMI = '/work/hk01/PMM/NASA/GPM.GMI/1C/V%s'%(verGMI)
    srcDirGMI  = baseDirGMI+ '/%04d/%02d/%02d'%(Year,Mon,Day)
    ssearchGMI = srcDirGMI + '/1C.GPM.GMI.*.HDF5'

    lsrcPathGMI = glob.glob(ssearchGMI)

    if len(lsrcPathGMI)==0:
        print 'No GMI file',Year,Mon,Day
        print ssearchGMI
        sys.exit()

    for srcPathGMI in lsrcPathGMI:
        oid = srcPathGMI.split('.')[-3]

        irank       = 1
        idxDir      = idxbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
        idxPathX    = idxDir + '/Xpy.%d.%s.npy'%(irank,oid)
        idxPathY    = idxDir + '/Ypy.%d.%s.npy'%(irank,oid)
    
        X    = np.load(idxPathX)
        Y    = np.load(idxPathY)
    
        a2x  = X[:,cx-w-ix0:cx+w+1-ix0]
        a2y  = Y[:,cx-w-ix0:cx+w+1-ix0]

        DatS2 = gmi.load_var_granule(srcPathGMI, 'S2/Tc')

        datout = f_match_fov.extract_3d(DatS2.T, a2x.T, a2y.T, -9999, -9999.).T

        datatype   = DatS2.dtype
        datout     = datout.astype(datatype)

        outbaseDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V%s/%s.ABp%03d-%03d.GMI.TcS2'%(fullverGMI, mwscan, cx-w, cx+w)
        outDir     = outbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
        outPath    = outDir + '/%s.%d.%s.npy'%('TcS2', irank, oid)

        util.mk_dir(outDir)
        np.save(outPath, datout)
        print outPath
