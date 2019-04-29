import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import myfunc.util as util
#import myfunc.IO.GPM.l2_dpr as l2_dpr
import h5py
import myfunc.IO.GPM.l1_gmi as l1_gmi
import glob
from datetime import datetime, timedelta
import numpy as np
import sys
from f_match_fov import *


gmi  = l1_gmi.L1_GMI()
#dpr  = l2_dpr.L2_DPR()
mwscan= 'S1'
radar = 'Ku'

#iDTime = datetime(2017,9,30)
#eDTime = datetime(2018,1,1)
iDTime = datetime(2017,1,1)
eDTime = datetime(2017,1,1)


#iDTime = datetime(2017,6,30)
#eDTime = datetime(2017,8,1)
#iDTime = datetime(2017,2,1)
#eDTime = datetime(2017,6,1)
#iDTime = datetime(2017,8,1)
#eDTime = datetime(2018,1,1)

#-- exclude missing files --
lDTimeTmp = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
lDTime = []
for DTime in lDTimeTmp:
    if ((datetime(2017,9,26)<=DTime) & (DTime <= datetime(2017,9,29))): continue

    lDTime.append(DTime)
#---------------------------

ix0 = 83   # in python indexing. GMI angle bins= 0, 1, 2, ..., 220 : in total=221
ex0 = 137  # in python indexing. GMI angle bins= 0, 1, 2, ..., 220 : in total=221
cx  = 110  # GMI center angle bin (py-idx)
cw  = 15    # extract this width around center
w   = int(cw/2)


verGMI = '05'
subverGMI = 'A'
fullverGMI = '%s%s'%(verGMI,subverGMI)

verDPR = '06'
subverDPR = 'A'
fullverDPR = '%s%s'%(verDPR,subverDPR)

baseDirGMI = '/work/hk01/PMM/NASA/GPM.GMI/1C/V%s'%(verGMI)
baseDirDPR = '/work/hk01/PMM/NASA/GPM.Ku/2A/V%s'%(verDPR)
idxbaseDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V%s/%s.ABp%03d-%03d.%s.V%s.IDX'%(fullverGMI, mwscan, ix0, ex0, radar, fullverDPR)

outrootDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V%s'%(fullverGMI)

lvar = ['/NS/SLV/precipRate']
#lvar = ['/NS/CSF/typePrecip']
#lvar = ['/NS/CSF/typePrecip','NS/PRE/heightStormTop','NS/CSF/flagAnvil','/NS/SLV/precipRate']
#lvar = ['/NS/CSF/typePrecip','NS/PRE/heightStormTop','NS/CSF/flagAnvil','/NS/VER/heightZeroDeg']
#lvar = ['/NS/Latitude','/NS/Longitude']
#lvar = ['/NS/VER/heightZeroDeg']


for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]

    srcDirDPR   = baseDirDPR + '/%04d/%02d/%02d'%(Year,Mon,Day)
    ssearch     = srcDirDPR  + '/2A.GPM.%s.*.V%s.HDF5'%(radar,fullverDPR)
    lsrcPathDPR = sort(glob.glob(ssearch))


    if len(lsrcPathDPR)==0:
        print 'No DPR file',Year,Mon,Day
        print ssearchDPR
        sys.exit()

    for srcPathDPR in lsrcPathDPR:
        oid = srcPathDPR.split('.')[-3]


        for irank in [1,2,3,4]:
            idxDir      = idxbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
            idxPathX    = idxDir + '/Xpy.%d.%s.npy'%(irank,oid)
            idxPathY    = idxDir + '/Ypy.%d.%s.npy'%(irank,oid)
    
            X    = np.load(idxPathX)
            Y    = np.load(idxPathY)
    
            a2x  = X[:,cx-w-ix0:cx+w+1-ix0]
            a2y  = Y[:,cx-w-ix0:cx+w+1-ix0]

            for var in lvar:
                #DatDPR = dpr.load_var_granule(srcPathDPR, var)
                with h5py.File(srcPathDPR) as hdpr:
                    DatDPR = hdpr[var][:]
 
                if   len(DatDPR.shape)==2:
                    datout = f_match_fov.extract_2d(DatDPR.T, a2x.T, a2y.T, -9999, -9999.).T
                elif len(DatDPR.shape)==3:
                    datout = f_match_fov.extract_3d(DatDPR.T, a2x.T, a2y.T, -9999, -9999.).T

                if var.split('/')[-1] in ['precipRate']:
                    datatype = 'int16'
                    miss     = -9999
                    datout   = (ma.masked_less(datout,0)*100).astype(datatype)
                    datout   = ma.masked_less(datout,0).filled(miss)
    
                else:
                    datatype   = DatDPR.dtype
                    datout     = datout.astype(datatype)


                varName    = var.split('/')[-1] 
                outbaseDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V%s/%s.ABp%03d-%03d.%s.V%s.%s'%(fullverGMI, mwscan, cx-w, cx+w, radar, fullverDPR, varName)
                outDir     = outbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
                outPath    = outDir + '/%s.%d.%s.npy'%(varName, irank, oid)


                util.mk_dir(outDir)
                np.save(outPath, datout)
                print outPath
