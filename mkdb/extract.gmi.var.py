import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import myfunc.util as util
import myfunc.IO.GPM.l2_dpr as l2_dpr
import myfunc.IO.GPM.l1_gmi as l1_gmi
import myfunc.IO.GPM.l2a_gprof_hdf5 as l2a_gprof_hdf5
import h5py

import glob
from datetime import datetime, timedelta
import numpy as np
import sys
from f_match_fov import *


gmi    = l1_gmi.L1_GMI()
gprof  = l2a_gprof_hdf5.L2A_GPROF_HDF5()

#iDTime = datetime(2017,6,30)
#eDTime = datetime(2017,8,1)
#iDTime = datetime(2016,12,31)
#eDTime = datetime(2018,1,1)

iDTime = datetime(2014,6,1)
eDTime = datetime(2015,5,31)


lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))

cx  = 110  # GMI center angle bin (py-idx)
cw  = 15    # extract this width around center
w   = int(cw/2)

verGMI = '05'
subverGMI = 'A'
fullverGMI = '%s%s'%(verGMI,subverGMI)
mwscan = 'S1'

outrootDir = '/work/hk02/utsumi/PMM/MATCH.GMI.V%s'%(fullverGMI)

#lvar = [['gmi','S1/Latitude'],['gmi','S1/Longitude'],['gmi','S1/SCstatus/SCorientation'],['gprof','S1/surfaceTypeIndex'],['gprof','S1/surfacePrecipitation'],['gmi','S1/Tc']
#lvar = [['gmi','S1/Latitude'],['gmi','S1/Longitude'],['gprof','S1/surfaceTypeIndex'],['gprof','S1/surfacePrecipitation'],['gmi','S1/Tc']]
lvar = [['gmi','S1/Latitude'],['gmi','S1/Longitude']]
#lvar = [['gprof','S1/surfacePrecipitation']]
#lvar = [['gprof','S1/surfaceTypeIndex']]
#lvar = [['gmi','S1/Tc']]
#lvar = [['gmi','S1/SCstatus/SCorientation']]

for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]

    for prod,var in lvar:

        if   prod=='gmi':
            baseDirGMI = '/work/hk02/PMM/NASA/GPM.GMI/1C/V%s'%(verGMI)
            srcDirGMI  = baseDirGMI+ '/%04d/%02d/%02d'%(Year,Mon,Day)
            ssearchGMI = srcDirGMI + '/1C.GPM.GMI.*.HDF5'

        elif prod=='gprof':
            baseDirGMI = '/work/hk02/PMM/NASA/GPM.GMI/2A/V%s'%(verGMI)
            srcDirGMI  = baseDirGMI+ '/%04d/%02d/%02d'%(Year,Mon,Day)
            ssearchGMI = srcDirGMI + '/2A.GPM.GMI.GPROF*.HDF5'

        else:
            print 'check prod',prod
            sys.exit()

        lsrcPathGMI = glob.glob(ssearchGMI)
    
        if len(lsrcPathGMI)==0:
            print 'No GMI file',Year,Mon,Day
            print ssearchGMI
            #sys.exit()
            continue
    
        for srcPathGMI in lsrcPathGMI:
            oid = srcPathGMI.split('.')[-3]


            if   prod=='gmi':
                with h5py.File(srcPathGMI) as h:
                    datoutTmp = h[var][:]
                #datoutTmp = gmi.load_var_granule(srcPathGMI,var)
            elif prod=='gprof':
                with h5py.File(srcPathGMI) as h:
                    datoutTmp = h[var][:]
                #datoutTmp = gprof.load_var_granule(srcPathGMI,var)
            else:
                print 'check prod',prod
                sys.exit()

            if len(datoutTmp.shape)==2:
                datout = datoutTmp[:,cx-w:cx+w+1] 
            elif len(datoutTmp.shape)==3:
                datout = datoutTmp[:,cx-w:cx+w+1,:] 
            elif len(datoutTmp.shape)==1:
                datout = datoutTmp
            else:
                print 'check dimension',datoutTmp.shape
                sys.exit() 


            varName    = var.split('/')[-1] 
            if len(datoutTmp.shape)==1:
                outbaseDir = '/tank/utsumi/PMM/MATCH.GMI.V%s/%s.GMI.%s'%(fullverGMI, mwscan, varName)

            else:
                outbaseDir = '/tank/utsumi/PMM/MATCH.GMI.V%s/%s.ABp%03d-%03d.GMI.%s'%(fullverGMI, mwscan, cx-w, cx+w, varName)


            outDir     = outbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
            outPath    = outDir + '/%s.%s.npy'%(varName, oid)

            util.mk_dir(outDir)
            np.save(outPath, datout)
            print outPath
    
    #sys.exit()
