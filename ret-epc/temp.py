from numpy import *
import numpy as np
import glob
import myfunc.util as util
iYM = [2017,2]
eYM = [2017,12]
lYM = util.ret_lYM(iYM, eYM)


for [Year,Mon] in lYM:
    print Year,Mon
    baseDir= '/work/hk01/utsumi/PMM/EPCDB/GMI.V05A.S1.ABp103-117'

    ssearch = baseDir + '/DPRGMI_MS_surfPrecipTotRate/%04d%02d/DPRGMI_MS_surfPrecipTotRate.?????.npy'%(Year,Mon)
    lmsPath = glob.glob(ssearch)
    for msPath in lmsPath:
        oid = int(msPath.split('.')[-2])
        msPath = baseDir + '/DPRGMI_MS_surfPrecipTotRate/%04d%02d/DPRGMI_MS_surfPrecipTotRate.%05d.npy'%(Year,Mon,oid)
        nsPath = baseDir + '/DPRGMI_NS_surfPrecipTotRate/%04d%02d/DPRGMI_NS_surfPrecipTotRate.%05d.npy'%(Year,Mon,oid)
    
        a2ms = np.load(msPath)
        a2ns = np.load(nsPath)
   
         
        if a2ms.shape != a2ns.shape:
            print Year,Mon,oid,a2ms.shape, a2ns.shape
