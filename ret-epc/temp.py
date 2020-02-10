import numpy as np
from numpy import *
import glob
import myfunc.util as util
from datetime import datetime, timedelta

iDTime = datetime(2014,6,8)
eDTime = datetime(2014,6,8)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))

#baseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/glb.relsurf01.minrec1000.maxrec10000'
baseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/glb.v03.minrec1000.maxrec10000'

for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]
    srcDir = baseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
    ssearch = srcDir + '/prwatprofNS.??????.y-9999--9999.nrec10000.npy'
    lsrcPath = glob.glob(ssearch)
    print ssearch
    print DTime, lsrcPath    
    for srcPath in lsrcPath:
        a = np.load(srcPath)
        print a.min(), a.max() 


#srcDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/EPCDB/samp.10000.GMI.V05A.S1.ABp103-117.01-12/DPRGMI_NS_precipTotWaterContRelSurf'
#ssearch = srcDir + '/?????.npy'
#lsrcPath = glob.glob(ssearch)
#for srcPath in lsrcPath:
#    a=np.load(srcPath)
#    print a.min(), a.max()
