from numpy import *
import numpy as np
import myfunc.util as util
from datetime import datetime, timedelta
import glob
iDTime = datetime(2014,6,1)
eDTime = datetime(2014,6,1)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))

baseDir1 = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/glb.minrec1000.maxrec10000'
#baseDir1 = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/glb.v03.minrec1000.maxrec10000'
baseDir2 = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/glb.v04.minrec1000.maxrec10000'


for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]
    ssearch = baseDir2 + '/%04d/%02d/%02d/nsurfNScmb.??????.y-9999--9999.nrec10000.npy'%(Year,Mon,Day)
    #lsrcPath = np.sort(glob.glob(ssearch))
    lsrcPath = np.sort(glob.glob(ssearch))[:6]

    s1 = 0.
    s2 = 0.
    for i,srcPath in enumerate(lsrcPath):
        oid = int(srcPath.split('.')[-4])
        #if oid not in [1454]: continue

        srcPath1 = baseDir1 + '/%04d/%02d/%02d/nsurfNScmb.%06d.y-9999--9999.nrec10000.npy'%(Year,Mon,Day,oid)
        srcPath2 = baseDir2 + '/%04d/%02d/%02d/nsurfNScmb.%06d.y-9999--9999.nrec10000.npy'%(Year,Mon,Day,oid)

        latPath1 = baseDir1 + '/%04d/%02d/%02d/lat.%06d.y-9999--9999.nrec10000.npy'%(Year,Mon,Day,oid)
        latPath2 = baseDir2 + '/%04d/%02d/%02d/lat.%06d.y-9999--9999.nrec10000.npy'%(Year,Mon,Day,oid)

        lonPath1 = baseDir1 + '/%04d/%02d/%02d/lon.%06d.y-9999--9999.nrec10000.npy'%(Year,Mon,Day,oid)
        lonPath2 = baseDir2 + '/%04d/%02d/%02d/lon.%06d.y-9999--9999.nrec10000.npy'%(Year,Mon,Day,oid)

        topidxdbPath1 = baseDir1 + '/%04d/%02d/%02d/top-idxdbNS.%06d.y-9999--9999.nrec10000.npy'%(Year,Mon,Day,oid)
        topidxdbPath2 = baseDir2 + '/%04d/%02d/%02d/top-idxdbNS.%06d.y-9999--9999.nrec10000.npy'%(Year,Mon,Day,oid)



        a2dat1 = ma.masked_less(np.load(srcPath1),0)
        a2dat2 = ma.masked_less(np.load(srcPath2),0)

        a2lat1 = np.load(latPath1)
        a2lat2 = np.load(latPath2)

        a2lon1 = np.load(lonPath1)
        a2lon2 = np.load(lonPath2)

        a2topidxdb1 = np.load(topidxdbPath1)
        a2topidxdb2 = np.load(topidxdbPath2)

        print ''
        #print srcPath2
        #print srcPath2
        print DTime,i,oid, a2dat1.mean(), a2dat2.mean()
        #print DTime,i,oid, a2dat1.count(), a2dat2.count()
        #print DTime,i,oid, s1, s2
        #print a2lat1.shape
        #print a2lat1.mean(), a2lat2.mean()
        #print a2topidxdb1[1000:1003,50:53]
        #print a2topidxdb2[1000:1003,50:53]
