import os, sys
import myfunc.util as util
import glob
from datetime import datetime, timedelta

iDTime = datetime(2014,10,1)
eDTime = datetime(2014,10,15)
lDTime = util.ret_lDTime(iDTime,eDTime, timedelta(days=1))

for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]
    srcDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/glb.relsurf01.minrec1000.maxrec10000/%04d/%02d/%02d'%(Year,Mon,Day)

    ssearch = srcDir + '/prwatprofNS.??????.y-9999--9999.nrec10000.npy'
    print ssearch
    lsrcPath = glob.glob(ssearch)

    for srcPath in lsrcPath:
        #fname = srcPath.split('/')[-1]
        ifname = os.path.basename(srcPath)
        ofname = 'prwatprofNS-rs' + '.' + '.'.join(ifname.split('.')[1:])

        #outPath = srcDir + '/' + ofname
        #os.rename(srcPath, outPath)
        #print outPath
