from numpy import *
import h5py
import numpy as np
import socket
import os, sys
import myfunc.util as util
from datetime import datetime, timedelta
import glob
import calendar
#*******************************
myhost = socket.gethostname()
if myhost =='shui':
    tankbaseDir = '/tank'
    workbaseDir = '/work'
    epcbaseDir  = '/tank/utsumi/PMM/retepc'
elif myhost == 'well':
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    workbaseDir = '/home/utsumi/mnt/lab_work'
    #epcbaseDir  = '/media/disk2/share/PMM/retepc/glb.wprof'
    epcbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc'

else:
    print 'check myhost'
    sys.exit()
#*******************************
iDTime = datetime(2014,9,1)
eDTime = datetime(2015,5,31)
#iDTime = datetime(2015,2,1)
#eDTime = datetime(2015,2,28)
lDTime = util.ret_lDTime(iDTime, eDTime, timedelta(days=1))
lYM_skip = [[2014,12],[2015,1],[2015,2]]

DB_MAXREC = 10000
DB_MINREC = 1000
expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)

miss_out= -9999.
lat0 = -60
lon0 = -180
dlatlon = 1.0
rettype = 'epcNScmb'
#rettype = 'epcNS'

ny,nx = 120,360
for DTime in lDTime:
    a2s = np.zeros([ny,nx],float32)
    a2n = np.zeros([ny,nx],int32)
    a2ss= np.zeros([ny,nx],float32)
    #-- EPC ---------
    Year,Mon,Day = DTime.timetuple()[:3]

    if [Year,Mon] in lYM_skip: continue

    srcDir = epcbaseDir + '/%s/%04d/%02d/%02d'%(expr,Year,Mon,Day)
    if rettype=='epcNS':
        ssearch  = srcDir + '/nsurfNS.??????.y-9999--9999.nrec%d.npy'%(DB_MAXREC)
    elif rettype=='epcNScmb':
        ssearch  = srcDir + '/nsurfNScmb.??????.y-9999--9999.nrec%d.npy'%(DB_MAXREC)
    else:
        print 'check rettype',rettype
        sys.exit()
    lsrcPath = sort(glob.glob(ssearch))
    
    for srcPath in lsrcPath:
        oid = int(srcPath.split('/')[-1].split('.')[1])
        if not os.path.exists(srcPath):
            continue
        print DTime,oid
    
        a2prec = ma.masked_less(np.load(srcPath)[:,83:137+1],0.01).filled(0.0)
        a2lat  = np.load(srcDir + '/lat.%06d.y-9999--9999.nrec%d.npy'%(oid, DB_MAXREC))[:,83:137+1]
        a2lon  = np.load(srcDir + '/lon.%06d.y-9999--9999.nrec%d.npy'%(oid, DB_MAXREC))[:,83:137+1]
    
        a1prec = a2prec.flatten()
        a1lat  = a2lat.flatten()
        a1lon  = a2lon.flatten()
        a1y    = np.floor(a1lat-lat0).astype(int32)
        a1x    = np.floor(a1lon-lon0).astype(int32)

        a1flagy = ma.masked_inside(a1y,0,ny-1).mask
        a1flagx = ma.masked_inside(a1x,0,nx-1).mask
        a1flag  = a1flagx * a1flagy
        a1prec  = a1prec[a1flag]
        a1y     = a1y[a1flag]
        a1x     = a1x[a1flag]
 
        for i in range(len(a1prec)):
            y,x = a1y[i], a1x[i]
            a2s[y,x]  = a2s[y,x] + a1prec[i]
            a2ss[y,x] = a2ss[y,x]+ a1prec[i]**2
            a2n[y,x]  += 1

    #-- Save -------
    outDir = tankbaseDir + '/utsumi/PMM/validprof/map-daily-uncond/%s.%s'%(rettype,expr)

    util.mk_dir(outDir)
    np.save(outDir+'/prec.sum.%04d%02d%02d.sp.one.npy' %(Year,Mon,Day), a2s)
    np.save(outDir+'/prec.num.%04d%02d%02d.sp.one.npy' %(Year,Mon,Day), a2n)
    np.save(outDir+'/prec.sum2.%04d%02d%02d.sp.one.npy'%(Year,Mon,Day), a2ss)
    print  outDir+'/prec.sum.%04d%02d%02d.sp.one.npy'  %(Year,Mon,Day)

