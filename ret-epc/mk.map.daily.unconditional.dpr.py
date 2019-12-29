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
elif myhost == 'well':
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    workbaseDir = '/home/utsumi/mnt/lab_work'

else:
    print 'check myhost'
    sys.exit()
#*******************************
#iDTime = datetime(2014,7,1)
#eDTime = datetime(2014,8,31)
iDTime = datetime(2014,9,1)
eDTime = datetime(2015,5,31)

lDTime= util.ret_lDTime(iDTime, eDTime, timedelta(days=1))
lYM_skip = [[2014,12],[2015,1],[2015,2]]


miss_out= -9999.
lat0 = -60
lon0 = -180
dlatlon = 1.0
rettype = 'dpr'
ny,nx = 120,360

for DTime in lDTime:
    a2s = np.zeros([ny,nx],float32)
    a2n = np.zeros([ny,nx],int32)
    a2ss= np.zeros([ny,nx],float32)
    #-----------
    Year,Mon,Day = DTime.timetuple()[:3]

    if [Year,Mon] in lYM_skip: continue

    srcDir = workbaseDir+ '/hk01/PMM/NASA/GPM.DPRGMI/2B/V06/%04d/%02d/%02d'%(Year,Mon,Day)
    ssearch  = srcDir + '/2B.GPM.DPRGMI.CORRA*.V06A.HDF5'
    lsrcPath = sort(glob.glob(ssearch))
    
    for srcPath in lsrcPath:
        oid = int(srcPath.split('/')[-1].split('.')[-3])
        if not os.path.exists(srcPath):
            continue
        print DTime,oid
   
        with h5py.File(srcPath,'r') as h: 
            a2prec = h['NS/surfPrecipTotRate'][:]
            a2lat  = h['NS/Latitude'][:]
            a2lon  = h['NS/Longitude'][:]


        a1prec = a2prec.flatten()
        a1lat  = a2lat.flatten()
        a1lon  = a2lon.flatten()
        a1y    = np.floor(a1lat-lat0).astype(int32)
        a1x    = np.floor(a1lon-lon0).astype(int32)

        a1flagp = ma.masked_greater_equal(a1prec,0).mask
        a1flagy = ma.masked_inside(a1y,0,ny-1).mask
        a1flagx = ma.masked_inside(a1x,0,nx-1).mask
        a1flag  = a1flagp * a1flagx * a1flagy
        a1prec  = a1prec[a1flag]
        a1y     = a1y[a1flag]
        a1x     = a1x[a1flag]
 
        for i in range(len(a1prec)):
            y,x = a1y[i], a1x[i]
            a2s[y,x]  = a2s[y,x] + a1prec[i]
            a2ss[y,x] = a2ss[y,x]+ a1prec[i]**2
            a2n[y,x]  += 1

    #-- Save -------
    outDir = tankbaseDir + '/utsumi/PMM/validprof/map-daily-uncond/%s'%(rettype)

    util.mk_dir(outDir)
    np.save(outDir+'/prec.sum.%04d%02d%02d.sp.one.npy' %(Year,Mon,Day), a2s)
    np.save(outDir+'/prec.num.%04d%02d%02d.sp.one.npy' %(Year,Mon,Day), a2n)
    np.save(outDir+'/prec.sum2.%04d%02d%02d.sp.one.npy'%(Year,Mon,Day), a2ss)
    print   outDir+'/prec.sum.%04d%02d%02d.sp.one.npy' %(Year,Mon,Day)

