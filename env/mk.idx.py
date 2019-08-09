from numpy import *
import h5py
import myfunc.util as util
from netCDF4 import *
from bisect import bisect_left
from datetime import datetime, timedelta
import glob
import socket

iDTime = datetime(2017,7,3)
eDTime = datetime(2017,7,3)
lDTimeDay = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
dprver = 'V06'
dprverfull='V06A'

lvar =['cape']

myhost = socket.gethostname()
if myhost =='shui':
    dprbaseDir = '/work/hk01/PMM/NASA/GPM.Ku/2A/%s'%(dprver)
    rabaseDir  = '/tank/utsumi/era5'
    outbaseDir = '/tank/utsumi/env/pair'
elif myhost == 'well':
    dprbaseDir = '/home/utsumi/mnt/lab_work/hk01/PMM/NASA/GPM.Ku/2A/%s'%(dprver)
    rabaseDir  = '/home/utsumi/mnt/lab_tank/utsumi/era5'
    outbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/env/pair'
else:
    print 'check myhost'
    sys.exit()

for DTimeDay in lDTimeDay:
    print DTimeDay
    YearDir,MonDir, DayDir = DTimeDay.timetuple()[:3]
    dprDir   = dprbaseDir + '/%04d/%02d/%02d'%(YearDir,MonDir,DayDir)
    #ssearch  = dprDir + '/2A.GPM.Ku.V8-20180723.20170101-S003624-E020857.016155.V06A.HDF5'
    ssearch  = dprDir + '/2A.GPM.Ku.*.%s.HDF5'%(dprverfull)
    lsrcPath = sort(glob.glob(ssearch))
    for srcPath in lsrcPath:
        oid  = srcPath.split('.')[-3]
        with h5py.File(srcPath,'r') as h5:
            a2prec= h5['NS/SLV/precipRateNearSurface'][:]
            a2lat = h5['/NS/Latitude'][:]
            a2lon = h5['/NS/Longitude'][:]
            a1year  = h5['/NS/ScanTime/Year'        ][:]
            a1mon   = h5['/NS/ScanTime/Month'       ][:]
            a1day   = h5['/NS/ScanTime/DayOfMonth'  ][:]
            a1hour  = h5['/NS/ScanTime/Hour'        ][:]
            a1mnt   = h5['/NS/ScanTime/Minute'      ][:]


        #*** Mask for precip=0 ***
        a2mask = ma.masked_less(a2prec, 0).mask
        #*** y,x for ERA ****
        latRA0= -90
        lonRA0= -180
        dlatRA = 0.25
        dlonRA = 0.25
        nyRA   = 721
        nxRA   = 1440

        a2pyera    = floor((a2lat[iy0:iy1] - latRA0 +dlatRA*0.5)/dlatRA).astype(int32)
        a2pxera    = floor((a2lon[iy0:iy1] - lonRA0 +dlonRA*0.5)/dlonRA).astype(int32)
        a2pxera    = ma.masked_equal(a2pxera, nxRA).filled(0)

        a1pyera = ma.masked_where(a2mask, a2pyera).compressed()
        a1pxera = ma.masked_where(a2mask, a2pxera).compressed()

        #*** y,x for DPR ****
        nydpr, nxdpr = a2lat.shape
        a2pxdpr, a2pydpr = np.meshgrid(nxdpr, nydpr)
        a1pydpr = ma.masked_where(a2mask, a2pydpr).compressed()
        a1pxdpr = ma.masked_where(a2mask, a2pxdpr).compressed()

        print a1pyera.shape, a1pxera.shape, a1pydpr.shape, a1pxdpr.shape


