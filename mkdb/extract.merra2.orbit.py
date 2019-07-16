from numpy import *
import h5py
import sys, os
import numpy as np
from netCDF4 import *
from bisect import bisect_left
from datetime import datetime, timedelta
import myfunc.util as util
import glob


iDTime = datetime(2017,12,31)
eDTime = datetime(2017,12,31)
lDTimeDay = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))

gmibaseDir = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05'

#raProd    = 'M2T1NXSLV'
raProd    = 'M2I1NXASM'
varName   = 't2m'
rabaseDir = '/work/hk01/utsumi/MERRA2'
draVar    = {'t2m':'T2M'}
raVar     = draVar[varName]

for DTimeDay in lDTimeDay:
    print DTimeDay
    YearDir,MonDir, DayDir = DTimeDay.timetuple()[:3] 
    gmiDir   = gmibaseDir + '/%04d/%02d/%02d'%(YearDir,MonDir,DayDir)

    ssearch  = gmiDir + '/1C.GPM.GMI.*.HDF5'
    lsrcPath = sort(glob.glob(ssearch))
    print ssearch
    print lsrcPath
    for srcPath in lsrcPath:
        #srcPath = '/home/utsumi/temp/1C.GPM.GMI.XCAL2016-C.20171206-S141617-E154850.021437.V05A.HDF5'
        gNum  = srcPath.split('.')[-3]
        with h5py.File(srcPath,'r') as h5:
            a2lat = h5['/S1/Latitude'][:]
            a2lon = h5['/S1/Longitude'][:]
            a1year  = h5['S1/ScanTime/Year'        ][:]
            a1mon   = h5['S1/ScanTime/Month'       ][:]
            a1day   = h5['S1/ScanTime/DayOfMonth'  ][:]
            a1hour  = h5['S1/ScanTime/Hour'        ][:]
            a1mnt   = h5['S1/ScanTime/Minute'      ][:]
 
        latRA0= -90
        lonRA0= -180
        dlatRA = 0.5
        dlonRA = 0.625
        nyRA   = 361
        nxRA   = 576
        
       
        lDTime   = []
        for y,m,d,H,M in map(None,a1year,a1mon,a1day,a1hour,a1mnt):
            lDTime.append( datetime(y,m,d,H,M) )
        
        year0,mon0,day0,hour0,mnt0 = lDTime[0].timetuple()[:5]
        year1,mon1,day1,hour1,mnt1 = lDTime[-1].timetuple()[:5]
        DTime0 = datetime(year0,mon0,day0,hour0)
        DTime1 = datetime(year1,mon1,day1,hour1)
        
        lDTimeHour = util.ret_lDTime(DTime0,DTime1,timedelta(hours=1))
        
        a2out = zeros(a2lat.shape,float32)
        for DTime in lDTimeHour:
            Year,Mon,Day,Hour,Mnt = DTime.timetuple()[:5]
            iy0 = bisect_left(lDTime,DTime-timedelta(minutes=30))
            iy1 = bisect_left(lDTime,DTime+timedelta(minutes=30))
        
            #--- MERRA2 ----
            raDir = rabaseDir + '/%s/%s/%04d%02d'%(raProd,raVar,Year,Mon)
            ssearch = raDir + '/%s.*.%04d%02d%02d.nc4'%(raVar,Year,Mon,Day)
            raPath = glob.glob(ssearch)[0]
            nc     = Dataset(raPath, 'r', format='NETCDF')
            #a1latRA= nc.variables['lat'][:]
            #a1lonRA= nc.variables['lon'][:]
            a2var  = nc.variables[raVar][Hour,:,:]
            #--- corresponding pixels --
            a2y    = floor((a2lat[iy0:iy1] - latRA0 +dlatRA*0.5)/dlatRA).astype(int32)
            a2x    = floor((a2lon[iy0:iy1] - lonRA0 +dlonRA*0.5)/dlonRA).astype(int32)
            a2x    = ma.masked_equal(a2x, nxRA).filled(0)
        
            a2out[iy0:iy1] = a2var[a2y.flatten(), a2x.flatten()].reshape(a2y.shape)
   
        #------------     
        outbaseDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.MERRA2.%s'%(varName) 
        outDir     = outbaseDir + '/%04d/%02d/%02d'%(YearDir,MonDir,DayDir)
        outPath    = outDir + '/%s.%s.npy'%(varName, gNum)
        util.mk_dir(outDir)
        np.save(outPath, a2out.astype(float32))
        print outPath
        
    
    
    
