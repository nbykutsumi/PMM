from numpy import *
import h5py
import sys, os
import numpy as np
from netCDF4 import *
from bisect import bisect_left
from datetime import datetime, timedelta
import myfunc.util as util
import glob


iDTime = datetime(2018,1,1)
eDTime = datetime(2018,1,1)
lDTimeDay = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))

#gpr_amsr2      = ["GCOMW1","AMSR2"]
#gpr_ssmis_f16  = ["F16","SSMIS"]
#gpr_atms_noaa20= ["NOAA20","ATMS"]
amsr2      = ["GCOMW1","AMSR2"]
ssmis_f16  = ["F16","SSMIS"]
ssmis_f17  = ["F17","SSMIS"]
ssmis_f18  = ["F18","SSMIS"]
atms_npp   = ["NPP","ATMS"]
atms_noaa20= ["NOAA20","ATMS"]
mhs_metopa = ["METOPA","MHS"]
mhs_metopb = ["METOPB","MHS"]

lspec = [amsr2, ssmis_f16, ssmis_f17, ssmis_f18, atms_npp, atms_noaa20, mhs_metopa, mhs_metopb]
#lspec = [atms_npp, atms_noaa20, mhs_metopa, mhs_metopb]

dmainscan = {'AMSR2':1, 'SSMIS':1, 'ATMS':1, 'MHS':1}


#raProd    = 'M2T1NXSLV'
raProd    = 'M2I1NXASM'
#lvarName  = ['tqv','t2m']
#lvarName  = ['tqv','t2m']
lvarName  = ['t2m']
#lvarName  = ['tqv']
#rabaseDir = '/work/hk01/utsumi/MERRA2'
rabaseDir = '/home/utsumi/mnt/lab_tank/utsumi/data/MERRA2'
draVar    = {'t2m':'T2M','tqv':'TQV'}

latRA0= -90
lonRA0= -180
dlatRA = 0.5
dlonRA = 0.625
nyRA   = 361
nxRA   = 576

miss = -9999.

for (sate,sensor) in lspec:
    pmwbaseDir = '/home/utsumi/mnt/lab_work/hk02/PMM/NASA/%s.%s/1C/V05'%(sate,sensor)
    mainscan = dmainscan[sensor]

    for DTimeDay in lDTimeDay:
        print DTimeDay
        YearDir,MonDir, DayDir = DTimeDay.timetuple()[:3] 
        pmwDir   = pmwbaseDir + '/%04d/%02d/%02d'%(YearDir,MonDir,DayDir)
    
        ssearch  = pmwDir + '/1C.%s.%s.*.HDF5'%(sate,sensor)
        lsrcPath = sort(glob.glob(ssearch))
        print ssearch
        print lsrcPath
        for srcPath in lsrcPath:
            #srcPath = '/home/utsumi/temp/1C.GPM.GMI.XCAL2016-C.20171206-S141617-E154850.021437.V05A.HDF5'
            gNum  = srcPath.split('.')[-3]
            with h5py.File(srcPath,'r') as h5:
                a2lat   = h5['/S1/Latitude'][:]
                a2lon   = h5['/S1/Longitude'][:]
                a1year  = h5['S1/ScanTime/Year'        ][:]
                a1mon   = h5['S1/ScanTime/Month'       ][:]
                a1day   = h5['S1/ScanTime/DayOfMonth'  ][:]
                a1hour  = h5['S1/ScanTime/Hour'        ][:]
                a1mnt   = h5['S1/ScanTime/Minute'      ][:]
    
            nypmw,nxpmw = a2lat.shape

            lDTime   = []
            for y,m,d,H,M in map(None,a1year,a1mon,a1day,a1hour,a1mnt):
                lDTime.append( datetime(y,m,d,H,M) )
            
            year0,mon0,day0,hour0,mnt0 = lDTime[0].timetuple()[:5]
            year1,mon1,day1,hour1,mnt1 = lDTime[-1].timetuple()[:5]
            DTime0 = datetime(year0,mon0,day0,hour0)
            DTime1 = datetime(year1,mon1,day1,hour1) + timedelta(hours=1)
            
            lDTimeHour = util.ret_lDTime(DTime0,DTime1,timedelta(hours=1))
    
    
            for varName in lvarName:     
                raVar = draVar[varName]
                a2out = zeros(a2lat.shape,float32)
                for DTime in lDTimeHour:
                    Year,Mon,Day,Hour,Mnt = DTime.timetuple()[:5]
                    iy0 = bisect_left(lDTime,DTime-timedelta(minutes=30))
                    iy1 = bisect_left(lDTime,DTime+timedelta(minutes=30))
               
                    ny = a2lat.shape[0]
                    if (iy1==0)or(iy0==ny):
                        continue 
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

                    #-- mask missing data in lat & lon --                
                    a2mask = ma.masked_outside(a2lat[iy0:iy1], -90,90).mask
                    a2mask = ma.masked_outside(a2lon[iy0:iy1], -180,180).mask + a2mask
                    a2y = ma.masked_where(a2mask, a2y).filled(0)
                    a2x = ma.masked_where(a2mask, a2x).filled(0)

                    #a2out[iy0:iy1] = a2var[a2y.flatten(), a2x.flatten()].reshape(a2y.shape)
                    a2out[iy0:iy1] = ma.masked_where(a2mask, a2var[a2y, a2x]).filled(miss)

           
                #------------     
                #outbaseDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.MERRA2.%s'%(varName) 
                outbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/MATCH.%s.%s.V05/S%d.ABp000-%03d.MERRA2.%s'%(sate,sensor,mainscan, nxpmw-1, varName) 
                outDir     = outbaseDir + '/%04d/%02d/%02d'%(YearDir,MonDir,DayDir)
                outPath    = outDir + '/%s.%s.npy'%(varName, gNum)
                util.mk_dir(outDir)
                np.save(outPath, a2out.astype(float32))
                print outPath
            
        
        
    
