from numpy import *
import h5py
import sys, os
import numpy as np
from netCDF4 import *
from bisect import bisect_left
from datetime import datetime, timedelta
import myfunc.util as util
import glob
import socket
import metpy.calc as mpcalc
from metpy.units import units

#noscreen = True
noscreen = False
ldtenv = [0]

iDTime = datetime(2017,7,1)
eDTime = datetime(2017,7,31)
lDTimeDay = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
dprver = 'V06'
dprverfull='V06A'

#lvar =['tp','cape','tcwv','mvimd']
#lvar =['cape']
#lvar =['r2m','q2m']
#lvar =['2t','ept2m','tv2m']
#lvar =['skt','r2m','q2m','tp','cape','tcwv','mvimd','2t','ept2m','tv2m']
#lvar = ['r_15','r_45','r_75']
#lvar = ['ept_15','ept_45','ept_75']
lvar = ['tv_15','tv_45','tv_75']
#lvar =['skt','tcwv','mvimd']

myhost = socket.gethostname()
if myhost =='shui':
    dprbaseDir = '/work/hk01/PMM/NASA/GPM.Ku/2A/%s'%(dprver)
    erabaseDir = '/tank/utsumi/era5'
    outbaseDir = '/tank/utsumi/env/pair'
elif myhost == 'well':
    dprbaseDir = '/home/utsumi/mnt/lab_work/hk01/PMM/NASA/GPM.Ku/2A/%s'%(dprver)
    erabaseDir = '/media/disk2/share/data/era5'
    outbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/env/pair'
    climbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/env/eramean'
else:
    print 'check myhost'
    sys.exit()

miss_out = -9999.
#**********************
def read_var_2d_hour(var,Year,Mon,Day,Hour):
    srcDir = erabaseDir + '/%s/%04d%02d'%(var,Year,Mon)
    srcPath= srcDir + '/%s.%04d.%02d.%02d.nc'%(var,Year,Mon,Day)
    ncvar = {'2t':'t2m', '2d':'d2m', 'sp':'sp'
            ,'cape':'cape','cin':'cin'
            ,'tp':'tp','ptype':'ptype'
            ,'tcwv':'tcwv','mvimd':'mvimd'
            ,'skt':'skt'
            }[var]

    with Dataset(srcPath) as np:
        a2var = np.variables[ncvar][Hour]
    return a2var
#**********************
def calc_var_2d_hour(var,Year,Mon,Day,Hour):
    if var =='r2m':
        a2t  = read_var_2d_hour('2t',Year,Mon,Day,Hour) * units.K
        a2dew= read_var_2d_hour('2d',Year,Mon,Day,Hour) * units.K
        a2var= mpcalc.relative_humidity_from_dewpoint(a2t, a2dew)
    
    elif var =='q2m':
        a2t  = read_var_2d_hour('2t',Year,Mon,Day,Hour) * units.K
        a2dew= read_var_2d_hour('2d',Year,Mon,Day,Hour) * units.K
        a2sp = read_var_2d_hour('sp',Year,Mon,Day,Hour) * units.Pa

        a2r  = mpcalc.relative_humidity_from_dewpoint(a2t, a2dew)
        a2w  = mpcalc.mixing_ratio_from_relative_humidity(a2r, a2t, a2sp)
        a2var= mpcalc.specific_humidity_from_mixing_ratio(a2w)

    elif var =='ept2m':
        a2t  = read_var_2d_hour('2t',Year,Mon,Day,Hour) * units.K
        a2dew= read_var_2d_hour('2d',Year,Mon,Day,Hour) * units.K
        a2sp = read_var_2d_hour('sp',Year,Mon,Day,Hour) * units.Pa
        a2var= mpcalc.equivalent_potential_temperature(a2sp, a2t, a2dew)

    elif var =='tv2m':
        a2t  = read_var_2d_hour('2t',Year,Mon,Day,Hour) * units.K
        a2dew= read_var_2d_hour('2d',Year,Mon,Day,Hour) * units.K
        a2sp = read_var_2d_hour('sp',Year,Mon,Day,Hour) * units.Pa

        a2r   = mpcalc.relative_humidity_from_dewpoint(a2t, a2dew)
        a2w   = mpcalc.mixing_ratio_from_relative_humidity(a2r, a2t, a2sp)
        a2var = mpcalc.virtual_temperature(a2t, a2w)  # Add @ 2019/9/29


    else:
        print 'check var in calc_var_2d_hour',var
        sys.exit()
    return a2var


def ret_a2eraclim(var, Year,Mon):
    if len(var.split('_')) == 2:
        varName, lev = var.split('_')
        lev = float(lev)/10
    else:
        varName = var
        lev = 0

    climDir = climbaseDir + '/%s'%(varName)
    sumPath= climDir + '/sum.%04.1fkm.%04d.%02d.npy'%(lev,Year,Mon)
    ssmPath= climDir + '/ssm.%04.1fkm.%04d.%02d.npy'%(lev,Year,Mon)
    numPath= climDir + '/num.%04.1fkm.%04d.%02d.npy'%(lev,Year,Mon)

    a2sum = np.load(sumPath)
    a2ssm = np.load(ssmPath)
    a2num = np.load(numPath)

    a2ave = ma.masked_where(a2num==0, a2sum) / a2num
    a2std = a2ssm - 2*a2ave*a2sum + a2num * np.square(a2ave)

    a2std = np.sqrt( a2std / a2num )
    a2std = ma.masked_invalid(a2std).filled(miss_out)

    return a2ave, a2std

#**********************
dYMpre = {}
for var in lvar:
    dYMpre[var] = (-9999,-9999)

d2ave = {}
d2std = {}
for DTimeDay in lDTimeDay:
    print DTimeDay
    YearDir,MonDir, DayDir = DTimeDay.timetuple()[:3] 
    dprDir   = dprbaseDir + '/%04d/%02d/%02d'%(YearDir,MonDir,DayDir)
    #ssearch  = dprDir + '/2A.GPM.Ku.V8-20180723.20170101-S003624-E020857.016155.V06A.HDF5'
    ssearch  = dprDir + '/2A.GPM.Ku.*.%s.HDF5'%(dprverfull)
    lsrcPath = sort(glob.glob(ssearch))
    if len(lsrcPath)==0:
        continue

    for srcPath in lsrcPath:
        oid  = srcPath.split('.')[-3]

        #if oid != '019015': continue   # test

        with h5py.File(srcPath,'r') as h5:
            a2prec= h5['NS/SLV/precipRateNearSurface'][:]
            a2lat = h5['/NS/Latitude'][:]
            a2lon = h5['/NS/Longitude'][:]
            a1year  = h5['/NS/ScanTime/Year'        ][:]
            a1mon   = h5['/NS/ScanTime/Month'       ][:]
            a1day   = h5['/NS/ScanTime/DayOfMonth'  ][:]
            a1hour  = h5['/NS/ScanTime/Hour'        ][:]
            a1mnt   = h5['/NS/ScanTime/Minute'      ][:]

 
      
        lDTime   = []
        for y,m,d,H,M in map(None,a1year,a1mon,a1day,a1hour,a1mnt):
            lDTime.append( datetime(y,m,d,H,M) )
        
        year0,mon0,day0,hour0,mnt0 = lDTime[0].timetuple()[:5]
        year1,mon1,day1,hour1,mnt1 = lDTime[-1].timetuple()[:5]
        DTime0 = datetime(year0,mon0,day0,hour0)
        DTime1 = datetime(year1,mon1,day1,hour1) + timedelta(hours=1)
        
        lDTimeHour = util.ret_lDTime(DTime0,DTime1,timedelta(hours=1))
      
        for dtenv in ldtenv: 
            for var in lvar:
                ny,nx = a2lat.shape 
    
                a1ave = array([])
                a1std = array([])
                for DTime in lDTimeHour:
                    iy0 = bisect_left(lDTime,DTime-timedelta(minutes=30))
                    iy1 = bisect_left(lDTime,DTime+timedelta(minutes=30))
                    if (iy1==0) or (iy0==ny):
                        continue
    
                    Year,Mon,Day,Hour,Mnt = (DTime + timedelta(hours=dtenv)).timetuple()[:5]
    
                    #--- Read ERA clim ------
                    if (Year,Mon) != dYMpre[var]:
                        d2ave[var], d2std[var] = ret_a2eraclim(var, Year,Mon) 
                        dYMpre[var] = (Year,Mon)

                        print 'Read',var
                    #--- corresponding pixels --
                    latRA0= 90   # from North to South
                    lonRA0= 0
                    dlatRA = 0.25
                    dlonRA = 0.25
                    nyRA   = 721
                    nxRA   = 1440
            
                    #a2y    = floor((a2lat[iy0:iy1] - latRA0 +dlatRA*0.5)/dlatRA).astype(int32)
                    a2y    = floor((latRA0 + dlatRA*0.5 - a2lat[iy0:iy1])/dlatRA).astype(int32)
                    a2x    = floor((a2lon[iy0:iy1] - lonRA0 +dlonRA*0.5)/dlonRA).astype(int32)
                    a2x    = ma.masked_equal(a2x, nxRA).filled(0)
    
    
                    if noscreen is True:          
                        a2prmask = False
                    else:
                        a2prmask= ma.masked_less_equal(a2prec[iy0:iy1],0).mask
                    
                    a1x    = ma.masked_where(a2prmask, a2x).compressed()
                    a1y    = ma.masked_where(a2prmask, a2y).compressed()
    
                    a1aveTmp  = d2ave[var][a1y,a1x]
                    a1stdTmp  = d2std[var][a1y,a1x]
                    a1ave = concatenate([a1ave, a1aveTmp])
                    a1std = concatenate([a1std, a1stdTmp])
    
                #**** If noscreen ********************
                if noscreen is True:
                    nydpr,nxdpr = a2prec.shape
                    a1var = a1var.reshape(nydpr,nxdpr)
    
                #**** Save interpolated variables ****
                #**** Save variable ****
                shead = ''
                if noscreen is True:
                    shead = shead + 'full.'
   
                if len(var.split('_'))==2:
                    varName, lev = var.split('_')
                    lev = float(lev)/10
                else:
                    varName = var
                    lev = 0


                outDir     = outbaseDir + '/%s/%04d/%02d/%02d'%(varName,YearDir,MonDir,DayDir)
                util.mk_dir(outDir)

                avePath = outDir + '/%save.%03dh.%04.1fkm.%s.npy'%(shead, dtenv, lev, oid)
                stdPath = outDir + '/%sstd.%03dh.%04.1fkm.%s.npy'%(shead, dtenv, lev, oid)

                if type(a1ave) is np.ma.core.MaskedArray:
                    a1ave = a1ave.filled(miss_out)
                if type(a1std) is np.ma.core.MaskedArray:
                    a1std = a1std.filled(miss_out)

                np.save(avePath, a1ave.astype(float32))
                np.save(stdPath, a1std.astype(float32))
                print avePath 
 
