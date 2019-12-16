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
import calendar

iYM = [2017,7]
eYM = [2017,7]
lYM = util.ret_lYM(iYM, eYM)
dprver = 'V06'
dprverfull='V06A'

lvar =['skt','r2m','q2m','tp','cape','tcwv','mvimd','2t','ept2m','tv2m']
#lvar =['cape']
#lvar =['r2m','q2m']
#lvar =['2t','ept2m','tv2m']
#lvar =['skt','tcwv','mvimd']
#lvar =['skt']
lev = 0  # km
myhost = socket.gethostname()
if myhost =='shui':
    dprbaseDir = '/work/hk01/PMM/NASA/GPM.Ku/2A/%s'%(dprver)
    erabaseDir = '/tank/utsumi/era5'
    outbaseDir = '/tank/utsumi/env/pair'
elif myhost == 'well':
    dprbaseDir = '/home/utsumi/mnt/lab_work/hk01/PMM/NASA/GPM.Ku/2A/%s'%(dprver)
    erabaseDir = '/media/disk2/share/data/era5'
    outbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/env/eramean'
else:
    print 'check myhost'
    sys.exit()

miss_out = -9999.
#--- corresponding pixels --
latRA0 = 90   # from North to South
lonRA0 = 0
dlatRA = 0.25
dlonRA = 0.25
nyRA   = 721
nxRA   = 1440

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

#**********************

a2one = np.ones([nyRA,nxRA]).astype(int32)

for var in lvar:
    for YM in lYM:
        a2sum = np.zeros([nyRA,nxRA])
        a2ssm = np.zeros([nyRA,nxRA])
        a2num = np.zeros([nyRA,nxRA]).astype(int32)

        Year,Mon = YM
        eDay = calendar.monthrange(Year,Mon)[1]
        #eDay = 2
        iDTime = datetime(Year,Mon,1,0)
        eDTime = datetime(Year,Mon,eDay,23)
        lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(seconds=60*60))
    
        for DTime in lDTime:
            Year,Mon,Day,Hour  = DTime.timetuple()[:4] 
            print var, DTime
            #--- Read ERA ******
            if var in ['q2m','r2m','ept2m','tv2m']:
                a2var = calc_var_2d_hour(var,Year,Mon,Day,Hour).data
            else:
                a2var = read_var_2d_hour(var,Year,Mon,Day,Hour).data
      
            a2var = ma.masked_invalid(ma.masked_less(a2var,0))
            a2sum = a2sum + a2var.filled(0.0)
            a2ssm = a2ssm + np.square(a2var.filled(0.0))
            a2num = a2num + ma.masked_where(a2var.mask, a2one).filled(0).astype(int32)

        #-- Save ---
        outDir = outbaseDir + '/%s'%(var)
        util.mk_dir(outDir)
        sumPath= outDir + '/sum.%04.1fkm.%04d.%02d.npy'%(lev,Year,Mon)
        ssmPath= outDir + '/ssm.%04.1fkm.%04d.%02d.npy'%(lev,Year,Mon)
        numPath= outDir + '/num.%04.1fkm.%04d.%02d.npy'%(lev,Year,Mon)

        np.save(sumPath, a2sum)
        np.save(ssmPath, a2ssm)
        np.save(numPath, a2num)


        print sumPath
