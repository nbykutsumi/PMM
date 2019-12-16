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

#lvar = ['r','ept','tv','w']
#lvar = ['r','ept']
#lvar = ['tv','w']
lvar = ['w']
#lvar = ['ept']

lev = 0  # km
lplev = [975,900,825,750,600,450,300,200]
a1p   = np.array(lplev).reshape(-1,1,1) * units.hPa
lzmeter = [1500,4500,7500]
nz = len(lzmeter)
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
def read_var_2d(var,Year,Mon,Day):
    srcDir = erabaseDir + '/%s/%04d%02d'%(var,Year,Mon)
    srcPath= srcDir + '/%s.%04d.%02d.%02d.nc'%(var,Year,Mon,Day)
    ncvar = {'2t':'t2m', '2d':'d2m', 'sp':'sp'
            ,'cape':'cape','cin':'cin'
            ,'tp':'tp','ptype':'ptype'
            ,'tcwv':'tcwv','mvimd':'mvimd'
            ,'skt':'skt'
            }[var]

    with Dataset(srcPath) as np:
        a3var = np.variables[ncvar][:]
    return a3var


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


def read_var_3d(var,Year,Mon,Day):
    srcDir = erabaseDir + '/%s/%04d%02d'%(var,Year,Mon)
    srcPath= srcDir + '/%s.%04d.%02d.%02d.nc'%(var,Year,Mon,Day)
    with Dataset(srcPath) as np:
        a3var = np.variables[var][:]
    return a3var


def read_var_3d_hour(var,Year,Mon,Day,Hour):
    srcDir = erabaseDir + '/%s/%04d%02d'%(var,Year,Mon)
    srcPath= srcDir + '/%s.%04d.%02d.%02d.nc'%(var,Year,Mon,Day)
    with Dataset(srcPath) as np:
        a3var = np.variables[var][Hour]
    return a3var


def interp_vertical(a3gph, a3dat,a2surfdat, a2orog, lh):
    flag2d = 0
    if len(a3dat.shape)==2:
        flag2d = 1
        nplev = a3gph.shape[0]
        a3gph = a3gph.reshape(nplev,1,-1)
        a3dat = a3dat.reshape(nplev,1,-1)
        a2surfdat=a2surfdat.reshape(1,-1)
        a2orog   =a2orog.reshape(1,-1)


    a3depth = a3gph - a2orog  # meter from suface
    nlayer,ny,nx = a3depth.shape
    X,Y = np.meshgrid(range(nx),range(ny))
    a3out = empty([len(lh),ny,nx],float32)
    for ihintp, hintp in enumerate(lh):
        #*** Fill under suface layer with surface variables ****
        a3datTmp = a3dat.copy()
        for ilayer in range(nlayer):
            a2flag = ma.masked_less(a3depth[ilayer], 0).mask
            if a2flag.sum()==0:
                continue
            a3depth[ilayer][a2flag] = 0
            a3datTmp[ilayer][a2flag] = a2surfdat[a2flag]
        #*** Find lower and upper layer ***
        a2ilow = ma.masked_less(a3depth, hintp).argmin(axis=0) -1
        a2iup  = a2ilow + 1
        #*** Interpolation ***
        a2mask = ma.masked_less(a3depth[-1], hintp).mask
        a2iup[a2mask] = a2ilow[a2mask]   # temporally replace

        #*** lower layer *****
        a2depthlow = a3depth[a2ilow,Y,X]-1
        a2datlow   = a3datTmp[a2ilow,Y,X]

        a2flag = ma.masked_equal(a2ilow,-1).mask
        if a2flag.sum() !=0:
            a2depthlow[a2flag] = 0
            a2datlow[a2flag] = a2surfdat[a2flag]

        #*** upper later *****
        a2depthup  = a3depth[a2iup, Y,X]+1
        a2datup    = a3datTmp[a2iup, Y,X]

        #*** Interpolation ***
        a2datintp = ((hintp-a2depthlow)*a2datup + (a2depthup-hintp)*a2datlow)/(a2depthup-a2depthlow)
        a3out[ihintp] = ma.masked_where(a2mask, a2datintp).filled(miss_out)

    if flag2d ==1:
        a3out = a3out.reshape(len(lh),-1)

    #x,y = 100,100
    #print 'lh'
    #print lh
    #print 'org dat'
    #print a3dat[:,x,y]
    #print 'depth'
    #print a3depth[:,x,y]
    #print a3out[:,x,y]
    #sys.exit() 

    return a3out



def read_zmeter(Year,Mon,Day):
    ''' geopotential height [m] '''
    g = 9.80665
    var = 'z'
    srcDir = erabaseDir + '/%s/%04d%02d'%(var,Year,Mon)
    srcPath= srcDir + '/%s.%04d.%02d.%02d.nc'%(var,Year,Mon,Day)
    with Dataset(srcPath) as np:
        a3var = np.variables[var][:]
    return a3var/g



def read_zmeter_hour(Year,Mon,Day,Hour):
    ''' geopotential height [m] '''
    g = 9.80665
    var = 'z'
    srcDir = erabaseDir + '/%s/%04d%02d'%(var,Year,Mon)
    srcPath= srcDir + '/%s.%04d.%02d.%02d.nc'%(var,Year,Mon,Day)
    with Dataset(srcPath) as np:
        a2var = np.variables[var][Hour]
    return a2var/g



def read_orogmeter():
    return np.load(erabaseDir + '/orog/orog.meter.na.npy')


#**********************

a3one = np.ones([nz,nyRA,nxRA]).astype(int32)

for var in lvar:
    for YM in lYM:
        a3sum = np.zeros([nz,nyRA,nxRA])
        a3ssm = np.zeros([nz,nyRA,nxRA])
        a3num = np.zeros([nz,nyRA,nxRA]).astype(int32)


        Year,Mon = YM
        eDay = calendar.monthrange(Year,Mon)[1]
        #eDay = 1
        for Day in range(1,eDay+1): 
            iDTime = datetime(Year,Mon,Day,0)
            eDTime = datetime(Year,Mon,Day,23)
            lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(seconds=60*60))

            #--- Read ERA ******
            a4zmeter = read_zmeter(Year,Mon,Day)[:,::-1,:,:] * units.m
            if var in ['ept','tv']:
                a4t = read_var_3d('t', Year,Mon,Day)[:,::-1,:,:] * units.K
                a4q = read_var_3d('q', Year,Mon,Day)[:,::-1,:,:] * units.g/units.g  # specific humidity
                a3tsurf = read_var_2d('2t',Year,Mon,Day) * units.K
                a3dsurf = read_var_2d('2d',Year,Mon,Day) * units.K
                a3sp    = read_var_2d('sp',Year,Mon,Day) /100.* units.hPa

            elif var in ['r']:
                a4r = read_var_3d(var, Year,Mon,Day)[:,::-1,:,:] * units.g/units.g  # specific humidity
                a3tsurf = read_var_2d('2t',Year,Mon,Day) * units.K
                a3dsurf = read_var_2d('2d',Year,Mon,Day) * units.K

            elif var=='w':
                a4w = read_var_3d('w', Year,Mon,Day)[:,::-1,:,:]


            #-------------------
            for DTime in lDTime:
                Year,Mon,Day,Hour  = DTime.timetuple()[:4] 
                print var, DTime


                if var in ['ept','tv']:
                    a3t = a4t[Hour]
                    a3q = a4q[Hour]
                    a2tsurf = a3tsurf[Hour]
                    a2dsurf = a3dsurf[Hour]
                    a2sp    = a3sp[Hour]
    
                elif var in ['r']:
                    a3r     = a4r[Hour]
                    a2tsurf = a3tsurf[Hour]
                    a2dsurf = a3dsurf[Hour]

                elif var=='w':
                    a3w     = a4w[Hour]
 
                #--- Make variables----
                if var =='ept':  # Equivalent potential temperature
                    #--- 3D ---
                    a3dew=mpcalc.dewpoint_from_specific_humidity(a3q, a3t, a1p)  # degC
                    a3var = mpcalc.equivalent_potential_temperature(a1p, a3t, a3dew)
    
                    #--- Surface -
                    a2surfvar = mpcalc.equivalent_potential_temperature(a2sp, a2tsurf, a2dsurf)
    
                elif var =='tv':  # Virtual temperature
                    a3m = mpcalc.mixing_ratio_from_specific_humidity(a3q)  # mixing ratio
                    a3var = mpcalc.virtual_temperature(a3t,  a3m)  # virtual temperature
    
                    #--- Surface -
                    #a2surfvar = mpcalc.equivalent_potential_temperature(a2sp, a2tsurf, a2dsurf) # comment out @ 2019/9/29
                    a2rsurf   = mpcalc.relative_humidity_from_dewpoint(a2tsurf, a2dsurf)
                    a2wsurf   = mpcalc.mixing_ratio_from_relative_humidity(a2rsurf, a2tsurf, a2sp)
                    a2surfvar = mpcalc.virtual_temperature(a2tsurf, a2wsurf)  # Add @ 2019/9/29
    
                elif var=='r':
                    #a3var = mpcalc.relative_humidity_from_specific_humidity(a3q, a3t, a1p)
                    a3var = a3r * units.g/units.g  # specific humidity
    
                    a2surfvar = mpcalc.relative_humidity_from_dewpoint(a2tsurf, a2dsurf)
    
                elif var=='w':
                    a3var = a3w
                    a2surfvar = np.zeros([nyRA,nxRA], float32)
    
                else:
                    print 'check var',var
                    sys.exit()
    
                #--- Read orography and zmeter ---
                a3zmeter = read_zmeter_hour(Year,Mon,Day,Hour)[::-1,:,:] 
                a2orog = read_orogmeter()
    
                #--- Interpolation ------
                a3intp = interp_vertical(a3zmeter, a3var, a2surfvar, a2orog, lzmeter)
    
                a3intp = ma.masked_equal(a3intp, miss_out) 
                a3sum  = a3sum + a3intp.filled(0.0)
                a3ssm  = a3ssm + np.square(a3intp).filled(0.0)
                a3num  = a3num + ma.masked_where(a3intp.mask, a3one).filled(0).astype(int32)
        #-- Save ---

        for ilev,lev in enumerate(lzmeter):
            a2sum = a3sum[ilev]
            a2ssm = a3ssm[ilev]
            a2num = a3num[ilev]
            lev = lev*0.001 
            outDir = outbaseDir + '/%s'%(var)
            util.mk_dir(outDir)
            sumPath= outDir + '/sum.%04.1fkm.%04d.%02d.npy'%(lev,Year,Mon)
            ssmPath= outDir + '/ssm.%04.1fkm.%04d.%02d.npy'%(lev,Year,Mon)
            numPath= outDir + '/num.%04.1fkm.%04d.%02d.npy'%(lev,Year,Mon)
    
            np.save(sumPath, a2sum)
            np.save(ssmPath, a2ssm)
            np.save(numPath, a2num)
    
    
            print sumPath
