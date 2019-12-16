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
#preenv   = True
#ldtenv = [-12,-6,-3,-1,0,1,3,6,12]
ldtenv = [-6,-1,0,6]

iDTime = datetime(2017,7,1)
eDTime = datetime(2017,7,1)
lDTimeDay = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
dprver = 'V06'
dprverfull='V06A'

lvar = ['r','ept','tv','w']
#lvar =['w','vo']
#lvar =['r','tv','ept']
#lvar =['tv','ept','w','r']
#lvar =['tv','ept','r']
lplev = [975,900,825,750,600,450,300,200]
a1p   = array(lplev).reshape(-1,1,1) * units.hPa
lzmeter = [1500,4500,7500]
#lzmeter = [1500]
nyRA   = 721
nxRA   = 1440


myhost = socket.gethostname()
if myhost =='shui':
    dprbaseDir = '/work/hk01/PMM/NASA/GPM.Ku/2A/%s'%(dprver)
    erabaseDir = '/tank/utsumi/era5'
    outbaseDir = '/tank/utsumi/env/pair'
elif myhost == 'well':
    dprbaseDir = '/home/utsumi/mnt/lab_work/hk01/PMM/NASA/GPM.Ku/2A/%s'%(dprver)
    #erabaseDir = '/home/utsumi/mnt/lab_tank/utsumi/era5'
    erabaseDir = '/media/disk2/share/data/era5'
    outbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/env/pair'
else:
    print 'check myhost'
    sys.exit()

miss_out = -9999.
#**********************
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
            a3depth[ilayer][a2flag] = 0
            a3datTmp[ilayer][a2flag] = a2surfdat[a2flag]
        #*** Find lower and upper layer ***
        a2ilow = ma.masked_less(a3depth, hintp).argmin(axis=0) -1
        a2iup  = a2ilow + 1

        #*** Interpolation ***
        a2mask = ma.masked_less(a3depth[-1], hintp).mask
        a2iup[a2mask] = a2ilow[a2mask]   # temporally replace

        #*** lower later *****
        a2depthlow = a3depth[a2ilow,Y,X]-1
        #a2datlow   = a3datTmp[a2ilow,Y,X]-1
        a2datlow   = a3datTmp[a2ilow,Y,X]   # 2019/11/26

        a2flag = ma.masked_equal(a2ilow,-1).mask
        a2depthlow[a2flag] = 0
        a2datlow[a2flag] = a2surfdat[a2flag]

        #*** upper later *****
        a2depthup  = a3depth[a2iup, Y,X]+1
        #a2datup    = a3datTmp[a2iup, Y,X]+1
        a2datup    = a3datTmp[a2iup, Y,X]   # 2019/11/26

        #*** Interpolation ***
        a2datintp = ((hintp-a2depthlow)*a2datup + (a2depthup-hintp)*a2datlow)/(a2depthup-a2depthlow)
        a3out[ihintp] = ma.masked_where(a2mask, a2datintp).filled(miss_out)

    if flag2d ==1:
        a3out = a3out.reshape(len(lh),-1)

    return a3out


def read_var_3d(var,Year,Mon,Day):
    srcDir = erabaseDir + '/%s/%04d%02d'%(var,Year,Mon)
    srcPath= srcDir + '/%s.%04d.%02d.%02d.nc'%(var,Year,Mon,Day)
    with Dataset(srcPath) as np:
        a3var = np.variables[var][:]
    return a3var


def read_var_surf(var,Year,Mon,Day):
    srcDir = erabaseDir + '/%s/%04d%02d'%(var,Year,Mon)
    srcPath= srcDir + '/%s.%04d.%02d.%02d.nc'%(var,Year,Mon,Day)
    ncvar = {'2t':'t2m', '2d':'d2m', 'sp':'sp'}[var]

    with Dataset(srcPath) as np:
        a2var = np.variables[ncvar][:]
    return a2var

def read_zmeter(Year,Mon,Day):
    ''' geopotential height [m] '''
    g = 9.80665
    var = 'z'
    srcDir = erabaseDir + '/%s/%04d%02d'%(var,Year,Mon)
    srcPath= srcDir + '/%s.%04d.%02d.%02d.nc'%(var,Year,Mon,Day)
    with Dataset(srcPath) as np:
        a2var = np.variables[var][:]
    return a2var/g


def read_var_3d_hour(var,Year,Mon,Day,Hour):
    srcDir = erabaseDir + '/%s/%04d%02d'%(var,Year,Mon)
    srcPath= srcDir + '/%s.%04d.%02d.%02d.nc'%(var,Year,Mon,Day)
    with Dataset(srcPath) as np:
        a3var = np.variables[var][Hour]
    return a3var


def read_var_2d_hour(var,Year,Mon,Day,Hour):
    srcDir = erabaseDir + '/%s/%04d%02d'%(var,Year,Mon)
    srcPath= srcDir + '/%s.%04d.%02d.%02d.nc'%(var,Year,Mon,Day)
    ncvar = {'2t':'t2m', '2d':'d2m', 'sp':'sp'}[var]

    with Dataset(srcPath) as np:
        a2var = np.variables[ncvar][Hour]
    return a2var

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
a2orog = read_orogmeter()

#**********************
for DTimeDay in lDTimeDay:
    print DTimeDay
    YearDir,MonDir, DayDir = DTimeDay.timetuple()[:3] 
    dprDir   = dprbaseDir + '/%04d/%02d/%02d'%(YearDir,MonDir,DayDir)
    #ssearch  = dprDir + '/2A.GPM.Ku.V8-20180723.20170101-S003624-E020857.016155.V06A.HDF5'
    ssearch  = dprDir + '/2A.GPM.Ku.*.%s.HDF5'%(dprverfull)
    lsrcPath = sort(glob.glob(ssearch))
    if len(lsrcPath)==0:
        continue
    ##*******************************
    for srcPath in lsrcPath:
        oid  = srcPath.split('.')[-3]

        #if oid < '019277': continue   # test

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

        for var in lvar:
            for dtenv in ldtenv:
                nplev = len(lplev)
                ny,nx = a2lat.shape 
    
                a2var = array([])
                a2zmeter = array([])
                for DTime in lDTimeHour:
    
                    iy0 = bisect_left(lDTime,DTime-timedelta(minutes=30))
                    iy1 = bisect_left(lDTime,DTime+timedelta(minutes=30))
                    if (iy1==0) or (iy0==ny):
                        continue
    
                    Year,Mon,Day,Hour,Mnt = (DTime + timedelta(hours=dtenv)).timetuple()[:5]
    
    
                    #--- Read ERA ******
                    a3zmeter = read_zmeter_hour(Year,Mon,Day,Hour)[::-1,:,:] * units.m
                    if var in ['ept','tv']:
                        a3t = read_var_3d_hour('t', Year,Mon,Day,Hour)[::-1,:,:] * units.K
                        a3q = read_var_3d_hour('q', Year,Mon,Day,Hour)[::-1,:,:] * units.g/units.g  # specific humidity
                        a2tsurf = read_var_2d_hour('2t',Year,Mon,Day,Hour) * units.K
                        a2dsurf = read_var_2d_hour('2d',Year,Mon,Day,Hour) * units.K
                        a2sp    = read_var_2d_hour('sp',Year,Mon,Day,Hour) /100.* units.hPa
    
                    elif var in ['r']:
                        a2tsurf = read_var_2d_hour('2t',Year,Mon,Day,Hour) * units.K
                        a2dsurf = read_var_2d_hour('2d',Year,Mon,Day,Hour) * units.K
    
                    #--- ERA5 ----
                    if var =='ept':  # Equivalent potential temperature
                        #--- 3D ---
                        a3dew=mpcalc.dewpoint_from_specific_humidity(a3q, a3t, a1p)  # degC
                        a3var = mpcalc.equivalent_potential_temperature(a1p, a3t, a3dew)
    
                        #--- Surface -
                        a2surfvar = mpcalc.equivalent_potential_temperature(a2sp, a2tsurf, a2dsurf)
    
                    elif var =='tv':  # Virtual temperature
                        a3w = mpcalc.mixing_ratio_from_specific_humidity(a3q)  # mixing ratio
                        a3var = mpcalc.virtual_temperature(a3t,  a3w)  # virtual temperature
    
                        #--- Surface -
                        #a2surfvar = mpcalc.equivalent_potential_temperature(a2sp, a2tsurf, a2dsurf) # comment out @ 2019/9/29
                        a2rsurf   = mpcalc.relative_humidity_from_dewpoint(a2tsurf, a2dsurf)
                        a2wsurf   = mpcalc.mixing_ratio_from_relative_humidity(a2rsurf, a2tsurf, a2sp)
                        a2surfvar = mpcalc.virtual_temperature(a2tsurf, a2wsurf)  # Add @ 2019/9/29
    
                    elif var=='r':
                        #a3var = mpcalc.relative_humidity_from_specific_humidity(a3q, a3t, a1p)
                        a3var = read_var_3d_hour(var, Year,Mon,Day,Hour)[::-1,:,:] * units.g/units.g  # specific humidity
                        
                        a2surfvar = mpcalc.relative_humidity_from_dewpoint(a2tsurf, a2dsurf)
    
                        ##-- test ---
                        #nall = 1038240.
                        #nmore= ma.masked_less_equal(a3var[2],1).count()
                        #nless= ma.masked_greater_equal(a3var[2],0).count()
                        #print a3var[2].min(),a3var[2].max(),'frac=',nless/nall, nmore/nall
    
                        #_,ny,nx = a3var.shape
                        #for y in range(ny):
                        #    for x in range(nx):
                        #        if a3var[2,y,x] > 1.2:
                        #            print y,x,a3var[2,y,x],a3q[2,y,x],a3t[2,y,x]
                        #sys.exit()
                        #-----------
                         
                    elif var=='w':
                        a3var = read_var_3d_hour('w', Year,Mon,Day,Hour)[::-1,:,:]
                        a2surfvar = np.zeros([nyRA,nxRA], float32)
    
                    else:
                        print 'check var',var
                        sys.exit()
    
                    
                    #--- corresponding pixels --
                    latRA0= 90   # from North to South
                    lonRA0= 0
                    dlatRA = 0.25
                    dlonRA = 0.25
            
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
    
                    a2varTmp  = a3var[:,a1y,a1x]
                    a2zmeterTmp = a3zmeter[:,a1y,a1x]
                    a1surfvarTmp= a2surfvar[a1y,a1x]
                    a1orogTmp   = a2orog[a1y,a1x]
    
                    if a2var.shape[0]==0: 
                        a2var    = a2varTmp
                        a2zmeter = a2zmeterTmp
                        a1surfvar= a1surfvarTmp
                        a1orog   = a1orogTmp
                    else:
                        a2var = concatenate([a2var, a2varTmp],axis=1)
                        a2zmeter=concatenate([a2zmeter, a2zmeterTmp],axis=1)
                        a1surfvar=concatenate([a1surfvar, a1surfvarTmp])
                        a1orog  =concatenate([a1orog, a1orogTmp])
    
                #*** Interpolation at zmeter ***
                a2intp = interp_vertical(a2zmeter, a2var, a1surfvar, a1orog, lzmeter)
    
                #**** If noscreen ********************
                if noscreen is True:
                    nydpr,nxdpr = a2prec.shape
                    a2intp = a2intp.reshape(len(lzmeter),nydpr,nxdpr)
                    a1surfvar=a1surfvar.reshape(nydpr,nxdpr)
    
    
                #**** Head for file name ***
                shead = ''
                if noscreen is True:
                    shead = shead+'full.'
    
                #**** Save interpolated variables ****
                outDir     = outbaseDir + '/%s/%04d/%02d/%02d'%(var,YearDir,MonDir,DayDir)
                util.mk_dir(outDir)
                for izmeter, zmeter in enumerate(lzmeter):
                    a1intp   = a2intp[izmeter]
    
                    outPath = outDir + '/%s%s.%03dh.%04.1fkm.%s.npy'%(shead, var, dtenv, zmeter*0.001, oid)
                    np.save(outPath, a1intp.astype(float32))
                print outPath
                #**** Save surface variable ****
    
                outPath = outDir + '/%s%s.%03dh.00.0km.%s.npy'%(shead, var, dtenv, oid)
                np.save(outPath, a1surfvar.astype(float32))
                print outPath 
        
