from numpy import *
import myfunc.util as util
from datetime import datetime, timedelta
import h5py
import os, sys, socket, glob
import numpy as np

iDTime = datetime(2017,7,1)
eDTime = datetime(2017,7,31)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
#lvar = ['Latitude','Longitude']
#lvar = ['convfrac','stratfrac','otherfrac','Latitude','Longitude']
#lvar = ['landSurfaceType']
lvar = ['heightStormTop']

noscreen = False
dprver = 'V06'
dprverfull='V06A'
myhost = socket.gethostname()
if myhost =='shui':
    dprbaseDir = '/work/hk01/PMM/NASA/GPM.Ku/2A/%s'%(dprver)
    erabaseDir = '/tank/utsumi/era5'
    outbaseDir = '/tank/utsumi/env/pair'
elif myhost == 'well':
    dprbaseDir = '/home/utsumi/mnt/lab_work/hk01/PMM/NASA/GPM.Ku/2A/%s'%(dprver)
    erabaseDir = '/media/disk2/share/data/era5'
    outbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/env/pair'
else:
    print 'check myhost'
    sys.exit()

miss_out = -9999.


#******************************
#def ret_convfrac(a2type):
#    a2conv = ma.masked_not_equal((a2type/10000000).astype(int32), 2)

def sum_9grids(a2in, miss_in=-9999., miss_out=-9999.):
    ny,nx = a2in.shape
    a2tmp = np.ones([ny+2,nx+2],float32)*miss_in
    a2tmp[1:-1,1:-1] = a2in

    a3tmp = np.ones([9,ny,nx],float32)*miss_in
    a3tmp[0] = a2tmp[:-2,:-2]
    a3tmp[1] = a2tmp[:-2,1:-1]
    a3tmp[2] = a2tmp[:-2,2:]
    a3tmp[3] = a2tmp[1:-1,:-2]
    a3tmp[4] = a2tmp[1:-1,1:-1]
    a3tmp[5] = a2tmp[1:-1,2:]
    a3tmp[6] = a2tmp[2:,:-2]
    a3tmp[7] = a2tmp[2:,1:-1]
    a3tmp[8] = a2tmp[2:,2:]

    a2sum = ma.masked_equal(a3tmp,miss_in).sum(axis=0)
    a2sum = a2sum.filled(miss_out)
    return a2sum

#******************************
for DTimeDay in lDTime:
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
            a2prec= h5['/NS/SLV/precipRateNearSurface'][:]
            a2mask = ma.masked_less_equal(a2prec,0).mask

        for var in lvar:
            with h5py.File(srcPath,'r') as h5:
                if var in ['Latitude','Longitude']:
                    hdfvar ='/NS/%s'%(var)
                    a2var = h5[hdfvar][:]

                elif var in ['landSurfaceType','heightStormTop']:
                    hdfvar = '/NS/PRE/%s'%(var)
                    a2var = h5[hdfvar][:]

                elif var in ['stratfrac','convfrac','otherfrac']:
                    hdfvar = '/NS/CSF/typePrecip'
                    a2type = h5[hdfvar][:]
                    a2type = (a2type/10000000).astype(int32)
                    ny,nx  = a2type.shape
                    a2wet  = ma.masked_where(a2type<0,   np.ones([ny,nx],float32)).filled(0.0)

                    if   var=='stratfrac':
                        a2bit = ma.masked_where(a2type !=1, np.ones([ny,nx],float32)).filled(0.0)

                    elif var=='convfrac':
                        a2bit = ma.masked_where(a2type !=2, np.ones([ny,nx],float32)).filled(0.0)

                    elif var=='otherfrac':
                        a2bit = ma.masked_where(a2type !=3, np.ones([ny,nx],float32)).filled(0.0)

                    a2sumbit = sum_9grids(a2bit)
                    a2sumwet = sum_9grids(a2wet)
                    a2var = (ma.masked_where(a2sumwet==0, a2sumbit) / a2sumwet)
                
                else:
                    print 'check var',var
                    sys.exit() 
   

            #-- Mask --- 
            if noscreen is True:
                a1var = a2var.flatten()
            elif noscreen is False:
                a1var = ma.masked_where(a2mask, a2var).compressed()
            else:
                print 'noscreen',noscreen
                sys.exit()

            #**** Save interpolated variables ****
            outDir     = outbaseDir + '/%s/%04d/%02d/%02d'%(var,YearDir,MonDir,DayDir)
            util.mk_dir(outDir)
            #**** Save variable ****
            if noscreen is True:
                outPath = outDir + '/full.%s.00.0km.%s.npy'%(var, oid)
            else:
                outPath = outDir + '/%s.00.0km.%s.npy'%(var, oid)
            np.save(outPath, a1var.astype(float32))
            print outPath

 




