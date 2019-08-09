from numpy import *
import numpy as np
import myfunc.util as util
import netCDF4
import socket

myhost = socket.gethostname()
if myhost =='shui':
    erabaseDir = '/tank/utsumi/era5'
    outDir = '/home/utsumi/temp/env'

elif myhost == 'well':
    erabaseDir = '/media/disk2/share/data/era5'
    outDir = '/home/utsumi/temp/env'

else:
    print 'check myhost'
    sys.exit()




def read_var_3d(var,Year,Mon,Day):
    srcDir = erabaseDir + '/%s/%04d%02d'%(var,Year,Mon)
    srcPath= srcDir + '/%s.%04d.%02d.%02d.nc'%(var,Year,Mon,Day)
    with netCDF4.Dataset(srcPath) as np:
        a3var = np.variables[var][:]
    return a3var


def read_var_surf(var,Year,Mon,Day):
    srcDir = erabaseDir + '/%s/%04d%02d'%(var,Year,Mon)
    srcPath= srcDir + '/%s.%04d.%02d.%02d.nc'%(var,Year,Mon,Day)
    ncvar = {'2t':'t2m', '2d':'d2m', 'sp':'sp'}[var]

    with netCDF4.Dataset(srcPath) as np:
        a2var = np.variables[ncvar][:]
    return a2var

def read_zmeter(Year,Mon,Day):
    ''' geopotential height [m] '''
    g = 9.80665
    var = 'z'
    srcDir = erabaseDir + '/%s/%04d%02d'%(var,Year,Mon)
    srcPath= srcDir + '/%s.%04d.%02d.%02d.nc'%(var,Year,Mon,Day)
    with netCDF4.Dataset(srcPath) as np:
        a2var = np.variables[var][:]
    return a2var/g

def read_orogmeter():
    return np.load(erabaseDir + '/orog/orog.meter.na.npy')

#***********************************************
Year,Mon,Day =2017,7,3
util.mk_dir(outDir)


a2orog  = read_orogmeter()
a4zmeter= read_zmeter(Year,Mon,Day)
a3zmeter= a4zmeter[0]

#**** orogmeter ***
outPath = outDir + '/orogmeter.npy'
np.save(outPath, a2orog.data)
print outPath

#**** zmeter ***
outPath = outDir + '/zmeter.npy'
np.save(outPath, a3zmeter.data)
print outPath

#**** 3D ****
for var in ['t','q']:
    a3var = read_var_3d(var,Year,Mon,Day)[0]
    outPath = outDir + '/%s.npy'%(var)
    np.save(outPath, a3var.data)
    print outPath
#**** 2D ****
for var in ['2t','2d','sp']:
    a2var = read_var_surf(var,Year,Mon,Day)[0]
    outPath = outDir + '/%s.npy'%(var)
    np.save(outPath, a2var.data)
    print outPath
