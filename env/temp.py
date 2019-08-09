import numpy as np
from  netCDF4 import Dataset

var = 't'
Year,Mon,Day = 2017,7,3
erabaseDir = '/home/utsumi/mnt/lab_tank/utsumi/era5'
srcDir = erabaseDir + '/%s/%04d%02d'%(var,Year,Mon)
srcPath= srcDir + '/%s.%04d.%02d.%02d.nc'%(var,Year,Mon,Day)


for i in range(5):
    print i
    with Dataset(srcPath,'r') as nc:
        dat= nc.variables['t']
        print dat.shape
    
    with Dataset(srcPath,'r') as nc:
        dat= nc.variables['t'][4]
        print dat.shape


