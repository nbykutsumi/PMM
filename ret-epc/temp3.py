from numpy import *
import netCDF4

srcDir = '/home/utsumi/bin/PMM/ret-epc'
srcPath= srcDir + '/test_netcdf.nc'
nc = netCDF4.Dataset(srcPath,'r')
print nc
