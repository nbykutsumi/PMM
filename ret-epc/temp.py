import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import h5py
#import netCDF4
srcPath='/home/utsumi/bin/PMM/ret-epc/GPM_EPC_012149_20160418_1228.NS_MS.nc'
#nc = netCDF4.Dataset(srcPath,'r')
#print nc.variables

with h5py.File(srcPath) as h:
    print h.items()
    a2lat0 = h['/latitude'][:]
    a2lon0 = h['/longitude'][:]
    a2lat1 = h['/DPR/latitude'][:]
    a2lon1 = h['/DPR/longitude'][:]


miss = -9999.
a2lat0 = ma.masked_equal(a2lat0,miss)
a2lat1 = ma.masked_equal(a2lat1,miss)
a2lon0 = ma.masked_equal(a2lon0,miss)
a2lon1 = ma.masked_equal(a2lon1,miss)

plt.scatter(a2lat0,a2lat1)
plt.savefig('/home/utsumi/temp/out/lat.png')
plt.clf()
plt.scatter(a2lon0,a2lon1)
plt.savefig('/home/utsumi/temp/out/lon.png')
plt.clf()
