import h5py
from numpy import *


srcPath = '/home/utsumi/temp/EPC/GPM_EPC_002421_20140802_0726.NS_MS.nc'
with h5py.File(srcPath) as h:
    a2precip = h['MS/precip'][:]
    a2lat    = h['latitude'][:]
    a2lon    = h['longitude'][:]

ny,nx = a2precip.shape
for iy in range(ny):
    for ix in range(nx):
        lat=a2lat[iy,ix]
        lon=a2lon[iy,ix]
        precip = a2precip[iy,ix]
        print iy,ix,lat,lon, '%.2f'%precip
