import matplotlib
matplotlib.use('Agg')
from numpy import *
import netCDF4
import numpy as np
from mpl_toolkits.basemap import Basemap
import h5py
import glob
import matplotlib.pyplot as plt

var = 'cape'
oid = 19015
Year,Mon,Day,Hour = 2017,7,3,22
zmeter = 0  # km
exvarPath = '/home/utsumi/mnt/lab_tank/utsumi/env/pair/%04d/%02d/%02d/full.%s.%04.1fkm.%06d.npy'%(Year,Mon,Day,var,zmeter*0.001,oid)

a2exvar = np.load(exvarPath)
a2exvar = ma.masked_less(a2exvar,0)
print exvarPath
print a2exvar
print a2exvar.min(),a2exvar.mean(),a2exvar.max()
dprbaseDir = '/home/utsumi/mnt/lab_work/hk01/PMM/NASA/GPM.Ku/2A/%s'%('V06')
dprDir   = dprbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
ssearch  = dprDir + '/2A.GPM.Ku.*.%06d.%s.HDF5'%(oid,'V06A')
dprPath  = glob.glob(ssearch)[0]

with h5py.File(dprPath,'r') as h5:
    a2prec= h5['NS/SLV/precipRateNearSurface'][:]
    a2lat = h5['/NS/Latitude'][:]
    a2lon = h5['/NS/Longitude'][:]
    a1year  = h5['/NS/ScanTime/Year'        ][:]
    a1mon   = h5['/NS/ScanTime/Month'       ][:]
    a1day   = h5['/NS/ScanTime/DayOfMonth'  ][:]
    a1hour  = h5['/NS/ScanTime/Hour'        ][:]
    a1mnt   = h5['/NS/ScanTime/Minute'      ][:]
    

#---- draw -------
lllat = 25
lllon = 120
urlat = 45
urlon = 140

fig = plt.figure(figsize=(8,4))

ax  = fig.add_axes([0.1,0.1,0.4,0.8])
M = Basemap(llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax, resolution='l')
im = M.scatter(a2lon, a2lat, c=a2exvar,vmin=0,vmax=200)
plt.colorbar(im, orientation='horizontal')

M.drawcoastlines()

ax  = fig.add_axes([0.5,0.1,0.4,0.8])
M = Basemap(llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax, resolution='l')
im = M.scatter(a2lon, a2lat, c=a2prec, vmin=0, vmax=10)
plt.colorbar(im, orientation='horizontal')
M.drawcoastlines()

figPath = '/home/utsumi/temp/env/orbit.era.%s.png'%(var)
plt.savefig(figPath)
print figPath


#eraDir = '/tank/utsumi/era5/%s/%04d%02d'%('2q',Year,Mon)
#print dprPath


