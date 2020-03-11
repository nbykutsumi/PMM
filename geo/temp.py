# %%
%matplotlib inline
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os, sys, glob
import h5py
from mpl_toolkits.basemap import Basemap

BBox   = [[-55,90],[55,200]] # Himawari
iabin = 103  # extracted GMI angle bin (start)
eabin = 117  # extracted GMI angle bin (end)
workbaseDir = '/home/utsumi/mnt/lab_work'
baseDir = workbaseDir + '/hk02/PMM/NASA/GPM.GMI/1C/V05'
Year,Mon,Day = 2017,9,1
oid = 19945
srcDir = baseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
srcPath = srcDir + '/1C.GPM.GMI.XCAL2016-C.20170901-S160609-E173843.019945.V05A.HDF5'
print srcPath
with h5py.File(srcPath, 'r') as h:
    a2lat = h['S1/Latitude'][:,iabin:eabin+1]   # 25 angle bins at the center (Same as Ka)
    a2lon = h['S1/Longitude'][:,iabin:eabin+1]  # Ku: 49 range bins. Ka: 25 range bins

iy,ey = 250, 1238
a2lat = a2lat[iy:ey+1]
a2lon = a2lon[iy:ey+1]


[[lllat,lllon],[urlat,urlon]] = BBox
mycm  = 'jet'
ssize = 2
fig = plt.figure(figsize=(6,6))
ax  = fig.add_axes([0.1,0.1,0.8,0.8])
M     = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon)

im    = M.scatter(a2lon, a2lat, c=a2lat, cmap=mycm, s=ssize, ax=ax)
M.drawcoastlines()
plt.show()
print a2lon.shape, a2lat.shape

plt.plot(np.arange(10))
plt.show()

# %%
