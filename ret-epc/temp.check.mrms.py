import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from numpy import *
import numpy as np
import h5py

mrmsPath = '/work/hk01/PMM/MRMS/match-GMI-orbit/GMI.MRMS.130W_55W_20N_55N.20140516.001204.00931-01219.npy'

gprofPath= '/work/hk01/PMM/NASA/GPM.GMI/2A/V05/2014/05/16/2A.GPM.GMI.GPROF2017v1.20140516-S005738-E023012.001204.V05A.HDF5'

a2mrms = np.load(mrmsPath)

iyey = mrmsPath.split('.')[-2]
iy,ey = map(int, iyey.split('-'))
sdate = mrmsPath.split('.')[-4]
oid   = int(mrmsPath.split('.')[-3])
scase = '%s %06d'%(sdate,oid)
#***************
# Read GPROF
#***************
with h5py.File(gprofPath) as h:
    a2gprof = h['S1/surfacePrecipitation'][:]
    a2lat   = h['S1/Latitude'][:]
    a2lon   = h['S1/Longitude'][:]

a2gprof = a2gprof[iy:ey+1]
a2lat   = a2lat[iy:ey+1]
a2lon   = a2lon[iy:ey+1]
#***************
# Read MRMS
#***************
a2mrms  = np.load(mrmsPath)

print a2gprof.shape, a2mrms.shape

#***************
# Figure
#***************
ssize = 1
[[lllat,lllon],[urlat,urlon]]=[[26,-100],[50,-70]]

fig = plt.figure(figsize=(10,10))
#-------------
# Scatter
#-------------
ax  = fig.add_axes([0.1,0.1,0.35,0.35])
xdat= ma.masked_less(a2mrms,0)
ydat= ma.masked_less(a2gprof,0)
ax.scatter(xdat,ydat,c='k',s=1)
plt.ylim([0,20])
plt.xlim([0,20])
plt.xlabel('MRMS [mm/h]')
plt.ylabel('GPROF [mm/h]')

#-------------
# MRMS
#-------------
ax  = fig.add_axes([0.1,0.5,0.35,0.35])
a2fig = ma.masked_less(a2mrms,0)
M = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax)
im = M.scatter(a2lon, a2lat, c=a2fig, cmap='jet', s=ssize, vmin=0, vmax=10)

dgrid      = 5
parallels  = arange(-90,90, dgrid)
meridians  = arange(-180,180,dgrid)
M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8, linewidth=0.5, fmt='%d')
M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8, linewidth=0.5, fmt='%d')
M.drawcoastlines()
plt.colorbar(im, orientation='horizontal')
plt.title('MRMS'+ ' ' + scase)

#-------------
# GPROF
#-------------
ax  = fig.add_axes([0.5,0.5,0.35,0.35])
a2fig = ma.masked_less(a2gprof,0)
M = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax)
im = M.scatter(a2lon, a2lat, c=a2fig, cmap='jet', s=ssize, vmin=0, vmax=10)

dgrid      = 5
parallels  = arange(-90,90, dgrid)
meridians  = arange(-180,180,dgrid)
M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8, linewidth=0.5, fmt='%d')
M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8, linewidth=0.5, fmt='%d')
M.drawcoastlines()
plt.colorbar(im, orientation='horizontal')
plt.title('GPROF' + ' ' + scase)

#-------------
# Save
#-------------

figPath = '/home/utsumi/temp/mrms.png'
plt.savefig(figPath)
print figPath
plt.clf()
