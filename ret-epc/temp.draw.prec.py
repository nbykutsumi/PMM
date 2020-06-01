# %%
%matplotlib inline
import numpy as np
from numpy import *
import glob, h5py
import myfunc.util as util
from datetime import datetime, timedelta
import EPCDB
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap



#BBox   = [[25,140],[45,160]]  # GCOMW1.AMSR2
#BBox   = [[-10,0],[10,30]]  # F16.SSMIS
#BBox   = [[-10,90],[20,120]]  # METOPA.MHS
BBox   = [[20,-130],[55,-55]]  # METOPA.MHS
sate   = 'NOAA20'
sensor = 'ATMS'
slabel = '000623.y1470-1569.nrec10000'
#srcdir = '/home/utsumi/temp/out/JPL/%s.%s'%(sensor,sate)
srcdir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/test/%s.%s/2018/01/01/'%(sensor,sate)
precpath= srcdir + '/nsurfNScmb.%s.npy'%(slabel)
latpath = srcdir + '/lat.%s.npy'%(slabel)
lonpath = srcdir + '/lon.%s.npy'%(slabel)
profpath= srcdir + '/prwatprofNS-rs.%s.npy'%(slabel)



a2prec = np.load(precpath)
a2lat  = np.load(latpath)
a2lon  = np.load(lonpath)
a3prof = np.load(profpath)

a2dat  = ma.masked_less(a2prec,0.1)
[[lllat,lllon],[urlat,urlon]] = BBox

fig = plt.figure(figsize=(8,8))
ax  = fig.add_axes([0.2,0.2, 0.7,0.7])
M   = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax)
im  = M.scatter(a2lon, a2lat, c=a2dat, s=2, vmax=10, cmap='jet')
plt.title('%s.%s %s'%(sate,sensor, slabel))
M.drawcoastlines()
M.drawparallels(np.arange(-90,90,2), labels=[1,0,0,0], linewidth=0.5, color='0.8')
M.drawmeridians(np.arange(-180,180,2), labels=[0,0,0,1], linewidth=0.5, color='0.8', rotation=30)
M.plot(-173.4, 8, 'x', color='k')
plt.colorbar(im)
plt.show()
# %%
