import matplotlib
matplotlib.use('Agg')
from numpy import *
import h5py
import numpy as np
import matplotlib.pyplot as plt


elevPath = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.gtopo.016155/2017/01/01/gtopo.016155.npy'

xPath = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.IDX/2017/01/01/Xpy.1.016155.npy'
yPath = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.IDX/2017/01/01/Ypy.1.016155.npy'


dprPath = '/work/hk01/PMM/NASA/GPM.Ku/2A/V06/2017/01/01/2A.GPM.Ku.V8-20180723.20170101-S003624-E020857.016155.V06A.HDF5'


a2elev1 = np.load(elevPath)[:,83:137+1]
a1elev1= a2elev1.flatten()

a1x     = np.load(xPath).flatten()
a1y     = np.load(yPath).flatten()
a1mask  = ma.masked_less(a1x,0).mask

a1elev1 = ma.masked_where(a1mask, a1elev1).compressed()
a1x     = ma.masked_where(a1mask, a1x).compressed()
a1y     = ma.masked_where(a1mask, a1y).compressed()

with h5py.File(dprPath) as h:
    a2elev2 = h['NS/PRE/elevation'][:]
    a2lat   = h['NS/Latitude'][:]
    a2lon   = h['NS/Longitude'][:]

a1elev2 = a2elev2[a1y,a1x]
a1lat   = a2lat[a1y,a1x]
a1lon   = a2lon[a1y,a1x]

#-- Figure 1 --
a1mask0 = ma.masked_outside(a1elev1, -5,5).mask
a1mask1 = ma.masked_less(a1elev2, 300).mask
a1mask  = a1mask0 + a1mask1
a1x     = ma.masked_where(a1mask, a1lon).compressed()
a1y     = ma.masked_where(a1mask, a1lat).compressed()
print a1x
print a1y
figPath = '/home/utsumi/temp/map.evel1.png'
plt.scatter(a1x, a1y)
plt.savefig(figPath)
print figPath
plt.clf()
#-- Figure 2 --
a1mask0 = ma.masked_greater(a1elev1, -1000).mask
a1mask1 = ma.masked_less(a1elev2, 300).mask
a1mask  = a1mask0 + a1mask1
a1x     = ma.masked_where(a1mask, a1lon).compressed()
a1y     = ma.masked_where(a1mask, a1lat).compressed()
print a1x
print a1y
figPath = '/home/utsumi/temp/map.evel2.png'
plt.scatter(a1x, a1y)
plt.savefig(figPath)
print figPath




#for i in range(len(a1elev1)):
#    d = abs(a1elev2[i] - a1elev1[i])
#    if d>200:
#        print a1lat[i],a1lon[i], a1elev2[i], a1elev1[i]
#    
#print a1elev1.shape, a1elev2.shape
#
#figPath = '/home/utsumi/temp/evel.png'
#plt.scatter(a1elev2, a1elev1)
#plt.savefig(figPath)
#print figPath
