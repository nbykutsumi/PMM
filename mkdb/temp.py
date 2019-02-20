import matplotlib
matplotlib.use('Agg')
from numpy import *
from f_match_fov import *
import h5py
import numpy as np
import matplotlib.pyplot as plt

dprPath = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.9ave.precipRateESurface/2017/01/01/precipRateESurface.016155.npy'

gmiPath = '/work/hk01/PMM/NASA/GPM.GMI/2A/V05/2017/01/01/2A.GPM.GMI.GPROF2017v1.20170101-S003624-E020857.016155.V05A.HDF5'

dprorgPath= '/work/hk01/PMM/NASA/GPM.Ku/2A/V06/2017/01/01/2A.GPM.Ku.V8-20180723.20170101-S003624-E020857.016155.V06A.HDF5'

a2dpr = np.load(dprPath)
h=h5py.File(gmiPath,'r')
a2gmi = h['S1/surfacePrecipitation'][:,83:137+1]

h=h5py.File(dprorgPath,'r')
a2dprorg=h['NS/SLV/precipRateESurface'][:]


print a2gmi.shape, a2dpr.shape

a2dpr = ma.masked_less(a2dpr,0)
a2gmi = ma.masked_less(a2gmi,0)

#plt.plot(a2dpr, a2gmi,'o')
#plt.savefig('/home/utsumi/temp/temp.png')

#ny,nx  = a2gmi.shape
#a2x, a2y = meshgrid(range(nx),range(ny))
#
#a1x = ma.masked_where(a2gmi<=0, a2x).compressed()
#a1y = ma.masked_where(a2gmi<=0, a2y).compressed()
#for i in range(len(a1x)):
#    print a1y[i], a1x[i]

xPath = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.IDX/2017/01/01/Xpy.1.016155.npy'
yPath = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.IDX/2017/01/01/Ypy.1.016155.npy'

a2x = np.load(xPath)
a2y = np.load(yPath)

y = 2834
x = 20

ydpr = a2y[y,x]
xdpr = a2x[y,x]

ave0 = a2dpr[y,x]

ave1 = 0
ldydx = [[dy,dx] for dy in [-1,0,1] for dx in [-1,0,1]]
for (dy,dx) in ldydx:
    dprorg = a2dprorg[ydpr+dy, xdpr+dx]
    print dprorg
    ave1 = ave1 + dprorg
ave1 = ave1/9.
print ''
print a2dpr[y,x], ave1
