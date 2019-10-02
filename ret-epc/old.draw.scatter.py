import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from numpy import *
import h5py
from bisect import bisect_left

clat    = 30.00
clon    = 269.0 -360  # -180 - +180
dlatlon = 6
dscan   = 30
BBox    = [[clat-dlatlon, clon-dlatlon],[clat+dlatlon,clon+dlatlon]]
[[lllat,lllon],[urlat,urlon]] = BBox


oid = 2421
iy, ey = 2034, 2084
#********************************
#-- My retrieval --
srcDir = '/home/utsumi/temp/out'
prcpPath = srcDir + '/esurf.%06d.y%04d-%04d.replace.npy'%(oid, iy, ey)
latPath  = srcDir + '/lat.%06d.y%04d-%04d.replace.npy'%(oid, iy, ey)
lonPath  = srcDir + '/lon.%06d.y%04d-%04d.replace.npy'%(oid, iy, ey)

a2epc = np.load(prcpPath)
a2epc = ma.masked_less(a2epc,0)

#------------
#- Read GPROF ----
gpDir = '/work/hk01/PMM/NASA/GPM.GMI/2A/V05/2017/01/01/2A.GPM.GMI.GPROF2017v1.20170101-S173441-E190714.016166.V05A.HDF5'
with h5py.File(gpDir) as h:
    a2gprof = h['S1/surfacePrecipitation'][:]

a2gprof = a2gprof[iy:ey+1]
a2gprof = ma.masked_less(a2gprof,0)

#------------
#- Read DPR 9grid mean ----
dprmeanPath = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.9ave.precipRateESurface/2017/01/01/precipRateESurface.016166.npy'
a2dpr   = np.load(dprmeanPath)
a2dpr   = a2dpr[iy:ey+1]
a2dpr  = ma.masked_less(a2dpr,0)


#*****************
#-- Draw figure ---
fig   = plt.figure(figsize=(10,4))
vmax  = 40
ssize = 1
#-- DPR(ave) vs EPC
xdat  = a2dpr
ydat  = ma.masked_less(a2epc[:,83:137+1], 0)
ax    = fig.add_axes([0.1,0.1,0.3,0.8])
im    = ax.scatter(xdat, ydat, color='k',s=ssize)
ax.set_xlabel('DPR(9-pixel ave.)')
ax.set_ylabel('EPC')
ax.set_ylim([0,vmax])
ax.set_xlim([0,vmax])

#-- GPROF vs EPC
xdat  = a2gprof
ydat  = a2epc

ax    = fig.add_axes([0.5,0.1,0.3,0.8])
im    = ax.scatter(xdat, ydat, color='k',s=ssize)
ax.set_xlabel('GPROF')
ax.set_ylabel('EPC')
ax.set_ylim([0,vmax])
ax.set_xlim([0,vmax])
#-- Save --
figPath = '/home/utsumi/temp/out/scatter.png'
plt.savefig(figPath)
print figPath
