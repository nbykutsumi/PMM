import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from numpy import *
import h5py
from bisect import bisect_left
import epcfunc
import sys

#clat    = 30.00
#clon    = 269.0 -360  # -180 - +180
clat    = 14
clon    = 2  # -180 - +180

dlatlon = 6
#dscan   = 55
dscan   = 25
BBox    = [[clat-dlatlon, clon-dlatlon],[clat+dlatlon,clon+dlatlon]]
[[lllat,lllon],[urlat,urlon]] = BBox
miss = -9999.

#oid = 16166
#iy, ey = 1952, 2012
oid = 2421
iy, ey = 2004, 2114
maxrec = 20000
srcDir = '/home/utsumi/temp/out'
#prcpPath = srcDir + '/esurf.%06d.y%04d-%04d.npy'%(oid, iy, ey)
#prcpPath = srcDir + '/esurf.%06d.y%04d-%04d.20000maxrec.npy'%(oid, iy, ey)
latPath  = srcDir + '/lat.%06d.y%04d-%04d.%dmaxrec.npy'%(oid, iy, ey, maxrec)
lonPath  = srcDir + '/lon.%06d.y%04d-%04d.%dmaxrec.npy'%(oid, iy, ey, maxrec)
idxPath  = srcDir + '/temp.idx.npy'

jplPath  = '/home/utsumi/bin/PMM/ret-epc/GPM_EPC_002421_20140802_0726.NS_MS.nc'  # HDF file

tbPath   = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05/2014/08/02/1C.GPM.GMI.XCAL2016-C.20140802-S062222-E075455.002421.V05A.HDF5'

#*****************
#- Read My data ----

#a2dat = np.load(prcpPath)
a2idx = np.load(idxPath)
a2lat = np.load(latPath)
a2lon = np.load(lonPath)

a2idx, a2lat, a2lon = epcfunc.extract_domain_2D(a2idx, a2lat, a2lon, clat, clon, dlatlon, dscan)

#a2dat = ma.masked_less_equal(a2dat,0)
a2idx = ma.masked_less_equal(a2idx,0)

#*****************
#- Read JPL data ----
with h5py.File(jplPath) as h:
    a2esurfjpl = h['MS/precip'][:]
    a2dpr      = h['DPR/precip_NS'][:]
    a2latjpl   = h['latitude'][:]
    a2lonjpl   = h['longitude'][:]
    a2idxjpl   = h['db_index'][:]

#a2esurfjpl, a2latjpl, a2lonjpl = epcfunc.extract_domain_2D(a2esurfjpl, a2latjpl, a2lonjpl, clat, clon, dlatlon, dscan)
#a2esurfjpl = ma.masked_less_equal(a2esurfjpl,0)
#a2dpr, a2latjpl, a2lonjpl = epcfunc.extract_domain_2D(a2dpr, a2latjpl, a2lonjpl, clat, clon, dlatlon, dscan)
#a2dpr = ma.masked_less_equal(a2dpr,0)
a2idxjpl, a2latjpl, a2lonjpl = epcfunc.extract_domain_2D(a2idxjpl, a2latjpl, a2lonjpl, clat, clon, dlatlon, dscan)

#-- test ---
np.save('/home/utsumi/temp/out/idx-jpl.npy', a2idxjpl)
ilat = a2latjpl[0,110]
elat = a2latjpl[-1,110]
iymy = bisect_left(a2lat[:,110], ilat)
eymy = bisect_left(a2lat[:,110], elat)
print ilat, a2lat[iymy,110]
print ilon, a2lon[eymy,110]
sys.exit()
#-----------
a2idxjpl = ma.masked_less_equal(a2idxjpl,0)


##-- test --------------------
#plt.scatter(a2idxjpl,a2idx)
#plt.xlabel('idx_db JPL')
#plt.ylabel('idx_db NU')
#plt.title('database index')
#plt.savefig('/home/utsumi/temp/out/scatter.idx.png')
#plt.clf()
#sys.exit()


#********************************
#-- Draw figure ---
fig   = plt.figure(figsize=(8,8))
ssize = 1

#-- My retrieval --
#ax    = fig.add_subplot(111)
ax    = fig.add_axes([0.1,0.5,0.3,0.3])
M     = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax)
im    = M.scatter(a2lon, a2lat, c=a2idx, cmap='jet', s=ssize, vmin=0, vmax=15624)
M.drawcoastlines()

dgrid      = 5
parallels  = arange(-90,90, dgrid)
meridians  = arange(-180,180,dgrid)
M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8, linewidth=0.5, fmt='%d')
M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8, linewidth=0.5, fmt='%d')
plt.colorbar(im, orientation='horizontal')

plt.title('db_index(NU)')
#*****************
#- Figure JPL ----
ax2 = fig.add_axes([0.1,0.1,0.3,0.3])
M   = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax2)
im    = M.scatter(a2lonjpl, a2latjpl, c=a2idxjpl, cmap='jet', s=ssize, vmin=0, vmax=15624)
M.drawcoastlines()
M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8, linewidth=0.5, fmt='%d')
M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8, linewidth=0.5, fmt='%d')
plt.colorbar(im, orientation='horizontal')

plt.title('db_index(JPL)')
#------------
#- Figure Difference ----
#a2dif = ma.masked_invalid(abs(a2idx-a2idxjpl)/a2idxjpl)
a2dif = ma.masked_invalid(abs(a2idx-a2idxjpl)/a2idxjpl.astype(float32))
a2dif = ma.masked_equal(a2dif,0)

ax2 = fig.add_axes([0.5,0.1,0.3,0.3])
M   = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax2)
im  = M.scatter(a2lonjpl, a2latjpl, c=a2dif, cmap='jet', s=ssize)
M.drawcoastlines()
M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8, linewidth=0.5, fmt='%d')
M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8, linewidth=0.5, fmt='%d')
plt.colorbar(im, orientation='horizontal')

plt.title('Difference abs(NU-JPL)/JPL')

#------------


outPath  = srcDir + '/db_index.png'
plt.savefig(outPath)
print outPath
