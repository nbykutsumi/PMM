import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from numpy import *
import h5py
from bisect import bisect_left
import epcfunc

#clat    = 30.00
#clon    = 269.0 -360  # -180 - +180

# QJRMS case, oid=012149, 2016/4/18
oid = 12149
iy, ey = 1038, 1098
clat    = 32. # QJRMS case, oid=012149
clon    = -94 # -180 - +180

# Africa case
#oid = 2421
#iy, ey = 2029, 2089
#clat    = 14 # Africa case
#clon    = 2  # -180 - +180

dlatlon = 6
#dscan   = 55
dscan   = 30
BBox    = [[clat-dlatlon, clon-dlatlon],[clat+dlatlon,clon+dlatlon]]
[[lllat,lllon],[urlat,urlon]] = BBox
miss = -9999.

srcDir = '/home/utsumi/temp/out'
prcpPath = srcDir + '/esurf.%06d.y%04d-%04d.nrec20000.npy'%(oid, iy, ey)
latPath  = srcDir + '/lat.%06d.y%04d-%04d.nrec20000.npy'%(oid, iy, ey)
lonPath  = srcDir + '/lon.%06d.y%04d-%04d.nrec20000.npy'%(oid, iy, ey)

#prcpPath = srcDir + '/esurf.%06d.y%04d-%04d.nrec20000.replace.npy'%(oid, iy, ey)
#latPath  = srcDir + '/lat.%06d.y%04d-%04d.nrec20000.replace.npy'%(oid, iy, ey)
#lonPath  = srcDir + '/lon.%06d.y%04d-%04d.nrec20000.replace.npy'%(oid, iy, ey)


#jplPath  = '/home/utsumi/bin/PMM/ret-epc/GPM_EPC_002421_20140802_0726.NS_MS.nc'  # HDF file
jplPath  = '/home/utsumi/bin/PMM/ret-epc/GPM_EPC_012149_20160418_1228.NS_MS.nc'  # HDF file

#tbPath   = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05/2014/08/02/1C.GPM.GMI.XCAL2016-C.20140802-S062222-E075455.002421.V05A.HDF5'

#*****************
#- Read My data ----

a2dat = np.load(prcpPath)
a2lat = np.load(latPath)
a2lon = np.load(lonPath)

a2dat, a2lattmp, a2lontmp = epcfunc.extract_domain_2D(a2dat, a2lat, a2lon, clat, clon, dlatlon, dscan)
a2dat = ma.masked_less_equal(a2dat,0)
#*****************
#- Read JPL data ----
with h5py.File(jplPath) as h:
    a2esurfjpl = h['MS/precip'][:]
    a2dpr      = h['DPR/precip_NS'][:]
    a2gprof    = h['GPROF/precip'][:]
    a2latjplOrg = h['latitude'][:]
    a2lonjplOrg = h['longitude'][:]

    

a2esurfjpl, a2lattmp, a2lontmp = epcfunc.extract_domain_2D(a2esurfjpl, a2latjplOrg, a2lonjplOrg, clat, clon, dlatlon, dscan)
a2esurfjpl = ma.masked_less_equal(a2esurfjpl,0)

a2dpr, a2lattmp, a2lontmp = epcfunc.extract_domain_2D(a2dpr, a2latjplOrg, a2lonjplOrg, clat, clon, dlatlon, dscan)
a2dpr = ma.masked_less_equal(a2dpr,0)

a2gprof, a2lattmp, a2lontmp = epcfunc.extract_domain_2D(a2gprof, a2latjplOrg, a2lonjplOrg, clat, clon, dlatlon, dscan)
a2gprof = ma.masked_less_equal(a2gprof,0)

#********************************
#-- Draw figure ---
fig   = plt.figure(figsize=(8,8))
ssize = 1

#-- JPL vs My --
#ax    = fig.add_subplot(111)
ax    = fig.add_axes([0.1,0.5,0.3,0.3])
ax.scatter(a2esurfjpl, a2dat)
ax.set_xlabel('JPL')
ax.set_ylabel('My')

#-- GPROF vs My --
#ax    = fig.add_subplot(111)
ax    = fig.add_axes([0.5,0.5,0.3,0.3])
ax.scatter(a2gprof, a2dat)
ax.set_xlabel('GPROF')
ax.set_ylabel('My')

#-- DPR vs My --
#ax    = fig.add_subplot(111)
ax    = fig.add_axes([0.1,0.1,0.3,0.3])
ax.scatter(a2dpr, a2dat)
ax.set_xlabel('DPR/NS')
ax.set_ylabel('My')

#-- DPR vs My --
#ax    = fig.add_subplot(111)
ax    = fig.add_axes([0.5,0.1,0.3,0.3])
ax.scatter(a2dpr, a2esurfjpl)
ax.set_xlabel('DPR/NS')
ax.set_ylabel('JPL')


#plt.suptitle('Replaced version')
plt.suptitle('oid=%06d (clat,clon)=(%.1f, %.1f)'%(oid, clat, clon))
#-- save ------
figPath = '/home/utsumi/temp/out/scatter.png'
#figPath = '/home/utsumi/temp/out/scatter.replace.png'
plt.savefig(figPath)
plt.clf()
print figPath
