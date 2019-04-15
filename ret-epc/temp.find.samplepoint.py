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
dscan   = 1
BBox    = [[clat-dlatlon, clon-dlatlon],[clat+dlatlon,clon+dlatlon]]
[[lllat,lllon],[urlat,urlon]] = BBox
miss = -9999.

#oid = 16166
#iy, ey = 1952, 2012
oid = 2421
iy, ey = 2034, 2084
maxrec = 20000
srcDir = '/home/utsumi/temp/out'
#prcpPath = srcDir + '/esurf.%06d.y%04d-%04d.npy'%(oid, iy, ey)
prcpPath = srcDir + '/esurf.%06d.y%04d-%04d.%dmaxrec.noelev.npy'%(oid, iy, ey, maxrec)
latPath  = srcDir + '/lat.%06d.y%04d-%04d.%dmaxrec.noelev.npy'%(oid, iy, ey, maxrec)
lonPath  = srcDir + '/lon.%06d.y%04d-%04d.%dmaxrec.noelev.npy'%(oid, iy, ey, maxrec)
idxPath  = srcDir + '/temp.idx.npy'

jplPath  = '/home/utsumi/bin/PMM/ret-epc/GPM_EPC_002421_20140802_0726.NS_MS.nc'  # HDF file

tbPath   = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05/2014/08/02/1C.GPM.GMI.XCAL2016-C.20140802-S062222-E075455.002421.V05A.HDF5'

#*****************
#- Read My data ----

a2dat = np.load(prcpPath)
a2idx = np.load(idxPath)
a2lat = np.load(latPath)
a2lon = np.load(lonPath)

#a2idx, a2lat, a2lon = epcfunc.extract_domain_2D(a2idx, a2lat, a2lon, clat, clon, dlatlon, dscan)
a2dat, a2lat, a2lon = epcfunc.extract_domain_2D(a2dat, a2lat, a2lon, clat, clon, dlatlon, dscan)

a2dat = ma.masked_less_equal(a2dat,0)
#a2idx = ma.masked_less_equal(a2idx,0)


#*****************
#- Read JPL data ----
with h5py.File(jplPath) as h:
    a2esurfjpl = h['MS/precip'][:]
    a2dpr      = h['DPR/precip_NS'][:]
    a2latjpl   = h['latitude'][:]
    a2lonjpl   = h['longitude'][:]
    a2idxjpl   = h['db_index'][:]


#a2idxjpl, a2latjpl, a2lonjpl, iyjpl, eyjpl = epcfunc.extract_domain_2D(a2idxjpl, a2latjpl, a2lonjpl, clat, clon, dlatlon, dscan, returnidx=True)
a2esurfjpl, a2latjpl, a2lonjpl, iyjpl, eyjpl = epcfunc.extract_domain_2D(a2esurfjpl, a2latjpl, a2lonjpl, clat, clon, dlatlon, dscan, returnidx=True)

a2esurfjpl = ma.masked_less(a2esurfjpl,0)

ny,nx = a2esurfjpl.shape
for iy in range(ny):
    for ix in range(nx):
        dat0 = a2dat[iy,ix]
        dat1 = a2esurfjpl[iy,ix]
        ddat = dat0-dat1
        if abs(ddat)>1:
            print iy,ix, dat0,dat1,dat0-dat1, a2lat[iy,ix],a2lon[iy,ix]
print ny,nx
##*****************
##- Read L1C Tb ---
#tbPath = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05/2014/08/02/1C.GPM.GMI.XCAL2016-C.20140802-S062222-E075455.002421.V05A.HDF5'
#with h5py.File(tbPath) as h:
#    a2latorg = h['S1/Latitude'][:]
#    a2lonorg = h['S1/Longitude'][:]
#
#
#
#a2dattmp, a2lattmp, a2lontmp, iyorg, eyorg = epcfunc.extract_domain_2D(a2latorg, a2latorg, a2lonorg, clat, clon, dlatlon, dscan)
#
#
