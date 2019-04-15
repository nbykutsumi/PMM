import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from numpy import *
import h5py
from bisect import bisect_left
import epcfunc

## QJRMS case, oid=012149, 2016/4/18
#oid = 12149
#iy, ey = 1038, 1098
#clat    = 32. # QJRMS case, oid=012149
#clon    = -94 # -180 - +180

# Africa case
oid = 2421
iy, ey = 2029, 2089
clat    = 14 # Africa case
clon    = 2  # -180 - +180

dlatlon = 5
#dscan   = 55
dscan   = 25
BBox    = [[clat-dlatlon, clon-dlatlon],[clat+dlatlon,clon+dlatlon]]
[[lllat,lllon],[urlat,urlon]] = BBox
miss = -9999.


srcDir = '/home/utsumi/temp/out'
#prcpPath = srcDir + '/esurf.%06d.y%04d-%04d.npy'%(oid, iy, ey)
prcpPath = srcDir + '/esurf.%06d.y%04d-%04d.nrec20000.npy'%(oid, iy, ey)
latPath  = srcDir + '/lat.%06d.y%04d-%04d.nrec20000.npy'%(oid, iy, ey)
lonPath  = srcDir + '/lon.%06d.y%04d-%04d.nrec20000.npy'%(oid, iy, ey)

#prcpPath = srcDir + '/esurf.%06d.y%04d-%04d.replace.npy'%(oid, iy, ey)
#latPath  = srcDir + '/lat.%06d.y%04d-%04d.replace.npy'%(oid, iy, ey)
#lonPath  = srcDir + '/lon.%06d.y%04d-%04d.replace.npy'%(oid, iy, ey)


#
jplPath  = '/home/utsumi/bin/PMM/ret-epc/GPM_EPC_002421_20140802_0726.NS_MS.nc'  # HDF file, Africa case
#jplPath  = '/home/utsumi/bin/PMM/ret-epc/GPM_EPC_012149_20160418_1228.NS_MS.nc'  # HDF file, QJRMS case


tbPath   = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05/2014/08/02/1C.GPM.GMI.XCAL2016-C.20140802-S062222-E075455.002421.V05A.HDF5'

#*****************
#- Read My data ----

a2dat = np.load(prcpPath)
a2lat = np.load(latPath)
a2lon = np.load(lonPath)


a2dat = ma.masked_less_equal(a2dat,0)

#*****************
#- Read JPL data ----
with h5py.File(jplPath) as h:
    a2esurfjpl = h['MS/precip'][:]
    #a2dpr      = h['DPR/precip_NS'][:]
    a2dpr      = h['CMB/precip_NS'][:]
    a2latjplOrg   = h['latitude'][:]
    a2lonjplOrg   = h['longitude'][:]


a2esurfjpl, a2latjpl, a2lonjpl = epcfunc.extract_domain_2D(a2esurfjpl, a2latjplOrg, a2lonjplOrg, clat, clon, dlatlon, dscan)
a2esurfjpl = ma.masked_less_equal(a2esurfjpl,0)

a2dpr, a2lattmp, a2lontmp = epcfunc.extract_domain_2D(a2dpr, a2latjplOrg, a2lonjplOrg, clat, clon, dlatlon, dscan)
a2dpr = ma.masked_less_equal(a2dpr,0)

##-- Read Tb data --
#with h5py.File(tbPath, 'r') as h:
#    a3tb1    = h['/S1/Tc'][:]
#    a3tb2org = h['/S2/Tc'][:]
#    a2lattb = h['/S1/Latitude'][:]
#    a2lontb = h['/S1/Longitude'][:]
#
#s2xPath= '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/2014/08/02/Xpy.1.002421.npy'
#s2yPath= '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/2014/08/02/Ypy.1.002421.npy'
#
#
##-- Matchup and Joint S1 and S2 Tb --
#a1x2  = np.load(s2xPath).flatten()
#a1y2  = np.load(s2yPath).flatten()
#
#a1mask= ma.masked_less(a1x2,0).mask
#a1x2  = ma.masked_less(a1x2,0).filled(0)
#a1y2  = ma.masked_less(a1y2,0).filled(0)
#
#nytmp, nxtmp, ztmp = a3tb2org.shape
#a2tb2 = a3tb2org[a1y2, a1x2]
#a2tb2[a1mask] = miss
#a3tb2 = a2tb2.reshape(nytmp,nxtmp,-1)
#a3tb = concatenate([a3tb1, a3tb2],axis=2)
#a2mask = np.any(ma.masked_outside(a3tb, 50, 350).mask, axis=2)
#
#print a2mask.shape, a2lattb.shape, a2lontb.shape
#a2mask, a2lattmp, a2lontmp = epcfunc.extract_domain_2D(a2mask, a2lattb, a2lontb, clat, clon, dlatlon, dscan)
#
##-- test ---
#plt.imshow(a2mask*1)
#plt.colorbar()
#plt.savefig('/home/utsumi/temp/out/temp.png')
#plt.clf()
#print 'check a2dat'
#print a2mask[50:55,105:110]




#********************************
#-- Draw figure ---
fig   = plt.figure(figsize=(8,8))
ssize = 1

#-- My retrieval --
#ax    = fig.add_subplot(111)
ax    = fig.add_axes([0.1,0.5,0.35,0.35])
M     = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax)
im    = M.scatter(a2lon, a2lat, c=a2dat, cmap='jet', s=ssize, vmin=0, vmax=20)
M.drawcoastlines()

dgrid      = 5
parallels  = arange(-90,90, dgrid)
meridians  = arange(-180,180,dgrid)
M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8, linewidth=0.5, fmt='%d')
M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8, linewidth=0.5, fmt='%d')
plt.colorbar(im, orientation='horizontal')

plt.title('NU')
#*****************
#- Figure JPL-EPC ----
ax2 = fig.add_axes([0.1,0.1,0.35,0.35])
M   = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax2)
im    = M.scatter(a2lonjpl, a2latjpl, c=a2esurfjpl, cmap='jet', s=ssize, vmin=0, vmax=20)
M.drawcoastlines()
M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8, linewidth=0.5, fmt='%d')
M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8, linewidth=0.5, fmt='%d')
plt.colorbar(im, orientation='horizontal')

plt.title('MS/precip')
#------------
#- Figure DPR ----
ax2 = fig.add_axes([0.5,0.1,0.35,0.35])
M   = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax2)
im  = M.scatter(a2lonjpl, a2latjpl, c=a2dpr, cmap='jet', s=ssize, vmin=0, vmax=20)
M.drawcoastlines()
M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8, linewidth=0.5, fmt='%d')
M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8, linewidth=0.5, fmt='%d')
plt.colorbar(im, orientation='horizontal')

plt.title('CMB/NS')

##------------
##- Read DPR 9grid mean ----
#dprmeanPath = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.9ave.precipRateESurface/2017/01/01/precipRateESurface.016166.npy'
#a2datfig   = np.load(dprmeanPath)
#a2datfig   = a2datfig[iy:ey+1]
#
#a2datfig = ma.masked_less_equal(a2datfig,0)
#
##- lat and lon from GPROF ----
#gpDir = '/work/hk01/PMM/NASA/GPM.GMI/2A/V05/2017/01/01/2A.GPM.GMI.GPROF2017v1.20170101-S173441-E190714.016166.V05A.HDF5'
#with h5py.File(gpDir) as h:
#    a2latgmi = h['S1/Latitude'][:]
#    a2longmi = h['S1/Longitude'][:]
#
#a2latfig = a2latgmi[iy:ey+1,83:137+1]
#a2lonfig = a2longmi[iy:ey+1,83:137+1]
#
#print a2datfig.shape, a2latfig.shape, a2lonfig.shape
##- fig --
#ax2 = fig.add_axes([0.5,0.5,0.3,0.3])
#M   = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax2)
#im    = M.scatter(a2lonfig, a2latfig, c=a2datfig, cmap='jet', s=ssize, vmin=0, vmax=20)
#print a2esurfgp
#M.drawcoastlines()
#M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8, linewidth=0.5, fmt='%d')
#M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8, linewidth=0.5, fmt='%d')
#plt.colorbar(im, orientation='horizontal')
#
#plt.title('DPR-9ave')
##------------


outPath  = srcDir + '/esurf.png'
plt.savefig(outPath)
print outPath
