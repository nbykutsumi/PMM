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

srcDir = '/home/utsumi/temp/out'
prcpPath = srcDir + '/esurf.npy'
latPath  = srcDir + '/lat.npy'
lonPath  = srcDir + '/lon.npy'

a2dat = np.load(prcpPath)
a2lat = np.load(latPath)
a2lon = np.load(lonPath)

a2dat = ma.masked_less_equal(a2dat,0)

ssize = 1
#********************************
# Functions
#--------------------------------
def extract_domain_2D(a2dat, a2lat, a2lon, clat, clon, dlatlon, dscan):

    nyTmp, nxTmp = a2lat.shape
    a1lat = a2lat[:,nxTmp/2]
    a1lon = a2lon[:,nxTmp/2]
    
    idx_latmax = np.argmax(a1lat)
    a1lat0 = a1lat[:idx_latmax+1]
    a1lat1 = a1lat[idx_latmax+1:]
    a1lon0 = a1lon[:idx_latmax+1]
    a1lon1 = a1lon[idx_latmax+1:]
    
    if (-180<=clat)and(clat <=180):
        #-- search first half: ascending --
        found_domain = 0
        idx_c  = bisect_left(a1lat0, clat)
        latTmp = a1lat0[idx_c]
        lonTmp = a1lon0[idx_c]
        if (clat-dlatlon<=latTmp)&(latTmp <=clat+dlatlon)&(clon-dlatlon<=lonTmp)&(lonTmp<=clon+dlatlon):
            found_domain = 1
        else:
            #-- search second half: descending --
            idx_c  = bisect_left(a1lat1[::-1], clat)
            idx_c  = len(a1lat) - idx_c -1
            latTmp = a1lat[idx_c]
            lonTmp = a1lon[idx_c]
    
            if (clat-dlatlon<=latTmp)&(latTmp <=clat+dlatlon)&(clon-dlatlon<=lonTmp)&(lonTmp<=clon+dlatlon):
                found_domain =1
    
        if found_domain==1:
            idx_first = idx_c - dscan
            idx_last  = idx_c + dscan
            a2odat  = a2dat[idx_first:idx_last+1,:]
            a2olat  = a2lat[idx_first:idx_last+1,:]    
            a2olon  = a2lon[idx_first:idx_last+1,:]
    
        else:
            print 'No matching scans in the target domain are found.'
            print 'Exit'
            sys.exit()
    
        print 'Extract target domain'
        print 'Extracted array size=', a2dat.shape
    
    else:
        pass
    
    return a2odat, a2olat, a2olon 

#********************************
#-- Draw figure ---
fig   = plt.figure(figsize=(8,8))

#-- My retrieval --
#ax    = fig.add_subplot(111)
ax    = fig.add_axes([0.1,0.5,0.3,0.3])
M     = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax)
im    = M.scatter(a2lon, a2lat, c=a2dat, cmap='jet', s=ssize, vmin=0, vmax=20)
M.drawcoastlines()

dgrid      = 5
parallels  = arange(-90,90, dgrid)
meridians  = arange(-180,180,dgrid)
M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8, linewidth=0.5, fmt='%d')
M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8, linewidth=0.5, fmt='%d')
plt.colorbar(im, orientation='horizontal')

plt.title('EPC')
#*****************
#- Read GPROF ----
gpDir = '/work/hk01/PMM/NASA/GPM.GMI/2A/V05/2017/01/01/2A.GPM.GMI.GPROF2017v1.20170101-S173441-E190714.016166.V05A.HDF5'
with h5py.File(gpDir) as h:
    a2esurfgp = h['S1/surfacePrecipitation'][:]
    a2latgp   = h['S1/Latitude'][:]
    a2longp   = h['S1/Longitude'][:]

print a2esurfgp.shape
print a2latgp.shape
a2esurfgp, a2latgp, a2longp = extract_domain_2D(a2esurfgp, a2latgp, a2longp, clat, clon, dlatlon, dscan)
a2esurfgp = ma.masked_less_equal(a2esurfgp,0)

#- fig --
ax2 = fig.add_axes([0.1,0.1,0.3,0.3])
M   = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax2)
im    = M.scatter(a2longp, a2latgp, c=a2esurfgp, cmap='jet', s=ssize, vmin=0, vmax=20)
print a2esurfgp
M.drawcoastlines()
M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8, linewidth=0.5, fmt='%d')
M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8, linewidth=0.5, fmt='%d')
plt.colorbar(im, orientation='horizontal')

plt.title('GPROF')
#------------
#- Read DPR ----
gpDir = '/work/hk01/PMM/NASA/GPM.Ku/2A/V06/2017/01/01/2A.GPM.Ku.V8-20180723.20170101-S173441-E190714.016166.V06A.HDF5'
with h5py.File(gpDir) as h:
    a2datfig  = h['NS/SLV/precipRateESurface'][:]
    a2latfig  = h['NS/Latitude'][:]
    a2lonfig  = h['NS/Longitude'][:]

a2datfig, a2latfig, a2lonfig = extract_domain_2D(a2datfig, a2latfig, a2lonfig, clat, clon, dlatlon, dscan*5)
a2datfig = ma.masked_less_equal(a2datfig,0)

#- fig --
ax2 = fig.add_axes([0.5,0.1,0.3,0.3])
M   = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax2)
im    = M.scatter(a2lonfig, a2latfig, c=a2datfig, cmap='jet', s=ssize, vmin=0, vmax=20)
print a2esurfgp
M.drawcoastlines()
M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8, linewidth=0.5, fmt='%d')
M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8, linewidth=0.5, fmt='%d')
plt.colorbar(im, orientation='horizontal')

plt.title('DPR')
#------------

outPath  = srcDir + '/esurf.png'
plt.savefig(outPath)
print outPath
