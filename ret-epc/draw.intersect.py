import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from numpy import *
import h5py
from bisect import bisect_left
import sys
import epcfunc

clat    = 30.00
clon    = 269.0 -360  # -180 - +180
dlatlon = 6
dscan   = 30
BBox    = [[clat-dlatlon, clon-dlatlon],[clat+dlatlon,clon+dlatlon]]
[[lllat,lllon],[urlat,urlon]] = BBox


oid = 16166
iy, ey = 1952, 2012
ssize = 1

#********************************
# Functions
#--------------------------------
def ret_domain_cy(a2lat, a2lon, clat, dlatlon):
    nyTmp, nxTmp = a2lat.shape
    a1lat = a2lat[:,nxTmp/2]
    a1lon = a2lon[:,nxTmp/2]

    idx_latmax = np.argmax(a1lat)
    a1lat0 = a1lat[:idx_latmax+1]
    a1lat1 = a1lat[idx_latmax+1:]
    a1lon0 = a1lon[:idx_latmax+1]
    a1lon1 = a1lon[idx_latmax+1:]

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
        return idx_c

    else:
        print 'No matching scans in the target domain are found.'
        print 'Exit'
        sys.exit()

#--------------------------------
def ave_9grids_3d(a3in, a1y, a1x, miss):
    '''
    returns 2-d array with the size of (nl,nz)
    a3in: (ny,nx,nz)
    nl = len(a1y)=len(a1x)
    output: (nl, nz)
    '''
    #-- Average 9 grids (over Linearlized Z)--
    nydpr,nxdpr,nzdpr= a3in.shape
    ldydx = [[dy,dx] for dy in [-1,0,1] for dx in [-1,0,1]]

    a1dprmask   = False

    a3datTmp    = empty([9,len(a1y),nzdpr], float32)

    for itmp, [dy,dx] in enumerate(ldydx):
        a1yTmp  = ma.masked_outside(a1y + dy, 0,nydpr-1)
        a1xTmp  = ma.masked_outside(a1x + dx, 0,nxdpr-1)
        a1dprmask= a1dprmask + a1yTmp.mask + a1xTmp.mask

        a2datTmp= a3in[a1yTmp.filled(0),a1xTmp.filled(0),:]

        a3datTmp[itmp,:] = a2datTmp


    a2datTmp = ma.masked_equal(a3datTmp,miss).mean(axis=0)
    a2datTmp[a1dprmask,:] = miss
    return a2datTmp


#--------------------------------
def average_2ranges_3d(a3in,miss=None,dtype=float32, fill=True):
    '''
    a3in: (ny, nx, nz) --> a2out: (ny, nx, nz/2)
    nz should be an even number
    '''
    ny,nx,nz = a3in.shape
    a4out = empty([2,ny,nx,nz/2], dtype)
    a4out[0] = a3in[:,:,::2]
    a4out[1] = a3in[:,:,1::2]
    if fill==True:
        a3out = ma.masked_equal(a4out,miss).mean(axis=0).filled(miss).astype(dtype)
    else:
        a3out = ma.masked_equal(a4out,miss).mean(axis=0).astype(dtype)
    return a3out


#********************************
# Read data
#--------------------------------

#-- Retrieval --
srcDir = '/home/utsumi/temp/out'
profPath = srcDir + '/prprof.%06d.y%04d-%04d.npy'%(oid, iy, ey)
latPath  = srcDir + '/lat.%06d.y%04d-%04d.npy'%(oid, iy, ey)
lonPath  = srcDir + '/lon.%06d.y%04d-%04d.npy'%(oid, iy, ey)

a3epc = np.load(profPath)
a2lat = np.load(latPath)
a2lon = np.load(lonPath)

##-- extract 83-137 --
#a2epc = np.load(prcpPath)[:,83:137+1]
#a2lat = np.load(latPath)[:,83:137+1]
#a2lon = np.load(lonPath)[:,83:137+1]

##-- extract center --
a2epc = np.load(profPath)[:,110,:]
a2lat = np.load(latPath)[:,110]
a2lon = np.load(lonPath)[:,110]


a2epc = ma.masked_less_equal(a2epc,0)*0.01

#-- Cooresponding DPR yx --------
srcDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.IDX/2017/01/01'
xPath= srcDir + '/Xpy.1.%06d.npy'%(oid)
yPath= srcDir + '/Ypy.1.%06d.npy'%(oid)

##-- extract center and domain --
a1xdpr = np.load(xPath)[iy:ey+1,110-83]
a1ydpr = np.load(yPath)[iy:ey+1,110-83]

#-- Read DPR ----
srcPath = '/work/hk01/PMM/NASA/GPM.Ku/2A/V06/2017/01/01/2A.GPM.Ku.V8-20180723.20170101-S173441-E190714.016166.V06A.HDF5'
with h5py.File(srcPath) as h:
    #a3dpr = h['NS/SLV/precipRate'][:]
    a3dpr = h['NS/SLV/precipRate'][:]
    a2latdpr = h['NS/Latitude'][:]
    a2londpr = h['NS/Longitude'][:]
    a2cfbBin = h['NS/PRE/binClutterFreeBottom'][:]
    a2surfBin= h['NS/PRE/binRealSurface'][:]


#- extract lower levels --
a3dpr = a3dpr[:,:,-44:]
#- average every 2-levels --
a3dpr = average_2ranges_3d(a3dpr, miss=-9999.9)
#- extract pixels matching with GMI --
a2dpr = ave_9grids_3d(a3dpr, a1ydpr, a1xdpr, miss=-9999.9)



#-- compair --
a2dprTmp = a3dpr[a1ydpr, a1xdpr,:]
a1cfbBin = a2cfbBin[a1ydpr,a1xdpr]
a1surfBin= a2surfBin[a1ydpr,a1xdpr]
print a2epc.shape, a2dpr.shape
for i in range(a2epc.shape[0]):
    cfbbin = a1cfbBin[i]
    surfbin= a1surfBin[i]
    print ''
    print i,'--------------------'
    print 'cfb-bin, surfbin=',cfbbin/2, surfbin/2 
    print a2dprTmp[i]
    print a2dpr[i]
    print a2epc[i]


#-- Draw figure ----------
vmin= 0
vmax= 10
cmap= 'jet'

fig = plt.figure(figsize=(8,5))

#- DPR --
ax  = fig.add_axes([0.1,0.6,0.8,0.3])
im  = ax.imshow(ma.masked_less(a2dpr.T,0.01), vmin=vmin, vmax=vmax, cmap=cmap)
plt.title('Ku')
#- EPC --
ax  = fig.add_axes([0.1,0.2,0.8,0.3])
im  = ax.imshow(ma.masked_less(a2epc.T,0.01), vmin=vmin, vmax=vmax, cmap=cmap)
plt.title('EPC')

#- options --
ax  = fig.add_axes([0.2,0.1,0.6,0.05])
plt.colorbar(im, cax=ax, orientation='horizontal')
#- Save --
figDir  = '/home/utsumi/temp/out'
figPath = figDir + '/intersect.png'
plt.savefig(figPath)
print figPath
