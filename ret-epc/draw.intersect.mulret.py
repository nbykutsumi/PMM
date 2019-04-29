import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from numpy import *
import h5py
from bisect import bisect_left
import epcfunc
import sys, os, glob

# QJRMS case, oid=012149, 2016/4/18
oid = 12149
iy, ey = 1038, 1098
clat    = 32. # QJRMS case, oid=012149
clon    = -94 # -180 - +180
DB_MAXREC = 20000
Year,Mon,Day = 2016,4,18

## Africa case
#oid = 2421
#iy, ey = 2029, 2089
#clat    = 14 # Africa case
#clon    = 2  # -180 - +180

dlatlon = 5
#dscan   = 55
dscan   = 30
BBox    = [[clat-dlatlon, clon-dlatlon],[clat+dlatlon,clon+dlatlon]]
[[lllat,lllon],[urlat,urlon]] = BBox
miss = -9999.


srcDir = '/home/utsumi/temp/out'
outDir = '/home/utsumi/temp/out'
#stamp = '%06d.y%04d-%04d.nrec%d'%(oid, idx_c-dscan, idx_c+dscan,DB_MAXREC)
stamp = '%06d.y%04d-%04d.nrec%d'%(oid, iy, ey, DB_MAXREC)

profNSPath   = srcDir + '/prprofNS.%s.npy'%(stamp)
profNScmbPath= srcDir + '/prprofNScmb.%s.npy'%(stamp)

latPath   = srcDir + '/lat.%s.npy'%(stamp)
lonPath   = srcDir + '/lon.%s.npy'%(stamp)
#
#jplPath  = '/home/utsumi/bin/PMM/ret-epc/GPM_EPC_002421_20140802_0726.NS_MS.nc'  # HDF file, Africa case
jplPath  = '/home/utsumi/bin/PMM/ret-epc/GPM_EPC_012149_20160418_1228.NS_MS.nc'  # HDF file, QJRMS case

#tbPath   = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05/2014/08/02/1C.GPM.GMI.XCAL2016-C.20140802-S062222-E075455.002421.V05A.HDF5'
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



#*****************
#- Read My data ----

a2NS    = np.load(profNSPath)[:,110,:]
a2NScmb = np.load(profNScmbPath)[:,110,:]
a1latMy = np.load(latPath)[:,110]
a1lonMy = np.load(lonPath)[:,110]

#***********************************
#- Read JPL data ----
with h5py.File(jplPath) as h:
    a3profjpl   = h['NS/precip_prof'][:]
    a2latjplOrg = h['latitude'][:]
    a2lonjplOrg = h['longitude'][:]

a2tmp, a2latjpl, a2lonjpl, iyjpl, eyjpl = epcfunc.extract_domain_2D(a2latjplOrg, a2latjplOrg, a2lonjplOrg, clat, clon, dlatlon, dscan, returnidx=True)
a3profjpl = a3profjpl[iyjpl:eyjpl+1,:,:]
a2profjpl = a3profjpl[:,110,:]
a2profjpl = ma.masked_less_equal(a2profjpl,0)

a2latjpl = a2latjplOrg[iyjpl:eyjpl+1,:]
a2lonjpl = a2lonjplOrg[iyjpl:eyjpl+1,:]

#***********************************
#- Read DPR data ----
#-- Cooresponding DPR yx --------
srcDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.IDX/%04d/%02d/%02d'%(Year,Mon,Day)
xPath= srcDir + '/Xpy.1.%06d.npy'%(oid)
yPath= srcDir + '/Ypy.1.%06d.npy'%(oid)

##-- extract center and domain --
a1xdpr = np.load(xPath)[iy:ey+1,110-83]
a1ydpr = np.load(yPath)[iy:ey+1,110-83]

#-- Read DPR ----
dprDir  = '/work/hk01/PMM/NASA/GPM.Ku/2A/V06/%04d/%02d/%02d'%(Year,Mon,Day)
dprPath = glob.glob(dprDir + '/2A.GPM.Ku.*.%06d.V06A.HDF5'%(oid))[0]
#srcPath = '/work/hk01/PMM/NASA/GPM.Ku/2A/V06/2017/01/01/2A.GPM.Ku.V8-20180723.20170101-S173441-E190714.016166.V06A.HDF5'
with h5py.File(dprPath) as h:
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

#********************************
#-- Draw figure ---
fig   = plt.figure(figsize=(8,8))
ssize = 1

#-- My retrieval --
for i in range(4):
    nx,nh = a2NS.shape
    a1y   = arange(nh)*0.25
    a1x   = a1latMy
    cmap  = 'jet'
    if i==0:
        ax = fig.add_axes([0.1,0.1,0.35,0.35])
        a2dat = ma.masked_less_equal(a2NS.T,0)*0.01
        im = ax.imshow(a2dat,cmap=cmap)
        stype = 'NS'

    elif i==1:
        ax = fig.add_axes([0.5,0.1,0.35,0.35])
        a2dat = ma.masked_less_equal(a2NScmb.T,0)*0.01
        im = ax.imshow(a2dat,cmap=cmap)
        stype = 'NScmb'
    elif i==2:
        ax = fig.add_axes([0.1,0.5,0.35,0.35])
        a2dat = ma.masked_less_equal(a2profjpl.T,0)*0.01
        im = ax.imshow(a2dat,cmap=cmap)
        stype = 'JPL-NScmb'

    elif i==3:
        ax = fig.add_axes([0.5,0.5,0.35,0.35])
        a2dat = ma.masked_less_equal(a2dpr.T,0)
        im = ax.imshow(a2dat,cmap=cmap)
        stype = 'DPR'

    print ''
    print ''
    print stype
    print a2dat.min(),a2dat.max()

    plt.colorbar(im, orientation='horizontal')
    plt.title(stype)
##------------
outPath  = outDir + '/prcp.prof.mulret.png'
plt.savefig(outPath)
print outPath
plt.clf()
