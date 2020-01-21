import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from numpy import *
import os, sys
import h5py


srcDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/himawari/temp/temp.201707010220.sir.01'
srcPath= srcDir + '/grid20.dat'
a=fromfile(srcPath,'float32').reshape(6000,6000)
plt.imshow(a)
plt.colorbar()
plt.savefig('/home/utsumi/temp/geo/temp.png')
plt.clf()

#srcDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/himawari/obt.ptype/2017/07/01/018971'
#srcPath= srcDir + '/tir.02.npy'
#ptypePath= srcDir + '/ptype.npy'
#air=np.load(srcPath)
#atype= np.load(ptypePath)
#
#
#alat   = np.load(srcDir + '/lat.npy')
#alon   = np.load(srcDir + '/lon.npy')
#adpry  = np.load(srcDir + '/dpry.npy')
#adprx  = np.load(srcDir + '/dprx.npy')
#
#
#dprDir = '/home/utsumi/mnt/lab_work/hk02/PMM/NASA/GPM.Ku/2A/V06/2017/07/01'
#dprPath= dprDir + '/2A.GPM.Ku.V8-20180723.20170701-S012000-E025234.018971.V06A.HDF5'
#with h5py.File(dprPath, 'r') as h:
#    a2type= h['NS/CSF/typePrecip'][:]
#    a2lat = h['NS/Latitude'][:]
#    a2lon = h['NS/Longitude'][:]
#
#
#i=50
#dpry = adpry[i]
#dprx = adprx[i]
#lat = alat[i]
#lon = alon[i]
#
#latman = a2lat[dpry,dprx]
#lonman = a2lon[dpry,dprx]
#
#print lat, latman
#print lon, lonman
#
#geoPath = '/home/utsumi/mnt/lab_tank/utsumi/PMM/himawari/temp/temp.201707010220.tir.02/grid20.dat'

#nygeo, nxgeo = 6000,6000
#lon0 = 85   # Himawari data boundary
#lat0 = -60  # Himawari data boundary
#dlon = 0.02
#dlat = 0.02
#a2geo = np.flipud(np.fromfile(geoPath,'float32').reshape(6000,6000))
#
#
#i=50
#a2ir = air[i]
#lat  = alat[i]
#lon  = alon[i]
#
#ix = int((lon-lon0)/dlon)
#iy = int((lat-lat0)/dlat)
#a2man = a2geo[iy-3:iy+3+1, ix-3:ix+3+1]
#
#print a2ir, a2man

#print atype
#atype = (atype/10000000).astype('int16')
#print atype
#
#nstra = ma.masked_not_equal(atype,1).count()
#nconv = ma.masked_not_equal(atype,2).count()
#nothr = ma.masked_not_equal(atype,3).count()
#print nstra, nconv, nothr
#
#for i,ptype in enumerate(atype):
#    if ptype==2:
#        plt.imshow(air[i])
#        plt.colorbar()
#        plt.savefig('/home/utsumi/temp/geo/temp.png')
#        plt.clf()
#        sys.exit()   
