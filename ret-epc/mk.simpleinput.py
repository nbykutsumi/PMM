import h5py
import glob
import numpy as np
import h5py
import itertools

#Year,Mon,Day = 2014,6,1
#oid     = 1453
#
#gmibaseDir  = '/home/utsumi/mnt/lab_work/hk02/PMM/NASA/GPM.GMI/1C/V05'
#matchbaseDir= '/media/disk2/share/PMM/MATCH.GMI.V05A'
#tankDir = '/home/utsumi/mnt/lab_tank'
#
#
#rnrPath = glob.glob(tankDir + '/utsumi/PMM/retepc/glb.v03.minrec1000.maxrec10000/%04d/%02d/%02d/nsurfNScmb.%06d.y-9999--9999.nrec10000.npy'%(Year,Mon,Day,oid))[0]
#srcPath = glob.glob(gmibaseDir + '/%04d/%02d/%02d/1C.GPM.GMI.XCAL2016-C.*.%06d.????.HDF5'%(Year,Mon,Day,oid))[0]
#s2xPath = glob.glob(matchbaseDir + '/S1.ABp000-220.GMI.S2.IDX/%04d/%02d/%02d/Xpy.1.%06d.npy'%(Year,Mon,Day,oid))[0]
#s2yPath = glob.glob(matchbaseDir + '/S1.ABp000-220.GMI.S2.IDX/%04d/%02d/%02d/Ypy.1.%06d.npy'%(Year,Mon,Day,oid))[0]
#
#
#a2rnr = np.load(rnrPath)
#a2x   = np.load(s2xPath)
#a2y   = np.load(s2yPath)
#
#a2test = np.load('/home/utsumi/mnt/lab_tank/utsumi/hometemp/share-epc/sampledata/oid-1453/nsurfNScmb.000000.y-9999--9999.nrec10000.npy')
#
#nytmp, nxtmp = a2rnr.shape
#lxtmp = range(nxtmp)
#lytmp = range(nytmp)
#for y,x in itertools.product(lytmp,lxtmp):
#    if a2rnr[y,x]>1:
#        print(y,x,a2rnr[y,x])
#        break
#
#print(y,x,a2rnr[y,x])
##-- Read TC ----
#with h5py.File(srcPath,'r') as h:
#    a3tb1 = h['/S1/Tc'][:]
#    a3tb2 = h['/S2/Tc'][:]
#    #a2lat.append( h['/S1/Latitude'][:])
#    #a2lon.append( h['/S1/Longitude'][:])
#
#
##y,x = np.unravel_index(np.argmax(a2rnr), a2rnr.shape)
#
#y2 = a2y[y,x]
#x2 = a2x[y,x]
#
#print(y,x)
#print(y2,x2)
#print(a2rnr[y,x],a2test[y,x])
#
##-- Make sample TC set ----
#atb1 = a3tb1[y,x,:].reshape(-1,1,9)
#atb2 = a3tb2[y2,x2,:].reshape(-1,1,4)
#print(atb1)
#print(atb2)
#
##-- Save sample 1D ------
#datdir = '/home/utsumi/temp/share-epc/sampledata'
#tbpath1  = datdir + '/tb1.npy'
#tbpath2  = datdir + '/tb2.npy'
#
#np.save(tbpath1, atb1)
#np.save(tbpath2, atb2)

#------------------------
# Read Sample TBs
#------------------------
datdir = '/home/utsumi/temp/share-epc/sampledata'
tbpath1  = datdir + '/tb1.npy'
tbpath2  = datdir + '/tb2.npy'

atb1   = np.load(tbpath1)
atb2   = np.load(tbpath2)

#------------------------
# Create and save HDF file
#------------------------

hdfpath = datdir + '/tb.HDF5'
ny,_,__ = atb1.shape
alat    = np.full([ny,1],0, 'int32')  # dummy
alon    = np.full([ny,1],0, 'int32')  # dummy
ainc    = np.full([ny,1,1],0, 'int32')  # dummy

with h5py.File(hdfpath, 'w') as h:
    h.create_group('S1')
    h.create_dataset('S1/Tc', data=atb1)

    h.create_group('S2')
    h.create_dataset('S2/Tc', data=atb2)

    h.create_dataset('S1/Latitude', data=alat)
    h.create_dataset('S1/Longitude',data=alon)
    h.create_dataset('S1/incidenceAngle',data=ainc)
    h.flush()

#------------------------
# Create and save dummy Xpy, Ypy file
#------------------------
oxpath = datdir + '/Xpy.npy'
oypath = datdir + '/Ypy.npy'

axdummy = np.zeros(ny).astype('int32').reshape([ny,1])
aydummy = np.arange(ny).astype('int32').reshape([ny,1])

np.save(oxpath, axdummy)
np.save(oypath, aydummy)

print(axdummy)
print(aydummy)
