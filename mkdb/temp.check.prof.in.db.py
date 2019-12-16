from numpy import *
import h5py
import numpy as np
import glob


#****************************************
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

#****************************************

idx_db = 4991

dbDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/EPCDB/GMI.V05A.S1.ABp103-117'

prwatPath= dbDir + '/DPRGMI_NS_precipTotWaterCont/201701/DPRGMI_NS_precipTotWaterCont.%05d.npy'%(idx_db)
aprofpmw= np.load(prwatPath)

apyxpmw = np.load(dbDir + '/pYXpmw/201701/pYXpmw.%05d.npy'%(idx_db))
aoid    = np.load(dbDir + '/gNum/201701/gNum.%05d.npy'%(idx_db))
amdhms  = np.load(dbDir + '/mdhms/201701/mdhms.%05d.npy'%(idx_db))

print apyxpmw
#i=3833
i=3765
oid = aoid[i]
Mon,Day = amdhms[i][:2]
ypmw,xpmw= apyxpmw[i]

dpridxDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.IDX/2017/%02d/%02d'%(Mon,Day)
a2dprx = np.load(dpridxDir + '/Xpy.1.%06d.npy'%(oid)).astype(int32)
a2dpry = np.load(dpridxDir + '/Ypy.1.%06d.npy'%(oid)).astype(int32)

dprx = a2dprx[ypmw,xpmw-83]
dpry = a2dpry[ypmw,xpmw-83]

a1profpmw = aprofpmw[i]


#--- HDF file ------
dprDir = '/home/utsumi/mnt/lab_work/hk01/PMM/NASA/GPM.DPRGMI/2B/V06/2017/%02d/%02d'%(Mon,Day)
ssearch= dprDir + '/2B.GPM.DPRGMI.CORRA2018.*.%06d.V06A.HDF5'%(oid)
dprPath= glob.glob(ssearch)[0]
print dprPath
with h5py.File(dprPath,'r') as h:
    a3profd = h['NS/precipTotWaterCont'][:,:,-60:]

a3profdpr = a3profd[dpry-1:dpry+2,dprx-1:dprx+2]

a1profdprfunc = ave_9grids_3d(a3profd, np.array([dpry]), np.array([dprx]), miss=-9999.9).reshape(-1,)
a1profdprman  = ma.masked_less(a3profdpr,0).mean(axis=(0,1)).reshape(-1,)

a1dfunc = a1profpmw - a1profdprfunc
a1dman  = a1profpmw - a1profdprman
#srcPath = '/work/hk01/PMM/NASA/GPM.Ku/2A/V06/2017/01/01/2A.GPM.Ku.V8-20180723.20170101-S003624-E020857.016155.V06A.HDF5'
#
#
#
#with h5py.File(srcPath) as h:
#    a3dat = h['/NS/SLV/precipRate'][:,15:49-15+1,:]
#    a2surf= h['/NS/SLV/precipRateESurface'][:,15:49-15+1]
#print a3dat.shape
#a2flag = ma.masked_greater(a3dat,0).mask.any(axis=2)
#print a2flag.sum()
#print ''
#
#a1flag = a2flag.flatten()
#a2dat  = a3dat.reshape(-1,176)[a1flag]
#a1surf = a2surf.flatten()[a1flag]
#for i in range(100):
#    print a1surf[i],'  ', a2dat[i,:]
#
#'''
#a3dat  = (ma.masked_less(a3dat,0)*100).astype(int16)
#a3dat  = a3dat.filled(-9999)
#a2flag = ma.masked_greater(a3dat,0).mask.any(axis=2)
#print a3dat.shape
#a2flag = ma.masked_greater(a3dat,0).mask.any(axis=2)
#print a2flag.sum()
#'''
