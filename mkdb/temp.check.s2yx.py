from numpy import *
import numpy as np

#srcDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.GMI.S2.IDX/2017/01/02'
srcDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.IDX/2017/01/02'
rev = 16171

xpath  = srcDir + '/Xpy.1.%06d.npy'%(rev)
ypath  = srcDir + '/Ypy.1.%06d.npy'%(rev)

a2x = np.load(xpath)
a2y = np.load(ypath)

ny,nx = a2x.shape
for iy in range(ny):
    for ix in range(nx):
        x = a2x[iy,ix]
        if x <0:
            print 'iy,ix=',iy,ix, 'loc=',x 

