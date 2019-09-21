from numpy import *
import numpy as np
import glob
import myfunc.util as util
import os, sys

#srcPath1 = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/glb.wprof.rnr/2014/06/01/top-tbNS.001453.y-9999--9999.nrec20000.npy'
#srcPath1 = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/test/2014/07/01/top-tbNS.001927.y1764-1784.nrec20000.npy'
#srcPath1 = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/test/2014/07/01/top-zmMS.001927.y1764-1784.nrec20000.npy'
srcPath1 = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/test/2014/07/01/top-prwatprofNS.001927.y1764-1784.nrec20000.npy'
a = np.load(srcPath1)

#srcPath2 = '/home/utsumi/temp/ret/top-tbNS.001453.y-9999--9999.nrec20000.npy'
#srcPath2 = '/home/utsumi/temp/ret/top-tbNS.001927.y1764-1784.nrec20000.npy'
#srcPath2 = '/home/utsumi/temp/ret/top-zmMS.001927.y1764-1784.nrec20000.npy'
srcPath2 = '/home/utsumi/temp/ret/top-prwatprofNS.001927.y1764-1784.nrec20000.npy'

b = np.load(srcPath2)

print 'shape'
print ma.masked_less(a,0).shape
print ma.masked_less(b,0).shape



print 'shape-flatten'
print ma.masked_less(a,0).flatten().shape
print ma.masked_less(b,0).flatten().shape


print 'count'
print ma.masked_less(a,0).count()
print ma.masked_less(b,0).count()


a=ma.masked_invalid(a)
b=ma.masked_invalid(b)
print 'sum'
print ma.masked_less(a,0).sum()
print ma.masked_less(b,0).sum()

ny,nx,nz = a.shape
icount =0
for y in range(ny):
    for x in range(nx):
        if a[y,x,nz-1] >0:
            #print map(a[y,x,:], b[y,x,:])
            if not (a[y,x,:] == b[y,x,:]).all():
                print 'y,x=',y,x
                print 'a=',a[y,x,:]
                print 'b=',b[y,x,:]
                icount = icount + 1
print icount    
sys.exit()


#X,Y = np.meshgrid(np.arange(nx), np.arange(ny))
#a1x = X.flatten()
#a1y = Y.flatten()


