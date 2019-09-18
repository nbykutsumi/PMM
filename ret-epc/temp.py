from numpy import *
import numpy as np
import glob
import myfunc.util as util

srcPath1 = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/test/2014/07/01/nsurfMScmb.001927.y1764-1784.nrec20000.npy'
srcPath2 = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/rnr/2014/07/01/nsurfMScmb.001927.y1764-1784.nrec20000.npy'

a= np.load(srcPath1)
b= np.load(srcPath2)

print ma.masked_less(a,0).count()
print ''
print ma.masked_less(b,0).count()
