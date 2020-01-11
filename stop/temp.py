import matplotlib
matplotlib.use('Agg')
from numpy import *
import numpy as np
import matplotlib.pyplot as plt


srcPath0 = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/glb.v03.minrec1000.maxrec10000/2014/10/14//nsurfNScmb.003556.y-9999--9999.nrec10000.npy'
srcPath1 = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/glb.stop-wgt-obs-01.minrec1000.maxrec10000/2014/10/14/nsurfNScmb.003556.y-9999--9999.nrec10000.npy'
srcPath2 = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/glb.stop-rng-obs-01.minrec1000.maxrec10000/2014/10/14/nsurfNScmb.003556.y-9999--9999.nrec10000.npy'

a0 = ma.masked_less(np.load(srcPath0)[:,83:137+1], 0.2)
a1 = ma.masked_less(np.load(srcPath1)[:,83:137+1], 0)
a2 = ma.masked_less(np.load(srcPath2)[:,83:137+1], 0)

print a0.sum(), a1.sum(), a2.sum()
print a0.shape, a1.shape, a2.shape
print a0.count(), a1.count(), a2.count()


fig = plt.figure()
plt.scatter(a0,a1)
plt.title('org vs weight')
plt.savefig('/home/utsumi/temp/ret/temp.org-wgt.png')

fig = plt.figure()
plt.scatter(a0,a2)
plt.title('org vs range')
plt.savefig('/home/utsumi/temp/ret/temp.org-rng.png')

fig = plt.figure()
plt.scatter(a1,a2)
plt.title('weight vs range')
plt.savefig('/home/utsumi/temp/ret/temp.wgt-rng.png')

