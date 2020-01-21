import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import numpy as np

#srcDir ='/home/utsumi/mnt/lab_tank/utsumi/PMM/himawari/temp/temp.201707010220.tir.01'
#srcPath= srcDir + '/grid20.dat'
#a=np.fromfile(srcPath, 'float32').reshape(6000,6000)
#print srcPath
#print a.min(), a.max()
#plt.imshow(a);plt.colorbar()
#plt.savefig('/home/utsumi/temp/geo/grid20.png')
#plt.clf()

srcDir ='/home/utsumi/mnt/lab_tank/utsumi/PMM/himawari/temp/temp.201707010220.vis.01'
srcPath= srcDir + '/grid10.dat'
a=np.fromfile(srcPath, 'float32').reshape(12000,12000)
print srcPath
print a.min(), a.max()
plt.imshow(a);plt.colorbar()
plt.savefig('/home/utsumi/temp/geo/grid10.png')
plt.clf()


srcDir ='/home/utsumi/mnt/lab_tank/utsumi/PMM/himawari/temp/temp.201707010220.ext.01'
srcPath= srcDir + '/grid05.dat'
a=np.fromfile(srcPath, 'float32').reshape(24000,24000)
print srcPath
print a.min(), a.max()
plt.imshow(a);plt.colorbar()
plt.savefig('/home/utsumi/temp/geo/grid05.png')
plt.clf()




