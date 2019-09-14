import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from numpy import *

rettype='gprof'
srcPath ='/home/utsumi/mnt/lab_tank/utsumi/validprof/map-uncond/%s/prec.sum.201406.sp.one.npy'%(rettype)
a = ma.masked_less_equal(np.load(srcPath), 0)
plt.imshow(a,origin='lower')
plt.colorbar()
figPath = '/home/utsumi/temp/validprof/temp1.%s.png'%(rettype)
plt.savefig(figPath)
plt.clf()
print figPath


srcPath ='/home/utsumi/mnt/lab_tank/utsumi/validprof/map-uncond/%s/prec.num.201406.sp.one.npy'%(rettype)
a = ma.masked_less_equal(np.load(srcPath), 0)
plt.imshow(a,origin='lower')
plt.colorbar()
figPath = '/home/utsumi/temp/validprof/temp2.%s.png'%(rettype)
plt.savefig(figPath)
plt.clf()
print figPath
