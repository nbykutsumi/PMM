import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import numpy as np


srcPath = '/tank/utsumi/validprof/const/orog.meter.sp.one.180x360.npy'

a2orog = np.load(srcPath)
a2orog = ma.masked_equal(a2orog,0)
plt.imshow(a2orog,origin='lower')
plt.colorbar()
figPath ='/tank/utsumi/validprof/const/orog.meter.sp.one.180x360.png'
plt.savefig(figPath)
print figPath





