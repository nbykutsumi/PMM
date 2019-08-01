import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from numpy import *

srcDir  = '/tank/utsumi/validprof/const'
srcPath = srcDir + '/landfrac.37SN.sa.1200x3600'

a2in = np.fromfile(srcPath,'float32').reshape(1200,3600)

nyout,nxout = 120,360
a2out= ones([nyout,nxout],float32)*(-9999.)
for iy in range(nyout):
    for ix in range(nxout):
        a2out[iy,ix] = a2in[iy*10:(iy+1)*10, ix*10:(ix+1)*10].mean()


a2outsa = a2out
a2outsp = np.concatenate([a2out[:,180:], a2out[:,:180]],axis=1)

#--- Save ----
outsaPath = srcDir + '/landfrac.sa.one.120x360.npy'
np.save(outsaPath, a2outsa)

outspPath = srcDir + '/landfrac.sp.one.120x360.npy'
np.save(outspPath, a2outsp)

#--- Figure ---
plt.imshow(a2outsa,origin='lower')
plt.colorbar(orientation='horizontal')
plt.savefig(outsaPath+'.png')
plt.clf()

plt.imshow(a2outsp,origin='lower')
plt.colorbar(orientation='horizontal')
plt.savefig(outspPath+'.png')
plt.clf()

print outspPath
