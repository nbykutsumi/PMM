from numpy import *
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


srcDir = '/tank/utsumi/validprof/pair.gprof/2017/01/01'
aprof1 = ma.masked_less(np.load(srcDir + '/profpmw.016155.npy'),0)
aprof2 = ma.masked_less(np.load(srcDir + '/profrad.016155.npy'),0)
asurf1 = ma.masked_less(np.load(srcDir + '/precpmw.016155.npy'),0)
asurf2 = ma.masked_less(np.load(srcDir + '/precrad.016155.npy'),0)

#i   = 1300  # interesting for oid=016155
i = 1350
a1y = range(36)

vmax1 = aprof1[i].max()
vmax2 = aprof2[i].max()
vmax  = max(vmax1, vmax2, asurf1[i], asurf2[i])

fig = plt.figure(figsize=(8,5))
ax  = fig.add_axes([0.1,0.1,0.3,0.8])
ax.plot(aprof1[i], a1y, 'o')
ax.plot(asurf1[i], 0, 'o')
plt.xlim([0,vmax])
plt.title('GPROF')

ax  = fig.add_axes([0.5,0.1,0.3,0.8])
ax.plot(aprof2[i], a1y, 'o')
ax.plot(asurf2[i], 0, 'o')
plt.xlim([0,vmax])
plt.title('CMB')

figPath= '/home/utsumi/temp/validprof/prof.png'
plt.savefig(figPath)
print figPath

print asurf1[i], asurf2[i]
print aprof2[i]
