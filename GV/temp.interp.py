import matplotlib
matplotlib.use('Agg')
from numpy import *
from gv_fsub import *
import matplotlib.pyplot as plt

Lat = arange(1,10)
Lon = arange(1,10)

a1lat = array([1.5, 5.1, 8.3])

a1lon = array([2.1, 7.9, 3.3])

a1dat = array([1, 3, 6])


a2out = gv_fsub.point2map(a1dat, a1lon, a1lat, Lon, Lat, 3).T

plt.imshow(a2out, origin='lower')
plt.colorbar()

outDir  = '/work/a01/utsumi/GPMGV/fig'
outPath = outDir + '/temp.png'
plt.savefig(outPath)
print a2out
print outPath
