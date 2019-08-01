import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import numpy as np


srcDir = '/tank/utsumi/validprof/pr.vs.metrics'
lprec = arange(0.25, 20+0.01, 0.5)

fig = plt.figure()
ax  = fig.add_axes([0.1,0.1,0.8,0.8])
for biaslev in [0,1,2,3]:
    for sgn in [-1,1]:
        s = np.load(srcDir + '/taylor.sum.bias=%dx%d.2017.01.npy'%(sgn,biaslev))
        n = np.load(srcDir + '/taylor.num.bias=%dx%d.2017.01.npy'%(sgn,biaslev))
        m = ma.masked_invalid(s/n)
    
        ax.plot(lprec,m, label=biaslev)

plt.legend()
plt.savefig('/home/utsumi/temp/validprof/temp.png')

