# %%
import matplotlib
matplotlib.use('Agg')
%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
from numpy import *
import glob
import os

srcdir = '/home/utsumi/bin/PMM/ret-epc'

precpath = srcdir + '/mrms.PRECIPRATE.GC.20180123.062000.30242.dat'
f=open(precpath,'r'); lines=f.readlines(); f.close()
l=[]
for line in lines[6:]:
    line = map(float, line.split())
    l.append(line)

a=np.array(l)
print a
a=ma.masked_less(a,0)

plt.imshow(a)
plt.colorbar()
plt.show()
# %%
