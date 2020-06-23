# %%
import matplotlib
matplotlib.use('Agg')
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
import os, sys
%matplotlib inline

epcpath = '/home/utsumi/mnt/lab_tank/utsumi/PMM/validprof/pair/gprof-shift/2014/06/01/epc-nsurfMScmb.001454.npy'
gprpath = '/home/utsumi/mnt/lab_tank/utsumi/PMM/validprof/pair/gprof-shift/2014/06/01/precpmw.001454.npy'

a=np.load(epcpath)
b=np.load(gprpath)

a=ma.masked_less(a,0.5)
b=ma.masked_less(b,0.5)

plt.scatter(b,a)
plt.show()

# %%
