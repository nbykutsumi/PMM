# %%
import matplotlib
matplotlib.use('Agg')
#%matplotlib inline
import numpy as np
from numpy import *
import glob, h5py
import myfunc.util as util
from datetime import datetime, timedelta
import JPLDB
import pandas as pd
import matplotlib.pyplot as plt
srcdir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/glb.relsurf01.minrec1000.maxrec10000/2014/06/01'
precpath = srcdir + '/nsurfNScmb.001454.y-9999--9999.nrec10000.npy'
aprec1 = ma.masked_less(np.load(precpath), 0)

srcdir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/glb.relsurf02.minrec1000.maxrec10000/2014/06/01'
precpath = srcdir + '/nsurfNScmb.001454.y-9999--9999.nrec10000.npy'
aprec2 = ma.masked_less(np.load(precpath), 0)

plt.scatter(aprec1,aprec2)
plt.savefig('/home/utsumi/temp/ret/temp.png')
# %%
