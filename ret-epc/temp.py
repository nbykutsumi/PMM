# %%
%matplotlib inline
import numpy as np
from numpy import *
import glob, h5py
import myfunc.util as util
from datetime import datetime, timedelta
import JPLDB
import pandas as pd
import matplotlib.pyplot as plt
# %%
srcdir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/test/glb.relsurf02.minrec1000.maxrec10000/2015/03/05'
convpath = srcdir + '/nsurfConvNScmb.005770.y-9999--9999.nrec10000.npy'
precpath = srcdir + '/nsurfNScmb.005770.y-9999--9999.nrec10000.npy'
aconv = np.load(convpath)
aprec = np.load(precpath)

#-- gprof --
gprpath = '/home/utsumi/mnt/lab_work/hk01/PMM/NASA/GPM.GMI/2A/V05/2015/03/05/2A.GPM.GMI.GPROF2017v1.20150305-S115403-E132637.005770.V05A.HDF5'
with h5py.File(gprpath,'r') as h:
    gconv = h['S1/convectivePrecipitation'][:]
    gprec = h['S1/surfacePrecipitation'][:]


aconv = ma.masked_less(ma.masked_where(aprec<0, aconv),0)
aprec = ma.masked_less(ma.masked_where(aprec<0, aprec),0)
gconv = ma.masked_less(ma.masked_where(aprec<0, gconv),0)
gprec = ma.masked_less(ma.masked_where(aprec<0, gprec),0)

afrac = ma.masked_where(aprec==0, aconv) / aprec
gfrac = ma.masked_where(gprec==0, gconv) / gprec
plt.clf()
plt.scatter(gfrac,afrac)

#

# %%
