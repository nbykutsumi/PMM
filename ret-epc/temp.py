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
srcdir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/test/2014/06/01/'
convpath = srcdir + '/nsurfConvNScmb.001453.y-9999--9999.nrec10000.npy'
precpath = srcdir + '/nsurfNScmb.001453.y-9999--9999.nrec10000.npy'
aconv = np.load(convpath)
aprec = np.load(precpath)

#-- gprof --
gprpath = '/home/utsumi/mnt/lab_work/hk01/PMM/NASA/GPM.GMI/2A/V05/2014/06/01/2A.GPM.GMI.GPROF2017v1.20140601-S010534-E023807.001453.V05A.HDF5'
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
plt.scatter(gprec,aprec)

#

# %%
