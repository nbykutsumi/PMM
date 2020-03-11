# %%
import matplotlib
import matplotlib.pyplot as plt
from numpy import *
import numpy as np
import myfunc.util as util
import h5py
import socket

hostname = socket.gethostname()
if hostname =="shui":
    workDir   = "/work"
elif hostname=="well":
    workDir   = '/home/utsumi/mnt/lab_work'

iYM = [1997,12]
eYM = [2015,3]
lYM = util.ret_lYM(iYM, eYM)

srcDir = workDir + '/hk02/PMM/NASA/TRMM.PR/3A.PR.M/06A'
srcPath = workDir + '/hk02/PMM/NASA/TRMM.PR/3A.PR.M/06A/2005/GPMTRM_KUR_0504_M_L3S_P3M_06A.h5'

with h5py.File(srcPath,'r'):
    acount = h['Grids/G2/heightStormTop/count'][:]
    amean  = h['Grids/G2/heightStormTop/mean'][:]
    astdev = h['Grids/G2/heightStormTop/stdev'][:]

print acount
# %%
