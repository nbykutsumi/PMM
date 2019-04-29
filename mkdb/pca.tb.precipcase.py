import numpy as np
import h5py
from numpy import *
from sklearn.decomposition import IncrementalPCA
import myfunc.util as util
from datetime import datetime, timedelta
import glob
import os,sys

iDTime = datetime(2017,1,1)
eDTime = datetime(2017,1,31)
dDTime = timedelta(days=1)
lDTime = util.ret_lDTime(iDTime,eDTime,dDTime)
baseDir = '/mnt/j/PMM/MATCH.GMI.V05A'
miss   = -9999.

#-- Read mean and std for normalization ---
coefDir = '/mnt/j/ret-ml/coef'
tc1meanPath = coefDir + '/mean.S1.ABp103-117.GMI.Tc.npy'
tc2meanPath = coefDir + '/mean.S1.ABp103-117.GMI.TcS2.npy'
tc1stdPath = coefDir + '/std.S1.ABp103-117.GMI.Tc.npy'
tc2stdPath = coefDir + '/std.S1.ABp103-117.GMI.TcS2.npy'

atc1mean = np.load(tc1meanPath)
atc2mean = np.load(tc2meanPath)
atc1std  = np.load(tc1stdPath)
atc2std  = np.load(tc2stdPath)
print atc1mean.shape
atcmean  = concatenate([atc1mean,atc2mean],axis=0)
atcstd   = concatenate([atc1std, atc2std],axis=0)
