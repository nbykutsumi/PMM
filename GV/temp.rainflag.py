from numpy import *
import numpy as np

baseDir = '/home/utsumi/mnt/wellshare'
srcPath = baseDir + '/GPMGV/L2A25/FLORIDA-STJ/200504/rainFlag.20050429153447-20050429153528.42483.npy'
aflag = np.load(srcPath)
print aflag

